"""Extended physics studies and numerical certification."""
from pathlib import Path
import numpy as np
from optics_portfolio.fiber import grid, sech, propagate, diagnostics, relative_error
from optics_portfolio.validation import certify_fiber, gaussian_exact, sampling_diagnostics
from optics_portfolio.output import result_dir, save_json, plt
from optics_portfolio.profiles import profile


def main(name="standard"):
    cfg = profile(name)
    out = result_dir(Path(__file__).parent)
    out.mkdir(parents=True, exist_ok=True)
    certificate = certify_fiber()
    t, w = grid(1024, 40)
    z, linear = propagate(np.exp(-t*t/2), t, 3, 120, g=0, save_every=4)
    width = diagnostics(linear, t, g=0)["rms_width"]
    _, spm = propagate(np.exp(-t*t/2), t, 3, 120, d=0, g=2, save_every=4)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), layout="constrained")
    axes[0, 0].plot(z, width/width[0], label="SSFM")
    axes[0, 0].plot(z, np.sqrt(1+z*z), "--", label="Analytical Gaussian")
    axes[0, 0].set(xlabel="Distance xi", ylabel="Width / input width", title="Dispersion-only broadening")
    axes[0, 0].legend()
    omega = np.fft.fftshift(w)
    for i in [0, 10, 30]:
        axes[0, 1].plot(t, abs(spm[i])**2, label=f"xi={z[i]:g}")
        spectral = abs(np.fft.fftshift(np.fft.fft(spm[i])))**2
        axes[1, 0].plot(omega, spectral/spectral.sum(), label=f"xi={z[i]:g}")
    axes[0, 1].set(xlim=(-4,4), xlabel="Retarded time", ylabel="Intensity", title="SPM keeps temporal intensity")
    axes[0, 1].legend()
    axes[1, 0].set(xlim=(-12,12), xlabel="Angular frequency", ylabel="Normalized spectral power", title="SPM changes the spectrum")
    axes[1, 0].legend()
    phase_mask=abs(t)<4
    # Unwrap only where the input is non-negligible; numerical tail phase is undefined.
    axes[1, 1].plot(t[phase_mask],np.unwrap(np.angle(spm[-1,phase_mask]/spm[0,phase_mask])),label="SSFM")
    axes[1, 1].plot(t[phase_mask],6*np.exp(-t[phase_mask]**2),"--",label="Exact nonlinear phase")
    axes[1, 1].set(xlim=(-4,4), xlabel="Retarded time", ylabel="Phase (rad)", title="Nonlinear phase imprint")
    axes[1, 1].legend()
    fig.savefig(out/"dispersion_spm.png")
    plt.close(fig)
    # Check higher-order recurrence at pi/2 with independent step and grid refinements.
    recurrence = []
    for order in [2, 3]:
        rows = []
        for n, steps in [(1024, cfg["fiber_steps"]), (2048, 2*cfg["fiber_steps"])]:
            tt, _ = grid(n, 40)
            _, uu = propagate(sech(tt), tt, np.pi/2, steps, g=order**2, save_every=steps)
            rows.append({"grid": n, "steps": steps,
                         "relative_amplitude_recurrence_error": relative_error(abs(uu[-1]), sech(tt)),
                         **sampling_diagnostics(uu)})
        recurrence.append({"N": order, "refinement": rows})
    # Separated moving sech pulses are an approximate initial two-soliton state.
    tc, _ = grid(1024, 64)
    left = sech(tc+6)*np.exp(1j*1.5*tc)
    right = sech(tc-6)*np.exp(-1j*1.5*tc)
    zc, collision = propagate(left+right, tc, 8, 2*cfg["fiber_steps"], save_every=max(1,cfg["fiber_steps"]//120))
    _, left_alone = propagate(left, tc, 8, 2*cfg["fiber_steps"], save_every=2*cfg["fiber_steps"])
    _, right_alone = propagate(right, tc, 8, 2*cfg["fiber_steps"], save_every=2*cfg["fiber_steps"])
    fig, axes = plt.subplots(1, 2, figsize=(11,4), layout="constrained")
    im = axes[0].pcolormesh(tc, zc, abs(collision)**2, shading="auto", cmap="magma")
    axes[0].set(xlim=(-12,12), xlabel="Retarded time", ylabel="Distance xi", title="Two moving pulses collide")
    fig.colorbar(im, ax=axes[0], label="Intensity")
    axes[1].plot(tc, abs(collision[-1])**2, label="Interacting pulses")
    axes[1].plot(tc, abs(left_alone[-1])**2+abs(right_alone[-1])**2, "--", label="Independent propagation")
    axes[1].set(xlim=(-12,12), xlabel="Retarded time", ylabel="Intensity", title="Outgoing position shifts")
    axes[1].legend()
    fig.savefig(out/"soliton_collision.png")
    plt.close(fig)
    # Central-window common samples distinguish resolution and domain effects.
    refinement = []
    ref_t, _ = grid(2048, 48)
    ref_initial = np.exp(-ref_t*ref_t/2-0.5j*ref_t*ref_t)
    _, ref = propagate(ref_initial, ref_t, 0.5, 2400, g=2, save_every=2400)
    for n, window, steps in [(256,24,300),(512,24,600),(1024,48,1200)]:
        tt, _ = grid(n, window)
        _, uu = propagate(np.exp(-tt*tt/2-0.5j*tt*tt), tt, 0.5, steps, g=2, save_every=steps)
        # Grids share sample positions; integer selection avoids interpolation error.
        idx = np.rint((tt-ref_t[0])/(ref_t[1]-ref_t[0])).astype(int)
        refinement.append({"grid": n, "window": window, "steps": steps,
                           "relative_error": relative_error(uu[-1], ref[-1,idx]), **sampling_diagnostics(uu)})
    invariant = diagnostics(collision, tc)
    np.savez_compressed(out/"extended_data.npz", t_linear=t, z_linear=z, linear=linear, spm=spm,
                        t_collision=tc, z_collision=zc, collision=collision)
    metrics = {"profile":name,"certificate":certificate,"recurrence":recurrence,
               "linear_field_error":relative_error(linear, gaussian_exact(t,z)),
               "spm_max_intensity_change":float(np.max(abs(abs(spm)**2-abs(spm[0])**2))),
               "compression_refinement":refinement,
               "collision": {**sampling_diagnostics(collision),
                   "relative_energy_drift":float(np.max(abs(invariant["energy"]/invariant["energy"][0]-1))),
                   "hamiltonian_drift_scaled":float(np.max(abs(invariant["hamiltonian"]-invariant["hamiltonian"][0]))/max(abs(invariant["hamiltonian"][0]),1e-12)),
                   "note":"approximate separated-pulse initial condition; no exact two-soliton solution claim"}}
    save_json(out/"advanced_metrics.json",metrics)
    if not certificate["passed"]:
        raise RuntimeError("Fiber reference certification failed")
    print("Extended fiber study complete", flush=True)


if __name__ == "__main__":
    import os
    main(os.environ.get("OPTICS_PROFILE","standard"))
