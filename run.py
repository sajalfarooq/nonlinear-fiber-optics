"""Run from any directory after pip install -e .."""
from pathlib import Path
import numpy as np
from optics_portfolio.fiber import grid, sech, propagate, diagnostics, relative_error
from optics_portfolio.output import setup, write_metrics, plt


def main():
    out = setup(Path(__file__).parent)
    t, w = grid()
    z, u = propagate(sech(t), t, zmax=3, steps=600, save_every=5)
    exact = np.exp(0.5j*z[:, None])*sech(t)
    diag = diagnostics(u, t)
    counts = np.array([40, 80, 160, 320])
    reference = propagate(np.exp(-t*t/2), t, 1.0, 5120)[1][-1]
    errors = [relative_error(propagate(np.exp(-t*t/2), t, 1.0, int(n))[1][-1], reference)
              for n in counts]
    slope = float(np.polyfit(np.log(1/counts), np.log(errors), 1)[0])
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), layout="constrained")
    for ax, n in zip(axes, [1, 2, 3]):
        zs, field = propagate(sech(t), t, np.pi, 1000, g=n*n, save_every=5)
        im = ax.pcolormesh(t, zs, abs(field)**2, shading="auto", cmap="magma")
        ax.set(xlim=(-6, 6), xlabel="Retarded time tau", ylabel="Distance xi", title=f"Soliton order N={n}")
        fig.colorbar(im, ax=ax, label="Normalized intensity")
    fig.savefig(out / "soliton_orders.png")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), layout="constrained")
    axes[0].loglog(1/counts, errors, "o-", label=f"Measured slope {slope:.3f}")
    axes[0].set(xlabel="Propagation step", ylabel="Relative complex-field error", title="Step-size convergence")
    axes[0].legend()
    axes[1].plot(z, diag["energy"]/diag["energy"][0]-1)
    axes[1].set(xlabel="Distance xi", ylabel="Relative energy drift", title="Unitary splitting diagnostic")
    fig.savefig(out / "validation.png")
    plt.close(fig)
    chirps = np.linspace(-2, 2, 25)
    minimum, distances = [], []
    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    for c in chirps:
        initial = np.exp(-t*t/2) * np.exp(0.5j*c*t*t)
        zs, field = propagate(initial, t, 2.0, 500, g=2, save_every=5)
        widths = diagnostics(field, t, g=2)["rms_width"]
        best = np.argmin(widths)
        minimum.append(float(widths[best]/widths[0]))
        distances.append(float(zs[best]))
        if np.isclose(c, [-2, 0, 2]).any():
            ax.plot(zs, widths/widths[0], label=f"Chirp C={c:g}")
    ax.set(xlabel="Distance xi", ylabel="RMS width / input width", title="Chirp and pulse compression")
    ax.legend()
    fig.savefig(out / "compression.png")
    plt.close(fig)
    np.savez_compressed(out / "study_data.npz", t=t, z=z, field=u, chirps=chirps,
                        minimum_width_ratio=minimum, compression_distance=distances)
    best = int(np.argmin(minimum))
    write_metrics(out, {"soliton_relative_field_error": relative_error(u, exact),
        "max_relative_energy_drift": float(np.max(abs(diag["energy"]/diag["energy"][0]-1))),
        "convergence_order": slope, "step_counts": counts.tolist(), "convergence_errors": errors,
        "chirp_scan": {"best_chirp": float(chirps[best]), "minimum_rms_width_ratio": minimum[best],
                       "distance": distances[best], "g": 2.0, "range": [-2, 2], "samples": 25}})
    print("Fiber study saved:", out)


if __name__ == "__main__":
    main()
