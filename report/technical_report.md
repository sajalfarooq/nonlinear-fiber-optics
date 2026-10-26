# Nonlinear fiber optics technical report

## Research question and findings

This study examines whether a spectral split-step solver reproduces the separate linear and nonlinear limits of pulse propagation, whether its second-order numerical accuracy is observable, and how chirp and soliton interactions change pulse evolution. The reference certification passes all declared checks. Higher-order recurrence improves under combined grid/step refinement, and the collision preserves energy while producing outgoing position differences relative to independently propagated pulses.

All data are synthetic. No novelty of the NLSE or split-step method is claimed. The contribution is a reproducible investigation that joins analytical checks, error diagnostics and a parameter study.

## Model and derivation

Under a slowly varying envelope approximation in a retarded frame, a common lossless fiber model is

$$\partial_z A=-i\frac{\beta_2}{2}\partial_{TT}A+i\gamma|A|^2A.$$

Set $T=T_0\tau$, $z=L_D\xi$, $A=\sqrt{P_0}u$, and $L_D=T_0^2/|\beta_2|$. Multiplying the resulting envelope equation by i yields

$$iu_\xi+\frac d2u_{\tau\tau}+g|u|^2u=0,\qquad d=-\mathrm{sign}(\beta_2),\quad g=\gamma P_0L_D.$$

The nonlinear length is $L_{NL}=1/(\gamma P_0)$, so g is also $N^2=L_D/L_{NL}$. Holding input amplitude fixed and changing g is a normalized comparison; physical peak power changes if material and width are held fixed.

Fourier transformation sends the second derivative to $-\omega^2$. Linear evolution over h multiplies each Fourier mode by $\exp(-id\omega^2h/2)$. Local nonlinear evolution multiplies the field by $\exp(ig|u|^2h)$. Strang splitting applies half of the linear step, a full nonlinear step and the remaining linear half. The global propagation-step error is second order when the field and temporal spectrum are resolved.

The energy $E=\int |u|^2d\tau$ is conserved by each substep. The Hamiltonian $H=\int[d|u_\tau|^2/2-g|u|^4/2]d\tau$ is conserved by the continuous model, but not exactly by splitting. Energy near machine precision therefore cannot substitute for phase error, step refinement or Hamiltonian monitoring.

## Experimental protocol

The baseline includes N=1,2,3 soliton evolution and a 25-point chirp scan. The extended study separates dispersion and SPM, checks N=2,3 recurrence at xi=pi/2 on two resolutions, follows counterpropagating moving pulses in retarded time, and compares compressed fields on common temporal sample coordinates. Reference certification uses a 40-unit window and 1024 temporal samples. Temporal tails and high-frequency energy are reported separately.

For the linear Gaussian, $u=q^{-1/2}\exp[-\tau^2/(2q)]$ with $q=1+id\xi$. This checks an infinite-domain analytical field, rather than only another FFT implementation. SPM-only intensity remains fixed while its spectrum changes. For d=g=1 the analytical soliton is $\mathrm{sech}(\tau)e^{i\xi/2}$.

## Measured certification

| Check | Measured value | Passed |
|---|---:|---|
| fundamental_soliton_complex_field | 5.11958e-07 | True |
| gaussian_infinite_domain_solution | 1.06283e-14 | True |
| strang_convergence_order | 1.9998 | True |

## Higher-order recurrence

| Soliton order | Temporal samples | Propagation steps | Relative amplitude error |
|---|---:|---:|---:|
| 2 | 1024 | 2400 | 7.19345e-06 |
| 2 | 2048 | 4800 | 1.79844e-06 |
| 3 | 1024 | 2400 | 0.000114602 |
| 3 | 2048 | 4800 | 2.86693e-05 |

Reducing the step while increasing temporal resolution reduces the recurrence error. These are combined refinement results; they do not separately estimate temporal-grid and propagation-step errors. The independent certification study isolates propagation-step order.

## Collision and chirp interpretation

The collision input is a sum of separated moving sech pulses, not an exact analytical two-soliton solution. Comparison with independent propagation illustrates an interaction-induced position difference. A quantitative integrability claim would additionally require fitted outgoing pulse parameters and independent control of initial overlap and radiation.

The chirp scan's best sampled width occurs at the chirp search boundary. It must be described as a constrained scan result. Compression depends on which metric is used: RMS width is sensitive to temporal tails, while FWHM and peak intensity can favour different pulse shapes.

## Reproduction and limitations

From the repository root run `python reproduce.py --profile standard`. Use `--run-tag experiment1` for an independent output directory. The five notebooks provide a conceptual route through dispersion, SPM, solitons, certification and compression. `results/advanced_metrics.json` is the authoritative machine-readable record.

The FFT assumes a periodic window. Edge-energy and high-frequency diagnostics are warnings, not rigorous error bounds. The model excludes loss, Raman response, higher-order dispersion and self-steepening. The N=3 baseline figure is less refined than the dedicated recurrence calculation. These ideal simulations should not be described as experimental fiber measurements.

Reference: G. P. Agrawal, *Nonlinear Fiber Optics*, 6th edition (2019), Chapters 2–5. Implementations are independently written with AI assistance.
