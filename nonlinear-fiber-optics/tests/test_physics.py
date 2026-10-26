import numpy as np
import pytest
from optics_portfolio.fiber import grid,sech,propagate,relative_error
from optics_portfolio.validation import gaussian_exact,sampling_diagnostics,certify_fiber

def test_soliton_and_second_order_convergence():
    # Larger window isolates splitting error from the periodic sech-tail floor.
    t, _ = grid(1024, 40)
    errors = []
    for steps in [40, 80, 160]:
        z, u = propagate(sech(t), t, 1, steps, save_every=steps)
        errors.append(relative_error(u[-1], sech(t)*np.exp(0.5j)))
        assert abs(np.sum(abs(u[-1])**2)/np.sum(sech(t)**2)-1) < 1e-12
    assert errors[-1] < 1e-4
    assert 3.5 < errors[0]/errors[1] < 4.5

def test_dispersion_matches_exact_fourier_propagator():
    t, w = grid(256, 20)
    initial = np.exp(-t*t/2)
    actual = propagate(initial, t, 1.7, 21, d=-1, g=0)[1][-1]
    expected = np.fft.ifft(np.fft.fft(initial)*np.exp(0.5j*w*w*1.7))
    assert relative_error(actual, expected) < 1e-12

def test_spm_exact_solution():
    t, _ = grid(128, 16)
    initial = sech(t)
    actual = propagate(initial, t, 0.7, 13, d=0, g=2)[1][-1]
    assert relative_error(actual, initial*np.exp(1.4j*initial**2)) < 1e-12

def test_infinite_domain_gaussian_both_dispersion_signs():
    t,_=grid(512,40)
    for d in [-1.,1.]:
        result=propagate(np.exp(-t*t/2),t,2,31,d=d,g=0)[1][-1]
        assert relative_error(result,gaussian_exact(t,2,d))<1e-10

def test_symmetric_split_time_reversal():
    t,_=grid(256,24)
    initial=np.exp(-t*t/2+0.2j*t*t)
    final=propagate(initial,t,1,120,g=1.4,save_every=120)[1][-1]
    back=propagate(final.conj(),t,1,120,g=1.4,save_every=120)[1][-1]
    assert relative_error(back.conj(),initial)<1e-12

def test_n2_recurrence():
    t,_=grid(1024,40)
    final=propagate(sech(t),t,np.pi/2,1600,g=4,save_every=1600)[1][-1]
    assert relative_error(abs(final),sech(t))<5e-4

def test_reference_certificate():
    assert certify_fiber()["passed"]

@pytest.mark.parametrize("bad",[np.nan,np.inf])
def test_nonfinite_fiber_rejected(bad):
    t,_=grid(32,12)
    with pytest.raises(ValueError): propagate(sech(t),t,zmax=bad)

@pytest.mark.parametrize("kwargs",[{"steps":1.5},{"save_every":0},{"g":np.nan}])
def test_invalid_fiber_controls(kwargs):
    t,_=grid(32,12)
    with pytest.raises(ValueError): propagate(sech(t),t,**kwargs)

def test_sampling_guard_detects_edge_energy():
    diagnostic=sampling_diagnostics(np.ones(128))
    assert diagnostic["edge_energy_fraction"]>.1
    assert diagnostic["high_frequency_energy_fraction"]==0

