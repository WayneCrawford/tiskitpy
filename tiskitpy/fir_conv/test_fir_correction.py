# test_fir_correction.py
import numpy as np
from fir_correction import FIRCorrection


def test_root_classification():
    lin_fir = [0.25, -0.5, 1, -0.5, 0.25]
    fc = FIRCorrection.from_zeros(lin_fir, sps=100, timetag=0)
    z = fc.z

    assert len(z) == 4
    assert len(fc.z_max) + len(fc.z_min) + len(fc.z_UC) == 4
    assert np.all(np.poly(np.roots(fc.firzeros)) == fc.firzeros)
    assert np.all(np.poly(fc.z) == fc.firzeros)


# def test_minimum_phase_filter_normalization():
#     fir = [1, -0.5, 0.25]
#     fc = FIRCorrection.from_zeros(fir, sps=100, timetag=0)
#     fir_mp = fc.compute_minimum_phase_filter()
# 
#     # Minimum-phase filter should sum to 1 after normalization
#     assert np.isclose(np.sum(fir_mp), 1.0)

def test_AR_MA_shapes():
    lin_fir = [0.25, -0.5, 1, -0.5, 0.25]
    fc = FIRCorrection.from_zeros(lin_fir, sps=100, timetag=0)
    a, b = fc._compute_AR_MA(fc.z)

    assert len(b) == len(a) + 1
    assert np.isfinite(a).all()
    assert np.isfinite(b).all()
