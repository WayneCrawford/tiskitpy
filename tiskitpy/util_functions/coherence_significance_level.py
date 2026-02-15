"""
coherence significance level
"""
import numpy as np


def coherence_significance_level(n_windows, prob=0.95):
    """
    Definition: L_1(alpha, q) = sqrt(1-alpha**(1/q))

    where alpha = 1-prob and 2(q+1) = nwinds (degree of freedom)

    For nwinds >> 1, L1 ~ sqrt(1-alpha**(2/nwinds))
    For a 95% signif level this comes out to
        sqrt(1-.05**(2/nwinds)) for nwinds >> 1.
    I previously used sqrt(2/nwinds) for the 95% signif level (alpha=0.05),
    but L1 is much closer to sqrt(6/nwinds).

    Args:
        n_windows (int): number of windows
        prob (float): significance level (between 0 and 1)
    """
    assert prob >= 0 and prob <= 1
    if n_windows < 1:
        raise ValueError(f"Can't have than one window: {n_windows=}")
    alpha = 1 - prob
    q = n_windows/2 - 1
    if q < 0:
        return np.nan
    elif q == 0:
        return 1.
    return np.sqrt(1 - alpha ** (1. / q))
