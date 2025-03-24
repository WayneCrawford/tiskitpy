import random

import numpy as np
from obspy.core.stream import Trace
from obspy.core import UTCDateTime


class TideCoefficients():
    def __init__(self, coeffs_dict=None):
        """
        Args:
            coeffs_dict (dict): keys = component names, values =
                (period (hours), amplitude (?), phase).  If phase is None,
                will randomly generate for each trace.
        """
        if coeffs_dict is None:
            coeffs_dict = self._earth_tides_default()
        assert isinstance(coeffs_dict, dict)
        for v in coeffs_dict.values():
            assert isinstance(v, list)
            assert len(v) == 3
        self.coeffs_dict = coeffs_dict

    @staticmethod
    def _earth_tides_default():
        return {"M2": [12.421, 384.83, None],
                "K1": [23.934, 191.78, None],
                "S2": [12, 179.05, None],
                "O1": [25.819, 158.11, None],
                "N2": [12.658, 73.69, None],
                "P1": [24.066, 70.88, None],
                "K2": [11.967, 48.72, None],
                "Mf": [13.661*12, 40.36, None]}

    def make_trace(self, ref_trace):
        tide_trace = ref_trace.copy()
        tide_trace.data = np.zeros(len(tide_trace.data))
        for k, v in self.coeffs_dict.items():
            period = v[0]
            amp = v[1]
            phase = v[2]
            if phase is None:
                phase = 360 * random.random()
            tide_trace.data += np.sin(2 * np.pi * tide_trace.times() / (3600 * period) + np.radians(phase)) * amp
        return tide_trace


if __name__ == "__main__":
    # Show an example
    tc = TideCoefficients()
    ref_trace = Trace(np.zeros(86400),
                      header={'sampling_rate': 0.01,
                              'starttime': UTCDateTime('2024-01-01T00:00:00')})
    tr = tc.make_trace(ref_trace)
    tr.plot()
