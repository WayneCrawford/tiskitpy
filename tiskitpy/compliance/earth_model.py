"""
1D Earth model
Author:  W Crawford
"""
import warnings

import numpy as np


class EarthModel1D():
    def __init__(self, prop_list):
        """
        Args:
            prop_list (list of lists): 1D model, each row is a
                layer with values [thick(m), rho(kg/m^3), vp(m/s), vs(m/s)]
        """
        if not isinstance(prop_list, (list, tuple)):
            raise ValueError(f'prop_list is not a list or tuple')
        for row in prop_list:
            if not isinstance(row, (list, tuple)):
                raise ValueError(f'{row=} of prop_list is not a list or tuple')
            if not len(row) == 4:
                raise ValueError(f'{row=} of prop_list does not have for elements')
            for c in row:
                if not isinstance(c, (float, int)):
                    raise ValueError(f'{element=} of {row=} of prop_list is not a number')
                if c < 0:
                    raise ValueError(f'{element=} of {row=} of prop_list is less than zero')
                elif c > 10000:
                    raise ValueError(f'{element=} of {row=} of prop_list is > 10000')
        self.thicks = np.array([x[0] for x in prop_list])
        self.rhos = np.array([x[1] for x in prop_list])
        self.vps = np.array([x[2] for x in prop_list])
        self.vss = np.array([x[3] for x in prop_list])
        for vs, vp in zip(self.vss, self.vps):
            if vs*np.sqrt(3) > vp:
                warnings.warn(f'{vs=} > {vp=} / sqrt(3)')
        
    def __str__(self):
        s = '<EarthModel1D>:\n'
        s += '        thickness (m) | rho (kg/m^3) | Vp (m/s) | Vs (m/s)\n'
        s += '        ------------- | ------------ | -------- | ----------\n'
        for t, r, vp, vs in zip(self.thicks.tolist(), self.rhos.tolist(),
                              self.vps.tolist(), self.vss.tolist()):
            s += f"         {t:12.0f} | {r:12.0f} | {vp:8.0f} | {vs:8.0f}\n"
        return s


if __name__ == '__main__':
    print('='*60)
    print('Example Earth Model')
    print('='*60)
    model = EarthModel1D([[1000, 3000, 3000, 1600],
                          [1000, 3000, 4000, 2300],
                          [1000, 3000, 5000, 2800],
                          [3000, 3000, 7500, 4300],
                          [3000, 3000, 8200, 4700]])
    print(model)
