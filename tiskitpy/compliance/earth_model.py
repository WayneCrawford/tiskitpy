"""
1D Earth model
Author:  W Crawford
"""
import warnings

import numpy as np
from matplotlib import pyplot as plt

from .compliance import Compliance


class EarthModel1D():
    def __init__(self, prop_list):
        """
        Args:
            prop_list (list of lists): 1D model, each row is a
                layer with values [thick(m), rho(kg/m^3), vp(m/s), vs(m/s)]
        """
        if not isinstance(prop_list, (list, tuple)):
            raise ValueError('prop_list is not a list or tuple')
        for row in prop_list:
            if not isinstance(row, (list, tuple)):
                raise ValueError(f'{row=} of prop_list is not a list or tuple')
            if not len(row) == 4:
                raise ValueError(f'{row=} of prop_list does not have for elements')
            for c in row:
                if not isinstance(c, (float, int)):
                    raise ValueError(f'element={c} of {row=} of prop_list is not a number')
                if c < 0:
                    raise ValueError(f'element={c} of {row=} of prop_list is less than zero')
                elif c > 10000:
                    raise ValueError(f'element={c} of {row=} of prop_list is > 10000')
        self.thicks = np.array([x[0] for x in prop_list])
        self.rhos = np.array([x[1] for x in prop_list])
        self.vps = np.array([x[2] for x in prop_list])
        self.vss = np.array([x[3] for x in prop_list])
        for vs, vp in zip(self.vss, self.vps):
            if vs*np.sqrt(3) > vp:
                warnings.warn(f'{vs=} > {vp=} / sqrt(3)')

    def __str__(self):
        s = '<EarthModel1D>:\n'
        s += '        thickness (m) | rho (kg/m^3) | Vp (m/s) | Vs (m/s) | Vs/Vp  | lambda (GPa) | mu (GPa)\n'
        s += '        ------------- | ------------ | -------- | -------- | ------ | ------------ | ----------\n'
        for t, r, vp, vs in zip(self.thicks.tolist(), self.rhos.tolist(),
                              self.vps.tolist(), self.vss.tolist()):
            s += f"         {t:12.0f} | {r:12.0f} | {vp:8.0f} | {vs:8.0f} | {vs/vp:6.3f} | {r*(vp*vp-2*vs*vs)/1e9:12.1f} | {r*vs*vs/1e9:8.1f}\n"
        return s

    def plot(self):
        prev_depth = 0
        depths, vps, vss, rhos, mus, comps = [], [], [], [], [], []
        for thick, rho, vp, vs in zip(self.thicks, self.rhos, self.vps, self.vss):
            for i in (1, 2):  # add tops and bottoms
                # Add layer tops
                depths.append(prev_depth)
                vps.append(vp)
                vss.append(vs)
                rhos.append(rho)
                mus.append(rho*vs*vs)
                comps.append(rho*vp*vp-2*rho*vs*vs)
                prev_depth += thick
            prev_depth -= thick
        fig, axs = plt.subplots(1, 2, sharey=True)
        axs[0].plot(vps, depths, 'r', label='vp')
        axs[0].plot(vss, depths, 'b', label='vs')
        axs[0].plot(rhos, depths, 'g', label='rho')
        axs[0].set_xlabel('velocity (m/s) or density (kg/m^3)')
        axs[0].set_ylabel('Depth (m))')
        axs[0].set_title('EarthModel1D')
        axs[0].invert_yaxis()
        axs[0].legend()
        axs[1].plot([x/1.e9 for x in comps], depths, 'r', label='lambda')
        axs[1].plot([x/1.e9 for x in mus], depths, 'b', label='mu')
        axs[1].set_xlabel('Lame parameter (GPa)')
        axs[1].legend()
        plt.show()

    def calc_ncompl(self, f, wdepth):
        """
        Return the normalized compliance for the given frequencies and water depth

        Wrapper for :meth:`tiskitpy.Compliance.calc_norm_compliance`
        """
        return Compliance.calc_norm_compliance(wdepth, f, self)


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
