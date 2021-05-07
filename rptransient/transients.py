#!env python3
""" Model and remove transients from BBOBS data"""

#################################################
# import obspy
from .periodic_transient import PeriodicTransient
from .utils import prep_filter as prep_filt

def_mag_limit = 5.85
def_days_per_magnitude = 1.5


class Transients():
    """
    Periodic Transients class
    """
    def __init__(self, transients=[]):
        """
        :param transients: list of PeriodicTransient
        """
        self.transients = transients
        for transient in transients:
            assert isinstance(transient, PeriodicTransient)

    def calc_timing(self, trace, eq_template, prep_filter=True):
        """
        Calculate transient time parameters

        :param trace: data
        :param eq_template: EQTemplate() object
        :param prep_filter: apply Wiedland prep filter before processing?
        """
        if prep_filter:
            inp = prep_filt(trace)
        else:
            inp = trace.copy()
        for transient in self.transients:
            self._print_announce(f'Finding {transient} times')
            transient.calc_timing(inp, eq_template)

    def calc_transients(self, trace, eq_template, plot=False,
                        prep_filter=True):
        """
        Calculate transient

        :param trace: data
        :param eq_template: EQTemplate() object
        :param plot: make information plots
        :param prep_filter: apply Wiedland prep filter before processing?
       """
        if prep_filter:
            inp = prep_filt(trace)
        else:
            inp = trace.copy()
        for transient in self.transients:
            self._print_announce(f'Fitting {transient}')
            transient.calc_transient(inp, eq_template, plots=plot)

    def remove_transients(self, trace, match=True, plot=False,
                          prep_filter=False):
        """
        Remove transient from data

        :param stream: data (trace?)
        :param match: individually match transients to data?
        :param plot: make information plots
        :param prep_filter: apply Wiedland prep filter before processing?
        """
        if self.transients[0].transient_model is None:
            print(f'{self.transient_model=}, did you run calc_transients?')
            return None
        if prep_filter:
            out = prep_filt(trace)
        else:
            out = trace.copy()
        for transient in self.transients:
            self._print_announce(f'Removing {transient}')
            out = transient.remove_transient(out, match, plots=plot)
        return out

    @staticmethod
    def _print_announce(text):
        print('='*75)
        print(text)
        print('='*75)


#########################################################################
# if __name__ == "__main__":
# 	sys.exit(main())
