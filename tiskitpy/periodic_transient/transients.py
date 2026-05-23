#!env python3
""" Model and remove transients from BBOBS data"""

#################################################
# import obspy
from .periodic_transient import PeriodicTransient
from .utils import prep_filter as prep_filt
from ..logger import init_logger

logger = init_logger()
def_mag_limit = 5.85
def_days_per_magnitude = 1.5


class Transients():
    """
    Periodic Transients class

    calcualates and stores a list of PeriodicTransents
    """
    def __init__(self, transients=[]):
        """
        Args:
            transients (list of PeriodicTransient): transients
        """
        self.transients = transients
        for transient in transients:
            assert isinstance(transient, PeriodicTransient)

    def calc_timing(self, trace, eq_remover, prep_filter=True):
        """
        Calculate transient time parameters

        Args:
            trace (:class ~obspy.core.Trace): data
            eq_remover (:class ~EQRemover): periods to remove because of EQs
            prep_filter (bool): apply Wiedland prep filter before processing?
        """
        if prep_filter:
            inp = prep_filt(trace)
        else:
            inp = trace.copy()
        for transient in self.transients:
            self._print_announce(f'Finding {transient} times')
            transient.calc_timing(inp, eq_remover)

    def calc_transients(self, trace, eq_remover, plot=False,
                        prep_filter=True):
        """
        Calculate transients

        Args:
            trace (~class `obspy.core.trace.Trace`): data
            eq_remover (~class `.EQRemover`): time spans to zero out
            plot (bool): make information plots
            prep_filter (bool): apply Wieland prep filter before processing?
       """
        if prep_filter:
            inp = prep_filt(trace)
        else:
            inp = trace.copy()
        for transient in self.transients:
            self._print_announce(f'Fitting {transient}')
            transient.calc_transient(inp, eq_remover, plots=plot)

    def remove_transients(self, trace, match=True, plot=False,
                          prep_filter=False):
        """
        Remove transient from data

        Args:
            trace (~class `obspy.core.trace.Trace`): data
            match (bool): individually match transients to data?
            plot (bool): make information plots
            prep_filter (bool): apply Wieland prep filter before processing?
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
