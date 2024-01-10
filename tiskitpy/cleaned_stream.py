from obspy.core import Stream

from .logger import init_logger
from .utils import CleanSequence

logger = init_logger()


class CleanedStream(Stream):
    """
    Stream subclass with added tag() & modified __str__() and plot() methods

    Simply injects Trace.stats.clean_sequence info into seed_id location codes
    before handing off to Stream's __str__() and plot() methods
    """
    def __str__(self, **kwargs):
        inp = Stream(CleanSequence.seedid_tag(self))
        return inp.__str__(**kwargs)

    def plot(self, **kwargs):
        inp = Stream(CleanSequence.seedid_tag(self))
        return inp.plot(**kwargs)

    def tag(self, the_tag, **kwargs):
        """Return object with tag added to clean_sequence
        
        Haven't figured out how to do it in place
        '"""
        return CleanSequence.tag(self, the_tag, **kwargs)
