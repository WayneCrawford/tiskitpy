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

    def select(self, **kwargs):
        """
        Adds selection by tiskitpy_id
        """
        out_stream =  super().select(**kwargs)
        if len(out_stream) == 0:
            logger.debug(f'Stream.select(**{kwargs}) returned no traces, '
                        f'tag with tiskitpy_ids and retry...')
            tagged = Stream(CleanSequence.seedid_tag(self))
            out_stream =  tagged.select(**kwargs)
            out_stream = CleanSequence.seedid_untag(out_stream)
        if len(out_stream) == 0:
            logger.debug(f'tagging with tiskit_py_ids returned no traces, '
                        f'try stripping tags from search terms')
            for key in kwargs.keys():
                if key == 'id':
                    subs = kwargs[key].split('.')  # individual id elements
                    kwargs[key] = '.'.join([x.split('-')[0] for x in subs])
                else:
                    kwargs[key] = kwargs[key].split('-')[0]
            out_stream =  super().select(**kwargs)
        return out_stream