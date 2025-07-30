"""
Copyright 2023 Wayne Crawford
"""
from .match_one_str import match_one_str
# from ..logger import init_logger
# logger = init_logger()

def get_full_id(match_str, stream):
    r"""
    Return stream trace's channel seed_id matching match_str
    
    Args:
        match_str (str): string to match (may have '*' and '?' wildcards)
        stream (:class:`obspy.core.Stream`): stream
    """
    return match_one_str(match_str, [x.get_id() for x in stream],
                         "match_str", "stream_ids")
    