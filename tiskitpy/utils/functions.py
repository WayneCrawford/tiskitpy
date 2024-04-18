"""
Copyright 2023 Wayne Crawford
"""
import fnmatch  # Allows Unix filename pattern matching

import numpy as np
from ..logger import init_logger

logger = init_logger()

def stream_unmask(stream):
    """
    Check if a stream is masked and, if so, unmask it
    Interpolates data in gaps
        """
    if np.any([np.ma.count_masked(tr.data) for tr in stream]):
        logger.warning('Unmasking masked data (usually a gap or overlap)')
        return stream.split().merge(fill_value='interpolate')
    return stream

def get_full_id(match_str, stream):
    """
    Return stream trace's channel seed_id matching match_str
    
    Args:
        match_str (str): string to match (may have \* and \? wildcards)
        stream (:class:`obspy.core.Stream`): stream
    """
    return match_one_str(match_str, [x.get_id() for x in stream],
                          "match_str", "stream_ids")
    

def match_one_str(one_str, str_list, one_str_name, str_list_name):
    """
    Match one_str, which may have wildcards, with str_list

    Args:
        one_str (str): string to match
        str_list (list): list of strings to be matched
        one_str_name (str): name to give one_str in error messages
        str_list_name (str): name to give str_list in error messages
    Returns:
        matched_str (str): str_list item that matches one_str

    Raises:
        TypeError of one_str is not a str
        ValueError if there is not exactly one match
    """

    if not isinstance(one_str, str):
        raise TypeError(f"Error: {one_str_name} is not a str")
    if not isinstance(str_list, list):
        raise TypeError(f"Error: {str_list_name} is not a list")
    for s in str_list:
        if not isinstance(s, str):
            raise TypeError(f"Error: {str_list_name} element is not a str")
    matches = fnmatch.filter(str_list, one_str)
    if len(matches) == 0:
        raise ValueError('No match for {}={} in {}={}'.format(
                         one_str_name, one_str, str_list_name, str_list))
    elif len(matches) > 1:
        raise ValueError('Multiple matches for {}={} in {}={}'.format(
                         one_str_name, one_str, str_list_name, str_list))
    return matches[0]
