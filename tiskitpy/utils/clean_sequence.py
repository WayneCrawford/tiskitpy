from obspy.core import Stream, Trace

from ..logger import init_logger
from .functions import get_full_id

logger = init_logger()
DOT_REPLACE_CHAR = '_'  # Character to replace id '.'s in cleaner strings


class CleanSequence:
    """
    Routines to store channel cleaning information

    In Trace and Stream objects, the information is stored in
    Trace.stats['clean_sequence']
    In SpectralDensity objects, it is stored in ...
    In DCRF objects, it is stored in...
    
    A string representation can be generated using the `string` method
    
    A string representation can be placed in / removed from a Trace/Stream's
    seed id using the  `seedid_tag()` and `seedid_untag()` methods
    """
    @staticmethod
    def tag(inp, clean_code, cleaned_ids=None, verbose=False):
        """
        Tag Trace stats with a clean_code

        Args:
            inp (:class:`obspy.core.Trace`, :class:`obspy.core.Stream`, a list,
                or None): the object to tag
            clean_code (str or list of str): transformation code ('ROT', ...)
                or the seed_id of the cleaned channel
                If a list, each str in the list will be appended
            cleaned_ids (list): a list of which Stream ids to apply to (can use
                wildcards). If None, applies to all of the Stream traces
        Returns:
            outp (:class:`obspy.core.Trace`, :class:`obspy.core.Stream`, or list):
                trace, stream, or list with modified clean_sequence
        Raises:
            ValueError if cleaned_ids is not None and trace_or_stream is a Trace
            TypeError if inp is not a Trace, Stream, list or NoneType
        """
        if inp is None:
            return CleanSequence.tag([], clean_code)
        if isinstance(inp, (Trace, list)):
            outp = inp.copy()
            if cleaned_ids is not None:
                raise ValueError('cleaned_ids is not None')
            if isinstance(clean_code, list):
                for x in clean_code:
                    outp = CleanSequence.tag(outp, x)
                return outp
            else:
                if not _is_seed_id(clean_code) and not _is_trans_code(clean_code):
                    raise ValueError('"{clean_code=}" is neither a seed_id nor a transformation code')
                if isinstance(inp, Trace):
                    if not 'clean_sequence' in outp.stats:
                        outp.stats['clean_sequence'] = [clean_code]
                    else:
                        outp.stats['clean_sequence'].append(clean_code)
                else:
                    outp.append(clean_code)
        elif isinstance(inp, Stream):
            outp = inp.copy()
            cleaned_ids = _expand_cleaned_ids(cleaned_ids, inp)
            for tr in outp:
                if tr.get_id() in cleaned_ids:
                    new_tr = CleanSequence.tag(tr, clean_code)
                    outp.remove(tr)
                    outp.append(new_tr)
        elif isinstance(inp, list):
            outp = inp.copy()
            if cleaned_ids is not None:
                raise ValueError('cleaned_ids is not None')
            for tr in outp:
                if tr.get_id() in cleaned_ids:
                    new_tr = CleanSequence.tag(tr, clean_code)
                    outp.remove(tr)
                    outp.append(new_tr)
        else:
            raise TypeError(f'inp ({type(inp)}) is neither Trace nor Stream')
        return outp

    @staticmethod
    def seedid_tag(tr_st, out_format='minimal'):
        """
        Add the CleanSequence's string representation to Trace(s) seed_id
        
        - This is useful for using obspy to plot each channel separately
        - This is problematic if you want to link to the inventory (to
          correct for the instrument response, for example)
        
        For example, if trace.stats.clean_sequence = [
            'XX.STA.00.BHX', 'XX.STA.00.BHY', 'XX.STA.00.BDH', 'XX.STA.01.BDH']
        and the input trace id is 'XX.STA.00.BHZ', the new trace id will be
        'XX.STA.00-X-Y-0_H-1_H.BHZ'
        """
        x = tr_st.copy()
        if isinstance(x, Stream):
            for tr in x:
                new_tr = CleanSequence.seedid_tag(tr, out_format)
                x.remove(tr)
                x.append(new_tr)
        elif isinstance(x, Trace):
            x = CleanSequence.seedid_untag(x)
            x.stats.location += CleanSequence.string(x, out_format)            
        else:
            raise TypeError('tr_st is neither a Trace nor a Stream')
        return x

    @staticmethod
    def seedid_untag(tr_st):
        """
        Remove the CleanSequence's string representation from the Trace(s) seed_id
        
        - This is useful for linking to the inventory (for example,
          for the instrument response)
        - This is problematic for using obspy.plot() (plots cleaned channels
          on same line as uncleaned)

        Args:
            tr_st (:class:`obspy.core.Trace` or `obspy.core.Trace`)
                Trace or Stream to remove clean_sequence info from the loc
        Returns:
            new_tr_st (:class:`obspy.core.Trace` or `obspy.core.Trace`)
        """
        x = tr_st.copy()
        if isinstance(x, Stream):
            for tr in x:
                new_tr = CleanSequence.seedid_untag(tr)
                x.remove(tr)
                x.append(new_tr)
        else:
            x.stats.location = x.stats.location.split('-')[0]            
        return x

    @staticmethod
    def string(inp, out_format='minimal'):
        """
        Returns the clean_sequence's string representation
        
        The form is "-" + id0 + "-" + id1 + ...
        
        By default, idx are minimized to contain the least characters
        uniquely identifying the seed id, from right to left

        Args:
            inp (:class:`obspy.core.Trace or `obspy.core.Stream):
                Trace or Stream from which to get the clean_sequence(s)
            out_format (str): how to format the output:
                'minimal': return the shortest string possible
                'min_code': use the shortest string containing full codes
                'full': return the full seedid string for each cleaned channel
        Returns:
            (str or list): depending on whether you input Trace or Stream
        """
        if isinstance(inp, Stream):
            return [CleanSequence.string(tr) for tr in inp]
        elif isinstance(inp, Trace):
            if not 'clean_sequence' in inp.stats:
                return ''
            return CleanSequence._clean_sequence_str(inp.stats.clean_sequence, out_format)
        else:
            raise TypeError('tr_st is neither a Trace nor a Stream')
                             
    @staticmethod
    def tiskitpy_id(seed_id, clean_sequence):
        """
        Returns a tiskitpy_id, based on a seed_id and a clean_sequence list
        
        Args:
            seed_id (str): a valid seed_id (has three '.'s)
            clean_sequence (list of str): ordered list of cleaning steps:
                either seed_ids or transformation codes
        Returns:
            (str): the tiskitpy_id
        """
        if not _is_seed_id(seed_id):
            try:
                seed_id = CleanSequence.seed_id(seed_id)
                logger.warning(f'{seed_id=} was not a valid seed_id, reduced to seed_id equivalent')
            except Exception:
                raise ValueError(f'{seed_id=} is neither a valid seed_id nor a valid tiskitpy_id')
        components = seed_id.split('.')
        components[2] += CleanSequence._clean_sequence_str(clean_sequence)
        return '.'.join(components)
                             
    @staticmethod
    def seed_id(tiskitpy_id):
        """
        Converts a tiskitpy_id to its base seed_id
        
        Args:
            tiskitpy_id (str): tiskitpy_id
        Returns:
            (str): the seed_id
        """
        if not _is_tiskitpy_id(tiskitpy_id):
            raise ValueError(f'{tiskitpy_id=} is not a valid tiskitpy_id')
        components = tiskitpy_id.split('.')
        components[2] = components[2].split('-')[0]
        return '.'.join(components)
                             
    def complete_seed_id(chan_id, verbose=False):
        """
        Change string to a seed_id
    
        Checks if there are 3 dots.  If not, adds "missing" dots to beginning
        """
        num_dots = chan_id.count('.')
        if num_dots == 3:
            new_seed_id = chan_id
        elif num_dots > 3:
            raise ValueError(f'More than 3 dots in input "{chan_id=}"')
        else:
            new_seed_id = (3-num_dots)*'.' + chan_id
            logger.info(f'Added {3 - num_dots} dots to start of seed_code: '
                        f'"{chan_id}"->"{new_seed_id}"')
        if not _is_seed_id(new_seed_id):
            raise ValueError('Created {new_seed_id=} is invalid')
        return new_seed_id

    @staticmethod
    def _clean_sequence_str(id_list, out_format='minimal'):
        """
        Return a string representing a clean_sequence list

        seed_ids must have three '.'s, transformation codes must have none
        neither can have '-'s or '_'s
        
        Args:
            id_list (list): list of cleaning steps (seed_id or transformation code)
            out_format (str): how to format the output:
                'minimal': return the shortest string possible (transformation
                    codes are kept in full)
                'min_code': use the shortest string containing full codes
                'min_level': use the shortest string containing all different levels
                'full': return the full seedid string for each cleaned channel

        Examples:
            >>> CleanSequence._clean_sequence_str(['NN.SSS.LL.BHZ', 'NN.SSS.LL.BH1', 'NN.SSS.LL.BDH'])
            '-Z-1-H'
            >>> CleanSequence._clean_sequence_str(['NN.SSS.LL.BHZ', 'NN.SSS.LL.BH1', 'NN.SSS.LL.BDH'], out_format='min_code')
            '-BHZ-BH1-BDH'
            >>> CleanSequence._clean_sequence_str(['NN.SSS.LL.BHZ', 'NN.SSS.LL.BH1', 'NN.SSS.LL.BDH'], out_format='full')
            '-NN_SSS_LL_BHZ-NN_SSS_LL_BH1-NN_SSS_LL_BDH'
            >>> CleanSequence._clean_sequence_str(['NN.SSS.LL.BHZ', 'OO.TTT.MM.BH1', 'PP.UUU.NN.BDH'])
            '-Z-1-H'
            >>> CleanSequence._clean_sequence_str(['NN.SSS.LL.BHZ', 'NN.SSS.LL.BH1', 'NN.SSS.LL.BLZ'])
            '-HZ-1-LZ'
            >>> CleanSequence._clean_sequence_str(['NN.SSS.LL.BHZ', 'NN.SSS.LL.BH1', 'NN.SSS.LL.BHZ'])
            ValueError: Could not create unique strs from id_list=['NN.SSS.LL.BHZ', 'NN.SSS.LL.BH1', 'NN.SSS.L2.BHZ']
            >>> CleanSequence._clean_sequence_str(['NN.SSS.LL.BHZ', 'NN.SSS.LL.BH1', 'NN.AAA.LL.BHZ'])
            '-S_L_Z-1-A_L_Z'
            >>> CleanSequence._clean_sequence_str(['NN.SSS.LL.BHZ', 'NN.SSS.LL.BH1', 'NN.AAA.LL.BHZ'], out_format='min_code')
            '-SSS_LL_BHZ-BH1-AAA_LL_BHZ'
            >>> CleanSequence._clean_sequence_str(['NN.SSS.LL.BHZ', 'NN.SSS.LL.BH1', 'NN.AAA.LL.BHZ'], out_format='min_level')
            '-SSS_LL_BHZ-SSS_LL_BH1-AAA_LL_BHZ'
        """
        # Error checking
        valid_out_formats = ['minimal', 'min_code', 'min_level', 'full']
        if not out_format in valid_out_formats:
            raise ValueError(f'{out_format=} not in {valid_out_formats}')
        if not isinstance(id_list, list):
            raise TypeError(f'id_list is type {type(id_list)}, not a list')
        if len(id_list) == 0:
            return ''
        for x in id_list:
            if not isinstance(x, str):
                raise TypeError('id_list element is not a str')
            if not ((len(x.split('.')) == 4) or ('.' not in x)):
                raise ValueError(f'"{x}" is neither a seed_id nor a transformation code')
        
        if out_format == 'full':    # Use full seed ids
            return '-' + '-'.join([DOT_REPLACE_CHAR.join(x.split('.')) for x in id_list])
        else:
            # Create shortened str describing each sequence element
            strs = [x if _is_trans_code(x) else None for x in id_list]
            unique_strs = [_is_trans_code(x) for x in id_list]
            fake_id_list = [x if _is_tiskitpy_id(x) else '...' for x in id_list]
            for i in range(3, -1, -1): # from channel to net, if needed
                strs, unique_strs = CleanSequence._prepend_unique_str(
                    [x.split('.')[i] for x in fake_id_list],
                    strs,
                    unique_strs,
                    out_format)
                if all(unique_strs):
                    return '-' + '-'.join(strs)
            raise ValueError(f'Could not create unique strs from {id_list=}')           
                    
    @staticmethod
    def _prepend_unique_str(codes, prevs, prev_uniques, out_format='minimal'):
        """
        Prepend a "minimum unique" string for each of the given codes.
        
        A "minimum unique" string is the shortest string, starting from the last
        character, which does not match an equivalent string from
        another code.
        
        Args:
            codes (list): code strings
            prevs (list of str): miniumum unique strings from the previous level.
                If None, then just return the present minimum string
            prev_uniques (list of bool): True if code at previous level was
                already unique.  In this case, prepend nothing.
            out_format (str): how to format the output:
                'minimal': return the shortest string possible, if previous level not unique
                'min_code': return the full subcode, if previous level not unique
                'min_level': return the full subcode
        Returns:
            short_strs (list of str): minimum strings from this level
            is_unique (list of bool): whether each this level's strings are unique
        Examples:
            >>> CleanSequence._prepend_unique_str(['BHZ', 'BH1', 'BDH'], [None, None, None], [False, False, False])
            (['Z', '1', 'H'], [True, True, True])
            >>> CleanSequence._prepend_unique_str(['BHZ', 'BH1', 'BDH'], [None, None, None], [False, False, False], 'min_code')
            (['BHZ', 'BH1', 'BDH'], [True, True, True])
            >>> CleanSequence._prepend_unique_str(['BHZ', 'BH1', 'BLZ'], [None, None, None], [False, False, False])
            (['HZ', '1', 'LZ'], [True, True, True])
            >>> CleanSequence._prepend_unique_str(['BHZ', 'BH1', 'BHZ'], [None, None, None], [False, False, False])
            (['Z', '1', 'Z'], , [False, True, False])
            >>> CleanSequence._prepend_unique_str(['00', '01', '02'], ['Z', '1', 'Z'], [False, True, False])
            (['0_Z', '1', '2_Z'], [True, True, True])
            >>> CleanSequence._prepend_unique_str(['00', '01', '02'], ['Z', '1', 'Z'], [False, True, False], 'min_code')
            (['00_Z', '1', '02_Z'], [True, True, True])
            >>> CleanSequence._prepend_unique_str(['00', '00', '00'], ['Z', '1', 'Z'], [False, True, False])
            (['0_Z', '1', '0_Z'], [False, True, False])           
            >>> CleanSequence._prepend_unique_str(['00', '00', '01'], ['Z', '1', 'Z'], [False, True, False])
            (['0_Z', '1', '1_Z'], [True, True, True])           
        """
        if out_format == 'minimal':
            base_codes = CleanSequence._min_codes(codes)
        else:
            base_codes = codes
        # Prepend minimum codes to previous minimum strings
        short_strs = []
        for code, prev, prev_uniq in zip(base_codes, prevs, prev_uniques):
            if prev is None:
                short_strs.append(code)
            elif prev_uniq is True and not out_format == 'min_level':
                short_strs.append(prev)
            else:
                short_strs.append(code + DOT_REPLACE_CHAR + prev)
        is_unique = [(short_strs.count(x) == 1) for x in short_strs] 
        return short_strs, is_unique

    @staticmethod
    def _min_codes(codes):
        """
        Return minimum unique back-loaded string for each code
        If there are two identical codes, returns the same value for each
        
        Examples:
            >>> CleanSequence._min_codes(['BHZ', 'BH1', 'BH2', 'BDH'])
            ['Z', '1', '2', 'H']
            >>> CleanSequence._min_codes(['BHZ', 'BH1', 'BH2', 'BLZ'])
            ['HZ', '1', '2', 'LZ']
            >>> CleanSequence._min_codes(['BHZ', 'BH1', 'BH2', 'BHZ'])
            ['Z', '1', '2', 'Z']
        """
        min_codes = []
        for code in codes:
            offset = 1
            other_codes = [x for x in codes if not x == code]
            while (code[-offset:]) in [x[-offset:] for x in other_codes]:
                offset += 1
                if offset > len(code):  # Don't run beyond length of code str
                    offset = len(code)
                    break
            min_codes.append(code[-offset:])
        return min_codes
    

def _expand_cleaned_ids(cleaned_ids, in_stream):
    """ Expand cleaned_ids """
    if cleaned_ids is None:
        return [tr.get_id() for tr in in_stream]
    return [get_full_id('*' + x, in_stream) for x in cleaned_ids]


def _is_tiskitpy_id(x):
    if (not len(x.split('.')) == 4) or  ('_' in x):
        return False
    return True

def _is_seed_id(x):
    if (not len(x.split('.')) == 4) or ('-' in x) or ('_' in x):
        return False
    return True

def _is_trans_code(x):
    if ('.' in x) or ('-' in x) or ('_' in x) or (len(x) > 5):
        return False
    return True