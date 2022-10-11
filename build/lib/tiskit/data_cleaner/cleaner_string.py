class CleanerString:
    """
    All routines for DataCleaner string creation and modification
    
    a CleanerString always starts with '-'
    """
    @staticmethod
    def make(in_chan, prev_remove_seq=""):
        """
        Generate a string describing the TransferFunction removal chain

        Args:
            in_chan: name of channel to remove
            prev_remove_seq (str): previous removal sequence string, used to
                avoid repeating an already-given channel code
        Returns:
            cleaner_str (str): string to add to end of removal sequence

        Example:
            >>> CleanerString.make('BHZ', '-X-Y')
            '-Z'
            >>> CleanerString.make('BHX', '-X-Y')
            '-HX'
            >>> CleanerString.make('X', '-X-Y')
            Traceback (most recent call last):
            ...
            ValueError: new_name('X') already in prev_remove_seq('-X-Y')
        """
        prev_components = prev_remove_seq.split("-")
        offset = 1
        while (new_name := in_chan[-offset:]) in prev_components:
            offset += 1
            if offset > len(in_chan):
                raise ValueError(
                    f"new_name('{new_name}') already in "
                    f"prev_remove_seq('{prev_remove_seq}')"
                )
        cleaner_str = "-" + new_name
        return cleaner_str

    @staticmethod
    def insert(cleaner_str, channel_name):
        """
        Insert cleaner string in current channel name

        Args:
            cleaner_str: cleaner string (made usine make_cleaner_str())
            channel_name (str): channel name to modify
        Returns:
            new_channel_name (str): modified channel name

        Example:
            >>> CleanerString.insert('-X-Y', 'BHZ'')
            '-X-Y.BHZ'
            >>> CleanerString.insert('-X-Y', 'XX.STA.00.BHZ'')
            'XX.STA.00-X-Y.BHZ'
        """
        if '.' in channel_name:
            elements = channel_name.split('.')
        else:
            elements = ['', channel_name]
        elements[-2] += cleaner_str
        
        return '.'.join(elements)

    @staticmethod
    def strip(in_str):
        """
        Strip transfer function removal string from a string

        Expects the cleaner string to be at the end of the next to last element,
        separated by '.'s (location code).  Also handles but complains
        about strings without '.', in which case it removes the cleaner_str
        from the end.  If there are only two elements and stripping leaves
        nothing in the first element, return without '.'

        Also handles lists of strings

        Example:
            >>> CleanerString.strip('00-1-2-H.BHZ')
            '00.BHZ'
            >>> CleanerString.strip('-1-2-H.BHZ')
            'BHZ'
            >>> CleanerString.strip('BHZ-1-2-H')
            'BHZ'
            >>> CleanerString.strip('BHZ')
            'BHZ'
        """
        if isinstance(in_str, list):
            return [CleanerString.strip(x) for x in in_str]
        if '.' in in_str:
            elements = in_str.split('.')
            had_loc = len(elements[-2])>0
            elements[-2] = elements[-2].split("-")[0]
            if len(elements)==2 and had_loc and len(elements[-2])==0:
                output = elements[-1]
            else:
                output = '.'.join(elements)
        else:
            output = in_str.split("-")[0]
            logger.warning("cleaner string was after channel, shouldn't be!")
        return output

#     @staticmethod
#     def strip_cleaner_str(in_str):
#         """
#         Strips transfer function removal string from a string
# 
#         Takes out everything after first '-'
# 
#         Also handles lists of strings
# 
#         Example:
#             >>> text = 'BHZ-1-2-H'
#             >>> strip_cleaner_str(text)
#             'BHZ'
#             >>> text = 'BHZ'
#             >>> strip_cleaner_str(text)
#             'BHZ'
#         """
#         print(f'in_str={in_str}')
#         if isinstance(in_str, list):
#             return [CleanerString.strip_cleaner_str(x) for x in in_str]
#         return in_str.split("-")[0]
# 
