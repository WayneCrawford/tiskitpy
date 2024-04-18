from ..logger import init_logger

logger = init_logger()


class CleanerString:
    """
    Routines to create, modify and strip data removal strings from a channel
    id string's location code.

    Each time a channel is cleaned, the name of the cleaning channel,
    preceded by "-", should be appended to the location code.
    
    For example, after removing the '1', '2' and 'H' channels, in that order,
    from a channel originally named XX.SSS..BHZ, the new name should be
    XX.SSS.-1-2-H.BHZ)
    """
    @staticmethod
    def make(removed_chan, prev_str="", full_str=False):
        """
        Generate a string describing the removed channel.
        Returns the shortest string, starting from the end of
        removed_chan, that does not repeat a previous removal string

        Args:
            removed_chan (str): name of channel that was removed (can be full
                or partial)
            prev_str (str): sequence of previous removals
            full_str (bool): don't try to shorten removed_chan
        Returns:
            cleaner_str (str): string to add to end of removal sequence

        Example:
            >>> CleanerString._make('BHZ', '-X-Y')
            '-Z'
            >>> CleanerString._make('BHZ', '-X-Y', full_str=True)
            '-BHZ'
            >>> CleanerString._make('BHX', '-X-Y')
            '-HX'
            >>> CleanerString._make('X', '-X-Y')
            Traceback (most recent call last):
            ...
            ValueError: new_name('X') already in prev_str('-X-Y')
        """
        if not len(prev_str) == 0:
            if not len(prev_str) > 1 and prev_str[0] == '-':
                raise ValueError(f'improper prev_str: "{prev_str}"')

        prev_components = prev_str.split("-")
        if full_str is True:
            offset = len(removed_chan)
        else:
            offset = 1
        while (new_name := removed_chan[-offset:]) in prev_components:
            offset += 1
            if offset > len(removed_chan):
                raise ValueError(
                    f"new_name('{new_name}') already in "
                    f"prev_str('{prev_str}')"
                )
        cleaner_str = "-" + new_name
        return cleaner_str

    @staticmethod
    def add(removed_chan, prev_str="", full_str=False):
        """
        CleanerString.make() + removed_chan

        Example:
            >>> CleanerString.make('BHZ', '-X-Y')
            '-X-Y-Z'
            >>> CleanerString.make('BHX', '-X-Y')
            '-X-Y-HX'
            >>> CleanerString.make('X', '-X-Y')
            Traceback (most recent call last):
            ...
            ValueError: new_name('X') already in prev_str('-X-Y')
        """
        return prev_str + CleanerString.make(removed_chan, prev_str, full_str)

    @staticmethod
    def insert(cleaner_str, orig_id):
        """
        Insert cleaner string into channel/id name

        Args:
            cleaner_str: cleaner string (must start with "-")
            orig_id (str): seed ID or channel code to modify
        Returns:
            new_seed_id (str): modified seed id

        Example:
            >>> CleanerString.insert('-X-Y', 'BHZ'')
            '-X-Y.BHZ'
            >>> CleanerString.insert('-X-Y', 'XX.STA.00.BHZ'')
            'XX.STA.00-X-Y.BHZ'
        """
        if len(cleaner_str) == 0:
            return orig_id
        if '.' in cleaner_str:
            raise ValueError('cleaner_str contains ".", did you misorder inputs?')
        if not cleaner_str[0] == '-':
            raise ValueError('cleaner_str does not start with "-"')
        if '.' in orig_id:
            elements = orig_id.split('.')
        else:
            elements = ['', orig_id]
        if cleaner_str in elements[-2]:
            raise ValueError(f"cleaner_str ('{cleaner_str}') already in "
                             f"orig_id.location ('{elements[-2]}')")
        elements[-2] += cleaner_str

        return '.'.join(elements)

    @staticmethod
    def insert_id(subtract_id, orig_id, full_str=False):
        """
        Return channel/id name with given channel subtracted

        Args:
            subtract_id (str): seed_id or channel code of channel that
                was subtracted
            orig_id (str): seed ID or channel code to modify
            full_str (bool): don't try to shorten subtract_id
        Returns:
            new_seed_id (str): modified seed id

        Example:
            >>> CleanerString.insert_id('BHY', 'BHZ')
            '-Y.BHZ'
            >>> CleanerString.insert_id('BHY', 'BHZ', full_str=True)
            '-BHY.BHZ'
            >>> CleanerString.insert_id('BHX', 'XX.STA.00.BHZ')
            'XX.STA.00-X.BHZ'
            >>> CleanerString.insert_id('BHY', 'XX.STA.00-X.BHZ')
            'XX.STA.00-X-Y.BHZ'
        """
        cleaner_str = CleanerString.extract(orig_id)
        new_cleaner_str = CleanerString.make(subtract_id, cleaner_str, full_str)
        return CleanerString.insert(new_cleaner_str, orig_id)

    @staticmethod
    def strip(in_str):
        """
        Strip frequency response function removal string from a string

        Expects the cleaner string to be at the end of the next to last
        element, separated by '.'s (location code).
        If there are only two elements and stripping leaves nothing in the
        first element, return without '.'

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
        return CleanerString._separate(in_str)[0]

    @staticmethod
    def extract(in_str):
        """
        Extract frequency response function removal string from a string

        Expects the cleaner string to be at the end of the next to last
        element, separated by '.'s (location code).

        Also handles lists of strings

        Example:
            >>> CleanerString.extract('00-1-2-H.BHZ')
            '-1-2-H'
            >>> CleanerString.extract('-1-2-H.BHZ')
            '-1-2-H'
            >>> CleanerString.extract('BHZ')
            ''
        """
        return CleanerString._separate(in_str)[1]

    @staticmethod
    def _separate(in_str):
        """
        Separates cleaner string and rest of seed id

        Expects the cleaner string to be at the end of the next to last
        element, separated by '.'s (location code).

        If there are only two elements and stripping leaves nothing in the
        first element, returns the last without '.' before it

        Also handles lists of strings
        
        Args:
            in_str (str): input seed id
        Returns:
            output (list):
                seed_id: seed_id without cleaner string
                cleaner_str: the cleaner string (starts with '-')
                

        Example:
            >>> CleanerString._separate('00-1-2-H.BHZ')
            ['00.BDH', '-1-2-H']
            >>> CleanerString._separate('-1-2-H.BHZ')
            ['BHZ', '-1-2-H']
            >>> CleanerString._separate('BHZ')
            ['BHZ', '']
            >>> CleanerString._separate(['XX.SSS.00-1-2-H.BHZ',
                                         'XX.SSS.00-1-2.BDH',
                                         'XX.SSS.00-1.BH2',
                                         'XX.SSS.00.BH1'])
            [['XX.SSS.00.BHZ', 'XX.SSS.00.BDH', 'XX.SSS.00.BH2', 'XX.SSS.00.BH1'],
             ['-1-2-H', '-1-2', '-1', '']]
        """
        if isinstance(in_str, list):
            temp =  [CleanerString._separate(x) for x in in_str]
            return [[x[0] for x in temp], [x[1] for x in temp]]
        if '.' in in_str:
            elements = in_str.split('.')
            had_loc = len(elements[-2]) > 0
            loc_elements = elements[-2].split("-")
            elements[-2] = loc_elements[0]
            if len(loc_elements) > 1:
                cleaner_str = '-' + '-'.join(loc_elements[1:])
            else:
                cleaner_str = ""
            if len(elements) == 2 and had_loc and len(elements[-2]) == 0:
                seed_id = elements[-1]
            else:
                seed_id = '.'.join(elements)
            return [seed_id, cleaner_str]
        else:
            assert "-" not in in_str
            seed_id = in_str.split("-")[0]
            return [seed_id, ""]
