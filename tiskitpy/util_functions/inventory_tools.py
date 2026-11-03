"""
Tools for adding processed (decimated, FIR-converted...) channels to an
obspy Inventory
"""
from pathlib import Path

from ..logger import init_logger

logger = init_logger()


def find_channel(inv, network, station, location, channel):
    """
    Return the Channel matching network, station, location and channel codes

    Args:
        inv (:class:`obspy.core.inventory.Inventory`): inventory
        network (str): network code
        station (str): station code
        location (str): location code
        channel (str): channel code
    Returns:
        (:class:`obspy.core.inventory.Channel`): the channel (not a copy)
    Raises:
        ValueError: if no channel or more than one channel matches
    """
    channels = [cha for net in inv if net.code == network
                for sta in net if sta.code == station
                for cha in sta if cha.location_code == location
                and cha.code == channel]
    seed_id = f'{network}.{station}.{location}.{channel}'
    if len(channels) == 0:
        raise ValueError(f'"{seed_id}" not found in inventory')
    if len(channels) > 1:
        raise ValueError(f'{len(channels)} channels match "{seed_id}"')
    return channels[0]


def insert_channel(inv, network, station, new_cha, source_location,
                   source_channel, force_overwrite=False):
    """
    Insert a new channel into a station of an inventory

    The channel is rejected if it has the same codes as its source channel.
    If another channel with the same codes exists, it is overwritten if
    force_overwrite is True, otherwise the new channel is rejected.

    Args:
        inv (:class:`obspy.core.inventory.Inventory`): inventory
        network (str): network code
        station (str): station code
        new_cha (:class:`obspy.core.inventory.Channel`): channel to insert
        source_location (str): location code of the channel new_cha was
            made from
        source_channel (str): channel code of the channel new_cha was
            made from
        force_overwrite (bool): overwrite an existing channel with the same
            codes
    Returns:
        (tuple):
            (:class:`obspy.core.inventory.Inventory`): inventory (a new
                object if a channel was overwritten)
            (bool): True if the channel was inserted
    """
    seed_id = f'{network}.{station}.{new_cha.location_code}.{new_cha.code}'
    if (new_cha.code == source_channel
            and new_cha.location_code == source_location):
        logger.error('New channel has same codes as source channel '
                     f'({seed_id}), rejected!')
        return inv, False
    if len(inv.select(network=network, station=station,
                      location=new_cha.location_code,
                      channel=new_cha.code)) > 0:
        boilerplate = 'New channel has same seed_id as an existing channel'
        if force_overwrite is True:
            logger.warning(f'{boilerplate} ({seed_id}), overwriting!')
            inv = inv.remove(network=network, station=station,
                             location=new_cha.location_code,
                             channel=new_cha.code)
        else:
            logger.warning(f'{boilerplate} ({seed_id}), rejected!')
            return inv, False
    for net in inv:
        if not net.code == network:
            continue
        for sta in net:
            if sta.code == station:
                sta.channels.append(new_cha)
                return inv, True
    raise ValueError(f'Did not find a home for {seed_id}')


def write_inventory_unique(inv, filename):
    """
    Write an inventory to a StationXML file, without overwriting a file

    If filename exists, a counter is added before '.xml'

    Args:
        inv (:class:`obspy.core.inventory.Inventory`): inventory
        filename (str): output filename
    Returns:
        (str): the name of the file written
    """
    of = str(filename)
    counter = 1
    while Path(of).exists():
        of = str(filename).replace('.xml', str(counter) + '.xml')
        counter += 1
    inv.write(of, format="STATIONXML")
    return of
