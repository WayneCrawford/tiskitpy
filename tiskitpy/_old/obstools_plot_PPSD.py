#!/usr/bin/env python3
"""
Plot station PPSD
"""
from obspy.clients.filesystem.sds import Client
from obspy import read_inventory, UTCDateTime
from obspy.signal import PPSD
import sys
import argparse

# inv_file='/Users/crawford/_Work/Parc_OBS/7_Missions/2019.Mayotte/obsinfo/1T.MOSO.STATION.xml'
# SDS_DIR='/Volumes/Wayne_Data/_NoBackup/MAYOBS_SDS'
# start_time=UTCDateTime("2019-03-01T00:00:00")
# end_time=UTCDateTime("2019-05-01T00:00:00")
# interval=86400


def main():
    arguments=_get_arguments()
    inv=read_inventory(arguments.sta_file)
    net,sta,loc=_get_default(arguments,inv)
    seed_id=f"{net}.{sta}.{loc}.{arguments.component}"
    
    now_time=arguments.start_time
    first_read=True
    client=Client(arguments.SDS_DIR)
    print(f'Calculating PPSDs for {seed_id}')
    while now_time < arguments.end_time - arguments.interval:
        print(f'\t{now_time}')   
        st=client.get_waveforms(net, sta, loc, arguments.component,
                                now_time, now_time + arguments.interval) 
        tr = st.select(id=seed_id)[0]
        if first_read:
            if arguments.component[1]=='D':
                ppsd=PPSD(tr.stats, metadata=inv,special_handling='hydrophone')
            else:
                ppsd = PPSD(tr.stats, metadata=inv)
            first_read=False
        ppsd.add(st)
        now_time += interval
    
    ppsd.save_npz(f'{seed_id}_PPSD.npz')
    ppsd.plot(f'{seed_id}_PPSD.png')
    #ppsd.plot_temporal([0.1,1,10])
    #ppsd.plot_spectrogram()
    return 0

#############################################################################
def _get_arguments():
    """
    Get command line arguments
    """
    # Set up
    p = argparse.ArgumentParser (description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('SDS_DIR',             help='SDS directory')
    p.add_argument('sta_file',            help='StationXML file')
    p.add_argument('component',           help='component to plot')
    p.add_argument('start_time',type=UTCDateTime,  help='PPSD start time')
    p.add_argument('end_time',  type=UTCDateTime,  help='PPSD end time')
    p.add_argument("-n","--network",  default='*', help="network code")
    p.add_argument("-s","--station",  default='*', help="station code")
    p.add_argument("-l","--location", default='*', help="location code")
    p.add_argument("-i","--interval", default=86400,
                   help="interval between PPSD calculations (seconds)")

    results = p.parse_args()
    return results

#############################################################################
def _get_default(arguments,inv):
    """Gets first net,sta,loc values from inventory, if not already set"""
    net=arguments.network
    sta=arguments.station
    loc=arguments.location
    if net=='*':
        net=inv[0].code
    if sta=='*':
        sta=inv.select(network=net)[0][0].code
    if loc=='*':
        loc=inv.select(network=net,station=sta)[0][0][0].location_code
    return net,sta,loc
    
#############################################################################
if __name__ == '__main__':
    sys.exit(main())
