#!/usr/bin/env python3
"""
Plot instrument response
"""
import sys
import argparse
from obspy import read_inventory, UTCDateTime
from obspy.signal.spectral_estimation import get_nlnm,get_nhnm
import numpy as np
import matplotlib.pyplot as plt

#############################################################################
def _get_arguments():
    """
    Get command line arguments
    """
    # Set up
    p = argparse.ArgumentParser (description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('sta_file',               help='StationXML file')
    p.add_argument('component',                help='component to plot')
    p.add_argument('-f','--min_freq',  type=float, default=None,
                    help='approx lower frequency limit')
    p.add_argument("-r","--stage_range", type=int, nargs=2, default=(None,None), 
                    help="plot only stages from MIN_STAGE to MAX_STAGE")
    p.add_argument('-e','--plot_each_stage', default=False, action='store_true',
                    help='plot each stage response (instead of product of all stages)')
    p.add_argument("-o","--output_type", default='VEL', choices=['VEL','ACC','DISP'], 
                    help="output type")
    p.add_argument("-n","--network",  default='*', help="network code")
    p.add_argument("-s","--station",  default='*', help="station code")
    p.add_argument("-l","--location", default='*', help="location code")

    # Parse 
    results = p.parse_args()
    return results

#############################################################################
def main():
    arguments=_get_arguments()
    
    inv=read_inventory(arguments.sta_file)
    net,sta,loc=_get_default(arguments,inv)
    SEED_id=f"{net}.{sta}.{loc}.{arguments.component}"
    try:
        channel=inv.select(network=net,station=sta,location=loc,channel=arguments.component)[0][0][0]
    except:
        print(f'{SEED_id} not found in {arguments.sta_file}')
    SEED_id=f"{net}.{sta}.{loc}.{channel.code}"
    if not arguments.min_freq:
        print('Calculating minimum frequency from lowest frequency pole:')
        paz=channel.response.get_paz()
        min_pole=np.min(np.abs(paz.poles))/(2*np.pi)
        print(f'    LF pole found at {min_pole:.2g} Hz ({1/min_pole:.2g} s)')
        arguments.min_freq=(min_pole)/10
        
    outf_base=f'{SEED_id}'
    if arguments.stage_range[0]:
        outf_base=outf_base+'_stages{:d}-{:d}'.format(
                            arguments.stage_range[0],arguments.stage_range[1])
    
    plt.figure(1)
    _plot_response(channel,SEED_id,arguments)
    plt.savefig(outf_base+'_response.png')
    
    return 0


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
def _plot_response(channel,seed_id,args):
    """ Plot instrument response """

    sr=channel.sample_rate
    t_samp=1./sr
    nfft=2**(np.ceil(np.log2(sr/args.min_freq)))
    resp=channel.response
    title=seed_id
    if args.plot_each_stage:
        if args.stage_range[0]==None:
            rng=range(1,len(resp.response_stages)+1)
        else:
            rng=range(args.stage_range[0],args.stage_range[1]+1)
        title=title+f', stages {rng[0]:d} to {rng[-1]:d}'
        for stage in rng:
            r,f=resp.get_evalresp_response(t_samp,nfft,
                                start_stage=stage,
                                end_stage=stage)
            plt.loglog(f[1:],abs(r[1:]),label=f'{stage:d}')
        plt.legend(fontsize='small',labelspacing=0.1)
    else:
        r,f=resp.get_evalresp_response(t_samp,nfft,
                                start_stage=args.stage_range[0],
                                end_stage=args.stage_range[1])
        plt.loglog(f[1:],abs(r[1:]),label='Min sensitivity')
    
    # Finish up plot
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Gain')
    plt.xlim([f[1],f[-1]])
    plt.title(title)
#############################################################################
if __name__ == '__main__':
    sys.exit(main())
