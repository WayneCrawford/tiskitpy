"""
Decimate seismological data using SAC FIR Filters

There are two main methods:

- Decimator.run(): decimates seismological data and updates channel/location codes
- Decimator.update_inventory(), which creates a new channel in the metadata, adding
in the new decimation stages.

Channel/location code rules:
    - If the change in sample rate also changes the band code, then the new
      band code is applied
    - otherwise, the location code is incremented by one.

Example:

    st = read('mydata.mseed')
    inv = read_inventory('mynetwork.xml')
    
    decim = Decimate([2,5])     # decimate by 10, in two steps
    st_dec = decim.run(st)      # decimates all the channels in the stream
    st_dec = decim.run(st, inv) # same as above, but makes sure new channel/loc
                                # code isn't already in the inventory
    inv_dec = decim.update_inventory(inv, st)  # creates decimated channels
                                               # for all channels in st

If the channel band code is the same as for the input, will change the
location code

examples:
decimator = Decimator([5,5,10])
# Output a decimated version of the input stream
outstream = decimator.run_stream(instream)

# Decimate all data with sampling rate 62.5 sps on the SDSClient,
# outputting the result to the same SDS client.  Passes other kwargs
# to client.get_waveforms()
chan_list = decimator.run_SDS(SDSClient, in_sps=62.5)

# Add new decimated channels' metatdata to inventory
inv_in = read_inventory('mynetwork.xml')
inv_out = decimate.add_stages(chan_list)
inv_out.write('mynetwork_decim.xml')
"""
from .decimator import Decimator
from .fir_filter import FIRFilter