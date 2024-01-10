from obspy.core.stream import read
from obspy.core.inventory import read_inventory
from tiskitpy import CleanRotator, SpectralDensity, CleanSequence

stream = read('data/XS.S11D.LH.2016.12.11.mseed', 'MSEED')
inv = read_inventory('data/XS.S11_decimated.station.xml', 'STATIONXML')
cr = CleanRotator(stream)
stream_rotated = cr.apply(stream)
    
z_compare = stream.select(channel='*Z') + stream_rotated.select(channel='*Z')
CleanSequence.stream_print(z_compare)
CleanSequence.stream_plot(z_compare, outfile='5_CleanRotator_z_compare.png')

sd_compare = SpectralDensity.from_stream(z_compare, inv=inv)
sd_compare.plot(overlay=True, outfile='5_CleanRotator_spect_compare.png')
