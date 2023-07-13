*******************************
TiSKit
*******************************

Time series data processing, based on [BP2010]_


Classes
=========================

:CleanRotator: Rotate data to minimize noise on vertical channel
:DataCleaner: Transfer_Function-based data cleaning
:Decimator: Decimate time series and update metadata with the decimator's
            response
:PeriodicTransient: Calculate and remove periodic transient (VERY manual!)
:SpectralDensity: Calculate and manipulate spectral density functions.
:TimeSpan: Specify time spans to be removed, kept, zeroed, etc.
:ResponseFunctions: Frequency response functions for a given input channel.
               
Functions
=========================

:fir2caus: Transform zero-phase data to minimum phase (only works for
           LCHEAPO loggers, need to update to calculate/work for any
           zero-phase filter)
:read_MSEED: Read MSEED data, even if the file is too big (> 2 GB)
             for obspy's read() function
:Peterson_noise_model: Return the Peterson High and Low Noise Models

Command-line programs
=========================

:tiskit_decimate_SDS: Decimate data in a SeisComp Data Structure database.
    Inserts the data into the same database and creates a new StationXML file
    (based on an existing StationXML file for the original database)

*I had to remove `SpectralDensity.from_ATACR()` (never tested
anyway) because readthedocs had a problem importing the `obstools` package
using `pip`*

.. [BP2010] Bendat J. S. and A. G. Piersol (1986), Random Data:
    Analysis and Measurement Procedures, 566 pp.