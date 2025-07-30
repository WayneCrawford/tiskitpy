# TiSKitPy

Routines for time series data processing

Uses the [obspy](https://docs.obspy.org) Stream (data) and Inventory (metadata)
classes


[Documentation](https://tiskitpy.readthedocs.io)


## Main Classes

- ``SpectralDensity``: Calculate and manipulate spectral density functions.
- ``Decimator``: Decimate time series and update metadata with the decimator's
  response
- ``CleanRotator``: rotate data to minimize noise on vertical channel
- ``DataCleaner``: Transfer_Function-based data cleaning
- ``ResponseFunctions``: Frequency response functions for a given input channel.
- ``Compliance``: Seafloor Compliance
- ``SeafloorSynthetic``: Generate synthetic seafloor data, including compliance signal
            
               
## Functions

- ``FIR_corr``: transform zero-phase data to minimum phase.  Based on
  [Scherbaum, Of Poles and Zeros](https://doi.org/10.1007/978-1-4020-6861-4).
  Only works for LCHEAPO loggers, need to update to calculate/work for any
  zero-phase filter)
- ``readMSEED``: read in MSEED data, including if the file is too big (> 2 GB)
               for obspy's read() function
- ``rptransient``: calculate and remove periodic transient (VERY manual!).  
 	Based on Matlab code by E Wielandt, used in
 	[Deen et al., 2017](https://doi.org/10.1002/2017GL074892) and
 	[Aminian et al., 2025](https://doi.org/10.1093/gji/ggaf253)
- ``PetersonNoiseModel``: return the Peterson High and Low Noise Models
  ([Peterson, 1993](https://doi.org/10.3133/ofr9332))
- ``plot_compliance_stack()``: plot, from top to bottom, Z PSD, P PSD,
  Z-P coherence and Z-P frequency response function.


## Installation

First, install `obspy` using the instructions on their webpage.
Then, in the pip/conda environment that contains obspy, type 
`pip install tiskitpy`