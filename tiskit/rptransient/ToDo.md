TO DO
======================

- Allow multiple startimes and an endtime
    - transient creation only within a startime block (program verifies only)
- Allow transients shorter than the period
    - shouldn't be difficult, as convolve will still work
- Allow transient to be applied to data outside of its "comb"
    - Just need transient and period?
    - Reconstruct a Dirac comb that best fits the data (should we have
      prep-filter parameters for other cases than ours?)
    - Should be able to take into account multiple starttimes, re-searching the
      best dirac comb afterwards
      - Allows processing of entire dataset, in order to create a continuous
        "cleaned" dataset
- Allow transient to be removed from other channels? (horizontala)
    - Would probably need to re-scale, at least
- Allow transient to work on data with different sampling rates
    - Useful for cleaning high-frequency data sets after calculating the 
      transient at low frequencies