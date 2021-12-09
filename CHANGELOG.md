v0: straight translation of Wielandt's Matlab routines
  
v1: Add features to handle earthquakes and noisy data:
- Remove "teeth" from the Dirac comb wherever the "eq_template" time series
  covering the same samples as the data contains a zero within the slice
  corresponding to that tooth
- Also remove teeth for slices that have significantly higher variance than
  the others
- Use this damaged comb and a clipped version of the data (using the clips
  variable) to generate the average transient.
- Apply this average transient and the full comb to clean the data
    
v2:
- In individual matching, allow further shifting of comb tooth:
    - In time
    - In height? (shouldn't allow too much liberty or may create errors)
- Recalculate transient based on this comb?
- Return a sparse dirac comb for the time series, which could be combined
  with others to create a dirac comb for the whole dataset
- Return the number of teeth used to calculate the master transient

v3:

- Remove parameters from code (now arguments to calls)
- Create classes
- Put on github/pypi

1.0a1
- Add automatic read of global earthquakes in order to remove them from noise/transient models
