The file `longtest_cleaning.py` is a longer, more involved test of data
cleaning.
It does not start with `test_` so that it is not automatically run by
`pytest` or by `python -m unittest`

# Things to change

- test_decimate reads in Scripps catalog whose FIR8 has sum[coeffs] = 0.9767,
  which creates a warning.