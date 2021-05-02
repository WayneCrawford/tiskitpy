## VERSIONS

v0.3:
- Handles spliced datafiles with zero-filled gaps:
  + Adds dates to headers in gap
  + Verifies that gap is the right size
    (no time tear afterwards)
- Handles incomplete datafiles
  + Warns the user that the data do not go as far as the
    directory claims
  + Changes the "directoryEntries" value in the header
