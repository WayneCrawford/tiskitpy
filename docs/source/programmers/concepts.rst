*******************************
Concepts
*******************************

`CleanSequence`s
=========================

One of the core functionalities of tiskitpy is noise removal.
`CleanSequence`s allow us to know what noise has been removed from a channel.
At heart, they are simply lists of the operations that have been performed,
the secret is where to store them and how to see them.

Where they are stored depends on the classtype:
|  :class:            |    location                             | value if no clean_sequence
| ------------------- | --------------------------------------- | ----------------------------
| `Trace`             | `stats.clean_sequence`                  | `stats.clean_sequence` does not exist
| `DataCleaner`       | derived property `clean_sequence`       | None
| `SpectralDensity`   | method `clean_sequence(channel)`        | None
| `ResponseFunctions` | derived property `input_clean_sequence` | None

How they are named depends on the cleaning operation:
  - removing coherent noise from another channel: **the cleaned channel's seed_id**
  - rotating the vertical channel to minimize noise: **'RRR'**