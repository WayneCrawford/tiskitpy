*******************************
clean_sequences
*******************************

``tiskitpy`` data cleaning can create multiple instances of a same data trace.
To keep track of each instance, cleaning steps are stored in ``clean_sequence``
list and an unique ``tiskitpy_id`` is generated for each instance

The ``clean_sequences`` list
====================================================

Every time a cleaning step is applied to a data object, a "tag" is added to
the object's ``clean_sequence`` list.  Tags are strings in one of two formats:

1. A ``seed_id`` code.  Indicates that coherent noise from that channel was subtracted
   from this channel.
2. A short, all caps text code, specifying a specific transformation.  Examples
   include ``'ROT'`` for the ``CleanRotator`` and ``'AVOID'`` and ``'SPANS'`` for specific
   spans specified to avoid or include when calling ``SpectralDensity.from_stream()``

Classes implementation
====================================================

Different classes store their clean_sequence lists in different places:

- ``obspy.Trace``: in ``self.stats.clean_sequences``
- ``SpectralDensity``: in ``self._clean_sequences``
- ``ReponseFunctions``: in ``self.input_clean_sequence``

``seed_id`` and ``tiskitpy_id``
====================================================

In order to distinguish between channels having the same ``seed_id`` but different
``clean_sequence``, while maintaining compatibility with
`obspy <https://docs.obspy.org>`, tiskitpy generates ``tiskitpy_id``,
which is the ``stream_id`` stuffed with minimal information about the 
``clean_sequence``.  For example, if you rotated the data on the channel with
``seed_id=XX.STA.00.BHZ``, then cleaned coherent noise from its
``BDH``, ``BD1`` and then ``BH2`` channels, the new channel's ``tiskitpy_id``
would be ``XX.STA.00-ROT-H-1-2.BHZ``.

The ``tiskitpy_id`` is shown when plotting spectra and coherencies and is used for
selecting SpectralDensity channels to plot.  

The ``tiskitpy_id`` is also shown and used when you ``plot()`` and ``print()``
using the ``CleanedStream`` subclass of :class:``obspy.Stream``.
