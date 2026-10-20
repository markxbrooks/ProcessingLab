The two files describe the same experiment with consistent geometry, but their intensities don't agree. Spot positions match closely, while the measured intensities barely correlate. So ProcessingLab's predictions are right, and the problem is how it measures and records intensities.

You asked for integrate.mtz. The most recent match is integrated.mtz, written by DIALS on Oct 18 in the thaumatin_9_1_master DIALS run folder from 2026-10-15. I compared that with process.mtz from this morning.

What agrees

Indexing and geometry. Both are P4 with a nearly identical cell, and the same reflections carry the same indices. Matched detector positions differ by only a few pixels, and rotation angles differ by about a quarter of a degree on average.
Coverage within the shared range. Almost every unique reflection ProcessingLab measured also appears in DIALS.
ProcessingLab	DIALS
Cell a, c (Å)	58.09, 150.33	57.94, 150.27
Images	999, 0 to 150°	750, 0.15 to 112.65°
High resolution limit (Å)	1.45	1.15
Rows	410,462	459,539
Unique reflections	82,926	133,685
Median I/σ(I)	0.86	1.53
What disagrees

Intensities hardly correlate. For the 59,010 reflections DIALS measures strongly, the intensity correlation is only 0.12. Two integrations of the same images would normally exceed 0.9. Even reflections that ProcessingLab records on a single image only reach 0.24.
The intensity ratio drifts with resolution and rotation. ProcessingLab reads about 2.4 times DIALS at low resolution and about 0.07 times at high resolution. By rotation, it starts at about 4.7 times over the first 15° and settles near 0.2 beyond 30°. A smooth scale factor cannot fix trends like these.
Multi-image reflections are recorded inconsistently. DIALS writes one summed row per reflection. ProcessingLab splits a reflection into a row per image. In 53,666 cases it marks more than one of those rows as fully recorded. It also has 111,867 rows with an intensity of exactly zero. Summing those rows back together does not fix the correlation.
The image headers in the file lack rotation ranges. Every ProcessingLab image entry has a start and end angle of zero. Scaling programs such as Aimless rely on these angles.
The Lorentz-polarization column follows a different convention. It is about ten times the DIALS value, and the two are anti-correlated. Some difference is expected, since the programs define this column differently, but the values don't line up even after inverting the DIALS ones.
ProcessingLab's uncertainties are about twice as large. This is why its I/σ(I) is lower.
Points 1 to 3 suggest that the per-image summation, or its background estimate, is biased. The partial-reflection bookkeeping makes things worse on top of that. A good next step is to take a few strong reflections that span several images and trace their pixels and background through ProcessingLab's summation.

I didn't change any code. The comparison scripts are in the scratchpad, for example /private/tmp/claude-501/-Users-brooks-projects-ElMo/65ab2b44-295b-43ee-a92d-fa0bddea3d07/scratchpad/mtzcmp.py, with three follow-up scripts named mtzcmp2.py to mtzcmp4.py in the same folder.