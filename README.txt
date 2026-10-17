MOSFLM 7.0.9
============

This program can be copied either as a pre-built executable (available
for some popular platforms) or in a build-it-yourself form (see the
appropriate directories).

We strongly recommend using the pre-built executables unless you have a 
properly installed copy of CCP4 version 5.0 or later. Even then, it is 
likely that you will gain very little by building your own copy.

Compiler notes
==============

(i) gfortran/gcc 4 is currently under active development; optimization of many
of the routines in Mosflm is not possible with this compiler, so the build is
performed with all optimization turned off. Surprisingly, some of the C
functionality seems to be rather fragile, as can be seen in the X11 display.

(ii) The Intel compilers ifort and icc for Mac OS X Intel-based machines 
give fast, stable executables and are well worth installing (as a time-limited
free download) as an alternative to gfortran/gcc 4. I include a method to use 
them, but cannot distribute executables built with them.


General notes
=============
   configure --with-ccp4-prefix=PFX
   make
   make install

You need to have at least CCP4 version 5.0 installed to build MOSFLM, 
including the X-windows libraries. I recommend installing the newest version 
you can find - currently 6.3.0.

The build-it-yourself version will not link with the 
now obsolete CCP4 version 4. 

The build-it-yourself version does not contain source code for the
REFIX autoindexing routine, but does have the DPS autoindexing
code. If you wish to build MOSFLM version 7.0.1 yourself, and want
REFIX, then e-mail Andrew Leslie (andrew@mrc-lmb.cam.ac.uk) or me 
(harry@mrc-lmb.cam.ac.uk) directly and ask for it. The file for
inclusion (refix_all.for) has changed so old versions will probably
not work.

The file 'mosflm.hlp' should be copied to the directory $CCP4_HELPDIR
so that the online help facility works.

Harry Powell 
harry@mrc-lmb.cam.ac.uk

20th August 2011
==============================================================================
