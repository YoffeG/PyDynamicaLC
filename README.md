# PyDynamicaLC
Python-compatible Light-curve with embedded planetary dynamics generator.

This code relies on several existing routines, the links to and installation guides of are listed below:

PyAstronomy
=======

This routine is used to generate a Mandel-Agol light-lightcurve

Information regarding installation can be found here: https://github.com/sczesla/PyAstronomy


TTVFast (C) - NOTE: the required modified version may only be downloaded here 
=======

This routine is one out of two possibilities to simulate planetary dynamics (the other being TTVFaster). TTVFast is a simplectic n-body integrator, whereas TTVFaster is a semi-analytic model accurate to first order in eccentricity which approximates TTVs using a series expansion.

The C version of the code is in the directory c_version, the Fortran version is in fortran_version. Both versions have specific README files.

All associated information regarding TTVFast (C) can be found here: https://github.com/kdeck/TTVFast/tree/master/c_version

TTVFaster
=======

This routine is one out of two possibilities to simulate planetary dynamics (the other being TTVFast). TTVFast is a simplectic n-body integrator, whereas TTVFaster is a semi-analytic model accurate to first order in eccentricity which approximates TTVs using a series expansion.

Information regarding installation can be found here: https://github.com/ericagol/TTVFaster

Citations
=======

If you use this code, please cite Yoffe, Ofir and Aharonson (2020), ApJ 

For TTVFast, please cite [Deck &amp; Agol (2014)](https://iopscience.iop.org/article/10.1088/0004-637X/787/2/132/pdf)

For TTVFaster, please cite [Agol &amp; Deck (2015)](http://arxiv.org/abs/1509.01623)

Please check back for updates to ensure that you are using the latest version.

-Gideon Yoffe, Aviv Ofir and Oded Aharonson
