# PyDynamicaLC
Pythonic photodynamical model generator, using three different approximations.

This code relies on several existing open-source routines, the links to and installation guides of are listed below:

PyAstronomy
=======

This routine is used to generate a Mandel-Agol light-curve

Information regarding installation can be found here: https://github.com/sczesla/PyAstronomy


TTVFast (C) - NOTE: the required modified version may only be downloaded here 
=======

This routine is one out of two possibilities to simulate planetary dynamics (the other being TTVFaster). TTVFast is a simplectic n-body integrator, whereas TTVFaster is a semi-analytic model accurate to first order in eccentricity which approximates TTVs using a series expansion.

We modified TTVFast to extract therefrom the instantaneous osculating Keplerian parameters at each time of mid-transit, with which each individual transit shape and timing will be determined in the "osculating" mode.

The C version of the code is in the directory c_version, the Fortran version is in fortran_version. Both versions have specific README files.

All associated information regarding TTVFast (C) can be found here: https://github.com/kdeck/TTVFast/tree/master/c_version

TTVFaster
=======

This routine is one out of two possibilities to simulate planetary dynamics (the other being TTVFast). TTVFast is a simplectic n-body integrator, whereas TTVFaster is a semi-analytic model accurate to first order in eccentricity which approximates TTVs using a series expansion.

Information regarding installation can be found here: https://github.com/ericagol/TTVFaster

INPUT
===

Each of the following is a dictionary which should contain the entries the exact entries listed here.

## Integration_Params:
    # t_min = t0 of integration [days]
    # t_max = integration limit [days]
    # dt = n-body sampling frequency [days^-1]

## Paths:
    # dyn_path: Path to TTVFast folder (should end with c_version/)
    # dyn_path_OS: Same path as dyn_path, except that it should be compatible with OS
    # dyn_file_name: Desired name of file
    # Coords_output_fileName: name of the desired cartesian coordinates file from TTVFast.

## Phot_Params: PHOTOMETRIC PARAMETERS
    # LC_times: list of times for which the lightcurve will be sampled [days]
    # LDcoeff1: linear limb-darkening coefficient
    # LDcoeff2: quadratic limb-darkening coefficient
    # r = numpy array of relative planet radii [R_star]
    # PlanetFlux: numpy array of additional fluxes from the planets (should be 0?)
    # transit_width: upper limit of the transit width relative to the orbital period

## Dyn_Params: DYNAMICAL PARAMETERS
    # LC_mode: string specifying which mode of lightcurve generation to use. Can be osc, ecc and circ (osculating, eccentric and quasi-circular, respectively)
    # m_star: stellar mass [m_sun]
    # r_star: stellar radius [r_sun]
    # masses: vector of planetary masses [m_earth]

    ### OSCULATING LIGHTCURVE - Keplerian parameters are initial conditions ###
    # dyn_coords: for LC_mode = osc, the initial conditions input can be input as either the Keplerian parameters or a Cartesian state vector. The input has to be specified as follows:
    
    # if dyn_coords == "keplerian":
        # p: vector of planetary orbital periods at t_min [days]
        # incs: vector of planetary inclinations at t_min [deg]
        # Omegas: vector of planetary inclinations at t_min [deg]
        # ecosomega: vector of planetary ecos(omega)s at t_min [deg]
        # esinomega: vector of planetary esin(omega)s at t_min [deg]
        # tmids: vector of first times of mid-transit for each planet [days]
        
    if dyn_coords == "cartesian":
        position_vec: vector of x, y, z position of each planet at t_min, in AU. Mind that the z coordinate should have its sign opposite from the convention (i.e. z -> -z)
        velocities_vec: vector of x_dot, y_dot, z_dot velocities of each planet at t_min, in AU/day. Mind that the z_dot coordinate should have its sign opposite from the convention (i.e. z -> -z)

    ### ECCENTRIC/CIRCULAR LIGHTCURVE - Keplerian parameters are AVERAGE values ###
    # p: vector of average planetary orbital periods [days]
    # incs: vector of average planetary inclinations [deg]
    # Omegas: vector of average planetary longitudes of ascending node [deg]
    # ecosomega: vector of average planetary ecos(omega)s [deg]
    # esinomega: vector of average planetary esin(omega)s [deg]
    # tmids: a mean reference epoch (i.e. the constant term in the linear fit to the times of transit)
    
 CLARIFICATION REGARDING t_mid FOR OSCULATING AND ECCENTRIC/QUASI-CIRCULAR CASES:
 
 Osculating: in this case, t_mid is simply the time of the first time of mid-transit with respect to t_min. The phase of the planet is then calculated by "reqinding" a Keplerian arc from that point to t_min (note: this is an approximation that may not be valid in very eccentric systems).
 
 Eccentric/quasi-Circular: In this case, the value t_min represents the LINEAR APPROXIMATION of the first time of mid-transit. That is, the time of mid-transit of the given epoch in a linear ephemeris regime (we extracted this value by fitting a line to the observed times of transit. The value t_mind would then be the value of the fitted line at the given epoch).

Citations
=======

If you use this code, please cite Yoffe, Ofir and Aharonson (2020), ApJ 

For TTVFast, please cite [Deck &amp; Agol (2014)](https://iopscience.iop.org/article/10.1088/0004-637X/787/2/132/pdf)

For TTVFaster, please cite [Agol &amp; Deck (2015)](http://arxiv.org/abs/1509.01623)

Please check back for updates to ensure that you are using the latest version.

-Gideon Yoffe, Aviv Ofir and Oded Aharonson
