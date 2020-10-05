from pylab import *
from PyDynamicaLC import *

####### INPUT #######
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
    # r = numpy array of relative planet radii [r/R_star]
    # PlanetFlux: numpy array of additional fluxes from the planets (should be 0?)
    # transit_width: transit width relative to the orbital period (0-1, where 1 means the transit duration is orbital period wide, and 0 means it's infinitely narrow. Our default value is 0.05)

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

    # if dyn_coords == "cartesian":
    #     position_vec: vector of x, y, z position of each planet at t_min, in AU. Mind that the z coordinate should have its sign opposite from the convention (i.e. z -> -z)
    #     velocities_vec: vector of x_dot, y_dot, z_dot velocities of each planet at t_min, in AU/day. Mind that the z_dot coordinate should have its sign opposite from the convention (i.e. z -> -z)

    ### ECCENTRIC/CIRCULAR LIGHTCURVE - Keplerian parameters are AVERAGE values ###
    # p: vector of average planetary orbital periods [days]
    # incs: vector of average planetary inclinations [deg]
    # Omegas: vector of average planetary longitudes of ascending node [deg]
    # ecosomega: vector of average planetary ecos(omega)s [deg]
    # esinomega: vector of average planetary esin(omega)s [deg]
    # tmids: the initial time of transit (according to the mean ephemeris)

 ## MultiNest_params: Optimization parameters for MultiNest
    # mode: Optimization mode. At this moment, this can only be "TTVFaster" (string)
    # data_LC: normalized flux of the multi-planet data light-curve. Should be the length of LC_times (np.array)
    # data_LC_err: errors of the normalized flux light-curve. Should be the length of data_LC/LC_times (np.array)
    # nPl: number of planets (integer)
    # mass_prior: list specifying mass prior parameters. Supposed to contain 3 components: "lin" or "log" strings specifying linear or logarithmic prior, and lower and upper limit numbers. In the case of "lin", the limit numbers are the absolute mass limits of the prior in [m_earth]. In the case of "log" the limit numbers are *powers* (of base ten) of the prior (the whole number in units of m_earth).
    # ex_prior: list specifying ecosomega prior. Similar to mass_prior with units of eccentricity, but in addition - a Rayleigh distribution of the prior can be chosen with a string "ray", then there is only one additional required number in the list - the scale-width of the distribution.
    # ey_prior: similar to ex_prior.
    # verbose: Boolean.
    # basename: a string containing the header of all MultiNest files to be generated
    # sampling_efficiency: MultiNest parameter (see https://johannesbuchner.github.io/PyMultiNest/). Recommended value: 0.5.
    # evidence_tolerance: MultiNest parameter (see link above). Recommended value: 0.01.

## Analyzer_params: Plotting and error-estimation of MultiNet output
    # err: can be either "percentile" or "chi2" (string). This performs a MultiNest-independent error estimation in the following manner:
    #     percentile: only the delta_loglike < 3sigma (relative to the best-fit) is considered. The best-fit is then the median with the ±1sigma uncertainties are the 16th and 84th percentiles.
    #     chi2: best-fit is unchanged, and the ±1sigma uncertainties are calculatesd a the absolute difference of the best-fit value and the minimal and maximal values in the 1sigma range of delta_chi2.
        # plot_posterior: Boolean. If true - plots the MultiNest posterior distribution for all parameters. 1-5 sigma ranges are color-coded.
        # plot_bestfit_LC: Boolean. If true - plots generates a light-curve with the best-fit values and plots it against the data.
        # plot_corner: uses corner.py (https://corner.readthedocs.io/en/latest/) to generate corner plots for all parameters within the delta_chi2 < 3sigma range. NOTE: this option requires corner.py to be installed!

generate_LC = True
run_new_MultiNest = False
run_MultiNest_output = False
fontsize = 12

LC_mode = "circ"
dyn_coords = "keplerian"

Dyn_path = ".../c_version_main/"
Dyn_path_OS = ".../c_version_main/"  # Path used by os to reach TTVFast C folder

LC = np.loadtxt(".../LC_inputTest.txt") # Uploading example noisy lightcurve and its uncertainties
LC_times, LC_noise = LC[0], LC[2]
verbose = True # Verbosity

Paths = {"dyn_path": Dyn_path, "dyn_path_OS": Dyn_path_OS, "dyn_fileName": "TestName", "Coords_output_fileName": "Coords_output_Test"}
Integration_Params = {"t_min": 120.53851, "t_max": 1591.001, "dt": 1./100.}
Phot_Params = {"LC_times": LC_times, "LDcoeff1": 0.3279, "LDcoeff2": 0.3061, "r": np.array([0.00538482, 0.0148575,  0.03481889]), "PlanetFlux": np.array([0., 0., 0.]), "transit_width": 0.05, "ReBin_N": 1., "ReBin_dt": 0.}

#### Generating light-curve based on desired mode ####
if LC_mode == "osc":
    if dyn_coords == "keplerian":
        ###### OSCULATING LIGHTCURVE, KEPLERIAN INITIAL CONDITIONS ######
        Dyn_Params = {"LC_mode": "osc", "dyn_coords": "keplerian", "m_star": 1.334, "r_star": 1.6, "masses": [1., 5., 9.], "p": [7.41, 15.41, 28.45],\
                     "incs": [90., 90., 90.], "tmids": [121.17, 121.94, 136.70], "Omegas": [- pi / 2, - pi / 2, - pi /2],\
                     "ecosomega": [1e-4, 1e-4, 1e-4], "esinomega": [1e-4, 1e-4, 1e-4]}
    if dyn_coords == "cartesian":
        ###### OSCULATING LIGHTfclose(Coords_output);CURVE, CARTESIAN INITIAL CONDITIONS ######
        positions_vec = np.array([[0.04176802491058396, -0.0011453810936292114, -0.07040784157434199], [0.07213230935274573, -0.001978043815797308, -0.11218928312952507], [-0.08330707409886152, 0.0022844831145995793, 0.1826915599942284]])
        velocities_vec = np.array([[-0.059693834798333026, 0.0016369505125177963, -0.03544160989100488], [-0.04573798568370711, 0.0012542470987065284, -0.02943137831461144], [0.040322512538516535, -0.0011057416195310865, 0.018398419955714666]])

        Dyn_Params = {"LC_mode": "osc", "dyn_coords": "cartesian", "m_star": 1.334, "r_star": 1.6, "masses": [1., 5., 9.], "xyz": positions_vec, "uvw": velocities_vec}
if LC_mode == "ecc":
    ###### ECCENTRIC LIGHTCURVE ######
    Dyn_Params = {"LC_mode": "circ", "m_star": 1.334, "r_star": 1.6, "masses": [1., 5., 9.], "p": [7.409776081165706, 15.40981173833136, 28.451153474441284],\
                 "incs": [88.429203673205, 88.42920367320512, 88.42920367320514], "tmids": [121.17024202356487, 121.93952937719575, 136.700951187845],\
                 "ecosomega": [-9.127376574448313e-05, 0.00034056562879647016, 4.214554692828396e-05], "esinomega": [1.9367434007565637e-05, -0.00016646251014989429, -7.494786805977495e-05], "a": np.array([11.00453233806568, 17.929493340940514, 26.984717192232818]), "r": np.array([0.00538482, 0.0148575,  0.03481889])}
if LC_mode == "circ":
    ###### CIRCULAR LIGHTCURVE ######
    Dyn_Params = {"LC_mode": "circ", "m_star": 1.334, "r_star": 1.6, "masses": [1., 5., 9.], "p": [7.409776081165706, 15.40981173833136, 28.451153474441284],\
                 "incs": [88.429203673205, 88.42920367320512, 88.42920367320514], "tmids": [121.17024202356487, 121.93952937719575, 136.700951187845],\
                 "ecosomega": [-9.127376574448313e-05, 0.00034056562879647016, 4.214554692828396e-05], "esinomega": [1.9367434007565637e-05, -0.00016646251014989429, -7.494786805977495e-05], "a": np.array([11.00453233806568, 17.929493340940514, 26.984717192232818]), "r": np.array([0.00538482, 0.0148575,  0.03481889])}

#### MultiNest optimization example ####
ex_prior = ["ray", 0.01] # Rayleigh distribution with scale-width of 0.01
ey_prior = ["ray", 0.01] # Rayleigh distribution with scale-width of 0.01
mass_prior = ["lin", 0.1, 50.] # Linear mass prior with the mass range being 0.1-50 m_earth

MultiNest_parms = {"mode": "TTVFaster", "data_LC": LC[1], "data_LC_err": LC_noise, "nPl": 3,  "mass_prior": mass_prior, "ex_prior": ex_prior, "ey_prior": ey_prior,\
                   "verbose": verbose, "basename": "PDLC_test", "sampling_efficiency": 0.5, "evidence_tolerance": 0.01}
analyzer_parms = {"plot_posterior": True, "plot_corner": True, "plot_bestfit_LC": True, "err": "percentile"}

if generate_LC == True:
    output = LightCurve_Gen(Dyn_Params, Phot_Params, Integration_Params, Paths, verbose)
    times, TTVs, lin_eph, multiplanet_lightcurve, lightcurve_individual = output["times"], output["ttvs"], output["lin_eph"], output["lightcurve_allPlanets"], output["lightcurve_singlePlanets"]

    #### Plotting TTVs ####
    for i in range(len(TTVs)):
        plt.subplot(len(TTVs), 1, i + 1)
        if i == 0:
            plt.title("TTVs", fontsize = fontsize)
        plt.plot(TTVs[i], 'o', label = "TTVs of planet {0}".format(i))
        plt.ylabel("TTVs [days]", fontsize = fontsize)
        plt.tick_params(labelsize = fontsize)
        if i == len(TTVs) - 1:
            plt.xlabel("Epoch", fontsize = fontsize)
    plt.show()

    #### Plotting individual light-curves ####
    for i in range(len(lightcurve_individual)):
        plt.subplot(len(lightcurve_individual), 1, i + 1)
        if i == 0:
            plt.title("Individual planet light-curves", fontsize = fontsize)
        plt.plot(LC_times, lightcurve_individual[i], label = "Light-curve planet {0}".format(i))
        plt.ylabel("Norm. flux", fontsize = fontsize)
        plt.tick_params(labelsize = fontsize)
        plt.legend(fontsize = fontsize)
        if i == len(lightcurve_individual) - 1:
            plt.xlabel("Time [days]", fontsize = fontsize)
    plt.show()

    #### Plotting multi-planet light-curve
    plt.figure()
    plt.title("Multi-planet light-curve", fontsize = fontsize)
    plt.plot(LC_times, multiplanet_lightcurve, color = 'black')
    plt.xlabel("Time [days]", fontsize = fontsize)
    plt.ylabel("Norm. Flux", fontsize = fontsize)
    plt.tick_params(labelsize = fontsize)
    plt.show()

if run_new_MultiNest == True:
    run_multinest(MultiNest_parms, Dyn_Params, Phot_Params, Integration_Params, Paths)

if run_MultiNest_output == True:
    multinest_analyzer(MultiNest_parms, analyzer_parms, Dyn_Params, Phot_Params, Integration_Params, Paths)

