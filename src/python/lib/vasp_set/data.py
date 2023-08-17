import os
from pathlib import Path
from typing import List

# Parameters that can be set in INCAR. The values which are None
# are not written and default parameters of VASP are used for them.

float_keys = [
    "aexx",  # Fraction of exact/DFT exchange
    "aggac",  # Fraction of gradient correction to correlation
    "aggax",  # Fraction of gradient correction to exchange
    "aldac",  # Fraction of LDA correlation energy
    "amin",  #
    "amix",  #
    "amix_mag",  #
    "bmix",  # tags for mixing
    "bmix_mag",  #
    "cshift",  # Complex shift for dielectric tensor calculation (LOPTICS)
    "deper",  # relative stopping criterion for optimization of eigenvalue
    "ebreak",  # absolute stopping criterion for optimization of eigenvalues
    # (EDIFF/N-BANDS/4)
    "efield",  # applied electrostatic field
    "emax",  # energy-range for DOSCAR file
    "emin",  #
    "enaug",  # Density cutoff
    "encut",  # Planewave cutoff
    "encutgw",  # energy cutoff for response function
    "encutfock",  # FFT grid in the HF related routines
    "hfscreen",  # attribute to change from PBE0 to HSE
    "kspacing",  # determines the number of k-points if the KPOINTS
    # file is not present. KSPACING is the smallest
    # allowed spacing between k-points in units of
    # $\AA$^{-1}$.
    "potim",  # time-step for ion-motion (fs)
    "nelect",  # total number of electrons
    "param1",  # Exchange parameter
    "param2",  # Exchange parameter
    "pomass",  # mass of ions in am
    "pstress",  # add this stress to the stress tensor, and energy E = V *
    # pstress
    "sigma",  # broadening in eV
    "smass",  # Nose mass-parameter (am)
    "spring",  # spring constant for NEB
    "time",  # special control tag
    "weimin",  # maximum weight for a band to be considered empty
    "zab_vdw",  # vdW-DF parameter
    "zval",  # ionic valence
    # The next keywords pertain to the VTST add-ons from Graeme Henkelman's
    # group at UT Austin
    "jacobian",  # Weight of lattice to atomic motion
    "ddr",  # (DdR) dimer separation
    "drotmax",  # (DRotMax) number of rotation steps per translation step
    "dfnmin",  # (DFNMin) rotational force below which dimer is not rotated
    "dfnmax",  # (DFNMax) rotational force below which dimer rotation stops
    "sltol",  # convergence ratio for minimum eigenvalue
    "sdr",  # finite difference for setting up Lanczos matrix and step
    # size when translating
    "maxmove",  # Max step for translation for IOPT > 0
    "invcurv",  # Initial curvature for LBFGS (IOPT = 1)
    "timestep",  # Dynamical timestep for IOPT = 3 and IOPT = 7
    "sdalpha",  # Ratio between force and step size for IOPT = 4
    # The next keywords pertain to IOPT = 7 (i.e. FIRE)
    "ftimemax",  # Max time step
    "ftimedec",  # Factor to dec. dt
    "ftimeinc",  # Factor to inc. dt
    "falpha",  # Parameter for velocity damping
    "falphadec",  # Factor to dec. alpha
    "clz",  # electron count for core level shift
    "vdw_radius",  # Cutoff radius for Grimme's DFT-D2 and DFT-D3 and
    # Tkatchenko and Scheffler's DFT-TS dispersion corrections
    "vdw_scaling",  # Global scaling parameter for Grimme's DFT-D2 dispersion
    # correction
    "vdw_d",  # Global damping parameter for Grimme's DFT-D2 and Tkatchenko
    # and Scheffler's DFT-TS dispersion corrections
    "vdw_cnradius",  # Cutoff radius for calculating coordination number in
    # Grimme's DFT-D3 dispersion correction
    "vdw_s6",  # Damping parameter for Grimme's DFT-D2 and DFT-D3 and
    # Tkatchenko and Scheffler's DFT-TS dispersion corrections
    "vdw_s8",  # Damping parameter for Grimme's DFT-D3 dispersion correction
    "vdw_sr",  # Scaling parameter for Grimme's DFT-D2 and DFT-D3 and
    # Tkatchenko and Scheffler's DFT-TS dispersion correction
    "vdw_a1",  # Damping parameter for Grimme's DFT-D3 dispersion correction
    "vdw_a2",  # Damping parameter for Grimme's DFT-D3 dispersion correction
    "eb_k",  # solvent permitivity in Vaspsol
    "tau",  # surface tension parameter in Vaspsol
    "langevin_gamma_l",  # Friction for lattice degrees of freedom
    "pmass",  # Mass for latice degrees of freedom
    "bparam",  # B parameter for nonlocal VV10 vdW functional
    "cparam",  # C parameter for nonlocal VV10 vdW functional
    "aldax",  # Fraction of LDA exchange (for hybrid calculations)
    "tebeg",  #
    "teend",  # temperature during run
    "andersen_prob",  # Probability of collision in Andersen thermostat
    "apaco",  # Distance cutoff for pair correlation function calc.
    "auger_ecblo",  # Undocumented parameter for Auger calculations
    "auger_edens",  # Density of electrons in conduction band
    "auger_hdens",  # Density of holes in valence band
    "auger_efermi",  # Fixed Fermi level for Auger calculations
    "auger_evbhi",  # Upper bound for valence band maximum
    "auger_ewidth",  # Half-width of energy window function
    "auger_occ_fac_eeh",  # Undocumented parameter for Auger calculations
    "auger_occ_fac_ehh",  # Undocumented parameter for Auger calculations
    "auger_temp",  # Temperature for Auger calculation
    "dq",  # Finite difference displacement magnitude (NMR)
    "avgap",  # Average gap (Model GW)
    "ch_sigma",  # Broadening of the core electron absorption spectrum
    "bpotim",  # Undocumented Bond-Boost parameter (GH patches)
    "qrr",  # Undocumented Bond-Boost parameter (GH patches)
    "prr",  # Undocumented Bond-Boost parameter (GH patches)
    "rcut",  # Undocumented Bond-Boost parameter (GH patches)
    "dvmax",  # Undocumented Bond-Boost parameter (GH patches)
    "bfgsinvcurv",  # Initial curvature for BFGS (GH patches)
    "damping",  # Damping parameter for LBFGS (GH patches)
    "efirst",  # Energy of first NEB image (GH patches)
    "elast",  # Energy of final NEB image (GH patches)
    "fmagval",  # Force magnitude convergence criterion (GH patches)
    "cmbj",  # Modified Becke-Johnson MetaGGA c-parameter
    "cmbja",  # Modified Becke-Johnson MetaGGA alpha-parameter
    "cmbjb",  # Modified Becke-Johnson MetaGGA beta-parameter
    "sigma_nc_k",  # Width of ion gaussians (VASPsol)
    "sigma_k",  # Width of dielectric cavidty (VASPsol)
    "nc_k",  # Cavity turn-on density (VASPsol)
    "lambda_d_k",  # Debye screening length (VASPsol)
    "ediffsol",  # Tolerance for solvation convergence (VASPsol)
    "deg_threshold",  # Degeneracy threshold
    "omegamin",  # Minimum frequency for dense freq. grid
    "omegamax",  # Maximum frequency for dense freq. grid
    "rtime",  # Undocumented parameter
    "wplasma",  # Undocumented parameter
    "wplasmai",  # Undocumented parameter
    "dfield",  # Undocumented parameter
    "omegatl",  # Maximum frequency for coarse freq. grid
    "encutgwsoft",  # Soft energy cutoff for response kernel
    "encutlf",  # Undocumented parameter
    "scissor",  # Scissor correction for GW/BSE calcs
    "dimer_dist",  # Distance between dimer images
    "step_size",  # Step size for finite difference in dimer calculation
    "step_max",  # Maximum step size for dimer calculation
    "minrot",  # Minimum rotation allowed in dimer calculation
    "dummy_mass",  # Mass of dummy atom(s?)
    "shaketol",  # Tolerance for SHAKE algorithm
    "shaketolsoft",  # Soft tolerance for SHAKE algorithm
    "shakesca",  # Scaling of each step taken in SHAKE algorithm
    "hills_stride",  # Undocumented metadynamics parameter
    "hills_h",  # Height (in eV) of gaussian bias for metadynamics
    "hills_w",  # Width of gaussian bias for metadynamics
    "hills_k",  # Force constant coupling dummy&real for metadynamics
    "hills_m",  # Mass of dummy particle for use in metadynamics
    "hills_temperature",  # Temp. of dummy particle for metadynamics
    "hills_andersen_prob",  # Probability of thermostat coll. for metadynamics
    "hills_sqq",  # Nose-hoover particle mass for metadynamics
    "dvvdelta0",  # Undocumented parameter
    "dvvvnorm0",  # Undocumented parameter
    "dvvminpotim",  # Undocumented parameter
    "dvvmaxpotim",  # Undocumented parameter
    "efermi",  # Undocumented parameter
    "enchg",  # Undocumented charge fitting parameter
    "tau0",  # Undocumented charge fitting parameter
    "encut4o",  # Cutoff energy for 4-center integrals (HF)
    "param3",  # Undocumented HF parameter
    "model_eps0",  # Undocumented HF parameter
    "model_alpha",  # Undocumented HF parameter
    "qmaxfockae",  # Undocumented HF parameter
    "hfscreenc",  # Range-separated screening length for correlations
    "hfrcut",  # Cutoff radius for HF potential kernel
    "encutae",  # Undocumented parameter for all-electron density calc.
    "encutsubrotscf",  # Undocumented subspace rotation SCF parameter
    "enini",  # Cutoff energy for wavefunctions (?)
    "wc",  # Undocumented mixing parameter
    "enmax",  # Cutoff energy for wavefunctions (?)
    "scalee",  # Undocumented parameter
    "eref",  # Reference energy
    "epsilon",  # Dielectric constant of bulk charged cells
    "rcmix",  # Mixing parameter for core density in rel. core calcs.
    "esemicore",  # Energetic lower bound for states considered "semicore"
    "external_pressure",  # Pressure for NPT calcs., equivalent to PSTRESS
    "lj_radius",  # Undocumented classical vdW parameter
    "lj_epsilon",  # Undocumented classical vdW parameter
    "lj_sigma",  # Undocumented classical vdW parameter
    "mbd_beta",  # TS MBD vdW correction damping parameter
    "scsrad",  # Cutoff radius for dipole-dipole interaction tensor in SCS
    "hitoler",  # Iterative Hirschfeld partitioning tolerance
    "lambda",  # "Spring constant" for magmom constraint calcs.
    "kproj_threshold",  # Threshold for k-point projection scheme
    "maxpwamp",  # Undocumented HF parameter
    "vcutoff",  # Undocumented parameter
    "mdtemp",  # Temperature for AIMD
    "mdgamma",  # Undocumented AIMD parameter
    "mdalpha",  # Undocumented AIMD parameter
    "ofield_kappa",  # Bias potential strength for interface pinning method
    "ofield_q6_near",  # Steinhardt-Nelson Q6 parameters for interface pinning
    "ofield_q6_far",  # Steinhardt-Nelson Q6 parameters for interface pinning
    "ofield_a",  # Target order parameter for interface pinning method
    "pthreshold",  # Don't print timings for routines faster than this value
    "qltol",  # Eigenvalue tolerance for Lanczos iteration (instanton)
    "qdr",  # Step size for building Lanczos matrix & CG (instanton)
    "qmaxmove",  # Max step size (instanton)
    "qdt",  # Timestep for quickmin minimization (instanton)
    "qtpz",  # Temperature (instanton)
    "qftol",  # Tolerance (instanton)
    "nupdown",  # fix spin moment to specified value
]

exp_keys = [
    "ediff",  # stopping-criterion for electronic upd.
    "ediffg",  # stopping-criterion for ionic upd.
    "symprec",  # precession in symmetry routines
    # The next keywords pertain to the VTST add-ons from Graeme Henkelman's
    # group at UT Austin
    "fdstep",  # Finite diference step for IOPT = 1 or 2
]

string_keys = [
    "algo",  # algorithm: Normal (Davidson) | Fast | Very_Fast (RMM-DIIS)
    "gga",  # xc-type: PW PB LM or 91 (LDA if not set)
    "metagga",  #
    "prec",  # Precission of calculation (Low, Normal, Accurate)
    "system",  # name of System
    "precfock",  # FFT grid in the HF related routines
    "radeq",  # Which type of radial equations to use for rel. core calcs.
    "localized_basis",  # Basis to use in CRPA
    "proutine",  # Select profiling routine
]

int_keys = [
    "ialgo",  # algorithm: use only 8 (CG) or 48 (RMM-DIIS)
    "ibrion",  # ionic relaxation: 0-MD 1-quasi-New 2-CG
    "icharg",  # charge: 0-WAVECAR 1-CHGCAR 2-atom 10-const
    "idipol",  # monopol/dipol and quadropole corrections
    "images",  # number of images for NEB calculation
    "iniwav",  # initial electr wf. : 0-lowe 1-rand
    "isif",  # calculate stress and what to relax
    "ismear",  # part. occupancies: -5 Blochl -4-tet -1-fermi 0-gaus >0 MP
    "ispin",  # spin-polarized calculation
    "istart",  # startjob: 0-new 1-cont 2-samecut
    "isym",  # symmetry: 0-nonsym 1-usesym 2-usePAWsym
    "iwavpr",  # prediction of wf.: 0-non 1-charg 2-wave 3-comb
    "kpar",  # k-point parallelization paramater
    "ldauprint",  # 0-silent, 1-occ. matrix written to OUTCAR, 2-1+pot. matrix
    # written
    "ldautype",  # L(S)DA+U: 1-Liechtenstein 2-Dudarev 4-Liechtenstein(LDAU)
    "lmaxmix",  #
    "lorbit",  # create PROOUT
    "maxmix",  #
    "ngx",  # FFT mesh for wavefunctions, x
    "ngxf",  # FFT mesh for charges x
    "ngy",  # FFT mesh for wavefunctions, y
    "ngyf",  # FFT mesh for charges y
    "ngz",  # FFT mesh for wavefunctions, z
    "ngzf",  # FFT mesh for charges z
    "nbands",  # Number of bands
    "nblk",  # blocking for some BLAS calls (Sec. 6.5)
    "nbmod",  # specifies mode for partial charge calculation
    "nelm",  # nr. of electronic steps (default 60)
    "nelmdl",  # nr. of initial electronic steps
    "nelmin",
    "nfree",  # number of steps per DOF when calculting Hessian using
    # finite differences
    "nkred",  # define sub grid of q-points for HF with
    # nkredx=nkredy=nkredz
    "nkredx",  # define sub grid of q-points in x direction for HF
    "nkredy",  # define sub grid of q-points in y direction for HF
    "nkredz",  # define sub grid of q-points in z direction for HF
    "nomega",  # number of frequency points
    "nomegar",  # number of frequency points on real axis
    "npar",  # parallelization over bands
    "nsim",  # evaluate NSIM bands simultaneously if using RMM-DIIS
    "nsw",  # number of steps for ionic upd.
    "nwrite",  # verbosity write-flag (how much is written)
    "vdwgr",  # extra keyword for Andris program
    "vdwrn",  # extra keyword for Andris program
    "voskown",  # use Vosko, Wilk, Nusair interpolation
    # The next keywords pertain to the VTST add-ons from Graeme Henkelman's
    # group at UT Austin
    "ichain",  # Flag for controlling which method is being used (0=NEB,
    # 1=DynMat, 2=Dimer, 3=Lanczos) if ichain > 3, then both
    # IBRION and POTIM are automatically set in the INCAR file
    "iopt",  # Controls which optimizer to use.  for iopt > 0, ibrion = 3
    # and potim = 0.0
    "snl",  # Maximum dimentionality of the Lanczos matrix
    "lbfgsmem",  # Steps saved for inverse Hessian for IOPT = 1 (LBFGS)
    "fnmin",  # Max iter. before adjusting dt and alpha for IOPT = 7 (FIRE)
    "icorelevel",  # core level shifts
    "clnt",  # species index
    "cln",  # main quantum number of excited core electron
    "cll",  # l quantum number of excited core electron
    "ivdw",  # Choose which dispersion correction method to use
    "nbandsgw",  # Number of bands for GW
    "nbandso",  # Number of occupied bands for electron-hole treatment
    "nbandsv",  # Number of virtual bands for electron-hole treatment
    "ncore",  # Number of cores per band, equal to number of cores divided
    # by npar
    "mdalgo",  # Determines which MD method of Tomas Bucko to use
    "nedos",  # Number of grid points in DOS
    "turbo",  # Ewald, 0 = Normal, 1 = PME
    "omegapar",  # Number of groups for response function calc.
    "taupar",  # (Possibly Depricated) Number of groups in real time for response function calc.
    "ntaupar",  # Number of groups in real time for response function calc.
    "antires",  # How to treat antiresonant part of response function
    "magatom",  # Index of atom at which to place magnetic field (NMR)
    "jatom",  # Index of atom at which magnetic moment is evaluated (NMR)
    "ichibare",  # chi_bare stencil size (NMR)
    "nbas",  # Undocumented Bond-Boost parameter (GH patches)
    "rmds",  # Undocumented Bond-Boost parameter (GH patches)
    "ilbfgsmem",  # Number of histories to store for LBFGS (GH patches)
    "vcaimages",  # Undocumented parameter (GH patches)
    "ntemper",  # Undocumented subspace diagonalization param. (GH patches)
    "ncshmem",  # Share memory between this many cores on each process
    "lmaxtau",  # Undocumented MetaGGA parameter (prob. max ang.mom. for tau)
    "kinter",  # Additional finer grid (?)
    "ibse",  # Type of BSE calculation
    "nbseeig",  # Number of BSE wfns to write
    "naturalo",  # Use NATURALO (?)
    "nbandsexact",  # Undocumented parameter
    "nbandsgwlow",  # Number of bands for which shifts are calculated
    "nbandslf",  # Number of bands included in local field effect calc.
    "omegagrid",  # Undocumented parameter
    "telescope",  # Undocumented parameter
    "maxmem",  # Amount of memory to allocate per core in MB
    "nelmhf",  # Number of iterations for HF part (GW)
    "dim",  # Undocumented parameter
    "nkredlf",  # Reduce k-points for local field effects
    "nkredlfx",  # Reduce k-points for local field effects in X
    "nkredlfy",  # Reduce k-points for local field effects in Y
    "nkredlfz",  # Reduce k-points for local field effects in Z
    "lmaxmp2",  # Undocumented parameter
    "switch",  # Undocumented dimer parameter
    "findiff",  # Use forward (1) or central (2) finite difference for dimer
    "engine",  # Undocumented dimer parameter
    "restartcg",  # Undocumented dimer parameter
    "thermostat",  # Deprecated parameter for selecting MD method (use MDALGO)
    "scaling",  # After how many steps velocities should be rescaled
    "shakemaxiter",  # Maximum # of iterations in SHAKE algorithm
    "equi_regime",  # Number of steps to equilibrate for
    "hills_bin",  # Update metadynamics bias after this many steps
    "hills_maxstride",  # Undocumented metadynamics parameter
    "dvvehistory",  # Undocumented parameter
    "ipead",  # Undocumented parameter
    "ngaus",  # Undocumented charge fitting parameter
    "exxoep",  # Undocumented HF parameter
    "fourorbit",  # Undocumented HF parameter
    "model_gw",  # Undocumented HF parameter
    "hflmax",  # Maximum L quantum number for HF calculation
    "lmaxfock",  # Maximum L quantum number for HF calc. (same as above)
    "lmaxfockae",  # Undocumented HF parameter
    "nmaxfockae",  # Undocumented HF parameter
    "nblock_fock",  # Undocumented HF parameter
    "idiot",  # Determines which warnings/errors to print
    "nrmm",  # Number of RMM-DIIS iterations
    "mremove",  # Undocumented mixing parameter
    "inimix",  # Undocumented mixing parameter
    "mixpre",  # Undocumented mixing parameter
    "nelmall",  # Undocumented parameter
    "nblock",  # How frequently to write data
    "kblock",  # How frequently to write data
    "npaco",  # Undocumented pair correlation function parameter
    "lmaxpaw",  # Max L quantum number for on-site charge expansion
    "irestart",  # Undocumented parameter
    "nreboot",  # Undocumented parameter
    "nmin",  # Undocumented parameter
    "nlspline",  # Undocumented parameter
    "ispecial",  # "Select undocumented and unsupported special features"
    "rcrep",  # Number of steps between printing relaxed core info
    "rcndl",  # Wait this many steps before updating core density
    "rcstrd",  # Relax core density after this many SCF steps
    "vdw_idampf",  # Select type of damping function for TS vdW
    "i_constrained_m",  # Select type of magmom. constraint to use
    "igpar",  # "G parallel" direction for Berry phase calculation
    "nppstr",  # Number of kpts in "igpar' direction for Berry phase calc.
    "nbands_out",  # Undocumented QP parameter
    "kpts_out",  # Undocumented QP parameter
    "isp_out",  # Undocumented QP parameter
    "nomega_out",  # Undocumented QP parameter
    "maxiter_ft",  # Max iterations for sloppy Remez algorithm
    "nmaxalt",  # Max sample points for alternant in Remez algorithms
    "itmaxlsq",  # Max iterations in LSQ search algorithm
    "ndatalsq",  # Number of sample points for LSQ search algorithm
    "ncore_in_image1",  # Undocumented parameter
    "kimages",  # Undocumented parameter
    "ncores_per_band",  # Undocumented parameter
    "maxlie",  # Max iterations in CRPA diagonalization routine
    "ncrpalow",  # Undocumented CRPA parameter
    "ncrpahigh",  # Undocumented CRPA parameter
    "nwlow",  # Undocumented parameter
    "nwhigh",  # Undocumented parameter
    "nkopt",  # Number of k-points to include in Optics calculation
    "nkoffopt",  # K-point "counter offset" for Optics
    "nbvalopt",  # Number of valence bands to write in OPTICS file
    "nbconopt",  # Number of conduction bands to write in OPTICS file
    "ch_nedos",  # Number dielectric function calculation grid points for XAS
    "plevel",  # No timings for routines with "level" higher than this
    "qnl",  # Lanczos matrix size (instanton)
]

bool_keys = [
    "addgrid",  # finer grid for augmentation charge density
    "kgamma",  # The generated kpoint grid (from KSPACING) is either
    # centred at the $\Gamma$
    # point (e.g. includes the $\Gamma$ point)
    # (KGAMMA=.TRUE.)
    "laechg",  # write AECCAR0/AECCAR1/AECCAR2
    "lasph",  # non-spherical contributions to XC energy (and pot for
    # VASP.5.X)
    "lasync",  # overlap communcation with calculations
    "lcharg",  #
    "lcorr",  # Harris-correction to forces
    "ldau",  # L(S)DA+U
    "ldiag",  # algorithm: perform sub space rotation
    "ldipol",  # potential correction mode
    "lelf",  # create ELFCAR
    "lepsilon",  # enables to calculate and to print the BEC tensors
    "lhfcalc",  # switch to turn on Hartree Fock calculations
    "loptics",  # calculate the frequency dependent dielectric matrix
    "lpard",  # evaluate partial (band and/or k-point) decomposed charge
    # density
    "lplane",  # parallelisation over the FFT grid
    "lscalapack",  # switch off scaLAPACK
    "lscalu",  # switch of LU decomposition
    "lsepb",  # write out partial charge of each band separately?
    "lsepk",  # write out partial charge of each k-point separately?
    "lthomas",  #
    "luse_vdw",  # Invoke vdW-DF implementation by Klimes et. al
    "lvdw",  # Invoke DFT-D2 method of Grimme
    "lvhar",  # write Hartree potential to LOCPOT (vasp 5.x)
    "lvtot",  # create WAVECAR/CHGCAR/LOCPOT
    "lwave",  #
    # The next keywords pertain to the VTST add-ons from Graeme Henkelman's
    # group at UT Austin
    "lclimb",  # Turn on CI-NEB
    "ltangentold",  # Old central difference tangent
    "ldneb",  # Turn on modified double nudging
    "lnebcell",  # Turn on SS-NEB
    "lglobal",  # Optmize NEB globally for LBFGS (IOPT = 1)
    "llineopt",  # Use force based line minimizer for translation (IOPT = 1)
    "lbeefens",  # Switch on print of BEE energy contributions in OUTCAR
    "lbeefbas",  # Switch off print of all BEEs in OUTCAR
    "lcalcpol",  # macroscopic polarization (vasp5.2). 'lcalceps'
    "lcalceps",  # Macroscopic dielectric properties and Born effective charge
    # tensors (vasp 5.2)
    "lvdw",  # Turns on dispersion correction
    "lvdw_ewald",  # Turns on Ewald summation for Grimme's DFT-D2 and
    # Tkatchenko and Scheffler's DFT-TS dispersion correction
    "lspectral",  # Use the spectral method to calculate independent particle
    # polarizability
    "lrpa",  # Include local field effects on the Hartree level only
    "lwannier90",  # Switches on the interface between VASP and WANNIER90
    "lsorbit",  # Enable spin-orbit coupling
    "lsol",  # turn on solvation for Vaspsol
    "lautoscale",  # automatically calculate inverse curvature for VTST LBFGS
    "interactive",  # Enables interactive calculation for VaspInteractive
    "lauger",  # Perform Auger calculation (Auger)
    "lauger_eeh",  # Calculate EEH processes (Auger)
    "lauger_ehh",  # Calculate EHH processes (Auger)
    "lauger_collect",  # Collect wfns before looping over k-points (Auger)
    "lauger_dhdk",  # Auto-determine E. window width from E. derivs. (Auger)
    "lauger_jit",  # Distribute wavefunctions for k1-k4 (Auger)
    "orbitalmag",  # Enable orbital magnetization (NMR)
    "lchimag",  # Use linear response for shielding tensor (NMR)
    "lwrtcur",  # Write response of current to mag. field to file (NMR)
    "lnmr_sym_red",  # Reduce symmetry for finite difference (NMR)
    "lzora",  # Use ZORA approximation in linear-response NMR (NMR)
    "lbone",  # Use B-component in AE one-center terms for LR NMR (NMR)
    "lmagbloch",  # Use Bloch summations to obtain orbital magnetization (NMR)
    "lgauge",  # Use gauge transformation for zero moment terms (NMR)
    "lbfconst",  # Use constant B-field with sawtooth vector potential (NMR)
    "nucind",  # Use nuclear independent calculation (NMR)
    "lnicsall",  # Use all grid points for 'nucind' calculation (NMR)
    "llraug",  # Use two-center corrections for induced B-field (NMR)
    "lbbm",  # Undocumented Bond-Boost parameter (GH patches)
    "lnoncollinear",  # Do non-collinear spin polarized calculation
    "bfgsdfp",  # Undocumented BFGS parameter (GH patches)
    "linemin",  # Use line minimization (GH patches)
    "ldneborg",  # Undocumented NEB parameter (GH patches)
    "dseed",  # Undocumented dimer parameter (GH patches)
    "linteract",  # Undocumented parameter (GH patches)
    "lmpmd",  # Undocumented parameter (GH patches)
    "ltwodim",  # Makes stress tensor two-dimensional (GH patches)
    "fmagflag",  # Use force magnitude as convergence criterion (GH patches)
    "ltemper",  # Use subspace diagonalization (?) (GH patches)
    "qmflag",  # Undocumented FIRE parameter (GH patches)
    "lmixtau",  # Undocumented MetaGGA parameter
    "ljdftx",  # Undocumented VASPsol parameter (VASPsol)
    "lrhob",  # Write the bound charge density (VASPsol)
    "lrhoion",  # Write the ionic charge density (VASPsol)
    "lnabla",  # Undocumented parameter
    "linterfast",  # Interpolate in K using linear response routines
    "lvel",  # Undocumented parameter
    "lrpaforce",  # Calculate RPA forces
    "lhartree",  # Use IP approx. in BSE (testing only)
    "ladder",  # Use ladder diagrams
    "lfxc",  # Use approximate ladder diagrams
    "lrsrpa",  # Undocumented parameter
    "lsingles",  # Calculate HF singles
    "lfermigw",  # Iterate Fermi level
    "ltcte",  # Undocumented parameter
    "ltete",  # Undocumented parameter
    "ltriplet",  # Undocumented parameter
    "lfxceps",  # Undocumented parameter
    "lfxheg",  # Undocumented parameter
    "l2order",  # Undocumented parameter
    "lmp2lt",  # Undocumented parameter
    "lgwlf",  # Undocumented parameter
    "lusew",  # Undocumented parameter
    "selfenergy",  # Undocumented parameter
    "oddonlygw",  # Avoid gamma point in response function calc.
    "evenonlygw",  # Avoid even points in response function calc.
    "lspectralgw",  # More accurate self-energy calculation
    "ch_lspec",  # Calculate matrix elements btw. core and conduction states
    "fletcher_reeves",  # Undocumented dimer parameter
    "lidm_selective",  # Undocumented dimer parameter
    "lblueout",  # Write output of blue-moon algorithm
    "hills_variable_w",  # Enable variable-width metadynamics bias
    "dvvminus",  # Undocumented parameter
    "lpead",  # Calculate cell-periodic orbital derivs. using finite diff.
    "skip_edotp",  # Skip updating elec. polarization during scf
    "skip_scf",  # Skip calculation w/ local field effects
    "lchgfit",  # Turn on charge fitting
    "lgausrc",  # Undocumented charge fitting parameter
    "lstockholder",  # Enable ISA charge fitting (?)
    "lsymgrad",  # Restore symmetry of gradient (HF)
    "lhfone",  # Calculate one-center terms (HF)
    "lrscor",  # Include long-range correlation (HF)
    "lrhfcalc",  # Include long-range HF (HF)
    "lmodelhf",  # Model HF calculation (HF)
    "shiftred",  # Undocumented HF parameter
    "hfkident",  # Undocumented HF parameter
    "oddonly",  # Undocumented HF parameter
    "evenonly",  # Undocumented HF parameter
    "lfockaedft",  # Undocumented HF parameter
    "lsubsrot",  # Enable subspace rotation diagonalization
    "mixfirst",  # Mix before diagonalization
    "lvcader",  # Calculate derivs. w.r.t. VCA parameters
    "lcompat",  # Enable "full compatibility"
    "lmusic",  # "Joke" parameter
    "ldownsample",  # Downsample WAVECAR to fewer k-points
    "lscaaware",  # Disable ScaLAPACK for some things but not all
    "lorbitalreal",  # Undocumented parameter
    "lmetagga",  # Undocumented parameter
    "lspiral",  # Undocumented parameter
    "lzeroz",  # Undocumented parameter
    "lmono",  # Enable "monopole" corrections
    "lrelcore",  # Perform relaxed core calculation
    "lmimicfc",  # Mimic frozen-core calcs. for relaxed core calcs.
    "lmatchrw",  # Match PS partial waves at RWIGS? (otherwise PAW cutoff)
    "ladaptelin",  # Linearize core state energies to avoid divergences
    "lonlysemicore",  # Only linearize semi-core state energies
    "gga_compat",  # Enable backwards-compatible symmetrization of GGA derivs.
    "lrelvol",  # Undocumented classical vdW parameter
    "lj_only",  # Undocumented classical vdW parameter
    "lvdwscs",  # Include self-consistent screening in TS vdW correction
    "lcfdm",  # Use coupled fluctuating dipoles model for TS vdW
    "lvdw_sametype",  # Include interactions between atoms of the same type
    "lrescaler0",  # Rescale damping parameters in SCS vdW correction
    "lscsgrad",  # Calculate gradients for TS+SCS vdW correction energies
    "lvdwexpansion",  # Write 2-6 body contribs. to MBD vdW correction energy
    "lvdw_relvolone",  # Undocumented classical vdW parameter
    "lberry",  # Enable Berry-phase calculation
    "lpade_fit",  # Undocumented QP parameter
    "lkproj",  # Enable projection onto k-points
    "l_wr_moments",  # Undocumented parameter
    "l_wr_density",  # Undocumented parameter
    "lkotani",  # Undocumented parameter
    "ldyson",  # Undocumented parameter
    "laddherm",  # Undocumented parameter
    "lcrpaplot",  # Plot bands used in CRPA response func. calc.
    "lplotdis",  # Plot disentangled bands in CRPA response func. calc.
    "ldisentangle",  # Disentangle bands in CRPA
    "lweighted",  # "Weighted" CRPA approach
    "luseorth_lcaos",  # Use orthogonalized LCAOs in CRPA
    "lfrpa",  # Use full RPA in CRPA
    "lregularize",  # Regularize projectors in CRPA
    "ldrude",  # Include Drude term in CRPA
    "ldmatrix",  # Undocumented parameter
    "lefg",  # Calculate electric field gradient at atomic nuclei
    "lhyperfine",  # Enable Hyperfine calculation
    "lwannier",  # Enable Wannier interface
    "localize",  # Undocumented Wannier parameter
    "lintpol_wpot",  # Interpolate WPOT for Wannier
    "lintpol_orb",  # Interpolate orbitals for Wannier
    "lintpol_kpath",  # Interpolate bandstructure on given kpath for Wannier
    "lintpol_kpath_orb",  # Interpolate orbitals on given kpath for Wannier
    "lread_eigenvalues",  # Use Eigenvalues from EIGENVALUES.INT file
    "lintpol_velocity",  # Interpolate electron velocity for Wannier
    "lintpol_conductivity",  # Interpolate conductivity for Wannier
    "lwannierinterpol",  # Undocumented Wannier parameter
    "wanproj",  # Undocumented Wannier parameter
    "lorbmom",  # Undocumented LDA+U parameter
    "lwannier90_run",  # Undocumented WANNIER90 parameter
    "lwrite_wanproj",  # Write UWAN files for WANNIER90
    "lwrite_unk",  # Write UNK files for WANNIER90
    "lwrite_mmn_amn",  # Write MMN and AMN files for WANNIER90
    "lread_amn",  # Read AMN files instead of recomputing (WANNIER90)
    "lrhfatm",  # Undocumented HF parameter
    "lvpot",  # Calculate unscreened potential
    "lwpot",  # Calculate screened potential
    "lwswq",  # Undocumented parameter
    "pflat",  # Only print "flat" timings to OUTCAR
    "qifcg",  # Use CG instead of quickmin (instanton)
    "qdo_ins",  # Find instanton
    "qdo_pre",  # Calculate prefactor (instanton)
    # The next keyword pertains to the periodic NBO code of JR Schmidt's group
    # at UW-Madison (https://github.com/jrschmidt2/periodic-NBO)
    "lnbo",  # Enable NBO analysis
]

list_int_keys = [
    "iband",  # bands to calculate partial charge for
    "kpuse",  # k-point to calculate partial charge for
    "ldaul",  # DFT+U parameters, overruled by dict key 'ldau_luj'
    "random_seed",  # List of ints used to seed RNG for advanced MD routines
    # (Bucko)
    "auger_bmin_eeh",  # 4 ints | Various undocumented parameters for Auger
    "auger_bmax_eeh",  # 4 ints | calculations
    "auger_bmin_ehh",  # 4 ints |
    "auger_bmax_ehh",  # 4 ints |
    "balist",  # nbas ints | Undocumented Bond-Boost parameter (GH patches)
    "kpoint_bse",  # 4 ints | Undocumented parameter
    "nsubsys",  # <=3 ints | Last atom # for each of up to 3 thermostats
    "vdw_refstate",  # ntyp ints | Undocumented classical vdW parameter
    "vdw_mbd_size",  # 3 ints | Supercell size for TS MBD vdW correction
    "nbands_index",  # nbands_out ints | Undocumented QP parameter
    "kpts_index",  # kpts_out ints | Undocumented QP parameter
    "isp_index",  # isp_out ints | Undocumented QP parameter
    "nomega_index",  # nomega_out ints | Undocumented QP parameter
    "ntarget_states",  # nbands ints | Undocumented CRPA parameter
    "wanproj_i",  # nions ints | Undocumented Wannier parameter
    "wanproj_l",  # ? ints | Undocumented Wannier parameter
]

list_bool_keys = [
    "lattice_constraints",  # 3 bools | Undocumented advanced MD parameter
    "lrctype",  # ntyp bools | Enable relaxed-core calc. for these atoms
    "lvdw_onecell",  # 3 bools | Enable periodicity in A, B, C vector for vdW
]

list_float_keys = [
    "dipol",  # center of cell for dipol
    "eint",  # energy range to calculate partial charge for
    "ferwe",  # Fixed band occupation (spin-paired)
    "ferdo",  # Fixed band occupation (spin-plarized)
    "magmom",  # initial magnetic moments e.g. List[float]
    "ropt",  # number of grid points for non-local proj in real space
    "rwigs",  # Wigner-Seitz radii
    "ldauu",  # ldau parameters, has potential to redundant w.r.t. dict
    "ldauj",  # key 'ldau_luj', but 'ldau_luj' can't be read direct from
    # the INCAR (since it needs to know information about atomic
    # species. In case of conflict 'ldau_luj' gets written out
    # when a calculation is set up
    "vdw_c6",  # List of floats of C6 parameters (J nm^6 mol^-1) for each
    # species (DFT-D2 and DFT-TS)
    "vdw_c6au",  # List of floats of C6 parameters (a.u.) for each species
    # (DFT-TS)
    "vdw_r0",  # List of floats of R0 parameters (angstroms) for each
    # species (DFT-D2 and DFT-TS)
    "vdw_r0au",  # List of floats of R0 parameters (a.u.) for each species
    # (DFT-TS)
    "vdw_alpha",  # List of floats of free-atomic polarizabilities for each
    # species (DFT-TS)
    "langevin_gamma",  # List of floats for langevin friction coefficients
    "auger_emin_eeh",  # 4 floats | Various undocumented parameters for Auger
    "auger_emax_eeh",  # 4 floats | calculations
    "auger_emin_ehh",  # 4 floats |
    "auger_emax_ehh",  # 4 floats |
    "avecconst",  # 3 floats | magnitude of magnetic moment (NMR)
    "magdipol",  # 3 floats | magnitude of magnetic dipole (NMR)
    "bconst",  # 3 floats | magnitude of constant magnetic field (NMR)
    "magpos",  # 3 floats | position for magnetic moment w/ 'nucind' (NMR)
    "bext",  # 3 floats | Undocumented (probably external magnetic field)
    "core_c",  # ntyp floats | pseudo-core charge magnitude (VASPsol)
    "sigma_rc_k",  # ntyp floats | width of pseudo-core gaussians (VASPsol)
    "darwinr",  # ntypd (?) floats | Undocumented parameter
    "darwinv",  # ntypd (?) floats | Undocumented parameter
    "dummy_k",  # ? floats | Force const. connecting dummy atoms to sys.
    "dummy_r0",  # ? floats | Minimum dist., ang., etc. for dummy atom DOFs
    "dummy_positions",  # 3 floats | Position of dummy atom(s?)
    "psubsys",  # <=3 floats | Coll. prob. for each of up to 3 thermostats
    "tsubsys",  # <=3 floats | Temp. for each of up to 3 thermostats
    "increm",  # ? floats | Undocumented advanced MD parameter
    "value_min",  # ? floats | Undocumented advanced MD parameter
    "value_max",  # ? floats | Undocumented advanced MD parameter
    "hills_position",  # ? floats | Dummy particle(s) pos. for metadynamics
    "hills_velocity",  # ? floats | Dummy particle(s) vel. for metadynamics
    "spring_k",  # ? floats | Spring constant for harmonic constraints
    "spring_r0",  # ? floats | Spring minima for harmonic constraints
    "spring_v0",  # ? floats | Initial velocity of harmonic constraints
    "hills_wall_lower",  # ? floats | Undocumented metadynamics parameter
    "hills_wall_upper",  # ? floats | Undocumented metadynamics parameter
    "efield_pead",  # 3 floats | homogeneous electric field for PEAD calc.
    "zct",  # ? floats | Undocumented charge fitting parameter
    "rgaus",  # ? floats | Undocumented charge fitting parameter
    "hfalpha",  # 10 floats | Undocumented HF parameter
    "mcalpha",  # 10 floats | Undocumented HF parameter
    "saxis",  # 3 floats | Coordinate for collinear spin calculations
    "vca",  # ? floats | Atom weight for VCA calculations
    "stm",  # 7 floats | "range for STM data"
    "qspiral",  # 3 floats | Undocumented parameter
    "external_stress",  # 6 floats | Target stress (adds w/ external_pressure)
    "m_constr",  # 3*nions floats | Local magmom assigned to each spin DOF
    "quad_efg",  # ntyp floats | Nuclear quadrupole moments
    "ngyromag",  # ntyp floats | Nuclear gyromagnetic ratios
    "rcrhocut",  # ntyp floats | Core density cutoff rad. for HF relcore calc
    "ofield_k",  # 3 floats | Undocumented parameter
    "paripot",  # ? floats | Undocumented parameter
    "smearings",  # ? floats | ismear,sigma smearing params to loop over
    "wanproj_e",  # 2 floats | Undocumented Wannier parameter
]

special_keys = [
    "lreal",  # non-local projectors in real space
]

dict_keys = [
    "ldau_luj",  # dictionary with L(S)DA+U parameters, e.g. {'Fe':{'L':2,
    # 'U':4.0, 'J':0.9}, ...}
]

keys: List[str] = [
    # 'NBLOCK' and KBLOCK       inner block; outer block
    # 'NPACO' and APACO         distance and nr. of slots for P.C.
    # 'WEIMIN, EBREAK, DEPER    special control tags
]

# Recomended POTCARS tags
potcar_path = Path(os.environ.get("POTCAR_PATH_PREFIX", ""))
reccomended_potcars = (
    {
        "H": "H",
        "He": "He",
        "Li": "Li_sv",
        "Be": "Be",
        "B": "B",
        "C": "C",
        "N": "N",
        "O": "O",
        "F": "F",
        "Ne": "Ne",
        "Na": "Na_pv",
        "Mg": "Mg",
        "Al": "Al",
        "Si": "Si",
        "P": "P",
        "S": "S",
        "Cl": "Cl",
        "Ar": "Ar",
        "K": "K_sv",
        "Ca": "Ca_sv",
        "Sc": "Sc_sv",
        "Ti": "Ti_sv",
        "V": "V_sv",
        "Cr": "Cr_pv",
        "Mn": "Mn_pv",
        "Fe": "Fe",
        "Co": "Co",
        "Ni": "Ni",
        "Cu": "Cu",
        "Zn": "Zn",
        "Ga": "Ga_d",
        "Ge": "Ge_d",
        "As": "As",
        "Se": "Se",
        "Br": "Br",
        "Kr": "Kr",
        "Rb": "Rb_sv",
        "Sr": "Sr_sv",
        "Y": "Y_sv",
        "Zr": "Zr_sv",
        "Nb": "Nb_sv",
        "Mo": "Mo_sv",
        "Tc": "Tc_pv",
        "Ru": "Ru_pv",
        "Rh": "Rh_pv",
        "Pd": "Pd",
        "Ag": "Ag",
        "Cd": "Cd",
        "In": "In_d",
        "Sn": "Sn_d",
        "Sb": "Sb",
        "Te": "Te",
        "I": "I",
        "Xe": "Xe",
        "Cs": "Cs_sv",
        "Ba": "Ba_sv",
        "La": "La",
        "Ce": "Ce",
        "Pr": "Pr_3",
        "Nd": "Nd_3",
        "Pm": "Pm_3",
        "Sm": "Sm_3",
        "Eu": "Eu_2",
        "Gd": "Gd_3",
        "Tb": "Tb_3",
        "Dy": "Dy_3",
        "Ho": "Ho_3",
        "Er": "Er_3",
        "Tm": "Tm_3",
        "Yb": "Yb_2",
        "Lu": "Lu_3",
        "Hf": "Hf_pv",
        "Ta": "Ta_pv",
        "W": "W_sv",
        "Re": "Re",
        "Os": "Os",
        "Ir": "Ir",
        "Pt": "Pt",
        "Au": "Au",
        "Hg": "Hg",
        "Tl": "Tl_d",
        "Pb": "Pb_d",
        "Bi": "Bi_d",
        "Po": "Po_d",
        "At": "At",
        "Rn": "Rn",
        "Fr": "Fr_sv",
        "Ra": "Ra_sv",
        "Ac": "Ac",
        "Th": "Th",
        "Pa": "Pa",
        "U": "U",
        "Np": "Np",
        "Pu": "Pu",
        "Am": "Am",
        "Cm": "Cm",
    },
)
mp_potcars = {
    "Ac": "Ac",
    "Ag": "Ag",
    "Al": "Al",
    "Ar": "Ar",
    "As": "As",
    "Au": "Au",
    "B": "B",
    "Ba": "Ba_sv",
    "Be": "Be_sv",
    "Bi": "Bi",
    "Br": "Br",
    "C": "C",
    "Ca": "Ca_sv",
    "Cd": "Cd",
    "Ce": "Ce",
    "Cl": "Cl",
    "Co": "Co",
    "Cr": "Cr_pv",
    "Cs": "Cs_sv",
    "Cu": "Cu_pv",
    "Dy": "Dy_3",
    "Er": "Er_3",
    "Eu": "Eu",
    "F": "F",
    "Fe": "Fe_pv",
    "Ga": "Ga_d",
    "Gd": "Gd",
    "Ge": "Ge_d",
    "H": "H",
    "He": "He",
    "Hf": "Hf_pv",
    "Hg": "Hg",
    "Ho": "Ho_3",
    "I": "I",
    "In": "In_d",
    "Ir": "Ir",
    "K": "K_sv",
    "Kr": "Kr",
    "La": "La",
    "Li": "Li_sv",
    "Lu": "Lu_3",
    "Mg": "Mg_pv",
    "Mn": "Mn_pv",
    "Mo": "Mo_pv",
    "N": "N",
    "Na": "Na_pv",
    "Nb": "Nb_pv",
    "Nd": "Nd_3",
    "Ne": "Ne",
    "Ni": "Ni_pv",
    "Np": "Np",
    "O": "O",
    "Os": "Os_pv",
    "P": "P",
    "Pa": "Pa",
    "Pb": "Pb_d",
    "Pd": "Pd",
    "Pm": "Pm_3",
    "Pr": "Pr_3",
    "Pt": "Pt",
    "Pu": "Pu",
    "Rb": "Rb_sv",
    "Re": "Re_pv",
    "Rh": "Rh_pv",
    "Ru": "Ru_pv",
    "S": "S",
    "Sb": "Sb",
    "Sc": "Sc_sv",
    "Se": "Se",
    "Si": "Si",
    "Sm": "Sm_3",
    "Sn": "Sn_d",
    "Sr": "Sr_sv",
    "Ta": "Ta_pv",
    "Tb": "Tb_3",
    "Tc": "Tc_pv",
    "Te": "Te",
    "Th": "Th",
    "Ti": "Ti_pv",
    "Tl": "Tl_d",
    "Tm": "Tm_3",
    "U": "U",
    "V": "V_pv",
    "W": "W",
    "Xe": "Xe",
    "Y": "Y_sv",
    "Yb": "Yb_2",
    "Zn": "Zn",
    "Zr": "Zr_sv",
}
