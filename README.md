# neath
NEATH is a framework for post-processing hydrodynamical simulations with a time-dependent astrochemical network, allowing the abundances of hundreds of molecules to be calculated at low computational cost. The code in this repository is an implementation based on a modified version of [UCLCHEM](https://uclchem.github.io/) ([Holdship et al. 2017](https://ui.adsabs.harvard.edu/abs/2017AJ....154...38H/abstract)), which has significantly diverged from the main branch.

The code treats each 'depth' point as an individual Lagrangian parcel of gas, which evolves completely independently of any external material. This has required some changes to how the code treats shielding from the external radiation field. The required input data are the evolution of gas properties with time for each gas parcel, presumed to be an SPH/tracer particle from a hydrodynamical simulation. An example is given in the 'data' folder. Note that while an unlimited numnber of particles or timesteps can be given, the timesteps must be the same for all particles. The chemical evolution is tied to the timesteps in the input data; the code will update the physical properties for each particle at the current time, evolve the chemistry up to the next timestep, then update the physical properties to their new values. This can lead to serious errors if the time resolution is too coarse compared to the rate of change of the system.

While the example data file includes the gas temperature and shielding column densities of total hydrogen, H2 and CO, in addition to properties not used by the code such as position and velocity, these can also be specified manually in postprocess.f90 (for example, enforcing an isothermal temperature of 10 K). The only absolute requirement is that the code is provided with a list of times, which it will attempt to evolve the chemistry through while updating the physical properties according to the routines in postprocess.f90. See [Priestley, Clark & Whitworth (2023)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.519.6392P/abstract) for an example where all required properties are given as a function of gas density.
