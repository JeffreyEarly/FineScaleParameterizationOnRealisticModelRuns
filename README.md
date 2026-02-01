# Pathways and Parameterizations Proposal Repo

The `GarrettMunkShearSpectrumResolutionTest.m` was used to produce some vertical shear spectrum and see how they dependended on resolution and the shape of the GM spectrum.

`TimeStepGarrettMunk.m` initialized a hydrostatic or Boussinesq transform, adds a GM spectrum, then time-steps for two days. This is supposed to serve as an adjusted IGW state for a particular resolution and transform.

`RestartWithCyclogeostrophicEddy.m` uses the output from `TimeStepGarrettMunk.m` and then perhaps removes the IGW field (or not), adds inertial oscillations at the surface (or not) and added a cyclogeostrophic cyclone or anticyclone (or not) and then time steps forward for 15 days.

`HovmullerWaveEnergyAndDissipation.m` creates a Hov-mueller diagram for the wave energy and dissipation of any of these simulations.

`MakeSpatialSpectralFluxMovie.m` is one silly attempt to make a movie showing both the physical and spectral evolution.

`ShowEnergyFlux2D.m` is a started script for showing the spectral fluxes.

`ShowVariousCyprusEddyFigures.m` shows cross-sections, vertical shear spectra, and other misc stuff.
