wvd = WVDiagnostics("fine-scale-hydrostatic-50km-128-222-anticyclone-7cms-cyclogeostrophic.nc");

%%
svv = wvt.forcingWithName("adaptive damping");
NoDampMask = (wvt.Kh < (svv.k_damp+svv.k_no_damp)/2) & (wvt.J < (svv.j_damp + svv.j_no_damp)/2);

i = 1;
flowComponents(i) = wvt.flowComponentWithName("wave");
flowComponents(i).name = "wave";
flowComponents(i).maskAp = flowComponents(i).maskAp & NoDampMask;
flowComponents(i).maskAm = flowComponents(i).maskAm & NoDampMask;

i = i+1;
flowComponents(i) = wvt.flowComponentWithName("inertial");
flowComponents(i).maskAp = flowComponents(i).maskAp & NoDampMask;
flowComponents(i).maskAm = flowComponents(i).maskAm & NoDampMask;

i = i+1;
flowComponents(i) = wvt.flowComponentWithName("geostrophic") + wvt.flowComponentWithName("mda");
flowComponents(i).name = "geostrophic";
flowComponents(i).maskA0 = flowComponents(i).maskA0 & NoDampMask;

i = i+1;
flowComponents(i) = wvt.flowComponentWithName("wave") + wvt.flowComponentWithName("inertial");
flowComponents(i).name = "damped wave";
flowComponents(i).maskAp = flowComponents(i).maskAp & ~NoDampMask;
flowComponents(i).maskAm = flowComponents(i).maskAm & ~NoDampMask;

i = i+1;
flowComponents(i) = wvt.flowComponentWithName("geostrophic") + wvt.flowComponentWithName("mda");
flowComponents(i).name = "damped geostrophic";
flowComponents(i).maskA0 = flowComponents(i).maskA0 & ~NoDampMask;

wvd.createReservoirGroup(name = "reservoir-damped-wave-geo-io", flowComponents=flowComponents);

%%
order = ["inertial_oscillation", "geostrophic", "wave", "damped_geostrophic", "damped_wave"];
% order = ["inertial_oscillation", "wave", "geostrophic", "damped_geostrophic", "damped_wave"];
wvd.plotSourcesSinksForReservoirGroup(timeIndices=1:61,name="reservoir-damped-wave-geo-io",customReservoirOrder=order);

%%
wvd.plotSourcesSinksForReservoirGroup(timeIndices=1:61,name="reservoir-damped-wave-geo");