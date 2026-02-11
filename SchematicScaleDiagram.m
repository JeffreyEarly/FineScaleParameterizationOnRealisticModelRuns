%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Making a schematic diagram of resolved energy flux and Ozmidov scales. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basedir = '/Volumes/SanDiskExtremePro/research/JeffreyEarly/FineScaleData/';
% basedir = '/Users/jearly/Dropbox/FineScaleData/';

% filename = '7cms/fine-scale-hydrostatic-50km-128-222-cyclone-7cms-cyclogeostrophic.nc';
% runname = '7cm/s, 50km, 128^2 x 222 cyclone';
% filename = '7cms/fine-scale-hydrostatic-50km-128-222-anticyclone-7cms-cyclogeostrophic.nc';
% runname = '7cm/s, 50km, 128^2 x 222 anticyclone';
filename = '25cms/fine-scale-hydrostatic-50km-256-443-anticyclone-cyclogeostrophic.nc';
runname = '25cm/s, 50km, 256^2 x 222 anticyclone';
% filename = 'fine-scale-hydrostatic-50km-256-443-one-half-dealias.nc';
% runname = '50km, 256^2 x 443 anticyclone';

% load the file
wvd = WVDiagnostics(fullfile(basedir,filename));
wvt = wvd.wvt;

% create figure directory
shouldExportFigures = false;
figureFolder = "./figures";
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder)
end

t = wvd.t_diag;

% set analysis time
wvd.iTime = length(t);

%% compute energy flux

% line colors
C = orderedcolors("gem"); 
colorDictionary = dictionary("vertical_scalar_diffusivity",{C(1,:)});
colorDictionary{"adaptive_damping"} = C(2,:);

timeIndices = (length(t)-4*4):length(t);
% timeIndices = length(t);

% load spectral energy fluxes
energy_fluxes = wvd.exactEnergyFluxesTemporalAverage(timeIndices=timeIndices);

% set forcing fluxes to plot
fluxesOfInterest = {"vertical_scalar_diffusivity","adaptive_damping"};

% setup forcing flux structure
clear forcing_fluxes;
for iFlux=1:length(fluxesOfInterest)
    forcing_fluxes(iFlux).color=colorDictionary{fluxesOfInterest{iFlux}};
    forcing_fluxes(iFlux).flux = energy_fluxes([energy_fluxes.name] == fluxesOfInterest{iFlux}).te/wvd.flux_scale;
    forcing_fluxes(iFlux).relativeAmplitude = 1.0;
    forcing_fluxes(iFlux).alpha = 1.0;
    forcing_fluxes(iFlux).fancyName = energy_fluxes([energy_fluxes.name] == fluxesOfInterest{iFlux}).fancyName;
end

% exact total flux
clear flux_advective
flux_advective.flux = energy_fluxes([energy_fluxes.name] == "nonlinear_advection").te/wvd.flux_scale;

% normalize fluxes
maxAmplitude = max(arrayfun( @(v) max(abs(v.flux(:))), forcing_fluxes));
for iFlux=1:length(forcing_fluxes)
    forcing_fluxes(iFlux).relativeAmplitude = max(abs((forcing_fluxes(iFlux).flux(:))))/maxAmplitude;
end


%% compute dissipation rate

int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);

svv = wvt.forcingWithName("adaptive damping");
Fp = sqrt(-wvt.uvMax * svv.damp) .* wvt. Ap;
Fm = sqrt(-wvt.uvMax * svv.damp) .* wvt. Am;
F0 = sqrt(-wvt.uvMax * svv.damp) .* wvt. A0;

Fu = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.phase.*(Fp) + wvt.UAm.*wvt.conjPhase.*(Fm),A0=wvt.UA0.*(F0));
Fv = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*wvt.phase.*(Fp) + wvt.VAm.*wvt.conjPhase.*(Fm),A0=wvt.VA0.*(F0));
Feta = wvt.transformToSpatialDomainWithG(Apm=wvt.NAp.*wvt.phase.*(Fp) + wvt.NAm.*wvt.conjPhase.*(Fm),A0=wvt.NA0.*(F0));

te_diss = Fu .* Fu + Fv .* Fv + shiftdim(wvt.N2,-2) .* Feta .* Feta;
te_diss_tot = int_vol( te_diss );

% units of te_diss are [m^2 s^{-3}]
% convert to GM/yr by dividing by wvd.flux_scale_units

% Ozmidov scale, m_O in Kunze, JPO, 2019
mOzmidov = sqrt( permute(wvt.N2.^(3/2), [2 3 1]) ./ te_diss );
LOzmidov = 2*pi./mOzmidov;

% rolloff wavenumber, m_c or m_f in Kunze, JPO, 2019
mCoriolis = sqrt( permute(wvt.f*wvt.N2, [2 3 1]) ./ te_diss );
LCoriolis = 2*pi./mCoriolis;

% horizontal Coriolis wavenumber
kCoriolis = sqrt( wvt.f^3 ./ te_diss );
LhCoriolis = 2*pi./kCoriolis;

fprintf("Maximum Ozmidov length scale %0.1f m.\n", max(LOzmidov(:)))
fprintf("Maximum vertical rolloff length scale %0.1f m.\n", max(LCoriolis(:)))
fprintf("Maximum horizontal rolloff length scale %0.1f m.\n", max(LhCoriolis(:)))


%% make figure

% plot
fig = figure('Units', 'points', 'Position', [50 50 400 400]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig = wvd.plotPoissonFlowOverContours(figureHandle=fig,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=3,jmax=2e-3,kmax=2e-3,inertialFlux=flux_advective,addFrequencyContours=true,addKEPEContours=false,yAxisLabel="deformation length");
% fig = wvd.plotPoissonFlowOverContours(figureHandle=fig,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=3,jmax=2e-3,kmax=2e-3,inertialFlux=flux_advective,addFrequencyContours=true,addKEPEContours=false,yAxisLabel="vertical mode");
ax = gca;
box on
title(runname,'Interpreter','none')

% save
if shouldExportFigures
    exportgraphics(fig,fullfile(figureFolder,sprintf("EnergyFluxExact_run%d_%d.png",runNumber,resolution)))
end


