%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Making a schematic diagram of resolved energy flux and Ozmidov scales. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basedir = '/Volumes/SanDiskExtremePro/research/JeffreyEarly/FineScaleData/';
% basedir = '/Users/jearly/Dropbox/FineScaleData/';

% filename = '7cms/fine-scale-hydrostatic-50km-128-222-cyclone-7cms-cyclogeostrophic.nc';
% shortname = '7cm/s cyclone';

% filename = '7cms/fine-scale-hydrostatic-50km-128-222-anticyclone-7cms-cyclogeostrophic.nc';
% shortname = '7cm/s anticyclone';

% filename = '25cms/fine-scale-hydrostatic-50km-256-443-anticyclone-cyclogeostrophic.nc';
% shortname = '25cm/s anticyclone';

% filename = '25cms/fine-scale-hydrostatic-50km-256-443-cyclone-cyclogeostrophic.nc';
% shortname = '25cm/s cyclone';

filename = 'fine-scale-hydrostatic-50km-256-443-one-half-dealias.nc';
shortname = 'GM IGW field';

% load the file
wvd = WVDiagnostics(fullfile(basedir,filename));
wvt = wvd.wvt;

% create figure directory
shouldExportFigures = true;
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

% timeIndices = (length(t)-4*4):length(t);
timeIndices = length(t);

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
% LOzmidov = 2*pi*sqrt( te_diss./permute(wvt.N2.^(3/2), [2 3 1]) );

% rolloff wavenumber, m_c or m_f in Kunze, JPO, 2019
mCoriolis = sqrt( permute(wvt.f*wvt.N2, [2 3 1]) ./ te_diss );
LCoriolis = 2*pi./mCoriolis;
% LCoriolis = 2*pi*sqrt( te_diss./permute(wvt.f*wvt.N2, [2 3 1]) );

% horizontal Coriolis wavenumber
kCoriolis = sqrt( wvt.f^3 ./ te_diss );
LhCoriolis = 2*pi./kCoriolis;
% LhCoriolis = 2*pi*sqrt( te_diss./wvt.f^3 );

fprintf("Maximum Ozmidov length scale %0.1f m.\n", max(LOzmidov(:)))
fprintf("Maximum vertical rolloff length scale %0.1f m.\n", max(LCoriolis(:)))
fprintf("Maximum horizontal rolloff length scale %0.1f m.\n", max(LhCoriolis(:)))

%% make figures

runname = sprintf("%s, $%d\\mathrm{km}^2$ x $%d\\mathrm{m}$, $%d^2$ x $%d$", shortname, wvt.Lx/1000, wvt.Lz, wvt.Nx, wvt.Nz);

%% histogram

fig = figure('Units', 'points', 'Position', [50 50 600 300]);
tl = tiledlayout(fig,2,3,"TileSpacing","compact","Padding","compact");

% line colors
C = orderedcolors("gem"); 
colorDictionary = dictionary("LO",{C(1,:)});
colorDictionary{"Lzf"} = C(2,:);
colorDictionary{"Lhf"} = C(3,:);

% histogram options
normalization = 'probability';
% normalization = 'cdf';

% peakOp = @(x) mean(x); peakOpName = 'mean value';
peakOp = @(x) median(x); peakOpName = 'median value';

% averaging indices
domainFrac = 1/3;
xInd = (floor(wvt.Nx/2)-floor(wvt.Nx*domainFrac/2)):(floor(wvt.Nx/2)+floor(wvt.Nx*domainFrac/2));
yInd = xInd;
zInd = 1:wvt.Nz;
zIndHist = 1:wvt.Nz;
% zIndHist = find(wvt.z>-1000);

tile = nexttile(tl,1,[2,1]);
hold on
LOzmidovProfile = squeeze( mean(mean(LOzmidov(xInd,yInd,zInd),2),1) );
LCoriolisProfile = squeeze( mean(mean(LCoriolis(xInd,yInd,zInd),2),1) );
LhCoriolisProfile = squeeze( mean(mean(LhCoriolis(xInd,yInd,zInd),2),1) );
plot(LOzmidovProfile,wvt.z,'Color',colorDictionary{'LO'},'linewidth',2,'DisplayName','L_{O}')
plot(LCoriolisProfile,wvt.z,'Color',colorDictionary{'Lzf'},'linewidth',2,'DisplayName','Lz_{f}')
plot(LhCoriolisProfile,wvt.z,'Color',colorDictionary{'Lhf'},'linewidth',2,'DisplayName','Lh_{f}')
% add grid spacing
plot(diff(wvt.z),wvt.z(1:end-1),'k--')
xline(wvt.Lx/wvt.Nx,'k:')
% fine tune plot
xlog
set(gca, 'XMinorTick', 'on')
xticks(logspace(-1, 3, 5));
xlabel('length scale (m)')
ylabel('depth (m)')
box on
% legend('Location','best')

LOzmidovHist = LOzmidov(:,:,zIndHist);
LCoriolisHist = LCoriolis(:,:,zIndHist);
LhCoriolisHist = LhCoriolis(:,:,zIndHist);

tile = nexttile(tl,2,[2,2]);
hold on
histLO = histogram(LOzmidovHist(:),'FaceColor',colorDictionary{'LO'},'EdgeColor','none','Normalization',normalization,'DisplayName','Ozmidov length');
histLzf = histogram(LCoriolisHist(:),'FaceColor',colorDictionary{'Lzf'},'EdgeColor','none','Normalization',normalization,'DisplayName','Vertical Coriolis length');
histLhf = histogram(LhCoriolisHist(:),'FaceColor',colorDictionary{'Lhf'},'EdgeColor','none','Normalization',normalization,'DisplayName','Horizontal Coriolis length');
xlog
xlim([0,1000])
xlabel('length scale (m)')
ylabel(normalization)
% yticklabels([])
box on
legend('Location','best')
% label peaks
% LO
% [maxCount, maxIdx] = max(histLO.Values);
% plot(histLO.BinEdges(maxIdx),maxCount,'k.','HandleVisibility','off')
% text(histLO.BinEdges(maxIdx),maxCount+3e-4,sprintf('%0.1f m',histLO.BinEdges(maxIdx)))
[~,idx] = min(abs(histLO.BinEdges-peakOp(LOzmidovHist(:))));
plot(peakOp(LOzmidovHist(:)),histLO.Values(idx),'k.','DisplayName',peakOpName)
text(peakOp(LOzmidovHist(:)),histLO.Values(idx),sprintf('%0.0f m',peakOp(LOzmidovHist(:))),'HorizontalAlignment','left','VerticalAlignment','bottom')
% Lzf
% [maxCount, maxIdx] = max(histLzf.Values);
% plot(histLzf.BinEdges(maxIdx),maxCount,'k.','HandleVisibility','off')
% text(histLzf.BinEdges(maxIdx),maxCount+3e-4,sprintf('%0.1f m',histLzf.BinEdges(maxIdx)))
[~,idx] = min(abs(histLzf.BinEdges-peakOp(LCoriolisHist(:))));
plot(peakOp(LCoriolisHist(:)),histLzf.Values(idx),'k.','HandleVisibility','off')
text(peakOp(LCoriolisHist(:)),histLzf.Values(idx),sprintf('%0.0f m',peakOp(LCoriolisHist(:))),'HorizontalAlignment','left','VerticalAlignment','bottom')
% Lhf
% [maxCount, maxIdx] = max(histLhf.Values);
% plot(histLhf.BinEdges(maxIdx),maxCount,'k.','HandleVisibility','off')
% text(histLhf.BinEdges(maxIdx),maxCount+3e-4,sprintf('%0.1f m',histLhf.BinEdges(maxIdx)))
[~,idx] = min(abs(histLhf.BinEdges-peakOp(LhCoriolisHist(:))));
plot(peakOp(LhCoriolisHist(:)),histLhf.Values(idx),'k.','HandleVisibility','off')
text(peakOp(LhCoriolisHist(:)),histLhf.Values(idx),sprintf('%0.0f m',peakOp(LhCoriolisHist(:))),'HorizontalAlignment','left','VerticalAlignment','bottom')

% model grid spacing
xline(min(diff(wvt.z)),'k--','DisplayName','\Deltaz')
xline(wvt.Lx/wvt.Nx,'k:','DisplayName','\Deltax')

% title
title(tl,runname,'Interpreter','latex')

% save
if shouldExportFigures
    exportgraphics(fig,fullfile(figureFolder,sprintf('TurbulenceScales_%s.png',replace(replace(shortname,'cm/s','cmps'),' ','-'))))
end


%% energy flux

% plot
fig = figure('Units', 'points', 'Position', [50 50 400 400]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
% fig = wvd.plotPoissonFlowOverContours(figureHandle=fig,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=3,jmax=2e-3,kmax=2e-3,inertialFlux=flux_advective,addFrequencyContours=true,addKEPEContours=false,yAxisLabel="deformation length");
fig = wvd.plotPoissonFlowOverContours(figureHandle=fig,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=3,jmax=2e-3,kmax=2e-3,inertialFlux=flux_advective,addFrequencyContours=true,addKEPEContours=false,yAxisLabel="vertical mode");
ax = gca;
box on
title(runname,'Interpreter','latex')

% save
if shouldExportFigures
    exportgraphics(fig,fullfile(figureFolder,sprintf('EnergyFluxExact_%s.png',replace(replace(shortname,'cm/s','cmps'),' ','-'))))
end


