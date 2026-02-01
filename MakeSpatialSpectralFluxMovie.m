%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Open the model output, create a new WVT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basedir = "/Users/jearly/Dropbox/FineScaleData/";
basedir = "";
% wvd = WVDiagnostics(basedir + "fine-scale-hydrostatic-50km-128-222-anticyclone-cyclogeostrophic.nc");
% wvd = WVDiagnostics("/Users/jearly/Dropbox/FineScaleData/fine-scale-hydrostatic-50km-128-222-cyclone-cyclogeostrophic.nc");
% wvd = WVDiagnostics("/Users/jearly/Dropbox/FineScaleData/fine-scale-hydrostatic-50km-128-222-anticyclone-cyclogeostrophic.nc");
wvd = WVDiagnostics("/Volumes/seattle_data1/jearly/FineScaleData/fine-scale-boussinesq-50km-128-222-no-eddy.nc");
wvt = wvd.wvt;

shouldShowEddyContours = false;
fluxArrowScale = 1.0; % I used 2 for the eddy+IO cases, and 0.5 for the IGW case.
fluxSmoothing = 5;

figureFolder = "./figures-spectral-movie";
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder)
end

t = wvd.t_diag;

xIndices = floor(wvt.Nx/2);
yIndices = 1:wvt.Ny;
zIndices = 1:wvt.Nz;
horzAxis = wvt.y/1e3;
vertAxis = wvt.z;

int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Define the phases
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure('Units', 'points', 'Position', [50 50 600 600]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

% tl = tiledlayout(2,2,TileSpacing="tight");
% tl = tiledlayout(2,2,TileSpacing="tight",Padding="tight");

wvd.iTime = 1;
t0 = wvt.t;

for iTime= 1:length(t)
    analysisIndices = (iTime-fluxSmoothing):(iTime + fluxSmoothing); %indices_phase{iPhase};
    analysisIndices(analysisIndices<1) = [];
    analysisIndices(analysisIndices>length(t)) = [];
    wvd.iTime = iTime;

    clf
    tl = tiledlayout(2,2,TileSpacing="tight",Padding="tight");

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Wave energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axTile = nexttile(tl,1);
    wave_energy = wvt.u_w.^2 + wvt.v_w.^2 + shiftdim(wvt.N2,-2) .* wvt.eta_w.^2;
    rv = wvt.diffX(wvt.v_g) - wvt.diffY(wvt.u_g);

    p3 = pcolor(horzAxis,vertAxis,squeeze(1e4*wave_energy(xIndices,yIndices,zIndices)).'); shading interp, hold on
    
    if shouldShowEddyContours
        contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[-0.01 -0.01],'k',linewidth=1.5)
        contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[0.01 0.01],'k',linewidth=1.5)
    end
    ylim([-750 0])
    clim([0 100])
    colormap(axTile,WVDiagnostics.cmocean('tempo'))
    title(axTile,"wave energy (" + wvt.waveEnergy/wvd.escale + " " + wvd.escale_units +")")
    cb_we = colorbar('south');
    cb_we.Label.String = "cm^2 s^{-2}";
    cb_we.Label.Position = [80 2.2750 0]; % positive above 80

    axTile.XTick = [];

    svv = wvt.forcingWithName("adaptive damping");
    Fp = sqrt(-wvt.uvMax * svv.damp) .* wvt. Ap;
    Fm = sqrt(-wvt.uvMax * svv.damp) .* wvt. Am;
    F0 = sqrt(-wvt.uvMax * svv.damp) .* wvt. A0;

    Fu = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.phase.*(Fp) + wvt.UAm.*wvt.conjPhase.*(Fm),A0=wvt.UA0.*(F0));
    Fv = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*wvt.phase.*(Fp) + wvt.VAm.*wvt.conjPhase.*(Fm),A0=wvt.VA0.*(F0));
    Feta = wvt.transformToSpatialDomainWithG(Apm=wvt.NAp.*wvt.phase.*(Fp) + wvt.NAm.*wvt.conjPhase.*(Fm),A0=wvt.NA0.*(F0));

    te_diss = Fu .* Fu + Fv .* Fv + shiftdim(wvt.N2,-2) .* Feta .* Feta;
    te_diss_tot = int_vol( te_diss );


    axTile = nexttile(tl,2);
    pcolor(wvt.x/1e3,wvt.z,log10(abs(squeeze(te_diss(xIndices,yIndices,zIndices)))).'), shading flat, hold on
    if shouldShowEddyContours
        contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[-0.01 -0.01],'k',linewidth=1.5)
        contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[0.01 0.01],'k',linewidth=1.5)
    end
    cb_diss = colorbar("south");
    cb_diss.Label.String = "log_{10} [m^2 s^{-3}]";
    cb_diss.Label.Position = [-8.5 2.2750 0]; % positive above -8.5
    % ylim([-2000 0])
    colormap(axTile,WVDiagnostics.cmocean('amp'))
    ylim([-750 0])
    clim([-10 -8])
    title("energy dissipation (" + te_diss_tot/wvd.flux_scale + " " + wvd.flux_scale_units +")")

     axTile.XTick = [];
     axTile.YTick = [];

    tExp = wvt.t - t0;
    days = floor(tExp/86400) ;
    hours = floor((tExp-86400*floor(tExp/86400))/3600);
    title(tl, "IGW field " + days + " days " + hours + " hours")

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fluxes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C = orderedcolors("gem");
    colorDictionary = dictionary("adaptive_damping",{C(2,:)});
    % colorDictionary{"M2_tidal_forcing"} = C(5,:);

    energy_fluxes = wvd.quadraticEnergyFluxesTemporalAverage(timeIndices=analysisIndices);
    [inertial_fluxes_g, inertial_fluxes_w, ks, js] = wvd.quadraticEnergyPrimaryTriadFluxesTemporalAverage2D(timeIndices=analysisIndices,outputGrid="full");

    clear ggg ggw www wwg
    ggg.flux = inertial_fluxes_g([inertial_fluxes_g.name] == "ggg").flux/wvd.flux_scale;
    ggw.flux = inertial_fluxes_g([inertial_fluxes_g.name] == "ggw").flux/wvd.flux_scale;
    www.flux = inertial_fluxes_w([inertial_fluxes_w.name] == "www").flux/wvd.flux_scale;
    wwg.flux = inertial_fluxes_w([inertial_fluxes_w.name] == "wwg").flux/wvd.flux_scale;

    clear forcing_fluxes;
    reservoirName = "te_gmda";
    fluxesOfInterest = {"adaptive_damping"};
    for i=1:length(fluxesOfInterest)
        forcing_fluxes(i).color=colorDictionary{fluxesOfInterest{i}};
        forcing_fluxes(i).flux = energy_fluxes([energy_fluxes.name] == fluxesOfInterest{i}).(reservoirName)/wvd.flux_scale;
        forcing_fluxes(i).relativeAmplitude = 1.0;
        forcing_fluxes(i).alpha = 1.0;
        forcing_fluxes(i).fancyName = energy_fluxes([energy_fluxes.name] == fluxesOfInterest{i}).fancyName;
    end

    maxAmplitude = max(arrayfun( @(v) max(abs(v.flux(:))), forcing_fluxes));
    for i=1:length(forcing_fluxes)
        forcing_fluxes(i).relativeAmplitude = max(abs((forcing_fluxes(i).flux(:))))/maxAmplitude;
    end
    
    axTile = nexttile(tl,3);
    wvd.plotPoissonFlowOverContours(figureHandle=axTile,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=fluxArrowScale,inertialFlux=wwg,addFrequencyContours=true);
    title("wwg")

    axTile = nexttile(tl,4);
    wvd.plotPoissonFlowOverContours(figureHandle=axTile,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=fluxArrowScale,inertialFlux=www,addFrequencyContours=true);
    title("www")
    axTile.YTick = [];
    axTile.YLabel.String = '';

    pause(0.5);
    exportgraphics(fig,figureFolder + "/" + "energy_flux_" + iTime + ".png",Resolution=300)
end