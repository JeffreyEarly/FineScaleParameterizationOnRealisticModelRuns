%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Open the model output, create a new WVT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basedir = "/Users/jearly/Documents/ProjectRepositories/FineScaleParameterizationOnRealisticModelRuns/";
wvdArray = cell(2,1);
%
wvdArray{1} = WVDiagnostics(basedir + "fine-scale-hydrostatic-50km-128-222-anticyclone-7cms-cyclogeostrophic.nc");
wvdArray{2} = WVDiagnostics(basedir + "fine-scale-hydrostatic-50km-128-222-cyclone-7cms-cyclogeostrophic.nc");


shouldShowEddyContours = true;
fluxArrowScale = 2.5; % I used 2 for the eddy+IO cases, and 0.5 for the IGW case.
fluxSmoothing = 5;
simTitle = "Anticyclone"; % "IGW field";

figureFolder = "./figures-proposal";
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder)
end

t = wvd.t_diag;

wvt = wvdArray{1}.wvt;
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

fig = figure('Units', 'points', 'Position', [50 50 1200 600]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

% tl = tiledlayout(2,2,TileSpacing="tight");
% tl = tiledlayout(2,2,TileSpacing="tight",Padding="tight");

wvd.iTime = 1;
t0 = wvt.t;

iTime = 50;

analysisIndices = (iTime-fluxSmoothing):(iTime + fluxSmoothing); %indices_phase{iPhase};
analysisIndices(analysisIndices<1) = [];
analysisIndices(analysisIndices>length(t)) = [];
shouldShowWaveEnergy = false;

clf
if shouldShowWaveEnergy
    tl = tiledlayout(2,4,TileSpacing="tight",Padding="tight");
    fudge = 0;
    nCols = 4;
else
    % tl = tiledlayout(2,3,TileSpacing="tight",Padding="tight");
    tl = tiledlayout(2,3,TileSpacing="compact",Padding="compact");
    fudge = -1;
    nCols = 3;
end
% tl = tiledlayout(2,4,TileSpacing="compact",Padding="compact");
% tl = tiledlayout(2,4);



for iWVD=1:2

    wvd = wvdArray{iWVD};
    wvd.iTime = iTime;
    wvt = wvd.wvt;

rv = wvt.diffX(wvt.v_g) - wvt.diffY(wvt.u_g);
wave_energy = wvt.u_w.^2 + wvt.v_w.^2 + shiftdim(wvt.N2,-2) .* wvt.eta_w.^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Wave energy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if shouldShowWaveEnergy
        axTile = nexttile(tl,(iWVD-1)*nCols + 1 + fudge);
        
        

        p3 = pcolor(horzAxis,vertAxis,squeeze(1e4*wave_energy(xIndices,yIndices,zIndices)).'); shading interp, hold on

        if shouldShowEddyContours
            contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[-0.01 -0.01],'k',linewidth=1.5)
            contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[0.01 0.01],'k',linewidth=1.5)
        end
        ylim([-750 0])
        clim([0 100])
        colormap(axTile,WVDiagnostics.cmocean('tempo'))
        % title(axTile,"wave energy (" + wvt.waveEnergy/wvd.escale + " " + wvd.escale_units +")")
        if iWVD == 1
            title(axTile,"wave energy")
            axTile.XTick = [];
        elseif iWVD == 2
            cb_we = colorbar('south');
            cb_we.Label.String = "cm^2 s^{-2}";
            cb_we.Label.Position = [80 2.2750 0]; % positive above 80
            xlabel("horizontal distance (km)")
        end
        if iWVD == 1
            axTile.YLabel.String = 'Anticyclone';
        elseif iWVD == 2
            axTile.YLabel.String = 'Cyclone';
        end



    end

    svv = wvt.forcingWithName("adaptive damping");
    Fp = sqrt(-wvt.uvMax * svv.damp) .* wvt. Ap;
    Fm = sqrt(-wvt.uvMax * svv.damp) .* wvt. Am;
    F0 = sqrt(-wvt.uvMax * svv.damp) .* wvt. A0;

    Fu = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.phase.*(Fp) + wvt.UAm.*wvt.conjPhase.*(Fm),A0=wvt.UA0.*(F0));
    Fv = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*wvt.phase.*(Fp) + wvt.VAm.*wvt.conjPhase.*(Fm),A0=wvt.VA0.*(F0));
    Feta = wvt.transformToSpatialDomainWithG(Apm=wvt.NAp.*wvt.phase.*(Fp) + wvt.NAm.*wvt.conjPhase.*(Fm),A0=wvt.NA0.*(F0));

    te_diss = Fu .* Fu + Fv .* Fv + shiftdim(wvt.N2,-2) .* Feta .* Feta;
    te_diss_tot = int_vol( te_diss );

    axTile = nexttile(tl,(iWVD-1)*nCols + 2+ fudge);
    pcolor(wvt.x/1e3,wvt.z,log10(abs(squeeze(te_diss(xIndices,yIndices,zIndices)))).'), shading flat, hold on
    % contour(wvt.x/1e3,wvt.z,squeeze(wave_energy(xIndices,yIndices,zIndices)).')
    if shouldShowEddyContours
        contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,1*[-0.01 -0.01],'k',linewidth=1.5)
        contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,1*[0.01 0.01],'k',linewidth=1.5)
    end
    if iWVD == 1
        title("energy dissipation")
        axTile.XTick = [];
    elseif iWVD == 2
        cb_diss = colorbar("south");
        cb_diss.Label.String = "log_{10} [m^2 s^{-3}]";
        cb_diss.Label.Position = [-8.5 2.2750 0]; % positive above -8.5
        xlabel("horizontal distance (km)")
    end
    % ylim([-2000 0])
    colormap(axTile,WVDiagnostics.cmocean('amp'))
    ylim([-750 0])
    clim([-10 -8])
    % title("energy dissipation (" + te_diss_tot/wvd.flux_scale + " " + wvd.flux_scale_units +")")

    axTile.YTick = [];

    tExp = wvt.t - t0;
    days = floor(tExp/86400) ;
    hours = floor((tExp-86400*floor(tExp/86400))/3600);
    % title(tl, simTitle + " " + days + " days " + hours + " hours")

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

    axTile = nexttile(tl,(iWVD-1)*nCols + 3+ fudge);
    wvd.plotPoissonFlowOverContours(figureHandle=axTile,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=fluxArrowScale,inertialFlux=www,addFrequencyContours=true,yAxisLabel = "vertical mode");
    % wvd.plotPoissonFlowOverContours(figureHandle=axTile,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=fluxArrowScale,inertialFlux=www,addFrequencyContours=true);
    if iWVD == 1
        title("www flux")
    end
    axTile.YTick = [];
    axTile.YLabel.String = '';
    % axTile.Legend.Visible = "off";
    if iWVD == 1
        axTile.XTick = [];
        axTile.XLabel.String = '';
    end
    pause(0.5);
    axTile = nexttile(tl,(iWVD-1)*nCols + 4+ fudge);
    wvd.plotPoissonFlowOverContours(figureHandle=axTile,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=fluxArrowScale,inertialFlux=wwg,addFrequencyContours=true,yAxisLabel = "vertical mode");
    % wvd.plotPoissonFlowOverContours(figureHandle=axTile,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=fluxArrowScale,inertialFlux=wwg,addFrequencyContours=true);
    if iWVD == 1
        title("wwg flux")
    end
    axTile.YAxisLocation = 'right';
    % axTile.YTick = [];
    % axTile.YLabel.String = '';
    % axTile.Legend.Visible = "off";
    if iWVD == 1
        axTile.XTick = [];
        axTile.XLabel.String = '';
    end
    pause(0.5);
end

% pause(0.5);
% exportgraphics(fig,figureFolder + "/" + "energy_flux_" + iTime + ".png",Resolution=300)