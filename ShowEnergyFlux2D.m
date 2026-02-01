%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Open the model output, create a new WVT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basedir = "/Users/jearly/Dropbox/FineScaleData/";
basedir = "";
% wvd = WVDiagnostics(basedir + "fine-scale-hydrostatic-50km-128-222-anticyclone-cyclogeostrophic.nc");
wvd = WVDiagnostics("/Users/jearly/Dropbox/FineScaleData/fine-scale-hydrostatic-50km-128-222-cyclone-cyclogeostrophic.nc");
wvt = wvd.wvt;

figureFolder = "./figures";
if ~exist(figureFolder, 'dir')
    mkdir(figureFolder)
end

t = wvd.t_diag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Define the phases
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indices_phase{1} = 1:4;
indices_phase{2} = 5:8;
indices_phase{3} = 4*(0:3) + 1;

fig_ggw = figure('Units', 'points', 'Position', [50 50 400 400]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

tl = tiledlayout(1,2);

for iPhase=3:3
    analysisIndices = 5; %indices_phase{iPhase};

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
    forcing_fluxes(i+1).flux = (inertial_fluxes_g([inertial_fluxes_g.name] == "tx-ggw").flux + inertial_fluxes_g([inertial_fluxes_g.name] == "tx-wwg").flux)/wvd.flux_scale;
    forcing_fluxes(i+1).color = 0.5*[0.5 0.5 1];
    forcing_fluxes(i+1).relativeAmplitude = 0.5;
    forcing_fluxes(i+1).alpha = 0.25;
    forcing_fluxes(i+1).fancyName = "wwg+ggw";

    % forcing_fluxes(i+1).flux = inertial_fluxes_g([inertial_fluxes_g.name] == "tx-ggw").flux/wvd.flux_scale;
    % forcing_fluxes(i+1).color = 0.5*[0.5 0.5 1];
    % forcing_fluxes(i+1).relativeAmplitude = 0.5;
    % forcing_fluxes(i+1).alpha = 0.25;
    % forcing_fluxes(i+1).fancyName = "ggw";
    %
    % i = i+1;
    % forcing_fluxes(i+1).flux = inertial_fluxes_g([inertial_fluxes_g.name] == "tx-wwg").flux/wvd.flux_scale;
    % forcing_fluxes(i+1).color = 0.5*[1 0.5 0.5];
    % forcing_fluxes(i+1).relativeAmplitude = 0.5;
    % forcing_fluxes(i+1).alpha = 0.25;
    % forcing_fluxes(i+1).fancyName = "wwg";

    maxAmplitude = max(arrayfun( @(v) max(abs(v.flux(:))), forcing_fluxes));
    for i=1:length(forcing_fluxes)
        forcing_fluxes(i).relativeAmplitude = max(abs((forcing_fluxes(i).flux(:))))/maxAmplitude;
    end

    % fig_ggg = figure('Units', 'points', 'Position', [50 50 400 400]);
    % set(gcf,'PaperPositionMode','auto')
    % set(gcf, 'Color', 'w');
    % 
    % fig_ggw = figure('Units', 'points', 'Position', [50 50 400 400]);
    % set(gcf,'PaperPositionMode','auto')
    % set(gcf, 'Color', 'w');
    
    tile = nexttile(tl,1);
    wvd.plotPoissonFlowOverContours(figureHandle=tile,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=2,inertialFlux=wwg,addFrequencyContours=true);
    title("wwg")
    % exportgraphics(fig_ggg,figureFolder + "/" + "energy_flux_2D_flow_geostrophic_ggg_phase" + iPhase + ".png",Resolution=300)

    tile = nexttile(tl,2);
    wvd.plotPoissonFlowOverContours(figureHandle=tile,vectorDensityLinearTransitionWavenumber=10^(-3.9),quiverScale=2,inertialFlux=www,addFrequencyContours=true);
    title("www")
    % exportgraphics(fig_ggw,figureFolder + "/" + "energy_flux_2D_flow_geostrophic_ggw_phase" + iPhase + ".png",Resolution=300)
end