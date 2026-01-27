% filename = "fine-scale-restart-with-eddy-no-igw-hydrostatic-200km-128-56.nc";
% filename = "fine-scale-restart-with-eddy-no-igw-hydrostatic-50km-64-111-one-half-dealias-take-2.nc";
% filename = "fine-scale-restart-with-eddy-no-igw-hydrostatic-200km-128-56-one-half-dealias.nc";
filename = "fine-scale-restart-with-eddy-no-igw-hydrostatic-50km-128-222-one-half-dealias-take-2.nc";
% filename = "fine-scale-restart-with-eddy-no-igw-hydrostatic-50km-64-111-one-half-dealias-take-2-cyclone.nc";
wvd = WVDiagnostics(filename);
wvt = wvd.wvt;
int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);

%%
wvd.plotEnergyOverTime(energyReservoirs = [EnergyReservoir.geostrophic_mda, EnergyReservoir.io, EnergyReservoir.igw, EnergyReservoir.total]);
title("medium-res cyclone")
exportgraphics(gcf,"energy-medium-res-cyclone.pdf")

%%
wvd.plotEnergyFluxOverTime();
title("medium-res cyclone")
ylim([-1.0 0.5])
exportgraphics(gcf,"energy-flux-medium-res-cyclone.pdf")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% y-z slice of (u,v,eta) through the entire domain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvd.iTime = 61;

xIndices = floor(wvt.Nx/2);
yIndices = 1:wvt.Ny;
zIndices = 1:wvt.Nz;
horzAxis = wvt.y/1e3;
vertAxis = wvt.z;

fig1 = figure('Units', 'points', 'Position', [50 50 750 800]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

tl = tiledlayout(2,1,TileSpacing="tight");

tile = nexttile;
eta = wvt.u_w;
eta = wvt.u_w.^2 + wvt.v_w.^2 + shiftdim(wvt.N2,-2) .* wvt.eta_w.^2;
rv = wvt.diffX(wvt.v_g) - wvt.diffY(wvt.u_g);

p3 = pcolor(horzAxis,vertAxis,squeeze(eta(xIndices,yIndices,zIndices)).'); shading interp, hold on

contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[-0.01 -0.01],'k',linewidth=1.5)
contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[0.01 0.01],'k',linewidth=1.5)
title(tile,'u of the wave field')
cb3 = colorbar('eastoutside');
% cb3.Label.String = 'cm';

% Non-transport (non-negative) dissipation
svv = wvt.forcingWithName("adaptive damping");
Fp = sqrt(-wvt.uvMax * svv.damp) .* wvt. Ap;
Fm = sqrt(-wvt.uvMax * svv.damp) .* wvt. Am;
F0 = sqrt(-wvt.uvMax * svv.damp) .* wvt. A0;

Fu = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.phase.*(Fp) + wvt.UAm.*wvt.conjPhase.*(Fm),A0=wvt.UA0.*(F0));
Fv = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*wvt.phase.*(Fp) + wvt.VAm.*wvt.conjPhase.*(Fm),A0=wvt.VA0.*(F0));
Feta = wvt.transformToSpatialDomainWithG(Apm=wvt.NAp.*wvt.phase.*(Fp) + wvt.NAm.*wvt.conjPhase.*(Fm),A0=wvt.NA0.*(F0));

te_diss = Fu .* Fu + Fv .* Fv + shiftdim(wvt.N2,-2) .* Feta .* Feta;
int_vol( te_diss )/wvd.flux_scale


tile = nexttile;
pcolor(wvt.x/1e3,wvt.z,log10(abs(squeeze(te_diss(xIndices,yIndices,zIndices)))).'), shading flat, hold on
contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[-0.01 -0.01],'k',linewidth=1.5)
contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[0.01 0.01],'k',linewidth=1.5)

cb = colorbar("eastoutside");
% ylim([-2000 0])
clim([-10 -8])
title("energy dissipation from adaptive damping")

title(tl,"day " + floor(wvt.t/86400))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Mode spectrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uz = wvt.diffZF(wvt.u);
vz = wvt.diffZF(wvt.v);
S2 = uz.*uz + vz.*vz;
Ri = shiftdim(wvt.N2,-2) ./ S2;
Ri(isinf(Ri)) = nan;
Ri(Ri>2) = nan;

figure
p3 = pcolor(horzAxis,vertAxis,squeeze(Ri(xIndices,yIndices,zIndices)).'); shading flat, hold on
clim([0 2])
colorbar("eastoutside");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% y-z slice of (u,v,eta) through the entire domain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt = 6;
dT = 10;



xIndices = floor(wvt.Nx/2);
yIndices = 1:wvt.Ny;
zIndices = 1:wvt.Nz;
horzAxis = wvt.y/1e3;
vertAxis = wvt.z;

fig1 = figure('Units', 'points', 'Position', [50 50 2000 800]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

tl = tiledlayout(2,Nt,TileSpacing="tight");

for iT=0:(Nt-1)
    wvd.iTime = 4*(iT*dT)+1;

    tile = nexttile(iT+1);
    eta = wvt.u_w;
    eta = wvt.u_w.^2 + wvt.v_w.^2 + shiftdim(wvt.N2,-2) .* wvt.eta_w.^2;
    rv = wvt.diffX(wvt.v_g) - wvt.diffY(wvt.u_g);

    p3 = pcolor(horzAxis,vertAxis,squeeze(eta(xIndices,yIndices,zIndices)).'); shading interp, hold on

    contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[-0.01 -0.01],'k',linewidth=1.5)
    contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[0.01 0.01],'k',linewidth=1.5)
    title(tile,"day " + floor(wvt.t/86400))
    % title(tile,'u of the wave field')
    % cb3 = colorbar('eastoutside');
    % cb3.Label.String = 'cm';

    % Non-transport (non-negative) dissipation
    svv = wvt.forcingWithName("adaptive damping");
    Fp = sqrt(-wvt.uvMax * svv.damp) .* wvt. Ap;
    Fm = sqrt(-wvt.uvMax * svv.damp) .* wvt. Am;
    F0 = sqrt(-wvt.uvMax * svv.damp) .* wvt. A0;

    Fu = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.phase.*(Fp) + wvt.UAm.*wvt.conjPhase.*(Fm),A0=wvt.UA0.*(F0));
    Fv = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*wvt.phase.*(Fp) + wvt.VAm.*wvt.conjPhase.*(Fm),A0=wvt.VA0.*(F0));
    Feta = wvt.transformToSpatialDomainWithG(Apm=wvt.NAp.*wvt.phase.*(Fp) + wvt.NAm.*wvt.conjPhase.*(Fm),A0=wvt.NA0.*(F0));

    te_diss = Fu .* Fu + Fv .* Fv + shiftdim(wvt.N2,-2) .* Feta .* Feta;
    int_vol( te_diss )/wvd.flux_scale


    tile = nexttile(iT + 1 + Nt);
    pcolor(wvt.x/1e3,wvt.z,log10(abs(squeeze(te_diss(xIndices,yIndices,zIndices)))).'), shading flat, hold on
    contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[-0.01 -0.01],'k',linewidth=1.5)
    contour(horzAxis,vertAxis,squeeze(rv(xIndices,yIndices,zIndices)).'/wvt.f,[0.01 0.01],'k',linewidth=1.5)

    % cb = colorbar("eastoutside");
    % % ylim([-2000 0])
    clim([-10 -8])
    % title("energy dissipation from adaptive damping")

    % title(tl,"day " + floor(wvt.t/86400))
end

tile = nexttile(iT+1);
cb3 = colorbar('eastoutside');
cb3.Label.String = 'wave energy (m^2/s^2)';

tile = nexttile(iT + 1 + Nt);
cb = colorbar("eastoutside");
cb.Label.String = 'energy dissipation log10(m^3/s^2)';