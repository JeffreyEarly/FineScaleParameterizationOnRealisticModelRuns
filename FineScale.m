filename = "hydrostatic_mean_flow_and_wave_forced_day3000.nc";
filename = "/Users/jearly/Documents/ProjectRepositories/FineScaleParameterizationOnRealisticModelRuns/fine-scale-hydrostatic-50km-128-222.nc";
wvt = WVTransform.waveVortexTransformFromFile(filename,shouldReadOnly=true,iTime=Inf);
wvt.addOperation(SpatialForcingOperation(wvt));
int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);

%%
figure
pcolor(wvt.x,wvt.y,wvt.zeta_z(:,:,end).'/wvt.f), shading flat

%%
uz = wvt.diffZF(wvt.u);
vz = wvt.diffZF(wvt.v);

S2 = uz.*uz + vz.*vz;

Egm = 6.3e-5;
j_star = 3;
b = 1300;
k_crit = 0.6;
N0 = 3*2*pi/3600;
S2_gm = (3*pi/2)*j_star*Egm*b*k_crit*N0*N0;

wkbScale = N0*N0./wvt.N2Function(wvt.z);

epsilon_0 = 6.7e-10; % W/kg, J/s/kg = N m/s/kg = kg m^2 / s^3 / kg = m^2 / s^3

epsilon_iw = epsilon_0*shiftdim(wkbScale,-2) .* (S2/S2_gm).^2;

%


tke_diss = squeeze(epsilon_iw(:,50,:));

figure
pcolor(wvt.x/1e3,wvt.z,log10(tke_diss).'), shading flat
cb = colorbar("eastoutside");
clim([-10 -8])

%% Flux from adaptive damping (this is positive and negative)
[Fu,Fv,Feta] = wvt.spatialFluxForForcingWithName("adaptive damping");

ke_diss = wvt.u .* Fu + wvt.v .* Fv;
pe_diss = shiftdim(wvt.N2,-2) .* wvt.eta .* Feta;
% p_flux = - (wvt.u .* wvt.diffX(wvt.p) + wvt.v .* wvt.diffY(wvt.p) + wvt.w .* wvt.diffZF(wvt.p))/wvt.rho0;
te_diss = ke_diss + pe_diss; % + p_flux;
int_vol( te_diss )



%% Non-transport (non-negative) dissipation
svv = wvt.forcingWithName("adaptive damping");
Fp = sqrt(-wvt.uvMax * svv.damp) .* wvt. Ap;
Fm = sqrt(-wvt.uvMax * svv.damp) .* wvt. Am;
F0 = sqrt(-wvt.uvMax * svv.damp) .* wvt. A0;

Fu_r = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.phase.*(Fp) + wvt.UAm.*wvt.conjPhase.*(Fm),A0=wvt.UA0.*(F0));
Fv_r = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*wvt.phase.*(Fp) + wvt.VAm.*wvt.conjPhase.*(Fm),A0=wvt.VA0.*(F0));
Feta_r = wvt.transformToSpatialDomainWithG(Apm=wvt.NAp.*wvt.phase.*(Fp) + wvt.NAm.*wvt.conjPhase.*(Fm),A0=wvt.NA0.*(F0));

ke_diss = Fu_r .* Fu_r + Fv_r .* Fv_r;
pe_diss = shiftdim(wvt.N2,-2) .* Feta_r .* Feta_r;
te_diss = ke_diss + pe_diss;
int_vol( te_diss )

%%
figure
pcolor(wvt.x/1e3,wvt.z,log10(abs(squeeze(te_diss(:,50,:)))).'), shading flat
cb = colorbar("eastoutside");
clim([-10 -8])

%% Richardson number
S2 = uz.*uz + vz.*vz;
Ri = shiftdim(wvt.N2,-2) ./ S2;
Ri0 = 1.0;

a = squeeze(Ri(:,50,:));
a(a>Ri0) = nan;
figure
pcolor(wvt.x/1e3,wvt.z,a.'), shading flat
cb = colorbar("eastoutside");
clim([min(a(:)) Ri0])

%% An attempt to manually mess with the adaptive damping flux. I don't think this is helpful
F = wvt.fluxForForcing;
F_damp = F{"adaptive damping"};
Fp = F_damp.Fp;
Fm = F_damp.Fm;
F0 = F_damp.F0;

Fu_r = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.phase.*real(Fp) + wvt.UAm.*wvt.conjPhase.*real(Fm),A0=wvt.UA0.*real(F0));
Fv_r = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*wvt.phase.*real(Fp) + wvt.VAm.*wvt.conjPhase.*real(Fm),A0=wvt.VA0.*real(F0));
Feta_r = wvt.transformToSpatialDomainWithG(Apm=wvt.NAp.*wvt.phase.*real(Fp) + wvt.NAm.*wvt.conjPhase.*real(Fm),A0=wvt.NA0.*real(F0));

ke_diss = wvt.u .* Fu_r + wvt.v .* Fv_r;
pe_diss = shiftdim(wvt.N2,-2) .* wvt.eta .* Feta_r;
te_diss = ke_diss + pe_diss;
int_vol( te_diss )

%%


Fp_sym = 0.5*(Fp + conj(Fp));
Fm_sym = 0.5*(Fm + conj(Fm));
F0_sym = 0.5*(F0 + conj(F0));

Fp_skew = 0.5*(Fp - conj(Fp));
Fm_skew = 0.5*(Fm - conj(Fm));
F0_skew = 0.5*(F0 - conj(F0));

[Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0);
[Ep_sym,Em_sym,E0_sym] = wvt.energyFluxFromNonlinearFlux(Fp_sym,Fm_sym,F0_sym);
[Ep_skew,Em_skew,E0_skew] = wvt.energyFluxFromNonlinearFlux(Fp_skew,Fm_skew,F0_skew);
sum(Ep(:) + Em(:) + E0(:))
sum(Ep_sym(:) + Em_sym(:) + E0_sym(:))
sum(Ep_skew(:) + Em_skew(:) + E0_skew(:))



%%
figure
pcolor(wvt.x/1e3,wvt.z,log10(abs(squeeze(ke_diss(:,50,:)))).'), shading flat
cb = colorbar("eastoutside");
clim([-11 -6])

%% Vertical shear spectrum (full water column)
Uz_quad = squeeze(uz(1,1,:));
Vz_quad = squeeze(vz(1,1,:));

figure
tiledlayout(1,2)
nexttile
plot(Uz_quad,wvt.z), hold on
plot(Vz_quad,wvt.z)

dz = min(abs(diff(wvt.z)));
Nz = floor(wvt.Lz/dz);
z = linspace(-wvt.Lz,0,Nz).';

z_quad = wvt.z;
z_quad(1) = -wvt.Lz;
z_quad(end) = 0;
Uz = interp1(z_quad,Uz_quad,z);
Vz = interp1(z_quad,Vz_quad,z);

nexttile
plot(Uz,z), hold on
plot(Vz,z)

Cz = Uz + sqrt(-1)*Vz;

[psi,lambda]=sleptap(size(Cz,1),3);
[omega,spp,snn,spn]=mspec(dz, Cz,psi);

figure
plot(omega/(2*pi),(2*pi)*spp + (2*pi)*snn), hold on
% plot(omega/(2*pi),)
xlabel("cpm")
ylabel("s^{-2}/cpm")
xscale('log')
yscale('log')