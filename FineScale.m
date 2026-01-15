filename = "hydrostatic_mean_flow_and_wave_forced_day3000.nc";
wvt = WVTransform.waveVortexTransformFromFile(filename,shouldReadOnly=true);

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

epsilon_0 = 6.7e-10; % W/kg

epsilon_iw = epsilon_0*shiftdim(wkbScale,-2) .* (S2/S2_gm).^2;

%

tke_diss = squeeze(epsilon_iw(:,50,:));

figure
pcolor(wvt.x/1e3,wvt.z,log10(tke_diss).'), shading flat
cb = colorbar("eastoutside");
clim([-11 -6])

%%
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