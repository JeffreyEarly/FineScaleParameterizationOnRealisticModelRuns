%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Setup the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% At 50 km horizontal domain size, 128 point resolution, isotropic scaling
% suggests 121 vertical grid points.
N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);
Lz = 4000;
L = 500e3;
Nxy = 512;
% Nz = WVStratification.verticalResolutionForHorizontalResolution(L,Lz,Nxy,N2=N2,latitude=27);
Nz = 700;

transform = 'hydrostatic';
% transform = 'boussinesq';

emptyFilename = "empty-"+transform+"-" + string(round(L)/1e3) + "km-" + string(Nxy) + "-" + string(Nz) + ".nc";
filename = "fine-scale-"+transform+"-" + string(round(L)/1e3) + "km-" + string(Nxy) + "-" + string(Nz) + ".nc";
if exist(emptyFilename,'file')
    wvt = WVTransform.waveVortexTransformFromFile(emptyFilename);
else
    if strcmp(transform,'boussinesq')
        wvt = WVTransformBoussinesq([L, L, Lz], [Nxy, Nxy, Nz], N2=N2,latitude=31);
    elseif strcmp(transform,'hydrostatic')
        wvt = WVTransformHydrostatic([L, L, Lz], [Nxy, Nxy, Nz], N2=N2,latitude=31);
    else
        error('unknown transform');
    end
    wvt.addForcing(WVAdaptiveDamping(wvt));
    wvt.writeToFile(emptyFilename,shouldOverwriteExisting=true);
end

%%

% wvt.initWithGMSpectrum();

wvt.removeAll;

GMAmplitude = 1.0;
j_star = 3;
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = GMAmplitude*L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;

Nmax = sqrt(max(wvt.N2));

j_slope = 1.0;

% Compute the proper vertical function normalization
H = (j_star.^2+(1:1024).^2).^(-j_slope);
H_norm = 1/sum(H);

% Do the same for the frequency function.
B_norm = 1/atan(sqrt(Nmax*Nmax/(wvt.f*wvt.f)-1));

H = @(j) H_norm*((j.^2+j_star.^2).^(-j_slope));
B = @(omega) B_norm*wvt.f./(omega.*sqrt(omega.*omega-wvt.f*wvt.f));
GM = @(omega,j) E*H(j) .* B(omega);

wvt.addWavesWithFrequencySpectrum(ApmSpectrum=GM,shouldOnlyRandomizeOrientations=true, shouldThrowErrorIfDensityViolation=false);


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

% model = WVModel(wvt);
% 
% model.createNetCDFFileForModelOutput("igw-simulation-day0-4.nc",outputInterval=wvt.inertialPeriod/4,shouldOverwriteExisting=true);
% 
% model.integrateToTime(wvt.t + 4*wvt.inertialPeriod);