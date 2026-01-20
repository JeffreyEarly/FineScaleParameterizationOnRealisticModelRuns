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

i = 1;
resolution(i).Nxy = 512;
resolution(i).Nz = 700;
resolution(i).L = 500e3;
i = i+1;
resolution(i).Nxy = 512;
resolution(i).Nz = 512;
resolution(i).L = 500e3;
i = i+1;
resolution(i).Nxy = 512;
resolution(i).Nz = 700;
resolution(i).L = 50e3;
i = i+1;
resolution(i).Nxy = 512;
resolution(i).Nz = 512;
resolution(i).L = 50e3;
i = i+1;
resolution(i).Nxy = 512;
resolution(i).Nz = 700;
resolution(i).L = 10e3;
i = i+1;
resolution(i).Nxy = 512;
resolution(i).Nz = 512;
resolution(i).L = 10e3;

transform = 'hydrostatic';

resolution(2:end) = [];

wvtArray = cell(length(resolution),1);
for iRes = 1:length(resolution)
    disp("Loading transforming...");

    Nxy = resolution(iRes).Nxy;
    Nz = resolution(iRes).Nz;
    L = resolution(iRes).L;

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
    wvtArray{iRes} = wvt;
end

%%
for iRes = 1:length(wvtArray)
    disp("Initializing spectrum...");
    wvt = wvtArray{iRes};

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
end

%%

figure

for iRes = 1:length(wvtArray)
wvt = wvtArray{iRes};

uz = wvt.diffZF(wvt.u);
vz = wvt.diffZF(wvt.v);
% 
% Uz_quad = squeeze(uz(1,1,:));
% Vz_quad = squeeze(vz(1,1,:));

Uz_quad = shiftdim(uz(1:10:end,1:10:end,:),2);
Vz_quad = shiftdim(vz(1:10:end,1:10:end,:),2);
Uz_quad = reshape(Uz_quad,wvt.Nz,[]);
Vz_quad = reshape(Vz_quad,wvt.Nz,[]);

dz = min(abs(diff(wvt.z)));
Nz = floor(wvt.Lz/dz);
z = linspace(-wvt.Lz,0,Nz).';

z_quad = wvt.z;
z_quad(1) = -wvt.Lz;
z_quad(end) = 0;
Uz = interp1(z_quad,Uz_quad,z);
Vz = interp1(z_quad,Vz_quad,z);

Cz = Uz + sqrt(-1)*Vz;

[psi,lambda]=sleptap(size(Cz,1),3);
[omega,spp,snn,spn]=mspec(dz, Cz,psi);

plot(omega/(2*pi),mean((2*pi)*spp + (2*pi)*snn,2)), hold on

end

xlabel("cpm")
ylabel("s^{-2}/cpm")
xscale('log')
yscale('log')

label = cell(length(resolution),1);
for iRes = 1:length(resolution)
    label{iRes} = "L=" + string(resolution(iRes).L/1e3) + "km, N=" + resolution(iRes).Nxy + ", Nz=" + resolution(iRes).Nz;
end
legend(label);

%exportgraphics(gca,"shear_spectrum.eps")