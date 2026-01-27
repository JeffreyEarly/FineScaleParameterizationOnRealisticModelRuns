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

L = 50e3;
Nxy = 256;
latitude = 31;

transform = 'hydrostatic';
% transform = 'boussinesq';

Nz = WVStratification.verticalResolutionForHorizontalResolution(L,Lz,Nxy,N2=N2,latitude=latitude);

emptyFilename = "empty-"+transform+"-" + string(round(L)/1e3) + "km-" + string(Nxy) + "-" + string(Nz) + "-one-half-dealias" + ".nc";
filename = "fine-scale-"+transform+"-" + string(round(L)/1e3) + "km-" + string(Nxy) + "-" + string(Nz) + "-one-half-dealias" + ".nc";
if exist(emptyFilename,'file')
    wvt = WVTransform.waveVortexTransformFromFile(emptyFilename);
else
    if strcmp(transform,'boussinesq')
        wvt = WVTransformBoussinesq([L, L, Lz], [Nxy, Nxy, Nz], N2=N2,latitude=latitude, Nj=floor((Nz-1)/2) );
    elseif strcmp(transform,'hydrostatic')
        wvt = WVTransformHydrostatic([L, L, Lz], [Nxy, Nxy, Nz], N2=N2,latitude=latitude, Nj=floor((Nz-1)/2) );
    else
        error('unknown transform');
    end
    wvt.addForcing(WVAdaptiveDamping(wvt));
    wvt.addForcing(WVVerticalDamping(wvt,nu=1.2e-6, kappa=1e-6));
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
model = WVModel(wvt);

model.createNetCDFFileForModelOutput(filename,outputInterval=3600,shouldOverwriteExisting=true);

model.integrateToTime(wvt.t + 2*86400);