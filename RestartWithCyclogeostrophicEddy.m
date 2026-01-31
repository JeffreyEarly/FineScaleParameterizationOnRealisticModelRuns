
resolution = 64;
eddysign = 1;

if eddysign == 1
    eddyname = "cyclone";
elseif eddysign == -1
    eddyname = "anticyclone";
end

if resolution == 64
    restartFile = "fine-scale-hydrostatic-50km-64-111-one-half-dealias.nc";
    filenamePrefix = "fine-scale-hydrostatic-50km-64-111-";
elseif resolution == 128
    restartFile = "fine-scale-hydrostatic-50km-128-222-one-half-dealias.nc";
    filenamePrefix = "fine-scale-hydrostatic-50km-128-222-";
elseif resolution == 256
    restartFile = "fine-scale-hydrostatic-50km-256-443-one-half-dealias.nc";
    filenamePrefix = "fine-scale-hydrostatic-50km-128-222-";
end

filename = filenamePrefix + eddyname + "-cyclogeostrophic.nc";

shouldAddInertialOscillations = true;
shouldAddEddy = true;

wvt = WVTransform.waveVortexTransformFromFile(restartFile,iTime=Inf);
wvt.removeAll;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Add some inertial oscillations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if shouldAddInertialOscillations
    U_io = 0.25; % m/s, amplitude of inertial oscillations
    Ld = 60; % m, Lh e-fold depth of initial NIOs

    u_NIO = @(z) U_io*exp((z/Ld));
    v_NIO = @(z) zeros(size(z));

    wvt.addInertialMotions(u_NIO,v_NIO);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Add an eddy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Consider a "shallow eddy" where the density anomaly sits close to the
% surface. This example was constructed in Early, Hernández-Dueñas, Smith,
% and Lelong (2024), https://arxiv.org/abs/2403.20269

if shouldAddEddy
    Le = wvt.Lx/8; % 7 is the minimum you'll want to go, 
    He = wvt.Lz/15;
    U = eddysign*0.25; % m/s
    x0 = (max(wvt.x)-min(wvt.x))/2;
    y0 = (max(wvt.y)-min(wvt.y))/2;

    H = @(z) exp(-(z/He).^2 );
    F = @(x,y) exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2);
    A = U*(sqrt(2)/Le)*exp(1/2);

    u = @(x,y,z) - A * (y - y0) .* F(x,y) .* H(z);
    v = @(x,y,z)   A * (x - x0) .* F(x,y) .* H(z);
    eta = @(x,y,z) - A * (Le*Le/He/He) * shiftdim(1 ./ wvt.N2,-2) .* z .* F(x,y) .* H(z) .* ( A*F(x,y) .* H(z) + wvt.f );

    wvt.addUVEta(u(wvt.X,wvt.Y,wvt.Z),v(wvt.X,wvt.Y,wvt.Z),eta(wvt.X,wvt.Y,wvt.Z))
end

%%

% N2 = wvt.N2 + squeeze( (N2(floor(wvt.Nx/2),floor(wvt.Ny/2),:)) );
model = WVModel(wvt);
model.createNetCDFFileForModelOutput(filename,outputInterval=6*3600);
model.integrateToTime(wvt.t + 15*86400);