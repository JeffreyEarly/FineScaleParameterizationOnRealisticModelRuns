
% restartFile = "fine-scale-hydrostatic-500km-512-512.nc";
restartFile = "fine-scale-hydrostatic-50km-256-443.nc";
filename = "fine-scale-restart-with-eddy-hydrostatic-50km-256-443.nc";
shouldAddInertialOscillations = true;
shouldAddEddy = true;

wvt = WVTransform.waveVortexTransformFromFile(restartFile,iTime=Inf);
% wvt.removeAll;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Add some inertial oscillations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if shouldAddInertialOscillations
    U_io = 0.25; % m/s, amplitude of inertial oscillations
    Ld = 120; % m, Lh e-fold depth of initial NIOs

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
    He = wvt.Lz/5;
    U = 0.25; % m/s
    x0 = (max(wvt.x)-min(wvt.x))/2;
    y0 = (max(wvt.y)-min(wvt.y))/2;

    H = @(z) exp(-(z/He/sqrt(2)).^2 );
    F = @(x,y) exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2);
    psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*H(z).*(F(x,y) - (pi*Le*Le/(wvt.Lx*wvt.Ly)));
    wvt.setGeostrophicStreamfunction(psi);
end

% N2 = wvt.N2 + squeeze( (N2(floor(wvt.Nx/2),floor(wvt.Ny/2),:)) );
model = WVModel(wvt);
model.createNetCDFFileForModelOutput(filename,outputInterval=3600,shouldOverwriteExisting=true);
model.integrateToTime(wvt.t + 4*86400);