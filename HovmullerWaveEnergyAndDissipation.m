filename = "fine-scale-restart-with-eddy-no-igw-hydrostatic-50km-128-222-one-half-dealias-take-2.nc";
% filename = "fine-scale-restart-with-eddy-no-igw-hydrostatic-50km-128-222-one-half-dealias-take-2-cyclone.nc";
[wvt,ncfile] = WVTransform.waveVortexTransformFromFile(filename,iTime=1);
t = ncfile.readVariables('wave-vortex/t');
%%
E_wave = zeros(length(t),wvt.Nz);
E_damp = zeros(length(t),wvt.Nz);
for i=1:length(t)
    % disp(i)
    wvt.initFromNetCDFFile(ncfile,iTime=i);
    E_w = wvt.u_w.^2 + wvt.v_w.^2 + shiftdim(wvt.N2,-2) .* wvt.eta_w.^2;
    E_wave(i,:) = squeeze( mean(mean(E_w,1),2) );

    svv = wvt.forcingWithName("adaptive damping");
    Fp = sqrt(-wvt.uvMax * svv.damp) .* wvt. Ap;
    Fm = sqrt(-wvt.uvMax * svv.damp) .* wvt. Am;
    F0 = 0*sqrt(-wvt.uvMax * svv.damp) .* wvt. A0;

    Fu = wvt.transformToSpatialDomainWithF(Apm=wvt.UAp.*wvt.phase.*(Fp) + wvt.UAm.*wvt.conjPhase.*(Fm),A0=wvt.UA0.*(F0));
    Fv = wvt.transformToSpatialDomainWithF(Apm=wvt.VAp.*wvt.phase.*(Fp) + wvt.VAm.*wvt.conjPhase.*(Fm),A0=wvt.VA0.*(F0));
    Feta = wvt.transformToSpatialDomainWithG(Apm=wvt.NAp.*wvt.phase.*(Fp) + wvt.NAm.*wvt.conjPhase.*(Fm),A0=wvt.NA0.*(F0));

    damp = Fu .* Fu + Fv .* Fv + shiftdim(wvt.N2,-2) .* Feta .* Feta;
    E_damp(i,:) = squeeze( sum(sum(damp,1),2) );
end

%%
figure
tl = tiledlayout(2,1);
% title(tl,"medium-res cyclone")
title(tl,"medium-res anticyclone")

nexttile
pcolor(t/86400,wvt.z,E_wave.'), shading flat
ylim([-1500 0])
cb = colorbar("eastoutside");
clim([0 2e-3])
cb.Label.String = 'wave energy (m^2/s^2)';

nexttile
pcolor(t/86400,wvt.z,log10(E_damp).'), shading flat
ylim([-1500 0])
clim([-10 -5])
cb = colorbar("eastoutside");
cb.Label.String = 'energy dissipation (m^2/s^3)';

exportgraphics(gcf,"hovmuller-medium-res-anticyclone.png")