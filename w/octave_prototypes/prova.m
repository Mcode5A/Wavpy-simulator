clear all
[kx, ky] = meshgrid([1e-3:1e-3:1e3], 0);
[B, S, PSI, PHI] = elfouhaily_spectrum(21, 0, 0.84, kx, ky);

figure(1);clf;
loglog(kx, B);
xlabel('wavenumber k [rd/m]');
ylabel('curvature spectrum B(k)');
grid on;
axis([1e-3 1e4 1e-4 1e0]);


figure(2);clf;
loglog(kx, S);
xlabel('wavenumber k [rd/m]');
ylabel('Omnidirectional elevation spectrum S(k)');
grid on;
axis([1e-3 1e4 1e-16 1e3]);


