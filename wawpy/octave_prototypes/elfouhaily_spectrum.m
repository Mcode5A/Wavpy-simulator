function [B, S, PSI, PHI] = elfouhaily_spectrum(U10, theta, Omega, kx, ky)
% [B, S, PSI, PHI] = elfouhaily_spectrum(U10, theta, Omega)
%
%  Output:
%     B     : Curvature Spectrum
%     S     : Omnidirectional spectrum
%     PSI   : Directional spectrum
%     PHI   : Angular spreading function
%
%  Input:
%     U10   : windspeed at 10m above surface [m/s]
%     theta : angle between wind and waves [deg]
%     Omega : inverse dimensionless wave age > 0.83
%     kx    : wavenumber in x direction [rad/m] (can be a matrix)
%     ky    : wavenumber in y direction [rad/m] (can be a matrix)
%

% ----- Checking validity range of input and derived parameters
if (Omega < 0.84)
  str = sprintf("Omega must be larger than 0.84, but Omega=%f.\n", Omega);
  error(str);
endif

Omegac = Omega*cos(theta*pi/180);
if (Omegac<0.84 | Omegac>5)
  str = sprintf("Omega_c must be larger than 0.83 and smaller than 5, but Omega_c=%f.\n", Omegac);
  error(str);
end;

k = sqrt(kx.^2 + ky.^2);
phi = atan2(ky, kx);

% ----- constants
g = 9.81; % [m/s^2] gravity acceleration

% ----- Computing the short wave part of the spectrum
cm = 0.23; % [m/s] Elfouhaily 1997.
km = (2*g/cm**2); % assuming deep waters => =370 rad/m
%%%c = sqrt(g./k);  % Phase velocity (deep waters assumption)
c = sqrt(g./k.*(1+(k/km).^2));

ufric = sqrt(1e-3 * (0.8 + 0.065*U10))*U10; %

if (ufric <= cm)
  alpham = (1+log(ufric/cm))*0.01;
elseif (ufric > cm)
  alpham = (1 + 3*log(ufric/cm))*0.01;
endif


Fm = exp(-0.25*((k./km -1 ).^2));  % -> matrix

Bh = 0.5 * alpham *cm./c.*Fm; % -> matrix

% ----- Computing long wave part of the spectrum
cp = U10/Omega;  % Phase velocity at the spectral peak
kp = g/(cp^2);  % Wave number of the spectral peak (assumtion of deep waters. TODO ref. needed)
alphap = 0.006*sqrt(Omega);
%c(find(isinf(c))) = 0;

if ((0.84<=Omegac) && (Omegac<=1))
  gmma = 1.7;
elseif ((1<Omegac) && (Omegac<=5))
  gmma = 1.7 + 6*log10(Omegac);
end


sigma = 0.08*(1+4*Omegac^(-3));
Gmma = exp(-((sqrt(k./kp)-1).^2)./(2*sigma^2));  % -> matrix

Jp = gmma.^Gmma;  % -> matrix
Lpm = exp(-5/4*(kp./k).^2);  % -> matrix
%Fp = Lpm.*Jp.*exp(-Omega/sqrt(10)*(sqrt(k./kp)-1));  % -> matrix
Fp = exp(-Omega/sqrt(10)*(sqrt(k./kp)-1));  % -> matrix
Bl = 0.5*alphap*cp./c.*Fp;  % -> matrix
%Bl(find(isnan(Bl))) = 0;

B = Lpm.*Jp.*(Bl + Bh); % -> matrix
S = k.^(-3).*B;  % -> matrix

% ----- Computing the Spreading function
a0 = 0.25*log(2); % Natural logarithm
ap = 4;
am = 0.13*ufric/cm;
Delta = tanh(a0 + ap*(c./cp).^(2.5) + am*(cm./c).^(2.5)); % -> matrix


% ----- Computing the directional spectrum
PHI = 1/(2*pi) * (1 + Delta .* cos(2*phi));  % -> matrix

PSI = S .* PHI ./ k;

PSI(find(isnan(PSI))) = 0;


endfunction
