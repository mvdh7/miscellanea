function [d13poc, ep] = rau1996(r, tempc, psal, talk, pco2, d13dic, ...
    si, phos)

% %% Test inputs
% r = 1e-5; % m
% tempc = 17; % degC
% psal = 35;
% talk = 2325; % umol / kg-sw
% pco2 = 400; % uatm
% d13dic = 1.7; % permille
% si = 0; % umol/kg-sw
% phos = 0; % umol/kg-sw

% Set base values from Table 1
ed = 0.7; % permille
ef = 25; % permille
% Qr = 4.06e-16 * 1; % mol / (cell * s)
P = 1e-4; % m / s
mu = 1.16e-5; % 1 / s

% DT = diffusivity of CO2 in m^2 / s
Ed = 19510; % J / mol
R = 8.3144598; % J / (K * mol)
tempk = tempc + 273.15; % K
DT = 5.019e-6 * exp(-Ed ./ (R * tempk)); % Eq. (8)
DT = DT .* (0.9508 - 7.389e-4 * tempc); % Eq. (9)
% DT = 1.45e-9; % Table 1 base value

% Use CO2SYS for [OH] and [CO2(aq)] and convert to volume units
co2s = CO2SYSv2(talk, pco2, 1, 4, psal, tempc, tempc, ...
    0, 0, si, phos, 3, 10, 3);
rho = cMP81(tempc, psal) * 1e-3; % kg / dm^3
oh = co2s(:, 10) * 1e-6 .* rho; % mol / dm^3
co2aq = co2s(:, 8) * 1e-6 .* rho; % mol / dm^3

% rk = reacto-diffusive length in m
Ek = 6.28e4; % J / mol
kExpRatio = exp(-Ek ./ (R * tempk)) / exp(-Ek / (R * 298.15));
k1 = 8500 * kExpRatio; % dm^3 / (mol * s)
k2 = 0.03* kExpRatio; % 1 / s
rk = sqrt(DT ./ (k1 .* oh + k2)); % Eqs. (6) & (7)
% rk = 2.06e-4; % Table 1 base value

% Qs = uptake rate per unit surface area in mol / (m^2 * s)
V = (r*1e6).^3 * 4*pi/3; % cell volume in um^3
gc = 3.154e-14 * V.^0.758; % Eq. (5) cell carbon content in mol-C
mui = mu * 2; % for 12:12 light:dark cycle
Qr = gc .* mui; % Eq. (4)
Qs = Qr ./ (4 * pi * r.^2); % mol / (m^2 * s)
% Qs = 3.23e-7; % Table 1 base value

% Rau et al. (1996) Eq. (2)
d13co2aq = d13dic + 23.644 - 9701.5 ./ tempk;

% Rau et al. (1996) Eq. (13)
d13poc = d13co2aq - ef + (ef - ed) * (Qs ./ (co2aq * 1e3)) .* ...
    (r ./ (DT .* (1 + r ./ rk)) + 1 / P);

% ep
ep = d13co2aq - d13poc;

end % function rau1996
