function [ Bo, Rs, Pb, Z, rho_G ] = calcStanding( T, Pav, gamma_g, API, Rp, rho_r_initial_guess  )
%CALCSTANDING PVT parameters calculation using Standing correlations
%   Bo, Rs, Bubble Point, Deviation Factor, and other gas properties are
%   obtained via the Standing correlations available in text books.

gamma_o = 141.5/(API + 131.5);
% CALC. PVT PARAMETERS (BLACK OIL, STANDING)
x = .0125*API - .00091*T; % T is in F
Rs = gamma_g*((Pav/18.2 + 1.4)*10^x)^1.2048; % scf/STB
Bo = .9759 + .000120*(Rs*(gamma_g/gamma_o)^.5 + 1.25*T)^1.2; % T is in F
Bw = 1; % water assumed incompressible
Pb = 18.2*(((Rp/gamma_g)^.83) * 10^(-x) - 1.4);

% CALC. GAS DEVIATION FACTOR
Tpc = 168 + 325*gamma_g - 12.5*gamma_g^2; % Rankine
Ppc = 677 + 15.0*gamma_g - 37.5*gamma_g^2; % psia
Tpr = (T+460)/Tpc;
Ppr = Pav/Ppc;
% standing equation for Z in terms of rho_r

A = [.3265 -1.0700 -.5339 .01569 -.05165 .5475 -.7361 .1844 .1056 .6134 .7210];
f = @(rho_r) (A(1) + A(2)/Tpr + A(3)/Tpr^3 + A(4)/Tpr^4 + A(5)/Tpr^5)*rho_r ...
    + (A(6) + A(7)/Tpr + A(8)/Tpr^2)*rho_r^2 - A(9)*(A(7)/Tpr + A(8)/Tpr^2)*rho_r^5 ...
    + A(10)*(1 + A(11)*rho_r^2)*(rho_r^2/Tpr^3)*exp(-A(11)*rho_r^2) + 1.0 ...
    - .27*Ppr/Tpr/rho_r;
rho_r = fsolve(f, rho_r_initial_guess, optimoptions('fsolve','Display','off'));

Z = .27*Ppr/Tpr/rho_r;
% rho_G & Bg ARE KNOW DERIVED
rho_G = 2.7*gamma_g*Pav/Z/(T+460); % PCF

end

