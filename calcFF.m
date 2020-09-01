function [ f_tp ] = calcFF( rho_n, v_m, D, n )
%CALCPCF calculates Poettmann & Carpenter or Fancher Brown friction factor
%   n switches between the two correlations, n = 1 for PC and n = 2 for FB

Re_numerator = rho_n*v_m*D; % lb/ft/sec
switch n
    case 1
        f_tp = 82.93*(Re_numerator)^(-2.47); % Poettmann & Carpenter
    case 2
        f_tp = 2.88*(Re_numerator)^(-1.25); % Fancher and Brown
end
end

