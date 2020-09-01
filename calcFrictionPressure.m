function [ f_tp ] = calcFrictionPressure( rho_m, v_m, mu_m, D, lambda_L, yl  )
%CALCFRICTIONPRESSURE returns the Beggs & Brill two-phase friction factor 

Re = 1488*rho_m*v_m*D/mu_m;
e = 0;

syms f_n
eqnf = 1/sqrt(f_n) == -4*log10(e/3.7065 - 5.0452/Re*log10((e^1.1098)/2.8257...
    + (7.149/Re)^.8981)); 
f_n = solve(eqnf, f_n);
f_n = double(f_n);

x = lambda_L/yl^2;

if x>1 && x<1.2
    S = log(2.2*x-1.2);
else
    S = log(x)/(-.0523+3.182*log(x)-.8725*(log(x))^2+.01853*(log(x))^4);    
end

f_tp = f_n*exp(S);


end

