function [ rho_bar, FlowPattern, yl ] = calcBB( lambda_L, v_m, v_sL, D, rho_L, rho_G, sigma, phi )
%CALCBB calculates the flow regime in which the two-phase fluid is flowing
%based on the Beggs and Brill method.
% The hydrostatic pressure drop is also calculated through this function.

phi = radtodeg(phi);
g = 32.17;
N_FR = v_m^2/g/D;
N_vL = 1.938*v_sL*(rho_L/sigma)^.25;

L_1 = 316*lambda_L^.302;
L_2 = .0009252*lambda_L^(-2.4684);
L_3 = .10*lambda_L^(-1.4516);
L_4 = .5*lambda_L^(-6.738);

S1 = 'Segregated';
S2 = 'Transition';
S3 = 'Intermittent';
S4 = 'Distributed';


% Determining Flow Regime
if lambda_L < .01 && N_FR < L_1 || lambda_L >= .01 && N_FR < L_2
    FlowPattern = S1;
elseif lambda_L >= .01 && L_2 < N_FR && N_FR <= L_3
    FlowPattern = S2;
elseif .01 <= lambda_L && lambda_L < .4 && L_3 < N_FR && N_FR <=L_1 || ...
        lambda_L >= .4 && L_3 < N_FR && N_FR <= L_4
    FlowPattern = S3;
elseif lambda_L < .4 && N_FR >= L_1 || lambda_L >= .4 && N_FR > L_4
    FlowPattern = S4;
end

% Determining Constants 
if strcmp(S1, FlowPattern) == 1
    a = .98; b = .4846; c = .0868; d = .011; e = -3.768; f = 3.539; g = -1.614;
elseif strcmp(S3, FlowPattern) == 1
    a = .845; b = .5351; c = .0173; d = 2.96; e = .305; f = -.4473; g = .0978;
elseif strcmp(S4, FlowPattern) == 1
    a = 1.065; b = .5824; c = .0609;
end

% Determining Caclulation Method
if strcmp(S1, FlowPattern) == 1 || strcmp(S3, FlowPattern) == 1 || ...
        strcmp(S4,  FlowPattern) == 1
    if strcmp (S4, FlowPattern) == 1
        C = 0;
    else
        C = (1-lambda_L)*log(d*lambda_L^e*N_vL^f*N_FR^g);
    end
    psi = 1 + C*(sin(degtorad(1.8*phi)) - .333*(sin(degtorad(1.8*phi))^3));
    yl0 = a*lambda_L^b/N_FR^c;
    yl = yl0*psi;
else % for the transition flow regime
    for u = 1:2
        if u == 1 % calculating yl for segregated flow regime
            a = .98; b = .4846; c = .0868; d = .011; e = -3.768; f = 3.539; g = -1.614;
            C = (1-lambda_L)*log(d*lambda_L^e*N_vL^f*N_FR^g);
            psi = 1 + C*(sin(degtorad(1.8*phi)) - .333*(sin(degtorad(1.8*phi))^3));
            yl0 = a*lambda_L^b/N_FR^c;
            yl = yl0*psi;
            yl_seg = yl;

        elseif u == 2 % calculating yl for intermittent flow regime
            a = .845; b = .5351; c = .0173; d = 2.96; e = .305; f = -.4473; g = .0978;
            C = (1-lambda_L)*log(d*lambda_L^e*N_vL^f*N_FR^g);
            psi = 1 + C*(sin(degtorad(1.8*phi)) - .333*(sin(degtorad(1.8*phi))^3));
            yl0 = a*lambda_L^b/N_FR^c;
            yl = yl0*psi;
            yl_int = yl;
        end
    end
    
    A = (L_3 - N_FR)/ (L_3 - L_2);
    B = 1 - A;
    yl = A*yl_seg + B*yl_int;
end

rho_bar = yl*rho_L + (1-yl)*rho_G; 


end

