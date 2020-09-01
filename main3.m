%% INITIALIZE MATLAB
clear all;
clc;
close all;
format long;

%% DASHBOARD

gc = 32.174; 

% %%%%%%%%%%%%%%%%%%%%% %
% 1. Entering Well Data %
% %%%%%%%%%%%%%%%%%%%%% %

Qo = 500; % STB/day
Qw = 200; % STB/day
Rp = 1000; % scf/STB
D = 1.995/12; % ft
depth = 4000; % ft
T = 130; % degrees Fahrenheit
e = .0006; 
Area = pi*D^2/4; % ft2
gamma_g = .75;
API = 30;
gamma_o = 141.5/(API + 131.5);
rho_w = 62.4; % PCF
rho_o = rho_w*gamma_o; % PCF
phi = degtorad(90); % tubing inclination angle 

% INCREMENTS
dL = 200; % ft
M = depth/dL + 1;
P = zeros(M, 1);
P(M) = 3000; % psia, pressure at depth
dP = zeros(M,1);
dP(M) = 50;
h = 0:dL:depth;
iteration = zeros(M,1);

%% MARCHING ALGORITHM 
for n = 1:2
    for i = M:-1:2
        crit = inf;
        while crit > 1e-12
            
            % 2. Guess the Pressure 
            P(i-1) = P(i) - dP(i);            
            
            % 3. Calculate Average P
            Pav(i) = (P(i) + P(i-1))/2;
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            % 4. Determine Fluid Properties @ Pav & T %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            
            % CALC. PVT PARAMETERS (BLACK OIL, STANDING)
            [ Bo(i), Rs, Pb, Z(i), rho_G(i) ] = calcStanding(T, Pav(i), gamma_g, API, Rp, 1);
            Bg = .0282793*Z(i)*(T+460)/Pav(i); % ft3/scf
            Bw = 1;
       
            % %%%%%%%%%%%%%%%%%%%%%%%%%% %
            % 5. Calc. Pressure Gradient %
            % %%%%%%%%%%%%%%%%%%%%%%%%%% %
            
            % TWO-PHASE PROPERTIES
            qo = Qo*Bo(i); % down-hole flow-rates, bbl/day
            qw = Qw*Bw; % bbl/day
            qG = (Rp-Rs)*Qo*Bg/5.615; % bbl/day

            fo = qo/(qo + qw); % Oil Fraction
            fw = qw/(qo + qw); % Water Fraction

            qL = fo*qo + fw*qw; % bbl/day

            lambda_L = qL/(qL+qG); % no-slip liquid hold-up
            lambda_G(i) = 1 - lambda_L; % gas void fraction

            rho_L(i) = rho_o*fo + rho_w*fw; % PCF
            rho_n = rho_L(i)*lambda_L + rho_G(i)*lambda_G(i); % PCF

            v_sL = qL/Area*5.615/24/60/60; % ft/s
            v_sG = qG/Area*5.615/24/60/60; % ft/s
            v_m = v_sG + v_sL; % ft/s

            % mu_L = mu_o*fo + mu_w*fw;
            % mu_n = mu_L*lambda_L + mu_G*lambda_G;
            % sigma_L = sigma_o*fo + sigma_w*fw;
            

            % PRESSURE GRADIENT        
            f_tp(i) = calcFF(rho_n, v_m, D, n); 
            dP_frictional(i) = (f_tp(i) * rho_n * v_m^2)/2/gc/D;
            dP_hydrostatic(i) = rho_n*sin(phi);
            dPdL = dP_hydrostatic(i) + dP_frictional(i); % PCF
            
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            % 6. Determining Error Criteria %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
            
            dP_guess(i) = P(i)-P(i-1); % psi
            dP(i) = dPdL*dL/144; % psi
            crit = abs(dP(i) - dP_guess(i));
            iteration(i) = iteration(i) + 1;
        end
    end
    if n == 1    
        P_PC = P;
    elseif n == 2
        P_FB = P;
    end
end

%% RESULTS
figure
plot(P_PC, -h, '+r' ,'LineWidth', 2)
grid on
title('Tubing Pressure Profile');
ylabel('Depth, ft');
xlabel('Pressure, psi');
hold on
plot(P_FB, -h, '*b' ,'LineWidth', 2)
legend('Poettmann & Carpenter', 'Fancher & Brown');

disp('Wellhead Pressure equals to: '); 
disp([num2str(P_PC(1)), ' using Poettman & Carpenter friction factor corr.']);
disp([num2str(P_FB(1)), ' using Fancher & Brown friction factor corr.']);

figure
h1 = area((P_FB - dP_hydrostatic'.*dL./144), -h, 0);
title('Tubing Pressure Profile, Pressure Drop Contributions');
ylabel('Depth, ft');
xlabel('Pressure, psi');
hold on
h2 = area((P_FB - dP_frictional'.*dL./144), -h, 0);
set(h1,'FaceColor',[0,0.25,0.25]);
set(h2,'FaceColor',[0,0.5,0.5]);
legend('Frictional Pressure Drop','Hydrostatic Pressure Drop');