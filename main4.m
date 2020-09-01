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


D = 6/12; % ft
depth = 4000; % ft
mu_L = .97; % cp
mu_G = .016; % cp
sigma_L = 8.41; % dynes/cm (or lbf/ft?)
T = 180; % degrees Fahrenheit
e = 0.00006; % ft 
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
P(M) = 1700; % psia, pressure at depth
dP = zeros(M,1);
dP(M) = 10;
h = 0:dL:depth;
iteration = zeros(M,1);

Qo = 10000; % STB/day
Qw = 200; % STB/day
Rp = 1000; % scf/STB


%% MARCHING ALGORITHM 
for i = M : -1 : 2
    crit(i) = inf;
    while crit(i) > 1e-6

        % 2. Guess the Pressure 
        P(i-1) = P(i) - dP(i);            

        % 3. Calculate Average P
        Pav(i) = (P(i) + P(i-1))/2;

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
        % 4. Determine Fluid Properties @ Pav & T %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        [ Bo(i), Rs, Pb, Z(i), rho_G ] = calcStanding(T, Pav(i), gamma_g, API, Rp, 10);
        Bg = .0282793*Z(i)*(T+460)/Pav(i); % ft3/scf
        Bw = 1;

        % %%%%%%%%%%%%%%%%%%%%%%%%%% %
        % 5. Calc. Pressure Gradient %
        % %%%%%%%%%%%%%%%%%%%%%%%%%% %

        % TWO-PHASE PROPERTIES % % % % % % % % % % % % % % % % % % % % % %
        % calculating average two-phase fluid properties 
        qo = Qo*Bo(i); % down-hole flow-rates, bbl/day
        qw = Qw*Bw; % bbl/day
        qG = (Rp-Rs)*Qo*Bg/5.615; % bbl/day

        fo = qo/(qo + qw); % Oil Fraction
        fw = qw/(qo + qw); % Water Fraction

        qL = fo*qo + fw*qw; % bbl/day

        lambda_L = qL/(qL+qG); % no-slip liquid hold-up
        lambda_G = 1 - lambda_L; % gas void fraction
        
        % average densities
        rho_L = rho_o*fo + rho_w*fw; % PCF
        rho_n = rho_L*lambda_L + rho_G*lambda_G; % PCF
        
        % superficial velocities
        v_sL = qL/Area*5.615/24/60/60; % ft/s
        v_sG = qG/Area*5.615/24/60/60; % ft/s
        v_m = v_sG + v_sL; % ft/s
        
        % average viscosity
        mu_n = mu_L*lambda_L + mu_G*lambda_G;

        % % PRESSURE GRADIENT % % % % % % % % % % % % % % % % % % % % % % %
        % calculating RhoBar, Holdup, and Flow Regime using the predefined 
        % Beggs & Brill function      
        [rho_bar(i), FlowRegime, yl(i)] = calcBB(lambda_L, v_m, v_sL, D,...
            rho_L, rho_G, sigma_L, phi);
        FR{i} = FlowRegime;
        
        % calculating Friction Factor using predefined function 
        f_tp(i) = calcFrictionPressure(rho_n, v_m, mu_n, D, lambda_L, yl(i));
        
        % now pressure gradients are derived based on rho_bar and friction
        % factor
        dP_frictional(i) = (f_tp(i) * rho_n * v_m^2)/2/gc/D;
        dP_hydrostatic(i) = rho_bar(i)*sin(phi);
        dPdL = dP_hydrostatic(i) + dP_frictional(i); % PCF

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
        % 6. Determining Error Criteria %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

        dP_guess = P(i) - P(i-1); % psi
        dP(i) = dPdL*dL/144; % psi
        crit(i) = abs(dP(i) - dP_guess);
        iteration(i) = iteration(i) + 1;
        if iteration(i) > 100
            error('Pressure did not converge');
            break
        end
    end
end



%% RESULTS
figure
plot(P, -h, '+r' ,'LineWidth', 2)
grid on
title('Tubing Pressure Profile');
ylabel('Depth, ft');
xlabel('Pressure, psi');

figure
I = imread('beggs.png');
if strcmp(FR{M}, 'Segregated') == 1
    I(170:270, 170:270, 1) = I(170:270, 170:270, 1) - 100;
elseif strcmp(FR{M}, 'Intermittent')
    I(105:155, 340:440, 1) = I(105:155, 340:440, 1) - 100;
elseif strcmp(FR{M}, 'Transition') == 1
    I(265:315, 365:465, 1) = I(265:315, 365:465, 1) - 120;
elseif strcmp(FR{M}, 'Distributed') == 1
    I(20:70, 220:320, 1) = I(20:70, 220:320, 1) - 100;   
end
imshow(I);

figure
plot(yl, -h, '+k', 'LineWidth', 2);
grid on
xlabel('Liquid Hold-up');
ylabel('Depth, ft');
title('Liquid Hold-up vs. Depth');