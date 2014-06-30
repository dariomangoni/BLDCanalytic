clc
clear all
close all

load('motor_data.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adattamento valori
Ri = Rs;
M0 = Bremanence/mi0;    % definizione estrapolata da 5.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Valore di test
r = (Ri - g + Ri)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EFFECT OF SLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(rotor_type,'internal')
    R0 = Ri - g;
    Rm = Ri - g - hm;
else
    error 'Motore non supportato'
end

% EFFECT OF SLOTTING ON THE EFFECTIVE FLUX
g1 = g + hm/mir;
tau_t = 2*pi*Ri/Qs;
gamma = 4/pi*( b0/(2*g1)*atan(b0/(2*g1)) - log(sqrt(1 + (b0/(2*g1))^2)));     %5.98

kCarter = tau_t/(tau_t - gamma*g1);

g_slotted = g + (kCarter -1)*g1;      % 5.95
Ri_slotted = Ri + (kCarter -1)*g1;     % 5.96

% EFFECT OF SLOTTING ON THE FIELD DISTRIBUTION - THE RELATIVE PERMEANCE CALCULATION
if strcmp(rotor_type,'internal')
    y_slotpar = r - Ri + g1;        % parametro y 5.113
else
    error 'Motore non supportato'
end

syms v_sym

a_slotpar_squared = 1 + (2*g1/b0)^2;  % parametro a 5.104

f = 0.5*log(...
    (sqrt(a_slotpar_squared + v_sym^2) + v_sym)/...
    (sqrt(a_slotpar_squared + v_sym^2) - v_sym)...
    )...
    + 2*g1/b0*atan(...
    2*g1/b0*v_sym/sqrt(a_slotpar_squared + v_sym^2)...
    ) - y_slotpar*pi/b0;

v_slotpar = single(solve(f, v_sym));   % 5.112

beta_slotpar = 0.5*(1 - (sqrt (1 + (b0/2/g1)^2 * (1 + v)^2))^(-1)); % 5.110

phi_mo = Bremanence/mi0/mir*hm; % 5.107

Bmax_slotpar = mi0*phi_mo/g1;   %5.108

lambda_0_tilde = 1/kCarter*(1-1.6*beta_slotpar*b0/tau_t);

n = 1:lambdaHarmonicNum;
rag1 = n*b0/tau_t;

lambda_comp = - beta_slotpar*4/(n*pi)*(0.5 + (rag1^2)/(0.78125 - 2*rag1^2))*sin(1.6*pi*rag1);

error








%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (RADIAL) MAGNETIC FIELD PRODUCED BY MAGNETS IN AIRSPACE
% Il vettore Br_theta contiene B(theta) radiale calcolato sulla
% circonferenza di raggio r.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_harmonicNum = 500;
theta = (0:0.005:pi)';  % The variable "theta" which is used in the 
                        % cylindrical model, is with reference to the
                        % centre of a magnet pole.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(stator_type,'slotted')
    Ri = Ri_slotted;
    g = g_slotted;
end


if strcmpi(rotor_type,'internal')
    R0 = Ri - g;
    Rm = Ri - g - hm;
else
    error 'Motore non supportato'
end

row = 1;
Br_comp = zeros(B_harmonicsnum, 2);
for n = 1:2:B_harmonicsnum

	np = n*p;
    
    Br_comp(row, 1) = n;
    
    Mn = 2*M0*alfa_p/(n*pi*alfa_p/2)*sin(n*pi*alfa_p/2);	%5.137
    
    if np ~= 1		
        Br_comp(row, 2) = (Mn*np/mir/((np^2) -1))*...
            ((np - 1) + 2*(Rm/R0)^(np + 1) - (np + 1)*(Rm/R0)^(2*np))/...
            ((mir + 1)/mir * (1 - (Rm/Ri)^(2*np)) - (mir - 1)/mir * ((R0/Ri)^(2*np) - (Rm/R0)^(2*np)))*...
            ((r/Ri)^(np -1)*(R0/Ri)^(np +1) + (R0/r)^(np +1));      % 5.62 mod
	else
		Br_comp(row, 2) = Mn/(2*mir)*...
			((R0/Ri)^2 - (Rm/Ri)^2 + (Rm/Ri)^2 * log((R0/Rm)^2))/...
			((mir +1)/mir*(1 - (Rm/Ri)^2) - (mir -1)/mir*((R0/Ri)^2 - (Rm/R0)^2))*...
			(1 + (Ri/r)^2);     % 5.57 mod (M1 = Mn(n = 1))
    end
    
    row = row +1;
end

Br_theta = zeros(length(theta), 2);

for row = 1:length(theta)   
    Br_theta(row, :) = [theta(row), sum(Br_comp(:, 2).*cos(Br_comp(:, 1)*p*theta(row)))];   
end

clear row

plot(Br_theta(:,1),Br_theta(:,2),'DisplayName','Br(1:32,1:2)','YDataSource','Br(1:32,1:2)');figure(gcf)
