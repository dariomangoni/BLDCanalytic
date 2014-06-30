clear all
clc

rotor_type = 'internal';
p = 2; %pole number (half of...)
Qs = 12; %slot number

alfa_p = 1; %Pole-arc/pole-pitch ratio 
g = 0.00075; %airgap
hm = 0.0045; %radial thickness of magnets
delta = g + hm; % effective airgap (ho sentito già ventisette definizioni di effective airgap)
Rs = 0.02975; %stator radius
b0 = 0.00205; %slot opening
Bremanence = 0.38; %magnet remanence
m = 3; %number of phases
mir = 1.05; %relative recoil permeability
mi0 = pi*4e-7; %absolute permability
magnetization = 'radial'; %magnetization direction
alfa_c = 0; %commutation angle; >0 for advanced commutation, <0 for retarded
l_ef = 0.05;

y_q = round(Qs/2/p);
alfa_y = 2*pi*y_q/Qs; %winding pitch (coil span)
W = 24; % number of series turns per phase
windoverlap = 'nonoverlapping';

if strcmp(rotor_type,'internal')
    Rm = Rs - g;
    Rr = Rm - hm;
elseif strcmp(rotor_type,'external')
    Rm = Rs + g;
    Rr = Rm + hm;
end

q = Qs/2/p/m; %stator slot per pole per phase
% theta0 = -alfa_c/p -5* pi/6/p; %initial angle (stator to rotor?)
theta0 = -alfa_c/p +1* pi/6/p; %initial angle (stator to rotor?)

save('motor_data.mat')