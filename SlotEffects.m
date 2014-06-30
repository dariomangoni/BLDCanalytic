clear all
close all
clc

syms v_sym

load('motor_data.mat');

r = 0.0405;
omega_r = 1500/60*2*pi;
alfa0 = b0/Rs;
lambda_ref = mi0/(g+hm/mir);



%% Total flux per pole reduction due to slotting
g1 = g + hm/mir;
tau_t = 2*pi*Rs/Qs; %statot tooth pitch
alfat = tau_t/Rs;
gamma = 4/pi*( b0/2/g1 * atan(b0/2/g1) - log( sqrt(1+(b0/2/g1)^2) ) );
Kc = tau_t/(tau_t-gamma*g1);

if strcmp(rotor_type,'internal')
    g_Carter = g + (Kc-1)*g1;
    Rs_Carter = Rs + (Kc-1)*g1;
else
    g_Carter = g + (Kc-1)*g1; %ATTENZIONE: è uguale a interni... dubbio
    Rs_Carter = Rs - (Kc-1)*g1;
end

%% Field distribution modification due to slotting
if strcmp(rotor_type,'internal')
    y = r-Rs-g1;
else
    y = Rs+g1-r;
end

asquared = 1 + (2*g1/b0)^2;
lambda0 = mi0/g1;

f(v_sym) = 0.5*log( (sqrt(asquared+v_sym^2) + v_sym)/(sqrt(asquared+v_sym^2) - v_sym) ) + 2*g1/b0*atan(2*g1/b0*v_sym/sqrt(asquared+v_sym^2))-y*pi/b0;
v = double(solve(f));

beta = 0.5*( 1 - 1/( sqrt(1+(b0/2/g1)^2)*(1+v^2) ) );

cont = 0;
alfa = 0:alfat/2/250:alfat/2;
lambda_tilde_array = zeros(length(alfa),3);
for alfa_ind = 1:length(alfa)
    if (alfa(alfa_ind)<0.8*alfa0)&&(alfa(alfa_ind)>=0)
        lambda = lambda0*(1-beta-beta*cos(pi/0.8/alfa0*alfa(alfa_ind)));
    elseif (alfa(alfa_ind)<=(alfat/2))&&(alfa(alfa_ind)>=0.8*alfa0)
        lambda = lambda0;
    else
        disp('Alfa error');
    end

    lambda_tilde = lambda/lambda0;
    
    cont = cont+1;     lambda_tilde_array(cont,:) = [alfa(alfa_ind), r, lambda_tilde];
end

figure('Name','Slot effects','NumberTitle','off');
plot([-flipud(lambda_tilde_array(2:end,1)*180/pi);lambda_tilde_array(:,1)*180/pi],[flipud(lambda_tilde_array(2:end,3));lambda_tilde_array(:,3)],'-x','MarkerSize',6);
set(gca,'XGrid','on')
set(get(gca,'XLabel'),'String','\alpha: angular position [°]');
set(get(gca,'YLabel'),'String','\lambda: relative permeance');
title('Relative permeance function','FontSize',13)
linex = get(gca,'XLim'); line(linex,[0 0],'LineStyle',':','Color','k');