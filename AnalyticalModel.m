clear all
close all
clc

% è stato modificato A3n che dava problemi. A3n ora è il vecchio A3n
% moltiplicato per Mn così da evitare una divisione per 0, le formule
% vedono Mn moltiplicato per le parentesi

load('motor_data.mat');

syms v_sym

r = Rs;
omega_r = 800/60*2*pi;
t = 0;
alfa_ma = omega_r*t + theta0;
alfa = (0:pi/500:2*pi)';

lambda_ref = mi0/(g+hm/mir);

data = zeros(length(alfa),8);
data(:,1) = alfa.*180/pi;
data(:,2) = alfa - ones(length(alfa),1).*alfa_ma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Slotting effects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Total flux per pole reduction
g1 = g + hm/mir;
tau_t = 2*pi*Rs/Qs; %stator tooth pitch
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

%% Lambda calc - Field distribution modification
if strcmp(rotor_type,'internal')
    y = r-Rs-g1;
else
    y = Rs+g1-r;
end

asquared = 1 + (2*g1/b0)^2;

f(v_sym) = 0.5*log( (sqrt(asquared+v_sym^2) + v_sym)/(sqrt(asquared+v_sym^2) - v_sym) ) + 2*g1/b0*atan(2*g1/b0*v_sym/sqrt(asquared+v_sym^2))-y*pi/b0;
if (f(0) == 0)
    v = 0;
else
    v = double(solve(f));
end

beta = 0.5*( 1 - 1/( sqrt(1+(b0/2/g1)^2)*(1+v^2) ) );

%% calc harmonic component
lambda_zero_tilde = (1-1.6*beta*b0/tau_t)/Kc;
alfa_sa = pi/Qs; %only if winding pitch (alfa_y ?) is an odd integer of slot pitch (almost always)
data(:,8) = ones(length(alfa),1).*lambda_zero_tilde;
for n=1:1500
    for alfa_ind = 1:length(alfa)
        
        lambda_mi_tilde = -beta*4/pi/n*( 0.5+(n*b0/tau_t)^2/(0.78125-2*(n*b0/tau_t)^2))*sin(1.6*pi*n*b0/tau_t);

        data(alfa_ind,8) = data(alfa_ind,8) + lambda_mi_tilde*cos(n*Qs*(alfa(alfa_ind)+alfa_sa));
    end
end

g = g_Carter;
Rs = Rs_Carter;

figure('Name','Slot effects','NumberTitle','off');
plot(data(:,1),data(:,8),'-x','MarkerSize',3);
set(gca,'XGrid','on')
set(gca,'XTick', data(1,1):30:data(end,1)+1);
set(get(gca,'XLabel'),'String','\alpha: angular position [°]');
set(get(gca,'YLabel'),'String','\lambda: relative permeance');
title('Relative permeance function','FontSize',13)
linex = get(gca,'XLim'); line(linex,[0 0],'LineStyle',':','Color','k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Magnets field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(rotor_type,'internal')
    Rm = Rs - g;
    Rr = Rm - hm;
elseif strcmp(rotor_type,'external')
    Rm = Rs + g;
    Rr = Rm + hm;
end

cont = 0;
for n=1:2:300
    np = n*p;
    
    if np == 1
        
        if strcmp(magnetization,'parallel')
            A1n = sin((np+1)*alfa_p*pi/2/p)/((np+1)*alfa_p*pi/2/p);
            A2n = 1;
            Mrn = Bremanence/mi0*alfa_p*(A1n+A2n);
            Mthetan = Bremanence/mi0*alfa_p*(A1n-A2n);
            Mn = Mrn + np*Mthetan;
            A3n = 2*Mrn-Mn;
        elseif strcmp(magnetization,'radial')
            Mrn = 2*Bremanence/mi0*alfa_p*sin(n*pi*alfa_p/2) / (n*pi*alfa_p/2);
            Mthetan = 0;
            Mn = Mrn;
            A3n = Mn;
        end
        
        KB = mi0/2/mir*( A3n*(Rm/Rs)^2 - A3n*(Rr/Rs)^2 + Mn*(Rr/Rs)*log((Rm/Rr)^2) )/...
            ((mir+1)/mir *(1-(Rr/Rs)^2) - (mir-1)/mir*( (Rm/Rs)^2 - (Rr/Rm)^2 ) );
        fBr = 1+(Rs/r)^2;
        fBtheta = -1+(Rs/r)^2;
        
    else
        
        if strcmp(magnetization,'parallel')
            A1n = ( sin((np+1)*alfa_p*pi/2/p) )/( (np+1)*alfa_p*pi/2/p );
            A2n = ( sin((np-1)*alfa_p*pi/2/p) )/( (np-1)*alfa_p*pi/2/p );
            Mrn = Bremanence/mi0*alfa_p*(A1n+A2n);
            Mthetan = Bremanence/mi0*alfa_p*(A1n-A2n);
            Mn = Mrn + np*Mthetan;
            A3n = (np-1/np)*Mrn + Mn/np;
        elseif strcmp(magnetization,'radial')
            Mrn = 2*Bremanence/mi0*alfa_p*sin(n*pi*alfa_p/2) / (n*pi*alfa_p/2);
            Mthetan = 0;
            Mn = Mrn;
            A3n = Mn*np;
        end
        
        
        if strcmp(rotor_type,'internal')
            KB = mi0*np/mir/(np^2-1) * ( A3n - Mn + 2*Mn*(Rr/Rm)^(np+1) - (A3n+Mn)*(Rr/Rm)^(2*np) )/...
                ( (mir+1)/mir * (1-(Rr/Rs)^(2*np)) - (mir-1)/mir * ( (Rm/Rs )^(2*np) - (Rr/Rm)^(2*np) ));
            fBr = (r/Rs)^(np-1)*(Rm/Rs)^(np+1) + (Rm/r)^(np+1);
            fBtheta = -(r/Rs)^(np-1)*(Rm/Rs)^(np+1) + (Rm/r)^(np+1);
        elseif strcmp(rotor_type,'external')
            KB = -mi0*np/mir/(np^2-1) * ( (A3n-Mn)*(Rm/Rr)^(2*np) +2*Mn*(Rm/Rr)^(np-1) - (A3n + Mn) ) /...
                ( (mir+1)/mir * (1-(Rs/Rr)^(2*np)) - (mir-1)/mir * ( (Rs/Rm)^(2*np) - (Rm/Rr)^(2*np) ) );
            fBr = (r/Rm)^(np-1) + (Rs/Rm)^(np-1)*(Rs/r)^(np+1);
            fBtheta = -(r/Rm)^(np-1) + (Rs/Rm)^(np-1)*(Rs/r)^(np+1);
        end
    end
    
    cont = cont+1;
    Bn(cont,:) = [n, KB * fBr ];
    
    for alfa_ind = 1:length(alfa)
        Br = KB * fBr * cos(np*(alfa(alfa_ind)-alfa_ma));
        Btheta = KB * fBtheta * sin(np*(alfa(alfa_ind)-alfa_ma));
        Mr = Mrn * cos(np*(alfa(alfa_ind)-alfa_ma));
        Mtheta = Mthetan * sin(np*(alfa(alfa_ind)-alfa_ma));

        data(alfa_ind,3) = data(alfa_ind,3) + Br;
        data(alfa_ind,4) = data(alfa_ind,4) + Btheta;
        data(alfa_ind,5) = data(alfa_ind,5) + Mr;
        data(alfa_ind,6) = data(alfa_ind,6) + Mtheta;
    end
    

end

load('motor_data.mat','Rs','g');

figure('Name','Magnets magnetization','NumberTitle','off'); hold on
plot(data(:,1),data(:,5),'r');
plot(data(:,1),data(:,6),'b');
legend('Mr','Mtheta');
set(gca,'XGrid','on')
set(gca,'XTick', data(1,1):30:data(end,1)+1);
set(get(gca,'XLabel'),'String','\theta: angular position [°]');
set(get(gca,'YLabel'),'String','M: magnetization [A/m]');
linex = get(gca,'XLim'); line(linex,[0 0],'LineStyle',':','Color','k');
title('Magnetization of magnets','FontSize',13)

figure('Name','Magnets field','NumberTitle','off'); hold on
plot(data(:,1),data(:,3),'r');
plot(data(:,1),data(:,4),'b');
legend('Br','Btheta');
set(gca,'XGrid','on')
set(gca,'XTick', data(1,1):30:data(end,1)+1);
set(get(gca,'XLabel'),'String','\theta: angular position [°]');
set(get(gca,'YLabel'),'String','B: magnetic field [T]');
linex = get(gca,'XLim'); line(linex,[0 0],'LineStyle',':','Color','k');
title('Magnetic field by magnets','FontSize',13)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Induced emf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = -pi/6/omega_r/p:0.00001:(2*pi-pi/6)/omega_r/p;
ea = zeros(length(t),1);
eb = ea;
ec = ea;
for n_ind = 1:length(Bn)
    n = Bn(n_ind,1);
    Phi = 2*Bn(n_ind,2)*Rs*l_ef*lambda_zero_tilde;
    Kdn = sin(q*n*p*pi/Qs) / (q*sin(n*p*pi/Qs)); %winding distribution factor (equation valid only if q is integer)
    Kpn = sin(n*p*alfa_y/2);
    Kdpn = Kdn*Kpn;
    for t_ind = 1:length(t)
        alfa_ma = omega_r*t(t_ind) + theta0;
        ea(t_ind) = ea(t_ind) + omega_r*Phi*W*Kdpn*sin(n*p*alfa_ma);
    end
end


figure('Name','Induced emf','NumberTitle','off'); hold on
plot(t,ea,'-x','MarkerSize',3);
legend('Phase A');
set(gca,'XGrid','on')
% set(gca,'XTick', data(1,1):30:data(end,1)+1);
set(get(gca,'XLabel'),'String','t: time [s]');
set(get(gca,'YLabel'),'String','ea: phase A emf [V]');
linex = get(gca,'XLim'); line(linex,[0 0],'LineStyle',':','Color','k');
title('Induced emf for phase A','FontSize',13)
save('prove\ea.mat','ea','t');
t = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Windings field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(windoverlap,'overlapping')
    fatt = 6;
elseif strcmp(windoverlap,'nonoverlapping')
    fatt = 3;
end

%% find nontriplen odd harmonics of current
cont2 = 0;
for cont=1:25
    if ( rem(cont,2) && rem(cont,3) ) % check if it's not multiple of 2 or 3
        cont2 = cont2+1;
        u_array(cont2) = cont;          % then add it in the array
    end
end

%% initialize matrixes
load CurrentHarmonics.mat
% I_array = 0:1/(length(u_array)-1):1;
% I_array = fliplr(I_array);
% I_array = [0.8 0.3 0.1 0.05 0 0 0 0 0];
theta_u_array = zeros(length(I_array),1);

for alfa_ind = 1:length(alfa)
    %% calc field generated by windings
    B_wind = 0;
    for u_index = 1:length(u_array) % current harmonics index

        u = u_array(u_index);
        theta_u = theta_u_array(u_index);
        I = I_array(u_array(u_index));

        sum_v = 0;
        for v = 0:1:50 % harmonics index
            if (~rem(v/p+u,fatt)) || (~rem(v/p+u,fatt)) || (v/p+u==0) || (v/p-u==0)
                Ksov = sin (v*b0/2/Rs) / (v*b0/2/Rs); % slot-opening factor
                Kdv = sin(q*v*pi/Qs) / (q*sin(v*pi/Qs)); %winding distribution factor (equation valid only if q is integer)
                Kpv = sin(v*alfa_y/2);
                Kdpv = Kdv*Kpv;
                Fv = delta*v/r*(r/Rs)^v * (1+(Rr/r)^(2*v))/(1-(Rr/Rs)^(2*v));
%                 sum_v = sum_v + Ksov*Kdpv*Fv/v * ( sin(u*p*omega_r*t+v*alfa+theta_u) + sin(u*p*omega_r*t-v*alfa+theta_u) );
                sum_v = sum_v + Ksov*Kdpv*Fv/v * ( sin(u*p*omega_r*t+v*alfa(alfa_ind)+theta_u) );
            end
        end
        B_wind = B_wind + I * sum_v;
    end
    B_wind = B_wind * mi0*3*W/pi/delta;
    data(alfa_ind,7) = B_wind;
end

figure('Name','Windings field','NumberTitle','off');
plot(data(:,1),data(:,7),'-x','MarkerSize',3);
set(gca,'XTick', data(1,1):30:data(end,1)+1);
set(gca,'XGrid','on')
set(get(gca,'XLabel'),'String','\alpha: angular position [°]');
set(get(gca,'YLabel'),'String','B: magnetic field [T]');
title('Magnetic field by windings','FontSize',13)
linex = get(gca,'XLim'); line(linex,[0 0],'LineStyle',':','Color','k');
