% PhaseCurrent2
clear all
close all
clc
load('ea.mat') %nella versione definiva sarà calcolato all'interno di Analytical Model


%% input arguments
Kem = 0.0987;
% Kem=max(ea);
p = 2;
R = 2.7;
L = 5.63e-3;
M = -2.15e-3;
Vps = 45;
omega_r = 800/60*2*pi;
% omega_m = omega_r/p;
omega_m = 800/60*2*pi;
DVt1 = 0.163;
DVt2 = 0.7;

%% Back-EMF
th = p*omega_r*t; %ottengo la forma d'onda emf in funzione della posizione angolare e non del tempo (corretto?)
emf = [ th', ea./max(ea), ea./max(ea), ea./max(ea) ]; %potrebbe essere il caso che non ci interessi proprio tutta la forma d'onda...

for th_ind = 1:length(emf)
    if emf(th_ind,1)>=pi/2
        shift = th_ind;
        break
    end
end
emf(:,3) = circshift(emf(:,3),[-shift,0]);
emf(:,4) = circshift(emf(:,4),[+shift,0]);
                                                                   

figure('Name','Back-Emf','NumberTitle','off');
plot(emf(:,1)/pi*180,emf(:,2),'b','LineWidth',2); hold on;
plot(emf(:,1)/pi*180,emf(:,3),'r');
plot(emf(:,1)/pi*180,emf(:,4),'k');
set(gca,'XTick', -30:30:330);
set(gca,'XGrid','on')
legend('Phase A','Phase B','Phase C');
set(get(gca,'XLabel'),'String','\theta: angular position [°]');
set(get(gca,'YLabel'),'String','e: ddp [V]');
title('Back-emf in phase A-B-C from -\pi/6 to +\pi/6','FontSize',13)
linex = get(gca,'XLim'); line(linex,[0 0],'Color','k');


C1 = R/(p*omega_m*(L-M));
C21 = (2*Vps -2*DVt1   +DVt2)/(3*p*omega_m*(L-M));
C22 = ( -Vps   +DVt1 -2*DVt2)/(3*p*omega_m*(L-M));
C23 = (  Vps   -DVt1        )/(2*p*omega_m*(L-M));
C3 = Kem/p/(L-M);

%% Theta1: conducting period
syms theta1_hyp;
flag = 1;
a = -pi/6;
b = +pi/6;
maxerr = 0.0000001;

f = @(theta1_hyp) conducting_angle(theta1_hyp, C1, C3, C21, C22, C23, emf);   % f: function
c = (a*f(b) - b*f(a))/(f(b) - f(a));
% disp('   Xn-1      f(Xn-1)      Xn      f(Xn)      Xn+1      f(Xn+1)');
% disp([a f(a) b f(b) c f(c)]);
while abs(f(c)) > maxerr
    a = b;
    b = c;
    c = (a*f(b) - b*f(a))/(f(b) - f(a));
%     disp([a f(a) b f(b) c f(c)]);
    flag = flag + 1;
    if(flag == 100)
        break;
    end
end

theta1 = c;
if theta1>pi/6
    theta1=pi/6;
    disp('Diode always conducting')
end
if theta1<-pi/6
    theta1=pi/6;
    disp('Soluzione trovata non corretta')
end
clear flag a b c f theta1_hyp maxerr
Dtheta = emf(2,1)-emf(1,1);

%% Phase A current
ia_history = [-pi/6,0];
%% calculating ia1
fun_integranda_old = (C21 - C3/3*(2*emf(1,2)-emf(1,3)-emf(1,4)))*exp(C1*emf(1,1));
for th_ind = 2:length(emf(:,1))
    if (emf(th_ind,1)>theta1) || (emf(th_ind,1)>pi/6)
        theta1_ind = th_ind-1;
        break
    end
    
    fun_integranda_new = (C21 - C3/3*(2*emf(th_ind,2)-emf(th_ind,3)-emf(th_ind,4)))*exp(C1*emf(th_ind,1));
    integral1 = (fun_integranda_new + fun_integranda_old)*Dtheta/2;
    ia_history(th_ind,:) = [emf(th_ind,1), exp(-C1*Dtheta)*( ia_history(th_ind-1,2) + exp(-C1*emf(th_ind-1,1))*integral1)];
    fun_integranda_old = fun_integranda_new;
end


%% calculating ia2
fun_integranda_old = (C23 - C3/2*(emf(th_ind-1,2)-emf(th_ind-1,3)))*exp(C1*emf(th_ind-1,1));
for th_ind = th_ind:length(emf(:,1))
    if emf(th_ind,1)>pi/6
        break
    end
    
    fun_integranda_new = (C23 - C3/2*(emf(th_ind,2)-emf(th_ind,3)))*exp(C1*emf(th_ind,1));
    integral1 = (fun_integranda_new + fun_integranda_old)*Dtheta/2;
    ia_history(th_ind,:) = [emf(th_ind,1), exp(-C1*Dtheta)*( ia_history(th_ind-1,2) + exp(-C1*emf(th_ind-1,1))*integral1)];
    fun_integranda_old = fun_integranda_new;
end

%% Phase C current
for th_ind = 1:length(emf(:,1))
    if emf(th_ind,1)>theta1
        break
    end
    Qc1 = C22 - C3/3*(2*emf(th_ind,4)-emf(th_ind,2)-emf(th_ind,3));
    fun_integranda1(th_ind,:) = [emf(th_ind,1), -Qc1*exp(C1*emf(th_ind,1))];
end
if ~exist('fun_integranda1','var')
    integral1 = 0;
    disp('Integrale zero')
else
    integral1 = simpson_mod(fun_integranda1); %Cc1
end
ic_history = [-pi/6, integral1*exp(C1*pi/6)];
% ic_history = [-pi/6, ia_history(end,2)];

%% calculating ic1
fun_integranda_old = (C22 - C3/3*(2*emf(th_ind,4)-emf(th_ind,2)-emf(th_ind,3)))*exp(C1*emf(1,1)); %7.23b
for th_ind = 2:length(emf(:,1))
    if (emf(th_ind,1)>theta1) || (emf(th_ind,1)>pi/6)
        break
    end
    
    fun_integranda_new = (C22 - C3/3*(2*emf(th_ind,4)-emf(th_ind,2)-emf(th_ind,3)))*exp(C1*emf(th_ind,1));
    integral1 = (fun_integranda_new + fun_integranda_old)*Dtheta/2;
    ic_history(th_ind,:) = [emf(th_ind,1), exp(-C1*Dtheta)*( ic_history(th_ind-1,2) + exp(-C1*emf(th_ind-1,1))*integral1)]; %7.23a
    fun_integranda_old = fun_integranda_new;
end

%% calculating ic2
for th_ind = th_ind:length(emf(:,1))
    if emf(th_ind,1)>pi/6
        break
    end
    
    ic_history(th_ind,:) = [emf(th_ind,1), 0]; %7.24
end

%% calculatin ib
ib_history = zeros(length(ia_history),2);
for th_ind = 1:length(ia_history)
    if ia_history(th_ind,1)<=theta1
        ib_history(th_ind,:) = [ia_history(th_ind,1), -ia_history(th_ind,2)-ic_history(th_ind,2)]; %7.12b
    elseif ia_history(th_ind,1)<=pi/6
        ib_history(th_ind,:) = [ia_history(th_ind,1), -ia_history(th_ind,2)]; %7.13b
    else
        break
    end
end


figure('Name','Current components','NumberTitle','off');
plot(ia_history(:,1)*180/pi,ia_history(:,2),'b'); hold on
plot(ib_history(:,1)*180/pi,ib_history(:,2),'r');
plot(ic_history(:,1)*180/pi,ic_history(:,2),'k');
set(gca,'XTick', sort([-30:30:330,theta1*180/pi]));
set(gca,'XGrid','on')
set(get(gca,'XLabel'),'String','\theta: angular position [°]');
set(get(gca,'YLabel'),'String','i: current [A]');
title('Currents in phase A-B-C from -\pi/6 to +\pi/6','FontSize',13)
linex = get(gca,'XLim'); line(linex,[0 0],'Color','k');
liney = get(gca,'YLim'); line([theta1,theta1]*180/pi,liney,'LineStyle','--','Color','k','LineWidth',2);

%% phase current
Ia = [ia_history(:,2);+ic_history(1:theta1_ind,2)+ia_history(1:theta1_ind,2);ia_history(theta1_ind:end,2);ic_history(:,2)];
I = [(-pi/6:pi/length(Ia):11*pi/6-0.0000001)',[Ia;-Ia],[Ia;-Ia],[Ia;-Ia]];
I(:,3) = circshift(I(:,3),[-round(size(I,1)/3),0]);
I(:,4) = circshift(I(:,4),[+round(size(I,1)/3),0]);

figure('Name','Currents','NumberTitle','off'); hold on
plot(I(:,1)*180/pi,I(:,2),'x','MarkerSize',3)
plot(I(:,1)*180/pi,I(:,3),'x','MarkerSize',2,'Color','r')
plot(I(:,1)*180/pi,I(:,4),'x','MarkerSize',2,'Color','k')
set(gca,'XGrid','on')
set(gca,'XTick', sort(unique([-30,30,90,theta1*180/pi,60+theta1*180/pi,120+theta1*180/pi,150,210,270,180+theta1*180/pi,180+60+theta1*180/pi,180+120+theta1*180/pi,180+150])));
set(gca,'YGrid','on')
% set(gca,'YTick', sort(unique([ia_history(1,2),ia_history(theta1_ind,2),ia_history(end,2),ia_history(1,2)+ic_history(1,2),...
%                         ia_history(theta1_ind,2)+ic_history(theta1_ind,2),ic_history(1,2),ic_history(end,2)])));
set(gca,'YTick', sort(unique([ia_history(1,2),ia_history(theta1_ind,2),ic_history(1,2),ic_history(end,2)])));
set(get(gca,'XLabel'),'String','\theta: angular position [°]');
set(get(gca,'YLabel'),'String','i: current [A]');
title('Currents in phase A -\pi/6 to +11/6\pi','FontSize',13)
linex = get(gca,'XLim'); line(linex,[0 0],'Color','k');
liney = get(gca,'YLim'); line([theta1,theta1]*180/pi,liney,'LineStyle','--','Color','k','LineWidth',2);


%% overall current
points = -pi/6:Dtheta:pi/6;
T = zeros(size(points,1),2);
Tave = 0;
Iave = 0;
for points_ind = 1:length(points)
    ia = interp1(I(:,1),I(:,2),points(points_ind));
    ea = interp1(emf(:,1),emf(:,2),points(points_ind));
    ib = interp1(I(:,1),I(:,3),points(points_ind));
    eb = interp1(emf(:,1),emf(:,3),points(points_ind));
    ic = interp1(I(:,1),I(:,4),points(points_ind));
    ec = interp1(emf(:,1),emf(:,4),points(points_ind));
    T(points_ind,:) = [points(points_ind), Kem*(ea*ia+eb*ib+ec*ic)];
    Tave = Tave + T(points_ind,2);
    Iave = Iave + ia;
end
Tave = Tave/length(points);
Iave = Iave/length(points);
eta_system = Tave*omega_m/Vps/Iave;
% eta_system = Tave*omega_m/Vwt/Iave;


line(linex,[Iave Iave],'LineStyle','--','Color','b','LineWidth',1.5);

figure('Name','Torque','NumberTitle','off'); hold on
plot(T(:,1)*180/pi,T(:,2),'x','MarkerSize',3)
set(gca,'XGrid','on')
set(gca,'XTick', -30:5:30);
set(get(gca,'XLabel'),'String','\theta: angular position [°]');
set(get(gca,'YLabel'),'String','T: torque [Nm]');
title('Torque -\pi/6 to +\pi/6','FontSize',13)
linex = get(gca,'XLim'); line(linex,[0 0],'Color','k');
line(linex,[Tave Tave],'Color','k','LineWidth',2);

clear fun_integranda1 fun_integranda_new fun_integranda_old integral1 linex liney shift th t ea th_ind eb ec ia ib ic points points_ind




