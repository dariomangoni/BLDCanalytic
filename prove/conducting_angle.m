function fun_theta1 = conducting_angle(theta1_hyp, C1, C3, C21, C22, C23, emf)

%% costruzione della funzione da passare all'integratore int((ea-eb)*exp(C1*th))
cont = 0;
for th_ind = 1:length(emf)
    if (emf(th_ind,1)>=(-pi/6))&&(emf(th_ind,1)<=(pi/6))
        cont = cont +1;
        fun_integranda1(cont,:) = [emf(th_ind,1),(emf(th_ind,2)-emf(th_ind,3))*exp(C1*emf(th_ind,1))];
    end
end
if  ~exist('fun_integranda1','var') || length(fun_integranda1) ==1
    disp('Integrale zero')
    integral1 = 0;
else
    integral1 = simpson_mod(fun_integranda1);
end


%% costruzione della funzione da passare all'integratore int((ea-eb)*exp(C1*th))
cont = 0;
for th_ind = 1:length(emf)
    if (emf(th_ind,1)>=(-pi/6))&&(emf(th_ind,1)<=theta1_hyp)
        cont = cont+1;
        fun_integranda2(cont,:) = [emf(th_ind,1),(2*emf(th_ind,4)-emf(th_ind,3)-emf(th_ind,2))*exp(C1*emf(th_ind,1))];
    end
end
if ~exist('fun_integranda2','var') || length(fun_integranda2) ==1
    integral2 = 0;
    disp('Integrale zero')
else
    integral2 = simpson_mod(fun_integranda2);
end


fun_theta1 = 1/C1*(C21+C22*exp(pi*C1/3))*(exp(C1*theta1_hyp)-exp(-C1*pi/6))+...
            C23/C1*(exp(C1*pi/6)-exp(C1*theta1_hyp))-...
            C3/2*integral1-...
            C3/6*(2*exp(C1*pi/3)-1)*integral2;
end