% gruppo_sin
clear all
close all
clc

cont2 = 0;
for cont = 1:50
    if ( rem(cont,2) && rem(cont,3) ) % check if it's not multiple of 2 or 3
        cont2 = cont2+1;
        u_array(cont2,1) = cont;          % then add it in the array
    end
end

p=2;
v=1:50;
for v_ind=1:length(v);
    for u_ind =1:length(u_array)
        array((v-1)*length(u_array)+u_ind,1) = v;
        array((v-1)*length(u_array)+u_ind,2) = u_array(u_ind);
        array((v-1)*length(u_array)+u_ind,3) = rem(v(v_ind)/p+u_array(u_ind),3);
        array((v-1)*length(u_array)+u_ind,4) = rem(v(v_ind)/p-u_array(u_ind),3);
    end
end


