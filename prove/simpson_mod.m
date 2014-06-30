function sol=simpson_mod(array)
h=2*(array(end,1)-array(1,1))/(length(array)-1);
fx = array(:,2);
sol=h/6*( fx(1) + 4*sum(fx(2:2:end-1)) + 2*sum(fx(3:2:end-2)) + fx(end));
return

