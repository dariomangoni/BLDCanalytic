function sol=simpson(fun,a,b,N)
h=(b-a)/N;
x=[a:h/2:b];
fx=feval(fun,x);
sol=h/6*(fx(1)+4*sum(fx(2:2:end-1))+2*sum(fx(3:2:end-2))+fx(end));
return