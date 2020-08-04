function [u,p]=sol_ex_C(rho,c0,u0,p0)
%Calculo de la solucion exacta del problema de Cauchy para las ecuaciones
%de la acustica lineal con condiciones iniciales u0, p0

%Solucion exacta 
u=@(x,t) u0(x-c0*t)/2+p0(x-c0*t)/(2*rho*c0)+u0(x+c0*t)/2-p0(x+c0*t)/(2*rho*c0);
p=@(x,t) rho*c0*(u0(x-c0*t)/2+p0(x-c0*t)/(2*rho*c0)-u0(x+c0*t)/2+p0(x+c0*t)/(2*rho*c0));
end