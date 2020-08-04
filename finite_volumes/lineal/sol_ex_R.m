function [u,p]=sol_ex_R(rho,c0,ul,pl,ur,pr)
%Calculo de la solucion exacta del problema de Riemann para las ecuaciones
%de la acustica lineal con condiciones iniciales ul, pl, ur, pr

%Coordenadas en la base de autovectores
alphal=-ul/2+pl/(2*rho*c0);
alphar=-ur/2+pr/(2*rho*c0);

%Solucion exacta para el problema de Riemann
u=@(x,t) ul.*(x+c0*t<0)+(ul-(alphar-alphal)).*((x+c0*t).*(x-c0*t)<=0)+...
        ur.*(x-c0*t>0);
p=@(x,t) pl.*(x+c0*t<0)+(pl+(alphar-alphal)*rho*c0).*((x+c0*t).*(x-c0*t)<=0)+...
        pr.*(x-c0*t>0);
end