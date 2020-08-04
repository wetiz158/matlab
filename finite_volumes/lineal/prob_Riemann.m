clear all
close all

%%Resolucion del problema de Riemann para las ecuaciones de la acustica
%%lineal con condiciones iniciales  ul. pl, ur, pr

% Definicion del dominio
a = -5;            % limite inferior del intervalo
b = 5;             % limite superior del intervalo

% Mallado del tiempo
dt = 0.00005;      % subintervalos de tiempo
t_max = 0.005;     % tiempo máximo de integracion

%Datos del problema
c0 = 343.2;         % velocidad del sonido en el agua
rho = 1.2754;         % densidad del agua

% Numero de Courant
cfl=1;

% Discretizacion a partir del numero de Courant
dx=c0*dt/cfl;
x = a:dx:b;
t = 0:dt:t_max;



%Condiciones iniciales
ul = 0;
ur = 0;
pl = 1e7;
pr = 1e5;

% Solucion exacta del problema de Riemann
[u,p]=sol_ex_R(rho,c0,ul,pl,ur,pr);



%Inicializacion
w0=zeros(2,length(x));
wnm1=w0;
f=zeros(2,length(x)-1);

%Instante inicial
for i=1:length(x)
    if x(i)<=0
        w0(:,i)=[ul;pl];
    else
        w0(:,i)=[ur;pr];
    end
end

wn=w0;

%Calculo  de la matriz valor absoluto
P=[1 -1; rho*c0 rho*c0];
A=[0 1/rho; rho*c0^2 0];
absA=P*diag([c0 c0])*P^(-1);

%Metodo de Godunov para un sistema lineal
for i=1:length(t)
    %Fujo numerico
    f(:,1:end)=0.5*(A*wn(:,1:end-1)+A*wn(:,2:end))-0.5*absA*(wn(:,2:end)-wn(:,1:end-1));
    
    % Calculo del instante de tiempo siguiente
    wnm1(:,2:end-1)=wn(:,2:end-1)-dt/dx*(f(:,2:end)-f(:,1:end-1));
    
    %Condiciones de contorno
    wnm1(:,1)=wn(:,2);
    wnm1(:,end)=wn(:,end-1);
    
    %Representacion grafica
    figure(2);
    tiledlayout(2, 1);
    nexttile
    plot(x, u(x,t(i)), 'ro');
    hold on
    plot(x, wn(1,:), 'b');
    xlabel('x'); ylabel('u(x,t)');
    title({'Velocidad del medio perturbado';['t = ', num2str(t(i))]});
    grid on;
    nexttile
    plot(x, p(x,t(i)), 'ro');
    hold on
    plot(x, wn(2,:), 'g');
    xlabel('x'); ylabel('p(x,t)');
    title({'Presión en el medio perturbado';['t = ', num2str(t(i))]});
    grid on;
    pause(0.5);
    
    %Actualizacion
    wn=wnm1;
end
    
    