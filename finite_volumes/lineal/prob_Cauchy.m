clear all
close all

%%Resolucion del problema de Cauchy para las ecuaciones de la acustica
%%lineal con condiciones iniciales  u0, p0

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
u0=@(x) exp(-2*x.^2);
p0=@(x) exp(-2*x.^2);

figure(1)
subplot(2,1,1)
plot(x, u0(x),'b')
xlabel('x'); ylabel('u0(x)');
title('Velocidad inicial del medio perturbado')
subplot(2,1,2)
plot(x, p0(x),'g')
xlabel('x'); ylabel('p0(x)');
title('Presion inicial en el medio perturbado')


% Solucion exacta del problema de Cauchy
[u,p]=sol_ex_C(rho,c0,u0,p0);
umax=max(u0(x));



%Inicializamos para utilizar el metodo de Godunov
w0=zeros(2,length(x));
wnm1=w0;
f=zeros(2,length(x)-1);

% Instante incial
w0(1,:)=u0(x);
w0(2,:)=p0(x);
wn=w0;

% Calculo de la matriz |A|
P=[1 -1; rho*c0 rho*c0];
A=[0 1/rho; rho*c0^2 0];
absA=P*diag([c0 c0])*P^(-1);

%Metodo de Godunov para un sistema lineal
for i=1:length(t)
    
    %Flujo numerico para el metodo de Godunov de un sistema lineal
    f(:,1:end)=0.5*(A*wn(:,1:end-1)+A*wn(:,2:end))-0.5*absA*(wn(:,2:end)-wn(:,1:end-1));
    
    %Calculo en el instante de tiempo siguiente
    wnm1(:,2:end-1)=wn(:,2:end-1)-dt/dx*(f(:,2:end)-f(:,1:end-1));
    
    %Condiciones de contorno
    wnm1(:,1)=wn(:,2);
    wnm1(:,end)=wn(:,end-1);
    
    %Representacion grafica 
    figure(2);
    tiledlayout(2, 1);
    nexttile
    plot(x, u(x,t(i)), 'ro'); % exacta
    hold on
    plot(x, wn(1,:), 'b'); % godunov
    legend('exacta','godunov')
    xlabel('x'); ylabel('u(x,t)');axis([a b 0 umax])
    title({'Velocidad del medio perturbado';['t = ', num2str(t(i))]});
    grid on;
    nexttile
    plot(x, p(x,t(i)), 'ro'); % exacta
    hold on
    plot(x, wn(2,:), 'b'); % godunov
    legend('exacta','godunov')
    xlabel('x'); ylabel('p(x,t)');%axis([a b -1e6 1e6])
    title({'Presión en el medio perturbado';['t = ', num2str(t(i))]});
    grid on;
    pause(0.1);
    
    % Actualizacion
    
    wn=wnm1;
end
   
    
    