%Analisis de precision del metodo de Godunov para el caso lineal 
clear all
close all

% Definicion del dominio
a = -5;            % limite inferior del intervalo
b = 5;             % limite superior del intervalo

% Mallado del tiempo
dt = 1e-6;       % subintervalo de tiempo
t_max = 0.001;     % tiempo máximo de integracion

% Datos
c0 = 343.2;         % velocidad del sonido en el agua
rho = 1.2754;         % densidad del agua
% Numero de Courant
cfl = [0.5, 1, 1.5];
%Inicializacion del vector de errores
error=zeros(size(cfl));

%Condiciones iniciales
u0=@(x) exp(-2*x.^2);
p0=@(x) exp(-2*x.^2);

% Solucion exacta del problema de Cauchy
[u,p]=sol_ex_C(rho,c0,u0,p0);

% Calculo de la matriz |A|
P=[1 -1; rho*c0 rho*c0];
A=[0 1/rho; rho*c0^2 0];
absA=P*diag([c0 c0])*P^(-1);


for j=1:length(cfl)
    % Calculo del mallado del dominio
    t = 0:dt:t_max;
    dx(j)=c0*dt/cfl(j);
    x = a:dx(j):b;
    dtdx=dt/dx(j);
    
    %Inicializacion
    w0=zeros(2,length(x));
    wnm1=w0;
    f=zeros(2,length(x)-1); 

    % Instante inicial
    w0(1,:)=u0(x);
    w0(2,:)=p0(x);

    wn=w0;
    %Metodo de Godunov para un sistema lineal
    for i=1:length(t)
        % Calculo del flujo numerico
        f(:,1:end)=0.5*(A*wn(:,1:end-1)+A*wn(:,2:end))-0.5*absA*(wn(:,2:end)-wn(:,1:end-1));
        
        % Calculo del instante siguiente de tiempo
        wnm1(:,2:end-1)=wn(:,2:end-1)-dtdx*(f(:,2:end)-f(:,1:end-1));
        %Condiciones de contorno
        wnm1(:,1)=wn(:,2);
        wnm1(:,end)=wn(:,end-1);  
    
        %Repesentacion grafica
        figure(2);
        tiledlayout(2, 1);
        nexttile
        umax=max(w0(1,:));
        plot(x, wn(1,:), 'b');
        xlabel('x'); ylabel('u(x,t)');axis([a b 0 umax])
        title({'Velocidad del medio perturbado';['t = ', num2str(t(i))]});
        grid on;
        nexttile
        pmax=max(w0(2,:));
        plot(x, wn(2,:), 'g');
        xlabel('x'); ylabel('p(x,t)');%axis([a b -1e6 1e6])
        title({'Presión en el medio perturbado';['t = ', num2str(t(i))]});
        grid on;
        pause(0.1);
        % Actualizacion
        wn=wnm1;
    end
    %Calculo del error en norma L-Infinito
    error(j)=max(max(abs(wnm1-[u(x,t(i));p(x,t(i))])))
end

    
    
