function BEGQC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Finite volume methods for Burgers equation
%
%       w + w w =0
%        t         x
%
%   Domain: [a,b]
%
%   Riemann Problem
%
%   w0(x)=wl if x<0; w0(x)=1-x if 0<=x<1  w0(x)=wr if x>1
%
%   Transmissive boundary conditions ??
%   Deltat is obtained from CFL number
%
%   Numerical methods:
%     * Godunov method
%     * Q-scheme van Leer / Roe
%     * Centered
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear all
 clc 
 disp('---------------------------------------------')
 disp('Finite volume methods for Burgers equation')
 disp('---------------------------------------------')
 disp('w + w w =0')
 disp(' t     x')
 disp('Domain: [a,b]')
 disp('Transmissive boundary conditions')
 disp('Riemann problem for Burgers equation')
 disp('w0(x)=wl if x<0; w0(x)=1-x if 0<=x<1  w0(x)=wr if x>1')
 disp('Exact solution')
 disp('---------------------------------------------')
 a=-3;
 disp(['Lower end of the interval a = ',num2str(a)]) 
 b=3;
 disp(['Upper end of the interval b = ',num2str(b)]) 
 m=200;
 disp(['Number of nodes m=', num2str(m)]) 
 deltax=(b-a)/(m);
 
 disp(['deltax =', num2str(deltax)])
 x=[a:deltax:b];
 cfl=0.9;
 disp(['Courant number = ', num2str(cfl)])
%
 wl=1;
 %wm=1-x;
 wr=0;
%
% Deltat from Courant number
 deltat=cfl*deltax/max(abs(wl),abs(wr));
 disp(['Time step from Courant number = ',...
       num2str(deltat)])
 tmax=5;
 disp(['Time end = ', num2str(tmax)])
 mt=tmax/deltat;
disp('---------------------------------------------')
disp(['Initial condition:'])
disp(['w0(x)=',num2str(wl),' if x<0'])
disp(['w0(x)=1-x',' if 0<=x<1'])
disp(['w0(x)=',num2str(wr),' if x>1'])
%
disp('---------------------------------------------')
%
if wl>wr
        s=(wl+wr)/2;
    disp(['The speed of the shock from'])
    disp(['the Rankine Hugoniot condition s=',...
           num2str(s)])
    disp(['The exact solution is a entropy shock']) 
    disp(['w(x,t)=',num2str(wl),...
         ' if (x)/t<=',num2str(s)])
    disp(['w(x,t)=',num2str(wr),...
         ' if (x)/t>',num2str(s)])
elseif wl<wr
        s=(wl+wr)/2;
    disp(['The exact solution is a rarefaction wave'])
    disp(['w(x,t)=',num2str(wl),...
         ' if (x)/t<=',num2str(wl)])
    disp(['w(x,t)= (x)/t',...
         ' if ',num2str(wl),'<= (x)/t<=',num2str(wr)])
    disp(['w(x,t)=',num2str(wr),...
        ' if (x)/t>',num2str(wr)])
elseif wl==wr
    disp('The exact solution is')
    disp('a constant solution w(x,t)=wl')
end
disp('---------------------------------------------')
disp(['Numerical methods:'])
disp([' * Godunov'])
disp([' * Q-scheme van Leer / Roe'])
disp([' * Centered'])
% 
disp('---------------------------------------------')
%
% Plot the initial condition at [a,b] 
%
 for i=1:m+1
     if(x(i)<0)
        w0(i)=wl;
     elseif(0<=x(i)&&x(i)<1)
         w0(i)=1-x(i);
     else 
        w0(i)=wr; 
     end
 end

 wmin=min(wl,wr)-0.2; %lower limit for the y axi
 wmax=max(wl,wr)+0.2; %upper limit for the y axi
 figure(1)
 %
 plot(x,w0,'ok')
 axis([a b wmin wmax])
 xlabel('x'); ylabel('w(x,0)');
title('Initial condition');
hold off
% 
% Inicialization
%
t=0;
%
we=w0;%Exact solution
wagod=w0;% Godunov method
waq=w0;% Q-scheme (Roe=van leer for Burgers equation) 
wacen=w0;% Centered scheme
went=w0; %entropy
%
dtdx=deltat/deltax;
%
for n=1:mt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Godunov method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wngod=god_btbc(wagod,dtdx,m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q-scheme 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wnq=qscheme_btbc(waq,dtdx,m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Centered scheme 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wncen=centered_scheme(wacen,dtdx,m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entropy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
went = entropia(wngod,t,x,m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Exact solution
%
t=t+deltat;
we=exact_burgers(x,t,m,wl,wr);
%    
%
figure(2)
    plot(x,wngod,'-x',x,wnq,'-+',x,wncen,'-*',x,we,'-ok')
    axis([a b wmin wmax])
    xlabel('x'); ylabel('w(x,t)');
    title(['Burgers equation, t = ',num2str(t)]);
    legend('Godunov','Q-scheme','Centered','Exact solution')
    pause(0.001)
% 
figure(3)
     plot(x,abs(we-wngod),'-x',...
          x,abs(we-wnq),'-+',x,abs(we-wncen),'-*')
    xlabel('x'); ylabel('|we(x,t)-w(x,t)|');
    axis([a b wmin wmax])
    title(['Burgers equation, t = ',num2str(t)]);
    legend('Error Godunov','Error Q-scheme','Error Centered')
    %
    pause(0.001)
    % graficar entropia
 figure(4)
     plot(x,abs(went),'-x',...
          x,abs(we-wnq),'-+',x,abs(we-wncen),'-*')
    xlabel('x'); ylabel('|we(x,t)-w(x,t)|');
    axis([a b wmin wmax])
    title(['entropia, t = ',num2str(t)]);
    legend('enttt')
    %
    pause(0.001)
% Update
wagod=wngod;
waq=wnq;
wacen=wncen;
end
disp('---------------------------------------------')
disp(['Error at time t=',num2str(t)])
disp(['* Godunov (m�x(abs(we-wngod)))=',...
                      num2str( max(abs(we-wngod)) )])
disp(['* Q-scheme (m�x(abs(we-wnq)))=',...
                       num2str( max(abs(we-wnq))   )])
disp(['* Centered (m�x(abs(we-wncen)))=',...
                       num2str( max(abs(we-wncen))  )])
disp('---------------------------------------------')
disp('---------------------------------------------')

disp('finish')
