function [wn]=centered_scheme(wa,dtdx,m)


wam1=zeros(1,m+1);
wap1=zeros(1,m+1);
f=zeros(1,m+1);

%
wap1(1:m)=wa(2:m+1);
% Transmissive boundary conditions 
wap1(m+1)=wa(m);
%
wam1(2:m+1)=wa(1:m);
% Transmissive boundary conditions
wam1(1)=wa(2);

%
% Flux
%
f(1:m+1)=0.5*wa.^2*dtdx.*(wap1(1:m+1)-wam1(1:m+1));
%

wn(1:m+1)=wa(1:m+1)-f(1:m+1);


end

