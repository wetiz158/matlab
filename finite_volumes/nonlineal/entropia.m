function [went]=entropia(wngod,t,x,m) % pasarlle w, t, x
%%%%%%%%%%%%%%%%%%%%%
fent=zeros(1,m+1);
%fentn=zeros(1,m+1);
went=zeros(1,m+1);

    for i=1:m
        
        if x<=t && t<1 % cond1
          fent(i)=1  ;
            went(i)=1;

        elseif t<x && x<=1 % cond2
              fent(i)=((1-x)/(1-t))^2  ;
            went(i)=(2/3)*((1-x)/(1-t))^3;
       %elseif x>1 && t<=1 % cond3
           
        end
    end
    %funcion
    went
    
    %flujo
    fent(2:m+1)=fw(1:m);

fwn(1)=fw(1);
%
went(1:m+1)=wa(1:m+1)-dtdx*(fw(1:m+1)-fwn(1:m+1));
