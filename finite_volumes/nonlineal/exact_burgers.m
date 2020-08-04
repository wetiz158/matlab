function [w] = exact_burgers(x,t,m,wl,wr)
%
% Exact Solution of Burgers Equation
%
 %w=zeros(1,m+1);
%   for i=1:m+1
%       if (x(i)<=t) && (t<1)
%           w(i) = wl ;
%       elseif (t<=x(i)) && (x(i)<=1)
%           w(i)=(1-x(i))/(1-t);
%       elseif (x(i)>1) && (t<=1)
%           w(i) = wr ;
%       end


if t<1
    w=wl.*(x<=t)+(1-x)/(1-t).*(t<=x).*(x<=1)+wr.*(x>1);
else
    w=wl.*(x<=1)+wr.*(x>1);
end
     
  
 
  
end

