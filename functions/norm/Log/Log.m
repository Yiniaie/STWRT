% define the hard threshold function, which is used above.
function y =Log(x,tau,lambda)
x0=sqrt(2*lambda)-tau;
y=sign(max(x-x0,0)).*(1/2).*((x-tau)+sqrt((x+tau).^2-2*lambda))+sign(max(-x0-x,0)).*(1/2).*((x+tau)+sqrt((x-tau).^2-2*lambda));
%y = sign(x).*max(abs(x)-tau,0);
