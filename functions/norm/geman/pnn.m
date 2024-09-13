

function X=pnn(Y,lambda,epsilon)
   X=zeros(size(Y));
%%the defination of p,q,delta,h,x1,x2,x3
if max(size(lambda,1),size(lambda,2))==1
    lambda=lambda*ones(size(Y));
end
p=-(epsilon+abs(Y)).^2/3;
q=epsilon*lambda-2/27*(epsilon+abs(Y)).^3;
delta=(q/2).^2+(p/3).^3;
indexh=find(delta<=0);
indexx1=find(delta>0);
x11=sign(-q(indexx1)/2+sqrt(delta(indexx1))).*abs(-q(indexx1)/2+sqrt(delta(indexx1))).^(1/3)+sign(-q(indexx1)/2-sqrt(delta(indexx1))).*abs(-q(indexx1)/2-sqrt(delta(indexx1))).^(1/3)-(2*epsilon-abs(Y(indexx1)))/3;
x1=real((-q(indexh)/2+sqrt(delta(indexh))).^(1/3)+(-q(indexh)/2-sqrt(delta(indexh))).^(1/3)-(2*epsilon-abs(Y(indexh)))/3);
x2=real((-1+sqrt(3)*1i)/2*(-q(indexh)/2+sqrt(delta(indexh))).^(1/3)+((-1+sqrt(3)*1i)/2).^2*(-q(indexh)/2-sqrt(delta(indexh))).^(1/3)-(2*epsilon-abs(Y(indexh)))/3);
x3=real(((-1+sqrt(3)*1i)/2).^2*(-q(indexh)/2+sqrt(delta(indexh))).^(1/3)+(-1+sqrt(3)*1i)/2*(-q(indexh)/2-sqrt(delta(indexh))).^(1/3)-(2*epsilon-abs(Y(indexh)))/3);
h=max(max(x1,x2),x3);
GX1=2*(gFun(abs(Y(indexx1)),max(x11,0),lambda(indexx1),epsilon)+lambda(indexx1));
Gh=2*(gFun(abs(Y(indexh)),max(h,0),lambda(indexh),epsilon)+lambda(indexh));



            X(indexh)=sign(max(abs(Y(indexh))-sqrt(Gh),0)).*sign(Y(indexh)).*max(h,0);
       
           X(indexx1)=sign(max(abs(Y(indexx1))-sqrt(GX1),0)).*sign(Y(indexx1)).*max(x11,0); 
 






end

function g=gFun(y,x,lambda,epsilon)
g=1/2*(abs(y)-x).^2-lambda*epsilon./(epsilon+x);
end





