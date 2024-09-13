function xstar = Logarithm( b1,lambda1,gamma)
if length(lambda1)==1
    lambda1=lambda1*ones(size(b1));
end
Sig=sign(b1);b=abs(b1(:));
lambda=lambda1(:);
%X 此处显示有关此函数的摘要
%gamma是范数自带参数
%lambda是惩罚系数;
n=length(b);
xstar=zeros(n,1);x0=zeros(n,1);
delta=(1+gamma.*b).^2-4*lambda*gamma^2/log(gamma+1);
index=find(delta>=0);
   b0=b(index);
   lambda0=lambda(index);
   x1=1/(2*gamma) *  (gamma.*b0-1+sqrt((1+gamma.*b0).^2-4*lambda0*gamma^2/log(gamma+1)));
   x2=1/(2*gamma) *  (gamma.*b0-1-sqrt((1+gamma.*b0).^2-4*lambda0*gamma^2/log(gamma+1)));
   
  x0(index)=max(max(x1,x2),0);
       f1=f(b,x0,gamma,lambda);
          f2=f(b,0,gamma,lambda);
   xstar(find(f1<=f2))=x0(find(f1<=f2));

    xstar=Sig.*reshape( xstar,size(b1));
   
end
 

function y = f(b,x,gamma,lambda)
%函数f
y=1/2*((x-b).^2) + lambda .* log(gamma.*x+1)./log(gamma+1);
end