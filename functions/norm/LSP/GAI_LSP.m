function [y] = GAI_LSP(x,lambda,gamma,Iter)

if nargin < 4
    Iter=10;
end
    
a0=max(sqrt(lambda)-gamma,0);
delta=a0+LSP_gradient(a0,gamma,lambda);
%===========================
tol=10.^(-4);%eps;%
  [m,n]=size(x);
 y=zeros(size(m,n));
 

 for i=1:m
     for j=1:n
         
         gradient_g= LSP_gradient(x(i,j),gamma,lambda);
         if gradient_g==0
             BARx_b=x(i,j);
             
         else
             if delta<x(i,j)
                 a=x(i,j);
                 sto=0;kkk=0;
                 while sto==0
                     kkk=kkk+1;
                     a1=x(i,j)-LSP_gradient(a,gamma,lambda);
                     a2=x(i,j)-LSP_gradient(a1,gamma,lambda);
                     % new_a=a1-(a2-a1)*(a1-a)/(a2-2*a1+a);
                     if abs(a2-2*a1+a)<tol || kkk>=Iter
                         BARx_b=a1-(a2-a1)*(a1-a)/(a2-2*a1+a);
                         break;
                     end
                     a=a1-(a2-a1)*(a1-a)/(a2-2*a1+a);  
                 end
             else
                 BARx_b=a0;
             end
             if LSP_fun(BARx_b,x(i,j),gamma,lambda)<=LSP_fun(0,x(i,j),gamma,lambda)
                 y(i,j)=BARx_b;   
             else
                 y(i,j)=0;
             end
         end
         
         
     end
 end
          






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Logarithm的gradient表达式
function gradient_g= LSP_gradient(x,gamma,lambda)
gradient_g=lambda./(x+gamma);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function g= LSP_fun(x,y,gamma,lambda)
g=1/2*(y-x).^2+lambda.*log(1+x./gamma);
return;

