function y= GAI_GST(x,lambda,p)
Iter=10;
a0=max((p*(1-p))^(1/(2-p)),0);
delta=a0+GST_gradient(a0,lambda,p);
%===========================
tol=10^(-4);
[m,n]=size(x);
y=zeros(size(m,n));
 
for i=1:m
    for j=1:n
        gradient_g= GST_gradient(x(i,j),lambda,p); 
        if gradient_g==0
            BARx_b=x(i,j);
        else
            if delta<x(i,j)
                a=x(i,j);
                sto=0;sss=0;
                while sto==0
                   sss=sss+1;
                    a1=x(i,j)-GST_gradient(a,lambda,p);
                    a2=x(i,j)-GST_gradient(a1,lambda,p);
                    % new_a=a1-(a2-a1)*(a1-a)/(a2-2*a1+a);
                    if abs(a2-2*a1+a)<tol || sss>=Iter
                        BARx_b=a1-(a2-a1)*(a1-a)/(a2-2*a1+a);
                        break;
                    end
                    a=a1-(a2-a1)*(a1-a)/(a2-2*a1+a);
                end
            else
                BARx_b=a0;
            end
            if GST_fun(BARx_b,x(i,j),lambda,p)<=GST_fun(0,x(i,j),lambda,p)
                y(i,j)=BARx_b;
            else
                y(i,j)=0;
            end
        end
        
    end
end

 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Lp的gradient表达式
function gradient_g= GST_gradient(x,lambda,p)
gradient_g=lambda.*p.*x.^(p-1);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function g= GST_fun(x,y,lambda,p)
g=1/2*(y-x).^2+lambda.*(abs(x).^p);
return;
