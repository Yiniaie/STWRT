% define the soft threshold function, which is used above.
function y = GAI_GST_col(x,tau,p)
%p=Par.p;
y=zeros(size(x));
n=size(x,2);
for i=1:n
    x_norm=GAI_GST(norm(x(:,i),'fro'),tau,p);
y(:,i)= x_norm/(norm(x(:,i),'fro')+eps)*x(:,i);
end
