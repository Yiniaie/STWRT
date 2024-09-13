% define the soft threshold function, which is used above.
function y = Log_col(x,lambda,normPar)
y=zeros(size(x));
n=size(x,2);
for i=1:n
    x_norm= Log(norm(x(:,i),'fro'),lambda,normPar);
y(:,i)= x_norm/(norm(x(:,i),'fro')+eps)*x(:,i);
end