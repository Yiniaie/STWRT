function [Delta]  =  Generalized_Soft_Thresholding(S, weight, p)
%p是范数中的参数
%weight是正则项系数，并不是权重
%p=Par.p;
J=5;
%sigmaX1=sqrt(max(S^2-4*c,0));
%weight=c./(sigmaX1.^(1/p)+eps);
%diagS=diag(S);
if length(weight)==1
weight=weight*ones(size(S));
end
weight=weight(:);
diagS=S(:);
Delta=diagS;
sigma0   =    abs(diagS);
tau_GST  =   (2*weight.*(1-p)).^(1/(2-p))   +   p*weight.*(2*(1-p)*weight).^((p-1)/(2-p));
%tau_GST= diag(tau_GST);
delta=[];
for i=1:length(diagS)
    
    if sigma0(i)>tau_GST(i)
        k=0;delta=sigma0(i);
        %weight=diag(weight); 
        for k=1:J 
            delta=sigma0(i)-weight(i)*p*(delta)^(p-1);           
            k=k+1;
        end
        Delta(i)   =   sign(diagS(i)).*(delta);
    else
        Delta(i)=0; 
        
    end
    
end
%Delta=diag(Delta);
Delta=reshape(Delta,size(S));