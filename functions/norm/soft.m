% define the soft threshold function, which is used above.
function y = soft(x,tau,Par)

y = sign(x).*max(abs(x)-tau,0);
