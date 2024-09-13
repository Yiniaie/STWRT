% 解决p=1 或者 0<p<1都可以

% 当style = ‘row’, 此函数求解优化问题
% min_A  0.5*||A-Z||_F^2 + lambda ||A||_{2,p}^p

% 当style = ‘col’, 此函数求解优化问题
% min_A  0.5*||A-Z||_F^2 + lambda ||A^T||_{2,p}^p
% 易知上述问题等价于
% min_A  0.5*||A^T-Z^T||_F^2 + lambda ||A^T||_{2,p}^p <=> min_B 0.5*||B-Z^T||_F^2 + lambda ||B||_{2,p}^p

function A = GroupSpecialF_2p(Z,lambda,p,style)
if nargin < 4
   style='col';
end
if strcmp(style,'col')
    Z = Z';
end
delta = lambda./sqrt(sum(Z.^2,2)).^(2-p);
y = ones(size(delta));

if p<1
    s_star  =  solve_Lp( y, delta, p );
else
    s_star = pos( y - delta );
end
A = repmat(s_star,1,size(Z,2)).*Z;
if strcmp(style,'col')
    A = A';
end