 function [tenB, tenT,tenW] = STWRT(tenD, opts)
% papers parameters batch=(100,30)lambda_L1L2Kappa=(0.5,100,1000)delta=(0.0001)lambda_tv=(0.1)rho=(1.1)mu=(1)T=20231114T000445
%%  parameters 

%% ours-mat
[m1, m2, m3] = size(tenD);
L1=0.5;
L2=100;
Kappa=1000;
lambda_tv =0.1;
delta= 0.0001;
nu  = 1;
mu = [nu, nu, nu];
eta = mu;

if isfield(opts, 'L1');              L1   = opts.L1;                        end
if isfield(opts, 'frame_num');       frame_num = opts.frame_num;            end
if isfield(opts, 'lambda_tv');       lambda_tv    = opts.lambda_tv;         end
if isfield(opts, 'delta');           delta    = opts.delta;                 end

alpha = [1, 1, delta]./(2+delta);%[1/3 1/3 1/3];%[2 2 0.0001]./4.0001;
beta = alpha;
lambda1 = L1*Kappa/(min(m1,m2)*frame_num);
lambda2 = L2*lambda1;
Punish_P=1;
Punish =1;
mu_bar = 1e7;
rho = 1.1;  %1.1  1.05    
tol = 1e-7;
maxIter = 100;
count=0;
count_T1 =0;
epsilon_num=0;


tenB = [];
tenT = [];
tenW = [];
%%  save parameters
opts.L1 = L1;
opts.L2 = L2;
opts.Kappa = Kappa;
opts.lambda1 = lambda1;
opts.lambda2 = lambda2;
opts.delta= delta;
opts.alpha = alpha;
opts.beta= beta;
opts.nu= nu;
opts.mu= mu;
opts.eta= eta;
opts.lambda_tv= lambda_tv;
opts.Punish_P= Punish_P;
opts.Punish= Punish;
opts.mu_bar= mu_bar;
opts.rho= rho;
opts.tol= tol;
opts.maxIter= maxIter;
opts.epsilon_total=[];
%%  adding 3DTV

Nway = [m1, m2,frame_num];

D1=zeros(Nway(1),Nway(1),Nway(3));
a=diag(-ones(Nway(1),1));
b=diag(ones(Nway(1)-1,1),1);
tempD=a+b;
tempD(end,1)=1;
D1(:,:,1) = tempD;

D2=zeros(Nway(2),Nway(2),Nway(3));
a=diag(-ones(Nway(2),1));
b=diag(ones(Nway(2)-1,1),-1);
tempD=a+b;
tempD(1,end)=1;
D2(:,:,1) = tempD;


D3=zeros(Nway(1),Nway(1),Nway(3));
D3(:,:,1) = -eye(Nway(1));
D3(:,:,end) = eye(Nway(1));

I=zeros(Nway(1),Nway(1),Nway(3));
I(:,:,1)=eye(Nway(1));

Fa = (1/sqrt(Nway(1)))*fft(eye(Nway(1)));
Fb = (1/sqrt(Nway(2)))*fft(eye(Nway(2)));

if rem(m3,frame_num)
    rem_frame = rem(m3,frame_num);
    tenD =cat(3,tenD,tenD(:,:,end-(frame_num-rem_frame)+1:end));
end

batch=ceil(m3/frame_num);

for b =1:batch
%     disp(opts.dataset_name+" batch:"+batch+'('+b+')');
    priorWeight=[];
    PW1=ones([m1, m2,frame_num]);
    epsilon=[];
    epsilon_num = 0;
    count = 0;
    B=[];
    T=[];
    D=[];
    nu  = opts.nu;
    mu = [nu, nu, nu];
    eta = mu;
    for b1=1:frame_num
        I_read0=tenD(:,:,frame_num*(b-1)+b1);  
        D(:,:,b1) = I_read0;
    end 
    %%  local prior structure
    ST=Prior_ST;
  
    priorWeight0 = ST.D3(D);
    if 1  %是否对结构张量进行高斯平滑
        priorWeight = preprocess_imdilate(priorWeight0);
        opts.Preprocessing_priorWeight_imdilate = 1;
    else
        opts.Preprocessing_priorWeight_imdilate = 0;
    end
    %%  初始化
    
    B = D; 
%     B = rand(size(D));
    if D==B
        opts.initialization_B ='D=B';
    elseif B==0
        opts.initialization_B ='B=0'; 
    else
        opts.initialization_B ='B=rand(size(B))'; 
    end

    weightTen = ones(size(D));

    O_1 = zeros(size(D));
    O_2 = zeros(size(D));
    O_3 = zeros(size(D));

    Q_1 = zeros(size(D));
    Q_2 = zeros(size(D));
    Q_3 = zeros(size(D));

    Y_1 = zeros(size(D));
    Y_2 = zeros(size(D));
    Y_3 = zeros(size(D));
    P =  zeros(size(D));
    T =  priorWeight0;
    N = zeros(size(D));
    iter = 0;
    converged = false;
    %%   updata X B Y T N O Q P
    while ~converged  && iter < maxIter
    iter = iter + 1;
    % update X

    temp = permute(B,[2 3 1]) -permute(O_1,[2 3 1])/mu(1);
    [X_1]=ipermute(ProTlogSum(temp,alpha(1)/mu(1),Punish_P),[2 3 1]);
    temp =permute(B,[3 1 2])-permute(O_2,[3 1 2])/mu(2);
    [X_2]=ipermute(ProTlogSum(temp,alpha(2)/mu(2),Punish_P),[3 1 2]);
    temp =  B-O_3/mu(3);
    [X_3] =ProTlogSum(temp,alpha(3)/mu(3),Punish_P);

    %      [X]=ProTlogSum(B-O/mu,alpha(1)/mu,Punish_P);
    %      [X,tnn_1,trank_1]=prox_tnn_Gfun(B-O/mu,alpha(1)/mu,fun,Punish_P);
    %       [X] = prox_tensorLaplace(B-O/mu, 1/mu, X,epsilon);
    % update B

    tempt1 = eta(1)*tprod(tran(D1),Y_1+Q_1./eta(1));
    tempt2 = eta(2)*tprod(Y_2+Q_2./eta(2),tran(D2));
    tempt3 = eta(3)*tprod(tran(D3),Y_3+Q_3./eta(3));
    tempt4 = mu(1)*(X_1+O_1./mu(1));
    tempt5 = mu(2)*(X_2+O_2./mu(2));
    tempt6 = mu(3)*(X_3+O_3./mu(3));
    tempt7 = nu*(D-T-N+P./nu);
    tempK3 = tempt1+tempt2+tempt3+tempt4+tempt5+tempt6+tempt7;

    tempK1 = nu*I+eta(1)*tprod(tran(D1),D1)+eta(3)*tprod(tran(D3),D3)+mu(1)*I+mu(2)*I+mu(3)*I;
    tempK2 = eta(2)*tprod(D2,tran(D2));

    tempK1f=fft(tempK1,[],3);
    tempK2f=fft(tempK2,[],3);
    tempK3f=fft(tempK3,[],3);

    Bf = zeros(Nway);
    for i=1:Nway(3)
        Ai=tempK1f(:,:,i);
        Bi=tempK2f(:,:,i);
        Ci=tempK3f(:,:,i);
        da=Ai(:,1); deigA=fft2(da);
        db=Bi(:,1); deigB=fft(db);        
        Sig=repmat(deigA,1,Nway(2))+repmat(deigB',Nway(1),1);       
        Sig=1./Sig;       
        temp=Sig.*(Fa*Ci*Fb');
        Bf(:,:,i)=Fa'*temp*Fb;
    end         
    B = real(ifft(Bf,[],3));
    % update Y_1 Y_2 Y_3
    Y_1 = soft( tprod(D1,B) - Q_1./eta(1),lambda_tv*beta(1)/eta(1), Punish);
    Y_2 = soft( tprod(B,D2) - Q_2./eta(2),lambda_tv*beta(2)/eta(2), Punish);
    Y_3 = soft( tprod(D3,B) - Q_3./eta(3),lambda_tv*beta(3)/eta(3), Punish);

    % updata T 

    w1 = weightTen.*lambda1./nu;
    T = soft(D - B -N+P./nu,w1,Punish);  
    sorted_T = sort(T(:), 'descend');
%     T(find(T<0))=0;
    count_T = length(find(sorted_T>0));
    if iter~=1
        count_T_W = compute_weighted_slices(T);
    else
        count_T_W=1;
    end
    weightTen = 1./ ((normalize_3D(abs(T))  +0.001).*priorWeight);% display_images(T)

    % updata N  

    N = (nu*(D -B-T)+P)/(2*lambda2+nu);

    % update Q_k P O_k
    Q_1 = Q_1 + eta(1)*(Y_1-tprod(D1,B));
    Q_2 = Q_2 + eta(2)*(Y_2-tprod(B,D2));
    Q_3 = Q_3 + eta(3)*(Y_3-tprod(D3,B));
    P = P     + nu*(D-B-T-N);
    O_1 = O_1 + mu(1)*(X_1-B);
    O_2 = O_2 + mu(2)*(X_2-B);
    O_3 = O_3 + mu(3)*(X_3-B);
    % update mu 

    mu = min(mu*rho, mu_bar);   
    eta = min(eta*rho, mu_bar); 
    nu = min(nu*rho, mu_bar); 
    % stop Criterion

    epsilon = norm3_fro(D-B-T-N)/norm3_fro(D);
    opts.epsilon_total(b,iter) = epsilon;
    % when the relative error  epsilon is less than tol
    if epsilon < tol 
        epsilon_num = epsilon_num+1;
        if epsilon_num>1
            converged = true;
            disp('the relative error reachs requirements');
        end
    end
    % the number of nonzero elements in target tensor no longer changes,iteration end
    if count_T1 == count_T
        count = count + 1;
    else
        count = 0;
    end
    if count == 3
        converged = true;
        disp('the number of nonzero elements in target tensor no longer changes');
    end
    count_T1=count_T;
    %  iteration information printf
    if mod( iter,1) == 0 || iter==1
        disp([' Iter: ' num2str(iter) ...
            ' |D-B-T-N|_F = ' num2str(norm3_fro(D-B-T-N))...
            ' epsilon = ' num2str(epsilon)...
            ' Number_T = ' num2str(count_T)]);
    end
    % Maximum iterations reached
    if ~converged && iter >= maxIter
            disp('Maximum iterations reached');
        converged = 1 ;
    end
    end
    tenB = cat(3,tenB,B);
    tenT = cat(3,tenT,T);
    tenW = cat(3,tenW,priorWeight0);

end
% figure 
%  subplot(1,3,1)
%  
%  r = tubalrank(X_1);
%  subplot(1,3,2)
%  r = tubalrank(X_2);
%  subplot(1,3,3)
%  r = tubalrank(X_3);
% arank = averagerank(M_3)
end

function weighted_slices = compute_weighted_slices(data)
    % 计算数据矩阵每个前切片的非零元素个数，并作为加权系数乘以每个前切片
    [n1, n2, n3] = size(data);
    % 初始化加权切片矩阵
    weighted_slices = ones([n1, n2, n3]);
    
    % 遍历每个前切片
    for i = 1:n3

        % 计算当前切片非零元素个数
        non_zero_count = length(find(data(:,:,i)>0));
%         non_zero_count1 = nnz(data(:,:,i));
        non_zero = non_zero_count/(n1*n2);
%         count_T_W = 850*non_zero;
        count_T_W = 1+tan(pi/2*non_zero);
        % 将非零元素个数作为加权系数乘以当前切片，并加到加权切片矩阵中
        weighted_slices(:, :, i) = count_T_W * weighted_slices(:,:,i);
    end
end