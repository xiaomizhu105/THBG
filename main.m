clc;clear all;
addpath([pwd, '/funs']);
data_ind = 6;
load('JAFFE_1024.mat'); % ‘ÿ»Î ˝æ›
X = X./max(X,[],2);
[num, dim] = size(X);
c = length(unique(Y));

%% obtain bipartite graph
tic
order = 2;
tau = 1;
m = fix(0.4*num);
k = 5;
rng(4);
style = 4; % 1-direct sample; 2-rand sample; 3-KNP; 4-kmeans sample
B = bipartite_gen(X,m,k,style);
toc

%% obtain high-order bipartite graphs
tic
Bho = high_order_bipartite_gen(B,order);
toc

%% Optimization
alpha = ones(order,1)/order;
tic
NITER = 30;
P = Bho{1};
lambda = 0.2;
f_value = [];
for iter = 1:NITER
    fprintf('The %d-th iteration...\n',iter);
    % update F
    [Fn,Fm] = solve_F(P,c);
    
    % update P 
    P = solve_P1(Bho,alpha,Fn,Fm,lambda);
    P = sparse(P);
    if ismember(Inf,P)==1
        P = Bho{1};
        break;
    end
    
    % update A_tensor
    B_tensor = cat(3,Bho{:});
    B_tensor = permute(B_tensor,[1,3,2]);
    [A_tensor,tnn,trank] = solve_A_tensor(B_tensor,tau);
    
    % update B
    A = permute(A_tensor,[1,3,2]);
    Bho = solve_B(A,P,alpha);
    B_tensor = cat(3,Bho{:});
    
    % update alpha
    alpha = solve_alpha(P,Bho); 
    
      
    Dn = sparse(diag(sum(P,2).^(-0.5)));
    Dm = sparse(diag(sum(P).^(-0.5)));
    
    X = Dn*P*Dm;
    
    [U,S,V] = svds(X,c+1);
    ev = diag(S);
    
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    
    f1 = 0;
    for o = 1:order
       f1 = f1 + 1/alpha(o)*norm(P-Bho{o},'fro');
    end
    f2 = 0;
    for o = 1:order
       f2 = f2 + tau/2*norm(A(:,:,o)-B_tensor(:,:,o),'fro');
    end
    f_value(iter) = tnn + f1 + f2;
   
   
    if fn1 < c-0.0000001
        lambda = 2*lambda;
    elseif fn2 > c+1-0.0000001
        lambda = lambda/2;
    else
        break;
    end
    
    
end
toc

%% obtain the result
S = sparse(num+m,num+m);
S(1:num,num+1:end) = P;
S(num+1:end,1:num) = P';
[clusternum, y] = graphconncomp(S);
y1 = y(1:num)';
y2 = y(num+1:end)';

result = ClusteringMeasure1(Y, y1)

f_new = [];
n = length(f_value);
f_new(1) = f_value(1);
for i = 2:n
    f_new(i) = f_new(i-1) - abs(f_value(i-1)-f_value(i));
end
f_new =f_new';