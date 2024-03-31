function P = solve_P(Bho,c,alpha)

if nargin < 3
   d = length(alpha);
   alpha = ones(d,1)/d;
end

[num,m] = size(Bho{1});

G = zeros(num,m);
for o = 1:length(Bho)
    G = G+1/alpha(o)*Bho{o};
end
G = sparse(G/sum(1./alpha));

% P = zeros(num,m);
% for o = 1:length(Bho)
%     P = P+alpha(o)*Bho{o};
% end

P = Bho{1};
P = sparse(P);
Dn = sparse(diag(sum(P,2).^(-0.5)));
Dm = sparse(diag(sum(P).^(-0.5)));

X = Dn*P*Dm;

[U,S,V] = svds(X,c);

Fn = sqrt(2)*Dn*U/2;
Fm = sqrt(2)*Dm*V/2;
lambda = 0.2;
islocal = 0; % only update the similarities of neighbors if islocal=1
NITER = 30;
%%
for iter = 1:NITER
    dist = L2_distance_1(Fn',Fm');
    P = zeros(num,m);
    for i = 1:num
        g0 = G(i,:);
        if islocal == 1
            idxg0 = find(g0>0);
        else
            idxg0 = 1:m;
        end
        gi = g0(idxg0);
        di = dist(i,idxg0);
        ad = gi-0.5*lambda*di; 
        P(i,idxg0) = EProjSimplex_new(ad);    
    end
    P = sparse(P);    
    Dn = sparse(diag(sum(P,2).^(-0.5)));
    Dm = sparse(diag(sum(P).^(-0.5)));
    
    X = Dn*P*Dm;
    if ismember(Inf,P)==1
        P = Bho{1};
        break;
    end
    
    [U,S,V] = svds(X,c+1);
    
    Fn = sqrt(2)*Dn*U(:,1:c)/2;
    Fm = sqrt(2)*Dm*V(:,1:c)/2;
    ev = diag(S);
    
    Fn_old = Fn;
    Fm_old = Fm;
    
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 < c-0.0000001
        lambda = 2*lambda;
    elseif fn2 > c+1-0.0000001
        lambda = lambda/2;   Fn = Fn_old; Fm = Fm_old;
    else
        break;
    end
end

% STemp = sparse(num+m,num+m);
% STemp(1:num,num+1:end) = P;
% STemp(num+1:end,1:num) = P';
% [clusternum, y] = graphconncomp(STemp);
% clusternum
% y1 = y(1:num)';
% y2 = y(num+1:end)';

% result = ClusteringMeasure1(Y, y1)


end
























