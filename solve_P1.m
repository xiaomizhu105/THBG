function P = solve_P1(Bho,alpha,Fn,Fm,lambda)


[num,m] = size(Bho{1});

G = zeros(num,m);
for o = 1:length(Bho)
    G = G+1/alpha(o)*Bho{o};
end
G = sparse(G/sum(1./alpha));

islocal = 0; % only update the similarities of neighbors if islocal=1
NITER = 1;
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
        ad = gi-0.5*lambda/sum(1./alpha)*di; 
        P(i,idxg0) = EProjSimplex_new(ad);    
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
























