function B = solve_B(A,P,alpha)

[num,m,d] = size(A);
B = cell(d,1);


for o = 1:d
    H = alpha(o)/(alpha(o)+1)*(1/alpha(o)*P+A(:,:,o));
    
    Btemp = zeros(num,m);
    for i = 1:num
        ad =  H(i,:);
        Btemp(i,:) = EProjSimplex_new(ad);
    end
    B{o} = Btemp;
    
end


end

%%
% for o = 1:d
%     H = alpha(o)/(alpha(o)+1)*(1/alpha(o)*P+A(:,:,o));
%     [U,S,V] = svd(H,'econ');
%     S = diag(S);
%     r = length(find(S>0));
%     if r>=1
%         S = S(1:r);
%         B{o} = U(:,1:r)*diag(S)*V(:,1:r)';
%     end
%     
% end