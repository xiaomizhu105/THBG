function Bho = high_order_bipartite_gen(B1,order)

disp('----------Generate High Order Bipartite Graphs----------');

Bho = cell(order,1);     

% 第一步归�?化步骤有待商�?
sigma = sparse(diag(sum(B1).^(-0.5)));
Bho{1} = B1*sigma;

[U,S,V] = svd(B1,'econ');

for d = 2:order
    Temp = U*S.^(2*d-1)*V';
    Temp = Temp./max(max(Temp));
    Temp(Temp<1e-5) = 0; 
    Bho{d} = Temp;
end


