function [Fn,Fm] = solve_F(P,c)

P = sparse(P);
Dn = sparse(diag(sum(P,2).^(-0.5)));
Dm = sparse(diag(sum(P).^(-0.5)));

X = Dn*P*Dm;

[U,S,V] = svds(X,c);

Fn = sqrt(2)*Dn*U/2;
Fm = sqrt(2)*Dm*V/2;
end
























