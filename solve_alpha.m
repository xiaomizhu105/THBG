function alpha = solve_alpha(P,Bho)

%%
order = length(Bho);

h = zeros(order,1);
for o = 1:order
    h(o) = norm(P-Bho{o}, 'fro');
end

alpha = h./sum(h);


end