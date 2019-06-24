function G0 = comp_G0(Geom,lmax,w_leb,Chi_e,basis_Y,phi0)

M = Geom.M;

G0 = zeros((lmax+1)^2,M); % G0(lm,j)
for j = 1:M
    vec = w_leb.*Chi_e(:,j);
    fac = vec*ones(1,(lmax+1)^2);
    G0(:,j) = -(fac.*basis_Y)'*phi0(:,j);
end


end