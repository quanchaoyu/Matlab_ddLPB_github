function C1 = comp_C1(e1,e2,Geom,lmax,w_leb,N_leb,basis_Y,Chi_e,Q)

M = Geom.M;
C1 = zeros((lmax+1)^2,(lmax+1)^2,M,M);

ind = 1:(lmax+1)^2;
list_l = floor(sqrt(ind-1));
for i = 1:M
    ri = Geom.R(i);
    for j = 1:M
        vec = w_leb.*Chi_e(:,j);
        fac_Y = vec*ones(1,(lmax+1)^2);
        
        fac_Q = 1/ri*(ones(N_leb,1)*list_l);
        C1(:,:,i,j) = e1/e2*(fac_Y.*basis_Y)'*(fac_Q.*Q(:,:,i,j));
    end
end

end