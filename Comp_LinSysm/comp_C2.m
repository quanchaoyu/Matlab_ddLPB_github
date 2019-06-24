function C2 = comp_C2(e1,e2,Geom,lmax,w_leb,N_leb,basis_Y,Chi_e,Q,coefi_der)

M = Geom.M;
C2 = zeros((lmax+1)^2,(lmax+1)^2,M,M);

ind = 1:(lmax+1)^2;
ind_l = floor(sqrt(ind-1))+1;
for i = 1:M
    for j = 1:M
        vec = w_leb.*Chi_e(:,j);
        fac_Y = vec*ones(1,(lmax+1)^2);
        
        fac_Q = ones(N_leb,1)*coefi_der(ind_l,i)';
        C2(:,:,i,j) = -(fac_Y.*basis_Y)'*(fac_Q.*Q(:,:,i,j));
    end
end

end