function A = comp_A(Geom,lmax,Inter,InterN,w_leb,N_leb,w,basis_Y,r_ijn,Y_s_ijn)

M = Geom.M;

A = zeros((lmax+1)^2,(lmax+1)^2,M,M);
ind = 1:(lmax+1)^2;
ind_l = floor(sqrt(ind-1))+1;
for j = 1:M
    % Case: i=j
    i = j;
    A(:,:,i,j) = eye((lmax+1)^2);

    % Case: i~=j
    for int_n = 1:InterN(j)
        i = Inter(j,int_n);
        ri = Geom.R(i);
        
        vec = w_leb.*w(:,i,j);
        fac_Y = vec*ones(1,(lmax+1)^2);
        
        fac_Y_s_ijn = ((r_ijn(:,i,j)/ri)*ones(1,(lmax+1))).^(ones(N_leb,1)*(0:lmax));
        A(:,:,i,j) = -(fac_Y.*basis_Y)'*(fac_Y_s_ijn(:,ind_l).*Y_s_ijn(:,:,i,j));
    end

end


end