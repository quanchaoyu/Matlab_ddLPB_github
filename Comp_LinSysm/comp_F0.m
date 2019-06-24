function F0 = comp_F0(e1,e2,Geom,lmax,w_leb,N_leb,Chi_e,basis_Y,phi0_der,Y_s_ijn,coef,coefk)

M = Geom.M;

%% c0
c0 = zeros((lmax+1)^2,M); % c0(lm,i)
for i = 1:M
    vec = w_leb.*Chi_e(:,i);
    fac = vec*ones(1,(lmax+1)^2);
    c0(:,i) = (fac.*basis_Y)'*phi0_der(:,i);
end

%% S
S = zeros(N_leb,M,M); % S(n,i,j)

ind = 1:(lmax+1)^2;
ind_l = floor(sqrt(ind-1))+1;
coef_Y = zeros(N_leb,(lmax+1)^2,M,M);
for i = 1:M
    for j = 1:M
        coef_Y(:,:,i,j)= (ones(N_leb,1)*coef(ind_l,i)').*coefk(:,ind_l,i,j);       
    end
end
for i = 1:M
    for j = 1:M
        coef_Yij = coef_Y(:,:,i,j).*Y_s_ijn(:,:,i,j);
        S(:,i,j) = coef_Yij*c0(:,i);
    end
end

%% F0
F0 = zeros((lmax+1)^2,M);%F0(lm,j)
for j = 1:M
    vec = w_leb.*Chi_e(:,j);
    fac = vec*ones(1,(lmax+1)^2);
    F0(:,j) = -e1/e2*(fac.*basis_Y)'*sum(S(:,:,j),2);
end

end