function Q = comp_Q(Geom,lmax,N_leb,Y_s_ijn,coefk,coef,P_Chi)
M = Geom.M;
Q = zeros(N_leb,(lmax+1)^2,M,M);

ind = 1:(lmax+1)^2;
ind_l = floor(sqrt(ind-1))+1;

coef_Y = zeros(N_leb,(lmax+1),M,M);
for i = 1:M
    for j = 1:M
        coef_Y(:,:,i,j)= (ones(N_leb,1)*coef(:,i)').*coefk(:,:,i,j);       
%         for l = 1:lmax
%             len = 2*l+1;
%             coef_Y(:,(l^2+1:l^2+len),i,j)= ((ones(N_leb,1)*coef(l,i)).*coefk(:,l,i,j))*ones(1,len);
%         end
    end
end

for i = 1:M
    for j = 1:M
        coef_Yij = coef_Y(:,ind_l,i,j).*Y_s_ijn(:,:,i,j);
        Q(:,:,i,j) = coef_Yij*P_Chi(:,:,i);
        
        if i == 2 && j == 1
            %coef_Yij(5,1)
        end
    end
end

end