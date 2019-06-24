function P_Chi = comp_P_Chi(Geom,lmax,x_leb,w_leb,N_leb,basis_Y,Chi_e)

M = Geom.M;
P_Chi = zeros((lmax+1)^2,(lmax+1)^2,M);

for i = 1:M
    vec = w_leb.*Chi_e(:,i);
    fac = vec*ones(1,(lmax+1)^2);
    P_Chi(:,:,i) = basis_Y'*(fac.*basis_Y);
end

end

% for k = 1:(lmax+1)^2
%     %k = l^2+(m+l)+1;
%     l = floor(sqrt(k-1));
%     m = k-l^2-l-1;
% end