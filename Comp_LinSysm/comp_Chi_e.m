function Chi_e = comp_Chi_e(Geom,Inter,InterN,x_leb,N_leb)
% Chi_e(n,j) is 1 or 0 representing if the Lebedeve point s_n is on the exterior part of Gamma_j

M = Geom.M;
C = Geom.centers;
R = Geom.R;
Chi_e = ones(N_leb,M);

for j = 1:M
    cj = C(j,:);
    rj = R(j);
    for int_n = 1:InterN(j)
        i = Inter(j,int_n);
        ci = C(i,:);
        ri = R(i);
        ind = sum((ones(N_leb,1)*(cj-ci)+rj*x_leb(:,:)).^2,2) < ri^2;
        Chi_e(ind,j) = 0;
%         for n = 1:N_leb
%             if norm(cj+rj*x_leb(n,:)-ci) < ri
%                 Chi_e(n,j) = 0;
%             end
%         end
    end
end


end