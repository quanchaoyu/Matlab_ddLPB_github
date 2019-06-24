function [w,Chi,N] = comp_w(Geom,Inter,InterN,x_leb,N_leb)

M = Geom.M;
C = Geom.centers;
R = Geom.R;

Chi = zeros(N_leb,M,M);
N = zeros(N_leb,M);
w = zeros(N_leb,M,M);

for j = 1:M
    cj = C(j,:);
    rj = R(j);
    for int_n = 1:InterN(j)
        i = Inter(j,int_n);
        ci = C(i,:);
        ri = R(i);
        Chi(:,i,j) = sum((ones(N_leb,1)*(cj-ci)+rj*x_leb(:,:)).^2,2) < ri^2;
        N(:,j) = N(:,j)+Chi(:,i,j);
    end
end

for j = 1:M
    for int_n = 1:InterN(j)
        i = Inter(j,int_n);
        ind  = find(N(:,j)~=0); % ind: list of overlapped Lebedev points
        w(ind,i,j) = Chi(ind,i,j)./N(ind,j);
    end
end

end