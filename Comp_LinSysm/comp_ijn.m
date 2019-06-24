function [r_ijn,s_ijn] = comp_ijn(Geom,x_leb,N_leb)
% compute the data r_ijn and s_ijn

M = Geom.M;
C = Geom.centers;
R = Geom.R;

r_ijn = zeros(N_leb,M,M); 
s_ijn = zeros(3,N_leb,M,M);

for i = 1:M
    for j = 1:M
        ci = C(i,:);
        cj = C(j,:);
        ri = R(i);
        rj = R(j);
        
        r_ijn(:,i,j) = sqrt(sum((ones(N_leb,1)*(cj-ci)+rj*x_leb(:,:)).^2,2));
        s_ijn(:,:,i,j) = (ones(N_leb,1)*(cj-ci)+rj*x_leb(:,:))'./(ones(3,1)*r_ijn(:,i,j)');
        
%         for n = 1:N_leb
%             s_ijn(:,n,i,j) = (cj-ci+rj*x_leb(n,:))'/r_ijn(n,i,j);
%         end
       
    end
end

end