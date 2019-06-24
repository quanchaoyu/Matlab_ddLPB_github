function Y_s_ijn = comp_Y_s_ijn(lmax,s_ijn,M,N_leb)

Y_s_ijn = zeros(N_leb,(lmax+1)^2,M,M);

vfac = Getvfac(lmax);
for i = 1:M
    for j = 1:M
        Y_s_ijn(:,:,i,j) = ylmbas(lmax,s_ijn(:,:,i,j)',vfac,0);
    end
end

end