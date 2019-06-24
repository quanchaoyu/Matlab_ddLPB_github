function Er = energy_reac(Geom,lmax,Xr,basis_Y)
M = Geom.M;
Er = 0;

if Geom.M == 1 % Gase 5, Kirkwood
    M = 6;
    C =  [0.2, 0.2, 0.2; 0.5, 0.5, 0.5; 0.8, 0.8, 0.8; -0.2, ...
              0.2, -0.2; 0.5, -0.5, 0.5; -0.8, -0.8, -0.8];
    charges = ones(6,1);
    vfac = Getvfac(lmax);
    
    for i = 1:M
        svec = C(i,:)/norm(C(i,:));
        basloc = ylmbas(lmax,svec,vfac,0);
        for ind = 1:(lmax+1)^2
            l = floor(sqrt(ind-1));
            Er = Er + 0.5*charges(i)*Xr(ind,1)*((norm(C(i,:))/Geom.R(1))^l)*basloc(ind);
        end
        
    end
    
    Er = Er*332.06364;
    disp(['Esolv : ', num2str(Er)])
    return
end

for j = 1:M
    Xr_j_0 = Xr((j-1)*(lmax+1)^2+1,1);
    psi_r_center = Xr_j_0*basis_Y(1,1); % 1/(2*sqrt(pi))
    Er = Er+0.5*Geom.charges(j)*psi_r_center;
end

end