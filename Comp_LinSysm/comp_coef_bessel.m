function [coefi,coefk,coef,coefi_der] = comp_coef_bessel(kappa,Geom,r_ijn,lmax,N_leb)
M = Geom.M;
coefi = zeros(N_leb,lmax+1,M,M);
coefk = zeros(N_leb,lmax+1,M,M);

V_Inf = 10^80;
V_zero = 10^-80;

for l = 0:lmax
    for i = 1:M
        ri = Geom.R(i);
        il_ri = besseli_ms(l,kappa*ri);
        kl_ri = besselk_ms(l,kappa*ri);
        
%         coefi(:,l,:,:) = besseli_ms(l,kappa*r_ijn)/il_ri;
%         coefk(:,l,:,:) = besseli_ms(l,kappa*r_ijn)/kl_ri;

        for j = 1:M
            x = r_ijn(:,i,j);
            coefi(:,l+1,i,j) = besseli_ms(l,kappa*x)/il_ri;
            coefk(:,l+1,i,j) = besselk_ms(l,kappa*x)/kl_ri;
            
            %% Overflow case: kappa is too BIG
            if il_ri >= V_Inf
                coefi(:,l+1,i,j) = zeros(size(x)); % 0
            end
            if kl_ri == 0
                coefk(:,l+1,i,j) = zeros(size(x)); % 0
            end
            
            %% Overflow case: kappa is 0
            if kappa <= V_zero
                coefi(:,l+1,i,j) = (x./ri).^l;
            end
            if kl_ri >= V_Inf
                coefk(:,l+1,i,j) = (ri./x).^(l+1); 
            end
            
        end

    end
end

coef = zeros(lmax+1,M);
coefi_der = zeros(lmax+1,M);
for l = 0:lmax
    for i = 1:M
        ri = Geom.R(i);
        il_ri = besseli_ms(l,kappa*ri);
        kl_ri = besselk_ms(l,kappa*ri);
        
        term_i = der_besseli_ms(l,kappa*ri)/il_ri;
        term_k = der_besselk_ms(l,kappa*ri)/kl_ri;
        
        %% Overflow case: kappa is too BIG
        if il_ri >= V_Inf
            term_i = 1;
        end
        if kl_ri <= V_zero
            term_k = -1;
        end
        
        %% Overflow case: kappa is 0, NOTE: coef is divied by (kappa*ri), but coefi_der is timed by (kappa*ri), which will be offseted
        if kappa <= V_zero
            kappa_term_i = l/ri + (l+1)*(kappa^2*ri)/((2*l+1)*(2*l+3)); 
        else
            kappa_term_i = kappa*term_i;
        end
        if kl_ri >= V_Inf
            kappa_term_k = - (l+1)/ri - l*(kappa^2*ri)/((2*l-1)*(2*l+1));
        else
            kappa_term_k = kappa*term_k;
        end
        
        coef(l+1,i) = (kappa_term_i - kappa_term_k)^(-1);
        coefi_der(l+1,i) = kappa_term_i;
        
        
    end
end

end