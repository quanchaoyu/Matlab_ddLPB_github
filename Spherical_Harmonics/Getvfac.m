function vfac = Getvfac(lmax)
sqrttwo = sqrt(2);

vfac = zeros(3,(lmax+1)^2);
for l = 0:lmax
    ind = l*l + l + 1;
    flm = sqrt((2*l+1)/(4*pi));
    vfac(1,ind) = flm;
    vfac(2,ind) = 0;
    vfac(4,ind) = 1;
    for m = 1:l 
        fac = sqrt(factorial(l-m)/factorial(l+m));
        fnorm = flm*fac;
        if(mod(m,2)==1) 
            fnorm = - fnorm;
            fac = - fac;
        end
        vfac(1,ind+m) = sqrttwo*fnorm;
        vfac(2,ind+m) = 0;
        vfac(1,ind-m) = 0;
        vfac(2,ind-m) = sqrttwo*fnorm;
        
        % for the complex SH's
        vfac(4,ind+m) = fac;
        vfac(4,ind-m) = 0;
    end
    if(l>0)
        vfac(3,ind) = factorial(l-1)/factorial(l+1);
    end
end
