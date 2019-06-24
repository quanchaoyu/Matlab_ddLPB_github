function plm = polleg(lmax,x,y)
% c
% c   computes the l,m associated legendre polynomial for -1 <= x <= 1
% c   using the recurrence formula
% c
% c   (l-m)p(l,m) = x(2l-1)p(l-1,m) - (l+m-1)p(l-2,m)
% c

Nx = size(x,1);
fact  = 1;
pmm   = ones(Nx,1);
somx2 = y;
plm = zeros(Nx,(lmax+1)^2);

for m = 0:lmax 
% c
% c  do pmm for each m:
% c      
    ind  = (m+1)*(m+1);
    plm(:,ind) = pmm;
    if(m ~= lmax) 
%        break;
%    else
        pmm1 = (2*m+1)*x.*pmm;
        ind2 = ind + 2*m + 2;
        plm(:,ind2) = pmm1;
        pmmo = pmm;
        for l = m+2:lmax
            pll   = ((2*l-1)*x.*pmm1 - (l+m-1)*pmm)/(l-m);
            ind = l*l + l + 1;
            plm(:,ind+m) = pll;
            pmm  = pmm1;
            pmm1 = pll;
        end
        pmm  = -fact*pmmo.*somx2;
        fact = fact + 2;
    end
end
