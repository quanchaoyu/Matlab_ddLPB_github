function [basis,grad,der_theta,der_phi] = ylmbas(lmax,x,vfac,der)
% c
% c   calculates a basis of real spherical harmonics up to order lmax.
% c   the routine computes as a first the generalized legendre polynomials
% c   and then weighs them with sine and cosine coefficients.
% c   this lasts are evaluated starting from cos(phi) = y/sqrt(1-z**2)
% c   using chebyshev polynomials.
% c
% c   input:  lmax  ... maximum angular momentum of the spherical harmonics 
% c                     basis
% c           x     ... unit vector (x \in s_2) where to evaluate the basis
% c           vfac  ... scratch array dimensioned 2*nbasis. Contains the
% c                     normalization factors for the spherical harmonics.
% c
% c   output: basis ... spherical harmonics at x. basis is dimensioned (lmax+1)^2
% c                     and the harmonics are ordered for increasing l and m 
% c                     (i.e.: 0,0; 1,-1; 1,0; 1,1; 2,-2; ...)
% c
% c                     grad: grad ylm in Cartesian coordinates
% c                     der_theta: derivative of ylm w.r.t. theta
% c                     der_phi: 1/sin(theta) times the derivative of ylm w.r.t. phi
% c
% c   scratch: vplm ... scratch array dimensioned (lmax + 1)^2. it hosts the
% c                     generalized legendre polynomials.
% c            vcos ... scratch array dimensioned lmax + 1. it hosts the cos(m*\phi)
% c                     values
% c            vsin ... scratch array dimensioned lmax + 1. it hosts the sin(m*\phi)
% c                     values.
% c
 
% c
% c   get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
% c   coordinates of x.
% c

if(nargin==3)
    der = 0;    % derivative is not computed
end

Nx = size(x,1);
if(der)
    grad = zeros(Nx,3,(lmax+1)^2);
else
    grad = 0;
end


basis = zeros(Nx,(lmax+1)^2);
vcos = zeros(Nx,lmax+1);
vsin = zeros(Nx,lmax+1);
 
cthe = x(:,3);
sthe = sqrt(1 - cthe.*cthe);
 
sind = sthe~=0;
sind_pos = size(x(sind,1),1) > 0;
if(sind_pos)
    cphi(sind,1) = x(sind,1)./sthe(sind);
    sphi(sind,1) = x(sind,2)./sthe(sind);
end
cphi(~sind,1) = 1;
sphi(~sind,1) = 0;
% c
% c   evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is 
% c   pointless if z = 1, as the only non vanishing terms will be the 
% c   ones with m=0.
% c
if(sind_pos)
    [vc,vs] = trgev(lmax,cphi(sind),sphi(sind));
    vcos(sind,:) = vc;
    vsin(sind,:) = vs;
end
vcos(~sind,1:lmax + 1) = 1;
vsin(~sind,1:lmax + 1) = 0;
VC = 0;
VS = cthe(~sind);


% c
% c   evaluate the generalized legendre polynomials.
% c
vplm = polleg(lmax,cthe,sthe);
% c
% c   now build the spherical harmonics. we will distinguish m=0,
% c   m>0 and m<0:
% c

if(der)
    ethe = [cthe.*cphi,cthe.*sphi,-sthe];
    if(sind_pos)
        ephi(sind,:) = [-sphi(sind)./sthe(sind), cphi(sind)./sthe(sind), zeros(size(sphi(sind)))];
    end
    Nss = size(sthe(~sind),1);
    ephi(~sind,:) = [zeros(Nss,1),ones(Nss,1),zeros(Nss,1)];
end


for l = 0:lmax
    ind = l^2 + l + 1;
    % c       m = 0
    basis(:,ind) = vfac(1,ind)*vplm(:,ind);
    if(der)
        if(l==0)
            grad(:,:,ind) = zeros(Nx,3);
            
            der_theta(:,ind) = zeros(Nx,1);
            der_phi(:,ind) = zeros(Nx,1);
        else
            
            dPlm = vplm(:,ind+1); %% the original codes is Correct!!
            
            grad(:,:,ind) = vfac(1,ind)*ethe.*(dPlm*ones(1,3));
            
            der_theta(:,ind) = vfac(1,ind)*dPlm;
            der_phi(:,ind) = zeros(Nx,1);
        end
    end
    for m = 1:l
        plm = vplm(:,ind+m);
        % c         m > 0
        basis(:,ind+m) = vfac(1,ind+m)*plm.*vcos(:,m+1);
        % c         m < 0
        basis(:,ind-m) = vfac(1,ind+m)*plm.*vsin(:,m+1);
        
        
        if(der)
            if(m<l)
                dPlm = 0.5*((l+m)*(l-m+1)*vplm(:,ind+m-1) - vplm(:,ind+m+1));
            else
                dPlm = 0.5*(l+m)*(l-m+1)*vplm(:,ind+m-1);
            end
            
            der_theta(sind,ind+m) = -vfac(1,ind+m).*dPlm(sind).*vcos(sind,m+1);
            der_theta(sind,ind-m) = -vfac(1,ind+m).*dPlm(sind).*vsin(sind,m+1);
            der_theta(~sind,ind+m) = -vfac(1,ind+m).*dPlm(~sind).*vcos(~sind,m+1);
            der_theta(~sind,ind-m) = -vfac(1,ind+m).*dPlm(~sind).*vsin(~sind,m+1);
            
            
            % compute the der w.r.t phi over sin(theta)
            der_phi(sind,ind+m) = - m*vfac(1,ind+m)*plm(sind).*vsin(sind,m+1)./sthe(sind);
            der_phi(sind,ind-m) =   m*vfac(1,ind+m)*plm(sind).*vcos(sind,m+1)./sthe(sind);
            % 1)after multiplying 1/sin(theta) 2)previous code, VC = 0, something wrong!
            der_phi(~sind,ind+m) = - m*vfac(1,ind+m)*(-dPlm(~sind)).*vsin(~sind,m+1);%- vfac(1,ind+m)*dPlm(~sind).*VC;%
            der_phi(~sind,ind-m) =   m*vfac(1,ind+m)*(-dPlm(~sind)).*vcos(~sind,m+1);%- vfac(1,ind+m)*dPlm(~sind).*VS;%

            grad(sind,:,ind+m) = vfac(1,ind+m)*(-ethe(sind,:).*(dPlm(sind).*vcos(sind,m+1)*ones(1,3)) - m*ephi(sind,:).*(plm(sind).*vsin(sind,m+1)*ones(1,3)));
            grad(sind,:,ind-m) = vfac(1,ind+m)*(-ethe(sind,:).*(dPlm(sind).*vsin(sind,m+1)*ones(1,3)) + m*ephi(sind,:).*(plm(sind).*vcos(sind,m+1)*ones(1,3)));
            grad(~sind,:,ind+m) = vfac(1,ind+m)*(-ethe(~sind,:).*(dPlm(~sind).*vcos(~sind,m+1)*ones(1,3)) - ephi(~sind,:).*(dPlm(~sind).*VC*ones(1,3)));
            grad(~sind,:,ind-m) = vfac(1,ind+m)*(-ethe(~sind,:).*(dPlm(~sind).*vsin(~sind,m+1)*ones(1,3)) - ephi(~sind,:).*(dPlm(~sind).*VS*ones(1,3)));

        end

    end
end

if der == 0
    grad = [];
    der_theta = [];
    der_phi = [];
end

end


