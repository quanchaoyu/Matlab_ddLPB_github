function kn = besselk_ms(n,x)
% modified spherical Bessel function of the second kind
% n is an integer, x is a set of points > 0
% return k_n

kn = sqrt(2./(pi.*x)).*besselk(n+1/2,x);

end