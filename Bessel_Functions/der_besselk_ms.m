function der_kn = der_besselk_ms(n,x)
% derivative of the modified spherical Bessel function of the first kind

der_kn = -1/(2*n+1)*(n*besselk_ms(n-1,x) + (n+1)*besselk_ms(n+1,x));

end