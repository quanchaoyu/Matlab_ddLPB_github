function der_in = der_besseli_ms(n,x)
% derivative of the modified spherical Bessel function of the first kind

der_in = 1/(2*n+1)*(n*besseli_ms(n-1,x) + (n+1)*besseli_ms(n+1,x));

end