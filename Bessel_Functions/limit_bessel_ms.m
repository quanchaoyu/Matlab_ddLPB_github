function limit_bessel_ms
% test the limit when kappa --> 0

kappa = 0.01;
r = 0:0.001:1;
x = kappa.*r;

for n = 0:10
    norm((besseli_ms(n,x)./besseli_ms(n,kappa) - r.^n))
end

end