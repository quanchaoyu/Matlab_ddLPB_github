function in = besseli_ms(n,x)
% modified spherical Bessel function of the first kind
% n is an integer, x is a set of points > 0
% return i_n

in = sqrt(pi./(2.*x)).*besseli(n+1/2,x);

index = isnan(in); % index of NaN 
if isempty(find(index)) == 0 % there exists x = 0
    if n == 0
        in(index) = 1;
    else
        in(index) = 0;
    end
end

end