function [leb_tmp,basis,vfac,grad,der_theta,der_phi] = spherical_harmonics(lmax)
% c
% c     creat the basis of spherical harmonics
% c     

% normalization factors for spherical harmonics
vfac = Getvfac(lmax); 

% lebedev quadrature rules
degree = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, ...
    350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, ...
    3470, 3890, 4334, 4802, 5294, 5810];
order = [3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,35,41,47,53,59,65,71,77, ...
    83,89,95,101,107,113,119,125,131];

ind = find(order>2*lmax); % at least twice order

global ind_degree_leb;
if isempty(ind_degree_leb)
    ind_degree_leb = 1;
end


global num_leb;
if isempty(num_leb)
    d = degree(ind(ind_degree_leb));
else
    d = num_leb;
end

%  create basis of spherical harmonics dimensioned the number of lebedev
%  quadratures (number of points) multiplying the number of spherical
%  harmonics (number of harmonics)

leb_tmp = getLebedevSphere(d); 
x = [leb_tmp.x,leb_tmp.y,leb_tmp.z]; 

global N_leb;
N_leb = size(x,1);

[basis,grad,der_theta,der_phi] = ylmbas(lmax,x,vfac,0);


% subplot(1,2,2)
% plot3(x(:,1)',x(:,2)',x(:,3)','r.')
% hold on;
% title(['Lebdev Quadratures: ',num2str(N_leb)])
% visusphere([0,0,0],1);
% hold on;

end