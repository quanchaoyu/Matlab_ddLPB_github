function [Er,Xr] = DD_PB(e1,e2,kappa,Rp,lmax,visu,filename)

%tic 


tokcal= 332.06364; 627.509469;

%% Data structure of the molecular geometry

Geom = Geom_ReadMolec(filename); %Geom = read_PDB(filename);

toang=0.52917721092;
Geom.R = toang*Geom.R;
Geom.centers = toang*Geom.centers;

Geom.R = Geom.R+Rp;

t = 0; % time (molecular dynamics)
Geom.centers = Geom.centers;

%% Lebedev quadrature points & weights, unit sphere

[leb_tmp,basis_Y] = spherical_harmonics(lmax); % x_leb: Lebdev nodes, w_leb: Lebdev weights
x_leb = [leb_tmp.x,leb_tmp.y,leb_tmp.z];
w_leb = leb_tmp.w;
N_leb = size(x_leb,1);

%phi0(n,j),phi0_der(n,j)
[phi0,phi0_der] = comp_phi0(Geom,x_leb,N_leb,e1);

%% Compute r_ijn & s_ijn, Y_s_ijn, depending on t!!!
%r_ijn(n,i,j),s_ijn(1:3,n,i,j)
[r_ijn,s_ijn] = comp_ijn(Geom,x_leb,N_leb);

%Y_s_ijn(n,lm,i,j)
Y_s_ijn = comp_Y_s_ijn(lmax,s_ijn,Geom.M,N_leb);

%% Characteristic function Chi, depending on t!!!
%intersection between balls, depending on t!!!
I = structure(Geom,Rp);
Inter = I.M_int; %intersection matrix
InterN = I.num_int; %intersection number for each SAS-ball

%Chi_e(n,j)
Chi_e = comp_Chi_e(Geom,Inter,InterN,x_leb,N_leb);

%w(n,i,j)
[w,Chi,N] = comp_w(Geom,Inter,InterN,x_leb,N_leb);

%% Compute coefi,coefk,coef; P_Chi; Q depending on t!!!
%coefi(n,l,i,j),coefk(n,l,i,j),coef(l,i),coefi_der(l,i): Bessel coefs
[coefi,coefk,coef,coefi_der] = comp_coef_bessel(kappa,Geom,r_ijn,lmax,N_leb);

%P_Chi(lm,l'm',M)
P_Chi = comp_P_Chi(Geom,lmax,x_leb,w_leb,N_leb,basis_Y,Chi_e);

%Q(n,lm,i,j)
Q = comp_Q(Geom,lmax,N_leb,Y_s_ijn,coefk,coef,P_Chi);

%% Compute C1(lm,l'm',i,j),C2(lm,l'm',i,j),F0(lm,j)
C1 = comp_C1(e1,e2,Geom,lmax,w_leb,N_leb,basis_Y,Chi_e,Q);
C2 = comp_C2(e1,e2,Geom,lmax,w_leb,N_leb,basis_Y,Chi_e,Q,coefi_der);
F0 = comp_F0(e1,e2,Geom,lmax,w_leb,N_leb,Chi_e,basis_Y,phi0_der,Y_s_ijn,coef,coefk);

%% Compute G0,A,B
%G0(lm,j), A(lm,l'm',i,j), B(lm,l'm',i,j)
G0 = comp_G0(Geom,lmax,w_leb,Chi_e,basis_Y,phi0);
A = comp_A(Geom,lmax,Inter,InterN,w_leb,N_leb,w,basis_Y,r_ijn,Y_s_ijn);
B = comp_B(Geom,lmax,Inter,InterN,w_leb,N_leb,w,basis_Y,coefi,Y_s_ijn);

%% Solve linear system
[Xr,Xe] = solver_lin(Geom,lmax,A,B,G0,C1,C2,F0);

%% Reaction energy
Er = tokcal*energy_reac(Geom,lmax,Xr,basis_Y);
disp('------------------- Electrostatic Reaction Energy -------------------');
disp(['Electrostatic reaction energy: ',num2str(Er)]);
%toc

%% Visu
if (visu)
    Visu_Psi(Geom,Xr,lmax)
end

%% Store data
global DATA;
DATA.basis_Y = basis_Y;
DATA.N_leb = N_leb;
DATA.M = Geom.M;
DATA.R = Geom.R;
DATA.psi0 = phi0;
DATA.Chi_e = Chi_e;
DATA.w_leb = w_leb;

end