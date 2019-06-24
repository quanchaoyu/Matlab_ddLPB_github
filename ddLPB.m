function ddLPB
clear all
InitCode_PB;

% input molecule name
filename = 'Benzene'; % HF, Formaldehyde, Benzene, Caffeine

% initial parameters
e1 = 1; 
e2 = 78.54; % water dielectric constant
kappa = 0.104;
global Rp; Rp = 0; % solvent radius
global lmax; lmax = 9; % maximum degree of spherical harmonics
global num_leb; num_leb = 194; % number of Lebedev points

% run ddLPB solver and visualize electrostatic potential
close all
visu = 1;
tic
DD_PB(e1,e2,kappa,Rp,lmax,visu,filename);
toc
view([90 0]);

end



