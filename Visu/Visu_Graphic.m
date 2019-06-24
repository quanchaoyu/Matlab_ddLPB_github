function Visu_Graphic

Xr = textread('/Users/quanchaoyu/Documents/Fortran/ddLPB/Xr.txt','%f \n');

data = textread('/Users/quanchaoyu/Documents/Fortran/ddLPB/files/1ajj.txt');
Geom.M = data(1,1);
Geom.charges = data(2:end,1);
Geom.centers = data(2:end,2:4);
Geom.R = data(2:end,5)';

%Geom.M = 10;
lmax = 7;11;

addpath('/Users/quanchaoyu/Documents/MATLAB/Matlab_PB');
InitCode_PB;

Visu_Psi(Geom,Xr,lmax);
%caxis([0.1 0.35])
light('position',20*[-1,-1,-1]);
view([-90 -90])

end