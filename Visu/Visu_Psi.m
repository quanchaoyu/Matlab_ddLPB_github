function Visu_Psi(Geom,Xr,lmax)
figure;

RES = 50;
for j = 1:Geom.M
    indj = (j-1)*(lmax+1)^2+1:j*(lmax+1)^2;
    Xr_j = Xr(indj);
    ShowPsi(Geom,j,Xr_j,lmax,RES);
end

%title('Reaction potential on VDW computed by ddLPB');
light('position',20*[1,1,1]);
lighting gouraud;

daxis = 3;
R = Geom.R';
C = Geom.centers;
xmin = min(C(:,1)-R)-daxis;xmax = max(C(:,1)+R)+daxis;
ymin = min(C(:,2)-R)-daxis;ymax = max(C(:,2)+R)+daxis;
zmin = min(C(:,3)-R)-daxis;zmax = max(C(:,3)+R)+daxis;
axis([xmin xmax ymin ymax zmin zmax]);
end