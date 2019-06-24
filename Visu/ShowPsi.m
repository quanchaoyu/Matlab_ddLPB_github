function ShowPsi(Geom,j,Xr_j,lmax,RES)

cj = Geom.centers(j,:);
rj = Geom.R(j);

%% Discrete points
% points on unit sphere
[X,Y,Z] = sphere(RES);

% coordinate matrix
x = [reshape(X,(RES+1)^2,1),reshape(Y,(RES+1)^2,1),reshape(Z,(RES+1)^2,1)];

% translate to physical sphere
X = rj*X + cj(1)*ones(size(X));
Y = rj*Y + cj(2)*ones(size(Y));
Z = rj*Z + cj(3)*ones(size(Z));

%% Psi
vfac = Getvfac(lmax);
% Ylm(x)
basis_Y_x = ylmbas(lmax,x,vfac,0);

% Psi_r(x) on Gamma_j
Psi_r_j = basis_Y_x(:,:)*Xr_j;

%% Visu

surf(X,Y,Z,reshape(Psi_r_j,RES+1,RES+1),'EdgeColor','none','FaceLighting','gouraud','FaceAlpha', 1);
hold on;
axis equal off; 
%set(gca,'visible','off','DataAspectRatio',[1 1 1])
%camzoom(1.0);
colorbar;
%colormap(hot)
%caxis([-0.01,0.02])

%set(gca,'CameraPosition',[1000,1000,0])
end