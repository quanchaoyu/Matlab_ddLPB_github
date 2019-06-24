function [Xr,Xe] = solver_lin(Geom,lmax,A,B,G0,C1,C2,F0)
M = Geom.M;

% U X = F
A_bf = modi_matr(A,M);
B_bf = modi_matr(B,M);
C1_bf = modi_matr(C1,M);
C2_bf = modi_matr(C2,M);
G0_col = modi_vec(G0,M);
F0_col = modi_vec(F0,M);

U = zeros(2*size(A_bf));
F = zeros(2*size(F0_col,1),1);
U = [A_bf+C1_bf,C2_bf;A_bf,-B_bf];
F = [G0_col+F0_col;G0_col];

U_org = [A_bf+C1_bf,C2_bf;C1_bf,B_bf+C2_bf];
F_org = [G0_col+F0_col;F0_col];

%% LU factorization solver
% tic
% %[L0,U0] = lu(U); % LU decomposition
% %X = U0\(L0\F);
% X = U\F;
% %X = linsolve(U,F);
% toc

% Xr = X(1:M*(lmax+1)^2,1);
% Xe = X(M*(lmax+1)^2+1:2*M*(lmax+1)^2,1);

t = cputime;
X = linsolve(U,F);
t_lu = cputime-t;

Xr = X(1:M*(lmax+1)^2,1);
Xe = X(M*(lmax+1)^2+1:2*M*(lmax+1)^2,1);

return;

%% GMRes solver
%U_pre = sparse(diag(diag(U)));          % diagonal preconditioner

% U_pre =@(X) [[A_bf MyZero];[MyZero B_bf]]\X;
% U_pre = [[A_bf MyZero];[MyZero B_bf]];

%AB = sparse([[A_bf+diag(diag(C1_bf)) MyZero];[MyZero -B_bf]]);
% % AB = sparse([[A_bf MyZero];[MyZero -B_bf]]);
% 
% ABd = sparse(diag(diag(AB)));
% % 
% U_pre =@(X) gmres(AB,X,10,10^-8,50,ABd);
% U_pre =@(X) gmres(AB,X,5,10^-3,50);

% tic
% X_gmres = gmres(U,F,10,10^-6,500,U_pre);
% %disp(['Difference between direct solver and GMRes solver: ',num2str(norm(X_gmres-X))])
% toc
% tic
% X_gmres = bicg(U,F,10^-6,800,U_pre);
% toc
% tic
% X_gmres = bicgstab(U,F,10^-6,800,U_pre);
% toc
% X = X_gmres;

% Xr = X(1:M*(lmax+1)^2,1);
% Xe = X(M*(lmax+1)^2+1:2*M*(lmax+1)^2,1);
% 

Tol = 10^-10;

t = cputime;
U_pre = sparse(diag(diag(U)));  
X_gmres = gmres(U,F,10,Tol,50,U_pre);
t_gmres = cputime-t;

% t = cputime;
% U_pre_org = [A_bf,zeros(size(C2_bf));zeros(size(C1_bf)),B_bf];  
% X_gmres_org = gmres(U_org,F_org,10,Tol,50,U_pre_org);
% t_gmres_org = cputime-t;

%% Like the global strategy
Nmax = 50;

Ad = sparse(diag(diag(A_bf))); % diagonal of A
Bd = sparse(diag(diag(B_bf)));
sA = sparse(A_bf); % sparse A
sB = sparse(B_bf); % sparse B

Fa = G0_col+F0_col;
Fb = F0_col;

t = cputime;
% Xr = bicgstab(A_bf,Fa,10^-1,50,Ad);
% Xe = bicgstab(B_bf,Fb,10^-1,50,Bd);
Xr = gmres(sA,Fa,10,10^-8,50,Ad);
Xe = gmres(sB,Fb,10,10^-8,50,Bd);

aux = C1_bf*Xr+C2_bf*Xe;

cnt = 0;
inc = 1;
relerr = 1;
while (relerr>Tol && cnt<Nmax)
    cnt = cnt + 1;
    %tol = max(inc,Tol);
    %tol = min(tol,0.1);
    tol = Tol;
    
    Xrold = Xr;
    Xeold = Xe;

    Fak = Fa-aux; 
    Fbk = Fb-aux; 

%     Xr = bicgstab(A_bf,Fak,tol,50,Ad,[],Xr);
%     Xe = bicgstab(B_bf,Fbk,tol,50,Bd,[],Xe);
     Xr = gmres(sA,Fak,10,tol,50,Ad,[],Xr);
     Xe = gmres(sB,Fbk,10,tol,50,Bd,[],Xe);

    aux = C1_bf*Xr+C2_bf*Xe;
    

%    res = sqrt(norm(A_bf*Xr + aux - Fa)^2 + norm(aux + B_bf*Xe - Fb)^2);

    inc = sqrt(norm(Xr-Xrold)^2 + norm(Xe-Xeold)^2);
    relerr = inc/sqrt(norm(Xr)^2 + norm(Xe)^2);
    
    %disp(['inc: ', num2str(inc)])
    %disp(['relerr: ', num2str(relerr)])
end
disp(['Number of iterations: ',num2str(cnt)])
disp(['Increment: ',num2str(inc)])

X_ite = [Xr;Xe];

t_ite = cputime-t;

res = sqrt(norm((A_bf+C1_bf)*Xr + C2_bf*Xe - Fa)^2 + norm(C1_bf*Xr + (B_bf+C2_bf)*Xe - Fb)^2);
disp(['Residual: ',num2str(res)])


%% Outputs
global err;
global Time;

err_gmres = norm(X_gmres-X)/norm(X);
err_ite = norm(X_ite-X)/norm(X);

err = [err_gmres,err_ite];
Time = [t_lu,t_gmres,t_ite];

end

function matr_bf = modi_matr(matr,M)

len = size(matr,1);
matr_bf = zeros(M*len,M*len);

for i = 1:M
    for j = 1:M
        indj = ((j-1)*len+1):(j*len);
        indi = ((i-1)*len+1):(i*len);
        matr_bf(indj,indi) = matr(:,:,i,j);
    end
end

end

function vec_col = modi_vec(vec,M)
len = size(vec,1);
vec_col = zeros(M*len,1);
for j = 1:M
    indj = ((j-1)*len+1):(j*len);
    vec_col(indj,1) = vec(:,j);
end
end