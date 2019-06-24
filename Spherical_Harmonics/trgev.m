function [cx,sx] = trgev(n,x,y)

Nx = size(x,1);
cx = zeros(Nx,n+1);
sx = zeros(Nx,n+1);

cx(:,1) = 1;
sx(:,1) = 0;
if(n>0)
    sx(:,2) = y;
    cx(:,2) = x;
    for i = 3:n+1
        cx(:,i) = 2*x.*cx(:,i-1) - cx(:,i-2);
        sx(:,i) = 2*x.*sx(:,i-1) - sx(:,i-2);
    end
end