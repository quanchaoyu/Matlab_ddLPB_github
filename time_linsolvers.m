function time_linsolvers


n = 10;

global err;
global Time;

t = zeros(1,3);

for i = 1:n
    TestsddPB;
    t = t+Time;
end

t = t/n; % average running time

t
%HF: err = [1.60554941879367e-06,1.32287227743264e-06];
T = [0.0440    0.1740    0.1850];
%Formaldehyde: err = [2.53660128772241e-06,2.19273433976886e-06];
T = [T;0.2570    0.4670    0.3720];
%Benzene: err = [2.95720135769766e-06,1.86573970315471e-06];
T = [T;6.0170   11.1620    4.6480];
%Caffein: err = [0.245205887354083,2.01005940118764e-06];
T = [T;50.5490  140.5670   10.5430];
end