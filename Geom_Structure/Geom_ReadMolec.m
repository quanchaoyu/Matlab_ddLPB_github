function [Geom,Data] = Geom_ReadMolec(name,s)

if(nargin==1) || isempty(s)
    s=0;
end

switch name
    case 'Kirkwood'
        A = [1   0.00000   0.00000   0   2];
        
    case '4Spheres'
        A = [[1   0.00000   0.00000   1.2   1];...
            [1   0.00000   0.00000   2.4   1];...
            [1   0.00000   0.6   1.8   1];...
            [1   0.00000   0.00000   3.6   1]];
    case '2Spheres'
        A = [[1   0.00000   0.00000   0   1];...
            [1   0.00000   0.00000   1.2   1]];
    case 'HF'
          A = [[0.33714   0.00000   0.00000  -1.58816   2.99956];...
             [-0.33714   0.00000   0.00000   0.17646+s   3.49637]];
%       %Rotation
%         if(s==0)
%             A = [[0.33714   0.00000   0.00000  0   2.99956];...
%             [-0.33714   0.00000   0.00000  1.76462   3.49637]];
%         else
%             A = [[0.33714   0.00000   0.00000  0   2.99956];...
%             [-0.33714   0.00000   1.76462*sin(s)  1.76462*cos(s)   3.49637]];
%         end
    case 'Formaldehyde'
        A = [[0.08130   0.00000   0.00000  -1.16691   4.00253];...
            [-0.20542   0.00000   0.00000   1.42202   3.63772];...
            [0.06206   0.00000   1.76689  -2.18736-s   2.99956];...
            [0.06206   0.00000  -1.76689  -2.18736-s   2.99956]];
%         A = [[0.08130   0.00000   0.00000  0            4.00253];...
%             [-0.20542   0.00000   0.00000   2.58893   3.63772];...
%             [0.06206   0.00000   1.76689  -1.02045   2.99956];...
%             [0.06206   0.00000  -1.76689  -1.02045   2.99956]];
    case 'Benzene'
        A =   [[-0.04192   0.00000   2.29035   1.32281   4.00253];...
              [-0.04192   0.00000   2.29035  -1.32281   4.00253];...
              [-0.04198   0.00000   0.00000  -2.64562   4.00253];...
              [-0.04192   0.00000  -2.29035  -1.32281   4.00253];...
              [-0.04192   0.00000  -2.29035   1.32281   4.00253];...
              [-0.04198   0.00000   0.00000   2.64562   4.00253];...
               [0.04193   0.00103   4.05914   2.34326   2.99956];...
               [0.04193   0.00103   4.05914  -2.34326   2.99956];...
               [0.04197   0.00000   0.00000  -4.68652   2.99956];...
               [0.04193  -0.00103  -4.05914  -2.34326   2.99956];...
               [0.04193  -0.00103  -4.05914   2.34326   2.99956];...
               [0.04197   0.00000   0.00000   4.68652   2.99956]];

    case 'Caffeine'
        A =    [[0.09581   1.64052   0.26811   0.00079   4.00253];...
               [0.17352   0.53040  -2.07460  -0.00051   4.00253];...
              [-0.25694  -1.93056  -2.33500  -0.00018   3.80402];...
               [0.33039  -3.37146  -0.33148   0.00013   4.00253];...
              [-0.27154  -5.88699  -0.46410  -0.00009   3.63772];...
              [-0.25299  -2.46950   2.04845   0.00073   3.80402];...
               [0.21982   0.08761   2.49865  -0.00032   4.00253];...
              [-0.28184   1.30432   4.81168  -0.00073   3.63772];...
              [-0.22781   4.22705  -0.09253   0.00089   3.80402];...
              [-0.20253   2.46931  -3.89829   0.00003   3.80402];...
               [0.12910   4.70148  -2.75142  -0.00017   4.00253];...
              [-0.09227   6.17186   1.91355  -0.00056   4.00253];...
              [-0.08846  -3.06677  -4.83502  -0.00048   4.00253];...
              [-0.08265  -3.86502   4.25693  -0.00006   4.00253];...
               [0.04393   6.41587  -3.85765  -0.00002   2.99956];...
               [0.09224   5.95218   3.08325   1.68217   2.99956];...
               [0.09228   5.95188   3.08333  -1.68106   2.99956];...
               [0.05299   8.04286   1.05689  -0.00040   2.99956];...
               [0.07429  -2.46586  -5.86225   1.68268   2.99956];...
               [0.07434  -2.46688  -5.86083  -1.68167   2.99956];...
               [0.10225  -5.11800  -4.65665  -0.00014   2.99956];...
               [0.08476  -5.05274   4.31652   1.68185   2.99956];...
               [0.08479  -5.05304   4.31660  -1.68138   2.99956];...
               [0.10651  -2.59079   5.87351   0.00130   2.99956]];

        
end

Geom.R = A(:,5)';
Geom.centers = A(:,2:4);
Geom.M = size(A,1);
Geom.charges = A(:,1);

Data.charges = A(:,1);
Data.x = Geom.centers;
