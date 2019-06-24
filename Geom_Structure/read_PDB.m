function [Geom,Data] = read_PDB(filename)
% c
% c   Read pdb-file and save data in Geom
% c

filename = [filename,'.pdb'];
fid = fopen(filename,'r');



%num_carbon = 0; % number of carbons
%charges_carbon = 0; % opposite total charges of all carbons
%{
Geom.centers = zeros(10^6,3);

DATA = textread(filename);
Geom.M = size(DATA,1);
Geom.centers = DATA(:,1:3);
Geom.R = DATA(:,4)';
return;
%}
line = fgetl(fid);
it = 0;

charge_negative = 0;
num_H = 0;
num_C = 0;
while ischar(line)
	if (strncmp('HETATM',line,6) || strncmp('ATOM',line,4)) 
        it = it + 1;
        switch(line(14))%14))
%%%%%%%%% UFF VdW table
            case  'H', color = [0.7 0.7 0.7]; r = 1.443; charge = 1;
            case  'C', color = [0.3 0.3 1.0]; r = 1.9255; charge = -0.2;
            case  'O', color = [0.3 1.0 0.3]; r = 1.75; charge = -0.1;
            case  'N', color = [1.0 0.3 1.0]; r = 1.83; charge = -0.15;
            case  'P', color = [1.0 0.0 1.0]; r = 2.0735; charge = -0.15;
            case  'S', color = [0.0 0.0 1.0]; r = 2.0175; charge = -0.1;
            otherwise, color = [1.0 0.0 0.0]; r = 1.5; charge = 0;
%%%%%%%%%% WIKI VdW table
%             case  'H', color = [0.7 0.7 0.7]; r = 1.2; charge = 1;
%             case  'C', color = [0.3 0.3 1.0]; r = 1.7; charge = 6;
%             case  'O', color = [0.3 1.0 0.3]; r = 1.52; charge = 8;
%             case  'N', color = [1.0 0.3 1.0]; r = 1.55; charge = 7;
%             otherwise, color = [1.0 0.0 0.0]; r = 1.0;
% %%%%%%%%%%% MSMS VdW table
%             case  'H', color = [0.7 0.7 0.7]; r = 1.20; 
%             case  'C', color = [0.3 0.3 1.0]; r = 1.74; 
%             case  'O', color = [0.3 1.0 0.3]; r = 1.4; 
%             case  'N', color = [1.0 0.3 1.0]; r = 1.7; 
%             case  'P', color = [1.0 0.0 1.0]; r = 1.80; 
%             %case  'F', color = [1.0 0.0 1.0]; r = 1.35; 
%             %case  'Cl', color = [1.0 0.0 1.0]; r = 1.8; 
%             case  'S', color = [0.0 0.0 1.0]; r = 1.80; 
%             otherwise, color = [1.0 0.0 0.0]; r = 2; 
        end
        
        c =  sscanf(line(31:54),'%f %f %f');
        
        Geom.R(it) = r;
        Geom.centers(it,:) = c;
        Geom.T(it) = line(14);
        Geom.charges(it) = charge;
        if charge < 0
            charge_negative = charge_negative+charge;
        elseif charge > 0
            num_H = num_H+1;            
        end
        if charge == -0.2
            num_C = num_C+1;
        end
%        Geom.color(it,:) = color;
        
    end
    
    line = fgetl(fid);
end
fclose(fid);

Geom.M = it;
Geom.centers = Geom.centers(1:it,:);

if num_H ~= 0
    charge_H = -charge_negative/num_H;
    for i = 1:Geom.M
        if Geom.charges(i) == 1
            Geom.charges(i) = charge_H;
        end
    end
else
    charge_C = -charge_negative/num_C-0.2;
    for i = 1:Geom.M
        if Geom.charges(i) == -0.2
            Geom.charges(i) = charge_C;
        end
    end
end

for i = 1:Geom.M
    if Geom.charges(i) == 1
        Geom.charges(i) = charge_H;
    end
end

Data.charges = Geom.charges;
Data.x = Geom.centers;

end