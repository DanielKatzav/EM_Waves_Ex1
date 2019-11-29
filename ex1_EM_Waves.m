%{

- Intro. to Electromagnetics and Waves class
- First Matlab excercise
- Created by Daniel Katzav, ID: 203389770
- Submission at 5/12/2019
- See PDF file attached

%}

clear all
close all
clc

mm = 10^-3;

w_out = 10*mm;
h_out = 4*mm;

w_in = 0.8*mm;
d_in = 2*mm;

delta = 0.1*mm;

HighPotential = 2;
LowPotential = 1;

Left_Box_x0 = (h_out/2) - (w_in/2);
Left_Box_y0 = (w_out/2) - (d_in/2) - w_in;

Right_Box_x0 = (h_out/2) - (w_in/2);
Right_Box_y0 = (w_out/2) + (d_in/2);

% Create a new Area
Area_S_Raw = CreateArea(h_out, w_out, delta, delta, HighPotential);
% Remove the Left rectangular area, as shown in the sketch 
Area_S_Left_Removed = Disregard_Area(Area_S_Raw, w_in/delta, w_in/delta, Left_Box_x0/delta, Left_Box_y0/delta, LowPotential);
% Remove the Right rectangular area, as shown in the sketch
Area_S = Disregard_Area(Area_S_Left_Removed, w_in/delta, w_in/delta, Right_Box_x0/delta, Right_Box_y0/delta, LowPotential);
phi_nm_0 = getPhiZero(Area_S);
[M,N] = size(Area_S);
Ns = M*N;


[NumOfIterations, finalerror]= finite_diff(Area_S, zeros(Ns,1), 10^-5, phi_nm_0);
disp(NumOfIterations);
disp(finalerror);



function NewArea = Disregard_Area(RawArea, height, width, start_y, start_x, LowPotential)
    % Remove the Inner parts of the disregarded area - i.e. set the values
    % there to 0
     RawArea(start_y:start_y+height, start_x:start_x+width) = 0;
     NewArea = RawArea;
    % Define the walls of the area to be with the Low Potential value
    for i = 0:width
        NewArea(start_y+i,start_x) = LowPotential;
        NewArea(start_y+i,start_x+width) = LowPotential;
    end
    for i = 0:height
        NewArea(start_y,start_x+i) = LowPotential;
        NewArea(start_y+height,start_x+i) = LowPotential;
    end
end

function mat = CreateArea(height, width, dx, dy, HigherPotential)
    N = ((width / dx)+1);
    M = ((height / dy)+1);
    
    % Define the Outer Walls of the Area to be with the Higher Potential
    % Value
    mat = HigherPotential*zeros(M,N);
    for i = 1:M
        mat(i,1) = HigherPotential;
        mat(i,N) = HigherPotential;
    end
    for i = 1:N
        mat(1,i) = HigherPotential;
        mat(M,i) = HigherPotential;
    end
    
end


function phi_zero = getPhiZero(Area)
    [M,N] = size(Area);
    phi_zero = zeros(M*N,1);

    for i = 1:M
        for j = 1:N
            Mean_Neighbors = 0;
             % Disregard the Area that we filterd 
             % whole Area and the outsides of the
            if (i+1 <= M && Area(i,j) ~= 0)
                Area(i+1,j)
                Mean_Neighbors = Mean_Neighbors+Area(i+1,j);
            end
            if (i-1 > 0 && Area(i,j) ~= 0)
                Area(i-1,j)
                Mean_Neighbors = Mean_Neighbors+Area(i-1,j);
            end
            if (j+1 <= N && Area(i,j) ~= 0)
                Area(i,j+1)
                Mean_Neighbors = Mean_Neighbors+Area(i,j+1);
            end 
            if (j-1 > 0 && Area(i,j) ~= 0)
                Area(i,j-1)
                Mean_Neighbors = Mean_Neighbors+Area(i,j-1);
            end
            
            phi_zero((i-1)*M + j,1) = 0.25 * Mean_Neighbors;
        end
        
        Area_S_potential = zeros(M,N);
        for j = 1:N 
            for k = 0:M-1
                Area_S_potential(k+1, j) = phi_zero(k*M+j,1);
            end
        end

        imagesc(Area_S_potential)
        colorbar
    end

end












