%{

- Intro. to Electromagnetics and Waves class
- First Matlab excercise
- Created by Daniel Katzav, ID: 203389770
- Submission at 5/12/2019
- See PDF file attached

%}
clc
clear all
clear variables
close all

% --------- Question 1 --------- %
mm = 10^-3;
w_out = 10*mm;
h_out = 4*mm;

w_in = 0.8*mm;
d_in = 2*mm;

delta = 0.1*mm;

HighPotential = 0.5;
LowPotential = 0;

Left_Box_x0 = (h_out/2) - (w_in/2);
Left_Box_y0 = (w_out/2) - (d_in/2) - w_in;

Right_Box_x0 = (h_out/2) - (w_in/2);
Right_Box_y0 = (w_out/2) + (d_in/2);

% Create a new Area
Area_S_Raw = CreateArea(h_out, w_out, delta, delta, LowPotential);
% Remove the Left rectangular area, as shown in the sketch 
Area_S_Left_Removed = Disregard_Area(Area_S_Raw, w_in/delta, w_in/delta, Left_Box_x0/delta, Left_Box_y0/delta, HighPotential);
% Remove the Right rectangular area, as shown in the sketch
Area_S = Disregard_Area(Area_S_Left_Removed, w_in/delta, w_in/delta, Right_Box_x0/delta, Right_Box_y0/delta, -HighPotential);

phi_nm_0 = getPhiZero(Area_S);
[M,N] = size(Area_S);
Ns = M*N;

% Looks like I have a problem here at Matrix of R, not sure where, can't
% really find it, but I have another box below the left box that is being
% disregarded - to be mentioned on the summarize in the PDF
% also, it seems like the iterations number is beyond 3 times of the actual
% number that supposed to be, so will need to state that as well
[NumOfIterations, gamma]= finite_diff(Area_S, zeros(Ns,1), 10^-5, phi_nm_0);



% --------- Question 2 --------- %

epsilon_0 = 8.8541878128e-12;   % [F/m]
epsilon_1 = 2.5*epsilon_0;
epsilon_2 = 3.9*epsilon_0;
miu_0 = 4*pi*(10^-7);                % [H/m]
miu_1 = miu_0;
miu_2 = miu_0;

C1 = GetCapacitance(gamma, epsilon_1);
C2 = GetCapacitance(gamma, epsilon_2);
Z_C1 = getImpedance(C1, getInductance(miu_1, gamma));
Z_C2 = getImpedance(C2, getInductance(miu_2, gamma));

v_1 = getWaveVelocity(miu_1, epsilon_1);
v_2 = getWaveVelocity(miu_2, epsilon_2);



ID = 203389770;
R_g = Z_C1*(ID /(ID +ID)); % Just Re-Stating here, I have made the excercise alone

Gamma_g = (Z_C1 - R_g)/(Z_C1 + R_g);
Gamma_21 = (Z_C2 - Z_C1)/(Z_C2 + Z_C1);
Gamma_12 = -Gamma_21;

T = 8;
el_1 = T*v_1;
el_2 = T*v_2;
T_p = T/10;

dt = T/50;
dz = (el_1 + el_2)/1000;

%Axis and Matrices Definition
t = 0:dt:8*T;
z = 0:dz:el_1+el_2;

[Z_cord,T_cord]=meshgrid(z,t);  %Defining a set of coords for V(t,z)
%V is defined by the sum of V1_p, V1_m, V2_p
V1_p=zeros(size(T_cord));   %Forward travelling wave on Z_C1
V1_m=zeros(size(T_cord));   %Backward travelling wave on Z_C1
V2_p=zeros(size(T_cord));   %Forward travelling wave on Z_C2

%This loop is the sum of 1 to 4 of the Gamma parameters multiplied by the
%corresponding wave
for k=1:4
    V1_p(:,1:el_1/dz) = V1_p(:,1:el_1/dz)+(Z_C1/(Z_C1+R_g))*(Gamma_12*Gamma_g)^(k-1)*(heaviside(T_cord(:,1:el_1/dz)-Z_cord(:,1:el_1/dz)./v_1-2*T*(k-1))-heaviside(T_cord(:,1:el_1/dz)-T/10-Z_cord(:,1:el_1/dz)./v_1-2*T*(k-1)));
    if k~=4
        V1_m(:,1:el_1/dz) = V1_m(:,1:el_1/dz)+(Z_C1/(Z_C1+R_g))*Gamma_12^k*Gamma_g^(k-1)*(heaviside(T_cord(:,1:el_1/dz)+z(:,1:el_1/dz)./v_1-2*T*k)-heaviside(T_cord(:,1:el_1/dz)-T/10+Z_cord(:,1:el_1/dz)./v_1-2*T*k));
        V2_p(:,el_1/dz+1:(el_1+el_2)/dz) = V2_p(:,el_1/dz+1:(el_1+el_2)/dz)+(Z_C1/(Z_C1+R_g))*(1*Gamma_12)*(Gamma_12*Gamma_g)^(k-1)*(heaviside(T_cord(:,el_1/dz+1:(el_1+el_2)/dz)-Z_cord(:,el_1/dz+1:(el_1+el_2)/dz)./v_2+el_1/v_2-T-2*T*(k-1))-heaviside(T_cord(:,el_1/dz+1:(el_1+el_2)/dz)-T/10-Z_cord(:,el_1/dz+1:(el_1+el_2)/dz)./v_2+el_1/v_2-T-2*T*(k-1)));
    else
        V1_m(:,1:el_1/dz) = V1_m(:,1:el_1/dz);
        V2_p(:,el_1/dz+1:(el_1+el_2)/dz) = V2_p(:,el_1/dz+1:(el_1+el_2)/dz);
    end
end

V = V1_p + V1_m + V2_p;
figure(4)
imagesc(20*log(abs(V)))
colorbar
title('Propagation of 20log|V(t,z)|')
xlabel('z [mm]')
ylabel('t [sec/50]')

function NewArea = Disregard_Area(RawArea, height, width, start_y, start_x, LowPotential)
    % Remove the Inner parts of the disregarded area - i.e. set the values
    % there to -1
     RawArea(start_y:start_y+height, start_x:start_x+width) = -1;
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

    for j = 1:N
        for i = 1:M
            Mean_Neighbors = 0;
            count = 0;
             % Disregard the Area that we filterd 
             % whole Area and the outsides of the
            if (i+1 <= M && Area(i,j) ~= -1)
                Area(i+1,j)
                Mean_Neighbors = Mean_Neighbors+Area(i+1,j);
                count = count +1;
            end
            if (i-1 > 0 && Area(i,j) ~= -1)
                Area(i-1,j)
                Mean_Neighbors = Mean_Neighbors+Area(i-1,j);
                count = count +1;
            end
            if (j+1 <= N && Area(i,j) ~= -1)
                Area(i,j+1)
                Mean_Neighbors = Mean_Neighbors+Area(i,j+1);
                count = count +1;
            end 
            if (j-1 > 0 && Area(i,j) ~= -1)
                Area(i,j-1)
                Mean_Neighbors = Mean_Neighbors+Area(i,j-1);
                count = count +1;
            end
            if (count ~= 0)
                phi_zero(i + (j-1)*M,1) = (1/count) * Mean_Neighbors;
            end
        end
        
        Area_S_potential = zeros(M,N);
        for n = 1:N 
            for m = 1:M
                Area_S_potential(m, n) = phi_zero(m+(n-1)*M,1);
            end
        end

        imagesc(Area_S_potential)
        colorbar
    end

end

function C = GetCapacitance(gamma, epsilon)
    C = gamma*epsilon;
end
function Z = getImpedance(C, L)
    Z = sqrt(C/L);
end
function L = getInductance(miu, gamma)
   L = miu/gamma;
end
function v = getWaveVelocity(miu, epsilon)
    v = 1/sqrt(miu*epsilon);
end









