%{

- Intro. to Electromagnetics and Waves class
- First Matlab excercise
- Created by Daniel Katzav, ID: 203389770

%}



function [iteration, gamma] = finite_diff(Area, f_0, halt_condition, b)
%{
Summary:
 FINITE_DIFF Iterative Algorithm to convolve the potential at a given
 cartisian point, by the mean of it's previous iteration 
 neighborspotential

Usage:
    Area - Matrix of MxN, that will represent the points of calulation
    f_0 - Starting condition of the requested iteration to calulate
    halt_condition - iteration solution convergance error rate
%}
[M,N] = size(Area);
Ns = M*N;

R_Matrix = CreateRMat(M, N);
iteration = 0;
ErrorStillBig = true;

f_n_1 = zeros(Ns,1);
f_n = f_0;

while(ErrorStillBig)
    f_n_1 = R_Matrix*f_n + b;
    [ErrorStillBig]= CheckError(f_n, f_n_1, halt_condition);
    f_n = f_n_1;
    iteration = iteration + 1;
end

figure(1)
Area_S_potential = reshape(f_n_1, [M,N]);
imagesc(Area_S_potential)
colorbar
title('Potential of Wave Conductor')
xlabel('y [mm]')
ylabel('x [mm]')

[x,y] = meshgrid(1:N, 1:M);
[ex,ey] = GetEModel(M, N, Area_S_potential);
figure(2)
quiver(x,y,ey,ex)
title('e of Wave Conductor')
xlabel('y')
ylabel('x')


gamma = GetGamma(M , N, ex, ey);
hy = ex/gamma;      % y-component of z*curl(grad(phi))/gamma
hx = -ey/gamma;     % x-component of z*curl(grad(phi))/gamma

figure(3)
quiver(x,y,hy,hx)
title('h of Wave Conductor')
xlabel('y')
ylabel('x')


end


function Proceed = CheckError(f_n, f_n_1, error_cond)
    
    error_rate = (norm(f_n_1 - f_n, 2) / norm(f_n, 2));
    if (error_rate < error_cond)
        Proceed = false;
    else
        Proceed = true;
    end
    
end
function R = CreateRMat(M, N)

    mm = 10^-3;
    w_out = 10*mm;
    h_out = 4*mm;
    w_in = 0.8*mm;
    d_in = 2*mm;
    delta = 0.1*mm;
    
    Left_Box_x0 = (h_out/2) - (w_in/2);
    Left_Box_y0 = (w_out/2) - (d_in/2) - w_in;

    Right_Box_x0 = (h_out/2) - (w_in/2);
    Right_Box_y0 = (w_out/2) + (d_in/2);
    % Left Inner square bottom and top limits
    Min_down_left = (Left_Box_x0)/(delta*2) + 1;     
    Min_up_left = (Left_Box_x0+w_in)/(delta*2) - 1; 
    
    Nin_down_left = (Left_Box_y0)/(delta*2) + 1;    
    Nin_up_left = (Left_Box_y0+w_in)/(delta*2) - 1;     
    
    % Right Inner square bottom and top limits
    Min_down_right = (Right_Box_x0)/(delta*2) + 1;     
    Min_up_right = (Right_Box_x0+w_in)/(delta*2) - 1; 
    
    Nin_down_right = (Right_Box_y0)/(delta*2) + 1;    
    Nin_up_right = (Right_Box_y0+w_in)/(delta*2) - 1;     
    
    R = zeros(M*N);
    Vect_MN = zeros(M*N ,1);      
    W = reshape(Vect_MN, [M,N]);
    
for j = 2 : N - 1
    for i = 2 : M - 1
        if j >= Nin_down_left-1 && j <= Nin_up_left+1 && i >= Min_down_left-1 && i <= Min_up_left+1
           continue; 
        end
        
        if j >= Nin_down_right-1 && j <= Nin_up_right+1 && i >= Min_down_right-1 && i <= Min_up_right+1
           continue; 
        end
        b(M*(j-1)+i) = 0.25*(W(i+1,j)+W(i-1,j)+W(i,j+1)+W(i,j-1));
        
        if i == 2 && j == 2
            R(M*(j-1)+i,M*(j-1)+i+1)=1/4;   %up
            R(M*(j-1)+i,M*j+i)=1/4;         %right
        
        
        elseif i == M-1 && j == 2
            R(M*(j-1)+i,M*(j-1)+i-1)=1/4;   %down
            R(M*(j-1)+i,M*j+i)=1/4;         %right
        
        
        elseif i == 2 && j == N-1
            R(M*(j-1)+i,M*(j-1)+i+1)=1/4;   %up
            R(M*(j-1)+i,M*(j-2)+i)=1/4;     %left
        
        
        elseif i == M-1 && j == N-1
            R(M*(j-1)+i,M*(j-1)+i-1)=1/4;   %down
            R(M*(j-1)+i,M*(j-2)+i)=1/4;     %left
        
        
        elseif (j == 2 && (i >= 3 && i <= M-2)) || (j == Nin_up_left+2 && (i >= Min_down_left-1 && i <= Min_up_left+1))
            %next to left edges
            R(M*(j-1)+i,M*(j-1)+i+1)=1/4;   %up
            R(M*(j-1)+i,M*(j-1)+i-1)=1/4;   %down
            R(M*(j-1)+i,M*j+i)=1/4;         %right
        
        elseif (j == N-1 && (i >= 3 && i <= M-2)) || (j == Nin_down_left-2 && (i >= Min_down_left-1 && i <= Min_up_left+1))
            %next to right edges
            R(M*(j-1)+i,M*(j-1)+i+1)=1/4;   %up
            R(M*(j-1)+i,M*(j-1)+i-1)=1/4;   %down
            R(M*(j-1)+i,M*(j-2)+i)=1/4;     %left
            
        elseif (i == 2 && (j >= 3 && j <= N-2)) || (i == Min_up_left+2 && (j >= Nin_down_left-1 && j <= Nin_up_left+1))
            %next to bottom edges
            R(M*(j-1)+i,M*(j-1)+i+1)=1/4;   %up
            R(M*(j-1)+i,M*(j-2)+i)=1/4;     %left
            R(M*(j-1)+i,M*j+i)=1/4;         %right
        
        
        elseif (i == M-1 && (j >= 3 && j <= N-2)) || (i == Min_down_left-2 && (j >= Nin_down_left-1 && j <= Nin_up_left+1))
            %next to top edges
            R(M*(j-1)+i,M*(j-1)+i-1)=1/4;   %down
            R(M*(j-1)+i,M*j+i)=1/4;         %right
            R(M*(j-1)+i,M*(j-2)+i)=1/4;     %left
            
                    elseif (j == 2 && (i >= 3 && i <= M-2)) || (j == Nin_up_right+2 && (i >= Min_down_right-1 && i <= Min_up_right+1))
            %next to left edges
            R(M*(j-1)+i,M*(j-1)+i+1)=1/4;   %up
            R(M*(j-1)+i,M*(j-1)+i-1)=1/4;   %down
            R(M*(j-1)+i,M*j+i)=1/4;         %right
        
        elseif (j == N-1 && (i >= 3 && i <= M-2)) || (j == Nin_down_right-2 && (i >= Min_down_right-1 && i <= Min_up_right+1))
            %next to right edges
            R(M*(j-1)+i,M*(j-1)+i+1)=1/4;   %up
            R(M*(j-1)+i,M*(j-1)+i-1)=1/4;   %down
            R(M*(j-1)+i,M*(j-2)+i)=1/4;     %left
            
        elseif (i == 2 && (j >= 3 && j <= N-2)) || (i == Min_up_right+2 && (j >= Nin_down_right-1 && j <= Nin_up_right+1))
            %next to bottom edges
            R(M*(j-1)+i,M*(j-1)+i+1)=1/4;   %up
            R(M*(j-1)+i,M*(j-2)+i)=1/4;     %left
            R(M*(j-1)+i,M*j+i)=1/4;         %right
        
        
        elseif (i == M-1 && (j >= 3 && j <= N-2)) || (i == Min_down_right-2 && (j >= Nin_down_right-1 && j <= Nin_up_right+1))
            %next to top edges
            R(M*(j-1)+i,M*(j-1)+i-1)=1/4;   %down
            R(M*(j-1)+i,M*j+i)=1/4;         %right
            R(M*(j-1)+i,M*(j-2)+i)=1/4;     %left
            
        else
            %everything else
            R(M*(j-1)+i,M*(j-1)+i+1)=1/4;   %up
            R(M*(j-1)+i,M*(j-1)+i-1)=1/4;   %down
            R(M*(j-1)+i,M*j+i)=1/4;         %right
            R(M*(j-1)+i,M*(j-2)+i)=1/4;     %left
        end
    end
end
end
function [ex,ey] = GetEModel(M,N, phi_i)

    mm = 10^-3;
    w_out = 10*mm;
    h_out = 4*mm;
    w_in = 0.8*mm;
    d_in = 2*mm;
    delta = 0.1*mm;
    
    Left_Box_x0 = (h_out/2) - (w_in/2);
    Left_Box_y0 = (w_out/2) - (d_in/2) - w_in;

    Right_Box_x0 = (h_out/2) - (w_in/2);
    Right_Box_y0 = (w_out/2) + (d_in/2);
    % Left Inner square bottom and top limits
    Min_down_left = (Left_Box_x0)/(delta*2) + 1;     
    Min_up_left = (Left_Box_x0+w_in)/(delta*2) - 1; 
    
    Nin_down_left = (Left_Box_y0)/(delta*2) + 1;    
    Nin_up_left = (Left_Box_y0+w_in)/(delta*2) - 1;     
    
    % Right Inner square bottom and top limits
    Min_down_right = (Right_Box_x0)/(delta*2) + 1;     
    Min_up_right = (Right_Box_x0+w_in)/(delta*2) - 1; 
    
    Nin_down_right = (Right_Box_y0)/(delta*2) + 1;    
    Nin_up_right = (Right_Box_y0+w_in)/(delta*2) - 1;     


ey = zeros(size(phi_i));     % y-component of grad(phi)
ex = zeros(size(phi_i));     % x-component of grad(phi)

for j = 2 : N - 1
    for i = 2 : M - 1
        if j >= Nin_down_left-1 && j <= Nin_up_left+1 && i >= Min_down_left-1 && i <= Min_up_left+1 
           ey(i,j) = 0;
           ex(i,j) = 0;
           continue; 
        end
        if j >= Nin_down_right-1 && j <= Nin_up_right+1 && i >= Min_down_right-1 && i <= Min_up_right+1 
           ey(i,j) = 0;
           ex(i,j) = 0;
           continue; 
        end
        
        ex(i,j) = -(phi_i(i+1,j) - phi_i(i-1,j))/(2*delta);      % Taking the negative of grad(phi) to get e
        ey(i,j) = -(phi_i(i,j+1) - phi_i(i,j-1))/(2*delta);
    end
end




end
function gamma = GetGamma(M,N, ex, ey)
    mm = 10^-3;
    w_out = 10*mm;
    h_out = 4*mm;
    delta = 0.1*mm;
    % Finding gamma using numerical the trapzoid integration method
    dx = delta*2:delta:w_out-delta;       %integration segment for x-axis with trapz
    dy = delta*2:delta:h_out-delta;       %integration segment for y-axis with trapz
    gamma = trapz(dy, ey(3:M-1,2)) + trapz(dx, -ex(M-1,3:N-1)) + trapz(dy, -ey(3:M-1,N-1)) + trapz(dx,ex(2,3:N-1));

end








