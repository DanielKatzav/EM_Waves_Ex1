%{

- Intro. to Electromagnetics and Waves class
- First Matlab excercise
- Created by Daniel Katzav, ID: 203389770

%}



function [iteration, error] = finite_diff(Area, f_0, halt_condition, b)
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
    [ErrorStillBig, error]= CheckError(f_n, f_n_1, halt_condition);
    f_n = f_n_1;
    iteration = iteration + 1;
end

Area_S_potential = zeros(M,N);
for n = 1:N 
    for m = 1:M
        Area_S_potential(m, n) = f_n_1(m+(n-1)*M,1);
    end
end

imagesc(Area_S_potential)
colorbar
end


function [Proceed, error_rate] = CheckError(f_n, f_n_1, error_cond)
    
    error_rate = (norm(f_n_1 - f_n, 2) / norm(f_n, 2))
    if (error_rate < error_cond)
        Proceed = false;
    else
        Proceed = true;
    end
    
end


function R = CreateRMat(M, N)
    R = zeros(M*N);
    for i = 0:M*N-1
        for j = 0:M*N-1
            if (j == mod(i+1,N*M-1) || j == mod(i-1,M*N-1) || j == mod((i+1)*M,N*M-1) || j == mod((i-1)*M,M*N-1))
                R(i+1,j+1) = 1/4;
            end
        end

    end
end

















