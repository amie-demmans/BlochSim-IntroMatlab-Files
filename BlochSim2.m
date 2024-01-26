%% Part a1: Matlab Vectors 
M = [1,0,0]';
disp('Matrix M')
disp(M)

%% Part a2: Transverse Relaxation
% A-2a) Magnetization M1 as only x-component 
T2 = 100;
t = 50;
Mx = exp(-t/T2);
M1 = M*Mx;
disp('Magnetization M1')
disp(M1)

% A-2b) M1 = A*M with A being 3x3 matrix
A = M1/M;
disp('Matrix A')
disp(A)

%% Part a3: Longitudinal Relaxation
% A-3a) Simulating only T1 relaxation 
%T1 = 600;
%M0 = 1;
%Mz = 1-exp(-t/T1);
%M2 = M*Mz;
%A2 = M2/M;
%disp('Matrix A for T1')
%disp(A2)
%B = M2-(A2*M);
%disp('Matrix B')
%disp(B))
