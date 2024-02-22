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

A = M1/M;
disp('Matrix A')
disp(A)