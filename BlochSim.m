%% Part A1: MATLAB Vectors  
M = [1,1,1]';    %Matrix here will be changing to account for different questions 
disp('Matrix M')
disp(M)

%% Part A2(a): Transverse Relaxation 
T2 = 100;
t = 50;
M1 = M*exp(-t/T2);
%A = [M1,0,0;0,M1,0;0,0,1];
disp ('Magnetization M1)')
disp(M1)
Mx = A*M;
disp('Matrix Mx')
disp(Mx)

%% Part A2(b): Transverse Relaxation (3x3 matrix) 
A = [M1,0,0;0,M1,0;0,0,1];
A = M1/M;
disp('Matrix A')
disp(A)

%% Part A3(a): Longitudinal Relaxation 
T1 = 600;
M0 = 1;
Mz = M0+(M-M0)*exp(-t/T1);
%columnIndex = 3;
%selectedColumn = A(:, columnIndex);
%B = Mz*selectedColumn;
A = (Mz-B)/M;
B = 
disp('Matrix B')
disp(B)
disp('Matrix A')
disp(A)

%% test
X = [1,0,0;0,1,0;0,0,1];
disp(X)



