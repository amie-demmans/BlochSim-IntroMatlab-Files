%% Part A-1: Matlab Vectors  
%  ------------------------- %
%% A-1): Defining equilibrium magnetization vector 

M = [1,0,0]';
disp('Matrix M')
disp(M)

%% Part A-2: Transverse Relaxation  
%  -------------------------------- %
%% A-2a) Magnetization vector due to T2 decay after 50ms (x-component only)
 
T2 = 100; % Transverse relaxation time (ms)
t = 50;   % time (ms)

Mx =(exp(-t/T2)); % x-component magnetization vector
M1 = M * Mx; % apply equilibrium magnetization vector  
disp('Magnetization M1')
disp(M1)

%% A-2b) Express M1 = A*M with A being 3x3 matrix

A = [Mx,0,0;0,Mx,0;0,0,1]; % State matrix A for transverse magnetization where Mx is that stated in A-2a
disp('Matrix A for T2')
disp(A)
M1_Transverse = A * M; % Apply equilibrium magentization vector to matrix A for transverse magnetization matrix 
disp ('Transverse Magnetization')
disp(M1_Transverse)

%% Part A-3: Longitudinal Relaxation  
%  ---------------------------------- %
%% A-3a) Expressing T1 relaxation in matrix form M1 = A*M+B

T1 = 600;        % Longitudinal relaxation time (ms)
M0 = 1;          % non-zero value
Mz = exp(-t/T1); % z-component magentization vector 

Mz1 = M0*Mz; % apply non-zero component to Mz to determine A
Mz2 = M0*(1-Mz); % apply non-zero component to 1-Mz to determine B 

A2 = [1,0,0;0,1,0;0,0, Mz1]; % State matrix A for longitudinal relaxation 
disp('Matrix A for T1') 
disp(A2)

B = [0,0,Mz2]'; % State 3x1 vector B 
disp('Matrix B')
disp(B)

M1_Longitudinal = (A2*M)+B; % Express T1 relaxation in matrix form 
disp('Longitudinal Magnitization')
disp(M1_Longitudinal)

%% A-3b) Combining T1 and T2 relaxation 

A3 = [Mx,0,0;0,Mx,0;0,0,Mz1]; % State matrix A including both transverse and longitudinal relaxation 
disp('Matrix A for Combination T1 and T2')
disp(A3)

M1_combined = (A3*M)+B; % Express T1 and T2 relaxation in matrix form 
disp('Combined Magnatization')
disp(M1_combined)

%% Part A-4: Rotations: Precession and Excitation  
%  ----------------------------------------------- %
%% A-4a) Simulating precession (rotation about z-axis)
% By adding 'd' at the end of sin,cos,tan, tells MATLAB to use degrees  

phi = 90; % flip angle (degrees)
Rz = [cosd(phi),-sind(phi),0;sind(phi),cosd(phi),0;0,0,1]; % State matrix Rz for simulating precession 
disp('Matrix Rz')
disp(Rz)

%% A-4b) using z-rot function to simulate precession

phi = pi/2;     % flip angle (radians)
Rz = zrot(phi); % call function zrot
disp('Matrix Rz using zrot function')
disp(Rz)

%% A-4c) Simulate excitation using xrot and yrot (rotation about x and y)

Rx = xrot(phi); % call function xrot
disp('Matrix Rx')
disp(Rx)

Ry = yrot(phi); % call function yrot 
disp('Matrix Ry')
disp(Ry)

%% A-4d) Rotation about transverse axis other than x and y 

phi2 = pi/4;  % flip angle (radians)
theta = pi/6; % axis rotation angle (radians) 

Rth = throt(phi2,theta); % call the function 
disp('Matrix Rth')
disp(Rth)

%% Part A-5: Free Precession
%  -------------------------- %
%% A-5a/b) Create function combining relaxation and precession effects then plot Mx,My,Mz vs t 

dt = 1;            % Time step for plotting (ms)
t = 1000;          % Response time (ms)
N = round(t/dt)+1; % number of time stops (round values up to nearest positive integer) --> could also use ceil to round up or down 
df = 10;           % off-resonance frequency (Hz)
T1 = 600;          % Longitudinal relaxation time (ms)
T2 = 100;          % Transverse relaxation time (ms)  

[Afp,Bfp]=freeprecess(dt,T1,T2,df); % call the function 

M = zeros(3,N);  % Generate 3xN matrix of zeros 
M(:,1)=[1,0,0]'; % Set first column of matrix

for k=2:N
	M(:,k) = Afp*M(:,k-1)+Bfp; % for each loop iteration, update the k-th column of M 
end

% Plot Mx, My and Mz 
figure ('Name','FID for a Single Species','NumberTitle', 'off')
time = (0:N-1)*dt;
plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;

%% Part B-1: Saturation Recovery 
%  ------------------------------- %
%% B-1a&b) Magnetization after first and second excitation assuming no off-resonance

df = 0;      % Off resonance frequency (Hz)
T1 = 600;    % Longitudinal relaxation time (ms)
T2 = 100;    % Transverse  relaxation time (ms)
TR = 500;    % Repetition time (ms)
TE = 1;      % Echo time (ms)
flip = pi/3; % Flip angle (radians) 

%Initialixing Matrices
M = [0,0,1]';       % Equilibrium magnetization vector 
Rflip = yrot(flip); % flip matrix with rotation in y-axis

[Ate,Bte] = freeprecess(TE,T1,T2,df); % free precession matrix for time = TE
[Atr,Btr] = freeprecess(TR,T1,T2,df); % free precession matrix for time = TR

% Now do the calculations
% Magnetization 1ms (TE) after first excitation 
M = Rflip*M;   %applying flip angle 
M = Ate*M+Bte; % applying free precession over 1ms (TE) 
disp('Magnetization after first excitation')
disp(M)

%Magnetization at time TE after second excitation (TR first then TE) 
M = [0,0,1]';  % restate initial magnetization
M = Rflip*M;   % applying flip angle 
M = Atr*M+Btr; % applying free precession over 500ms (TR)
M = Rflip*M;   % apply flip angle again (second pulse) 
M = Ate*M+Bte; % apply free precession 1ms (TE) after second excitation
disp('Magnetization after second excitation')
disp(M)

%% B-1c) Plotting Mx, My, and Mz vs time every 1ms for 5s starting at equilibrium

df = 0;             % Off resonance frequency (Hz)
dt = 1;             % Time step for plotting (ms)
T1 = 600;           % Longitudinal relaxation time (ms)
T2 = 100;           % Transverse relaxation time (ms)
TR = 500;           % Repetition time (ms)
flip = pi/3;        % Flip angle (radians)
NTR = round(TR/dt); % Number of repitition times
excitations = 10;   % Number of excitations

M = [0,0,1]';       % Equilibrium magnetization
Rflip = yrot(flip); % Flip matrix with rotation around y-axis

[Afp,Bfp]=freeprecess(dt,T1,T2,df); %free precession matrix for plotting with time step dt

count = 1; % Initializing count to 1

for n = 2 : excitations  
    M(:,count) = Rflip*M(:,count); % for each loop iteration, the column vector M(:,count) is multiplied by flip angle  

    for k = 1: NTR       
        count = count+1;           % increment count by 1 for each loop iteration
        M(:,count) = Afp*M(:,count-1)+Bfp; % for each loop iteration, the column vector M(:,count) is updated using Afp*M(:,count-1)+Bfp
    end

end

% Plot Mx, My and Mz starting at equilibrium
figure ('Name','Magnetization over first 10 excitations','NumberTitle', 'off')
time = (0:count-1)*dt;
plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;

%% B-1d) Steady state magnetization at the point just prior to the exciation

T1 = 600;     % Longitudinal Relaxation time (ms)
T2 = 100;     % Transverse Relaxation time (ms) 
TR = 500;     % Repition time (ms) 
flip = pi/3;  % Flip angle (radians)

Rflip = yrot(flip); % flip matrix with rotation along y-axis
[Atr,Btr]=freeprecess(TR,T1,T2,df); % free precession matrix for t = TR 

% steady state magnetization at point prior to excitation 
% --> apply the inverse of a 3x3 identity matrix subtracted by the Atr*Rflip
% --> multiply that inverse by Btr to get a new column vector 

M = inv(eye(3)-Atr*Rflip)*Btr; 
disp('Steady State Magnetization at the point prior to excitation')
disp(M)

%% B-1e) Steady state magnetization at echo time (TE = 1ms)
df = 0;      % Off resonance frequency (Hz)
flip = pi/3; % Flip angle (radians)
T1 = 600;    % Longitudinal relaxation time (ms)
T2 = 100;    % Transverse relaxation time (ms)
TE = 1;      % Echo time (ms)
TR = 500;    % Repitition time (ms) 

%Call the function sssignal to display steady state magnetization @ TE=1ms
Mss = sssignal(flip,T1,T2,TE,TR,df);
disp('Steady state magnetization at TE = 1ms')
disp(Mss)

%% B-1f) Calculate the result of B-1d on paper by neglecting the transverse magnetization prior to excitation 
  % Calculations all done on paper but the result is...

Mss = [0,0,0.7224]';
disp('Steady State Magnetization Neglecting transverse magentization prior to excitation')
disp(Mss) 

%% B-1g) Multiply matrix Atr by what is given to get same answer as B-1f

df = 0;      % Off resonance frequency (Hz)
flip = pi/3; % Flip angle (radians)
T1 = 600;    % Longitudinal relaxation time (ms)
T2 = 100;    % Transverse relaxation time (ms)
TR = 500;    % Repitition time (ms) 

%Call the function srsignal to display steady state magnetization
% --> Neglects transverse magnetization

Mss = srsignal(flip,T1,T2,TR,df);
disp('Steady state magnetization keeping transverse term')
disp(Mss)

%% Part B-2: Spin-Echo Sequences  
%  ------------------------------- %
%% B-2a) Plotting magnetization components as a function of time

df = 10;  % Off resonance frequency (Hz) 
t = 1000; % Response time 
T1 = 600; % Longitudinal relaxation time (ms) 
T2 = 100; % Transverse relaxation time (ms)
TR = 500; % Repitition time (ms) 
TE = 50;  % echo time (ms) 
dt = 1;   % time step (ms)

NTE = round((TE/2)/dt);  %number of echo times
NTR = round(TR-TE/2/dt); %number of repitition times 

[Afp,Bfp]=freeprecess(dt,T1,T2,df); %free precession matrix for plotting with time step dt

M = zeros(3,NTE+NTR); % generate a matrix with 3 rows and NTE+NTR columns 
M(:,1)=[0,0,1]'; %initial mangetization vector 

Rflip=yrot(pi/2); % create your flip matrix with rotation around y
Refocus=xrot(pi); % create your refocusing pulse matrix with rotation around x 

M(:,2) = Afp*Rflip*M(:,1)+Bfp; % computes magnetization after first echo which involves 90 degree pulse 

for k=3:(NTE+1) 
	M(:,k) = Afp*M(:,k-1)+Bfp; % updates magnetization at TE using freeprecess 
end

M(:,NTE+2) = Afp*Refocus*M(:,NTE+1)+Bfp; % computes magnetization after refocusing pulse (180 degree pulse) 

for k=2:(NTR-1) 
    M(:,k+NTE+1) = Afp*M(:,k+NTE)+Bfp; % updates magnetization at TR using freeprecess 
end

% Plot magnetization for all spin-echo signals 
figure ('Name','Magnetization as a Function of Time for Spin-Echo','NumberTitle', 'off')
time = (0:NTE+NTR-1)*dt;
plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
axis([min(time) max(time) -1 1]);
grid on;

%% B-2b) Plotting mangitude, phase and complex mean as a function of time 

t = 1000; % Response time 
T1 = 600; % Longitudinal relaxation time (ms) 
T2 = 100; % Transverse relaxation time (ms)
TR = 500; % Repitition time (ms) 
TE = 50;  % Echo time (ms) 
dt = 1;   % Time step (ms)

NTE = round((TE/2)/dt); %number of echo times
NTR = round(TR/dt);     %number of repitition times 

spins = 5; % number of spins 

% Store magnitude and phase information for each spin in the simulation 
Mag = zeros(spins,NTE+NTR); % matrix sotring magnitude of complex numbers
phase = zeros(spins,NTE+NTR); 

for freq = 1:spins
    df = rand*100-50; % randmonizing off-resonance frequency between -50 to 50 Hz

[Afp,Bfp]=freeprecess(dt,T1,T2,df); %free precession matrix for plotting with time step dt

M = zeros(3,NTE+NTR); % generate a matrix with 3 rows and NTE+NTR columns 
M(:,1)=[0,0,1]'; %initial mangetization vector 

Rflip=yrot(pi/2); % create your flip matrix with rotation around y
Refocus=xrot(pi); % create your refocusing pulse matrix with rotation around x 

M(:,2) = Afp*Rflip*M(:,1)+Bfp; % computes magnetization after first echo which involves 90 degree pulse 

for k=3:(NTE+1) 
	M(:,k) = Afp*M(:,k-1)+Bfp; % updates magnetization at TE using freeprecess 
end

M(:,NTE+2) = Afp*Refocus*M(:,NTE+1)+Bfp; % computes magnetization after refocusing pulse (180 degree pulse) 

for k=2:(NTR-1) 
    M(:,k+NTE+1) = Afp*M(:,k+NTE)+Bfp; % updates magnetization at TR using freeprecess 
end

Magnitude(freq,:) = (M(1,:)+1i*M(2,:));


end 

% Plot Magnitude and Phase of all Signals as a function of time 
figure ('Name','Magnitude and Phase of the Spins as a Function of time','NumberTitle', 'off')
time = (0:NTE+NTR-1)*dt;
subplot(2, 1, 1);
plot(time, abs(Magnitude'));
title('Magnitude of Spins');
xlabel('Time (ms)');
ylabel('Magnitude');
grid on;

subplot(2, 1, 2);
plot(time, angle(Magnitude'));
title('Phase of Spins');
xlabel('Time (ms)');
ylabel('Phase (radians)');
grid on;
legend('Spin 1', 'Spin 2', 'Spin 3', 'Spin 4', 'Spin 5', 'Spin 6', 'Spin 7', 'Spin 8', 'Spin 9', 'Spin 10');

% Plot Magnitude of the complex mean of all signals
figure ('Name','Complex Mean of Signals as a Function of time','NumberTitle', 'off')
time = (0:NTE+NTR-1)*dt;
plot(time, abs(mean(Magnitude)));
title('Magnitude of Complex Mean');
xlabel('Time (ms)');
ylabel('Magnitude');
grid on;

%% B-2c) Create function similar to sssignal but includes 180 degree pulse 

df = rand*100-50; % Off-resonance frequency (Hz)
T1 = 600;         % Longitudinal relaxation time (ms)
T2 = 100;         % Transverse relaxation time (ms) 
TE = 50;          % Echo time (ms) 
TR = 1000;        % Repitition time (ms) 

Mss = sesignal(T1,T2,TE,TR,df); % call the function
disp('Mx for spin echo pulse')
disp(Mss)

%% B-2d) Fast spin echo (FSE) sequence signal amplitudes 

T1 = 600;  % Longitudinal relaxation time (ms) 
T2 = 100;  % Transverse relaxation time (ms)
TE = 50;   % Echo time (ms) 
TR = 1000; % Repitition time (ms) 
ETL = 8;   % Echo-train length / total number of spin echoes 
df = 0;    % Off-resonance frequency (Hz) 

[Msig] = fsesignal(T1,T2,TE,TR,df,ETL);
disp('Fast Spin Echo Signal Amplitudes')
disp(Msig)

%% Part B-3) Gradient-Spoiled Sequences 
% -------------------------------------- %
%% B-3a) Simulate magnetizations in each voxel for GE sequences 

flip = pi/3; % Flip angle (radians)
T1 = 600;    % Longitudinal relaxation time (ms) 
T2 = 100;    % Transverse relaxation time (ms) 
TE = 2;      % Echo time (ms) 
TR = 10;     % Repitition time (ms) 
df = 0;      % Off resonance frequency (Hz) 
phi = pi/2;  % Dephasing angle (radians)

 Mss = gssignal(flip,T1,T2,TE,TR,df,phi); % call the function 
disp('Magnetizations in each voxel for GE sequences')
disp(Mss)

%% B-3b) Average magnetization when a gradient is placed at the end of TR

flip = pi/3; % Flip angle (radians)
T1 = 600;    % Longitudinal relaxation time (ms) 
T2 = 100;    % Transverse relaxation time (ms) 
TE = 2;      % Echo time (ms) 
TR = 10;     % Repitition time (ms) 
df = 0;      % Off resonance frequency (Hz) 

Mss = gresignal(flip,T1,T2,TE,TR,df); % Call the function
disp('Average Magnetization for GE')
disp(Mss)

%% B-3c) Changing input variable for srsignal 

flip = pi/3; % Flip angle (radians)
T1 = 600;    % Longitudinal relaxation time (ms) 
T2 = 100;    % Transverse relaxation time (ms) 
TE = 2;      % Echo time (ms) 
TR = 10;     % Repitition time (ms) 
df = 0;      % Off resonance frequency (Hz) 

Mss = srsignal(flip,T1,T2,TE,TR,df);
disp('Magnetization using srsignal for GE sequences')
disp(Mss)

%% Part B-4) Steady-State Free-Precession 
% ---------------------------------------- %
%% B-4a) Plotting Signal magnitude and phase for refocused SSFP as a function of frequency 

flip = pi/3;     % Flip angle 
T1 = 600;        % Longitudinal relaxation time (ms)
T2 = 100;        % Transverse relaxation time (ms) 
TR = 10;         % Repetition times (ms) 
TE = 0:2.5:10;   % Echo time (ms) --> Repeat TE from 0-10 in increments on 2.5
df = (-100:100); % Off resonance frequency (Hz) 


Sig = zeros(length(df),length(TE)); % Create matrix of zeros df by TE

for n = 1:length(TE)  

    for k = 1:length(df) 
    Mss = sssignal(flip,T1,T2,TE(n),TR,df(k)); % Call sssignal but with reference to TE(n) and df(k) as they are ranges of values
    Sig(k,n) = Mss(1)+1i*Mss(2); % Apply complex number representation to retain both magntiude and phase info in one variable 
    end 

end 

figure ('Name','Refocused SSFP','NumberTitle', 'off')

% Plot the magnitude
subplot(2, 1, 1);
plot(df, abs(Sig));
title('Magnitude of Spins');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
legend('TE = 0', 'TE = 2.5', 'TE = 5', 'TE = 7.5', 'TE = 10');

% Plot the phase
subplot(2, 1, 2);
plot(df, angle(Sig));
title('Phase of Spins');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
axis ([min(df) max(df) -pi pi]);
grid on;
legend('TE = 0', 'TE = 2.5', 'TE = 5', 'TE = 7.5', 'TE = 10');

%% B-4b) Plotting signal magnitude and phase for gradient spoiled sequence

flip = pi/3;     %Flip angle 
T1 = 500;        % Longitudinal relaxation time (ms)
T2 = 100;        % Transverse relaxation time (ms) 
TR = 10;         % Repetition times (ms) 
TE = 0:2.5:10;   % Echo time (ms) --> Repeat TE from 0-10 in increments on 2.5
df = (-100:100); % Off resonance frequency (Hz) 


Sig = zeros(length(df),length(TE)); % Create matrix of zeros df by TE

for n = 1:length(TE)

    for k = 1:length(df) 
    Mss = gresignal(flip,T1,T2,TE(n),TR,df(k)); % Call sssignal but with reference to TE(n) and df(k) as they are ranges of values
    Sig(k,n) = Mss(1)+1i*Mss(2); % Apply complex number representation to retain both magntiude and phase info in one variable 
    end 

end 


figure ('Name','Gradient Spoiled Sequence','NumberTitle', 'off')

% Plot the magnitude
subplot(2, 1, 1);
plot(df, abs(Sig));
title('Magnitude of Spins');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
legend('TE = 0', 'TE = 2.5', 'TE = 5', 'TE = 7.5', 'TE = 10');

% Plot the phase
subplot(2, 1, 2);
plot(df, angle(Sig));
title('Phase of Spins');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
axis ([min(df) max(df) -pi pi]);
grid on;
legend('TE = 0', 'TE = 2.5', 'TE = 5', 'TE = 7.5', 'TE = 10');

%% B-4c) Repeating B-4a but with TE = TR/2 

flip = pi/3;     % Flip angle 
T1 = 600;        % Longitudinal relaxation time (ms)
T2 = 100;        % Transverse relaxation time (ms) 
TR = [2, 6, 10]; % Repetition times (ms) --> Repeat TR from 0-10 in increments of 2
df = (-500:500); % Off resonance frequency (Hz) 


Sig = zeros(length(df),length(TE)); % Create matrix of zeros df by TE

for n = 1:length(TR) 

    for k = 1:length(df) 
    Mss = sssignal(flip,T1,T2,TR(n)/2,TR(n),df(k)); % Call sssignal but with reference to TE(n) and df(k) as they are ranges of values
    Sig(k,n) = Mss(1)+1i*Mss(2); % Apply complex number representation to retain both magntiude and phase info in one variable 
    end 

end 

figure ('Name','Refocused SSFP with TE = TR/2','NumberTitle', 'off')

% Plot the magnitude
subplot(2, 1, 1);
plot(df, abs(Sig));
title('Magnitude of Spins');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
legend ('TR = 2', 'TR = 6',  'TR = 10');

% Plot the phase
subplot(2, 1, 2);
plot(df, angle(Sig));
title('Phase of Spins');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
axis ([min(df) max(df) -pi pi]);
grid on;
legend( 'TR = 2',  'TR = 6',  'TR = 10');

%% B-4d) Plotting Magnitude and phase using phase increment on RF pulse

flip = pi/3;              % Flip angle 
T1 = 600;                 % Longitudinal relaxation time (ms)
T2 = 100;                 % Transverse relaxation time (ms) 
TR = 5;                   % Repetition times (ms) 
TE = 2.5;                 % Echo time (ms) --> Repeat TE from 0-10 in increments on 2.5
df = (-500:500);          % Off resonance frequency (Hz) 
phi = [0,pi/2,pi,1.5*pi]; % Dephasing angle increment

Sig = zeros(length(df),length(TE)); % Create matrix of zeros df by TE

for n = 1:length(phi) 

    for k = 1:length(df) 
    Mss = ssfp(flip,T1,T2,TE,TR,df(k),phi(n)); % Call sssignal but with reference to TE(n) and df(k) as they are ranges of values 
    Sig(k,n) = Mss(1)+1i*Mss(2); % Apply complex number representation to retain both magntiude and phase info in one variable 
    end

end 

figure ('Name','GE Sequence with "on-resonance" magnetization','NumberTitle', 'off')

% Plot the magnitude
subplot(2, 1, 1);
plot(df, abs(Sig));
title('Magnitude of Spins');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
legend('phi = 0', 'phi = pi/2', 'phi = pi', 'phi = 3pi/2');

% Plot the phase 
subplot(2, 1, 2);
plot(df, angle(Sig));
title('Phase of Spins');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');
axis ([min(df) max(df) -pi pi]);
grid on;
legend('phi = 0', 'phi = pi/2', 'phi = pi', 'phi = 3pi/2');

%% Part B-5) RF-Spoiled Sequences  (Confused by this whole section)
% ---------------------------------------- %
%% B-5a) Signal magnitude and phase as a function of excitation number (my way)
%T1 = 600; % T1 relaxation time (ms) 
%T2 = 100; % T2 relaxation time (ms) 
%TR = 10; % repitition time (ms) 
%TE = 2; % echo time (ms) 
%flip = pi/6; % flip angle (radians) 
%df = 0; % off-resonance frequency (Hz) 

%increment = 117/180*pi; % Increment of excitations (switched to radians by dividing ny 180degrees*pi)
%excitations = 100; % Number of excitations 

%magnitude = zeros(1,excitations); % Create matrix of zeroes 1 x excitations 
%phase = zeros(1, excitations); % Create matrix of zeroes 1 x excitations 

%RFphase = 0; %initial excitation pulse



% Simulate the signal for each pulse 
%for N = 1:excitations 
 %   RFphase = mod(phi+(N-1)*increment,360);
 %   [Mss] = sssignal(flip,T1,T2,TE,TR,df);
 %   magnitude = abs(Mss);
 %   phase = angle(Mss);
 %   phi = phi-RFphase;
%end

%time = (0:excitations)*TR+TE;
%figure ('Name','GE Sequence with "on-resonance" magnetization','NumberTitle', 'off')
%subplot(2, 1, 1);
%plot(1:excitations, magnitude);
%title('Magnitude of Spins');
%xlabel('Excitation');
%ylabel('Magnitude');
%grid on;
%legend('phi = 0', 'phi = pi/2', 'phi = pi', 'phi = 3pi/2');

%subplot(2, 1, 2);
%plot(1:excitations,phase);
%title('Phase of Spins');
%xlabel('Frequency (Hz)');
%ylabel('Phase (radians)');
%axis ([min(df) max(df) -pi pi]);
%grid on;

%% B-5a) Signal magntiude and phase as a function of excitation number (bloch sim answer)  

df = 0;		      % Off-resonance frequency (Hz)
T1 = 600;	      % Longitudinal relaxation time (ms)
T2 = 100;	      % Transverse relaxation time (ms) 
TE = 2;		      % Echo time (ms) 
TR = 10;          % Repitition time (ms)
flip = pi/6;	  % Flip angle (radians) --> listed as alpha on pulse sequence diagram
inc = 117/180*pi; % Increment of excitations (switched to radians by dividing ny 180degrees*pi)

Nex = 100; % Number of excitations 

Nf = 100;	% Simulate 50 different gradient-spoiled spins.
phi = (1:Nf)/Nf*2*pi;  % phase of excitation

M=zeros(3,Nf,Nex+1); % Generate a 3 x Nf array of zeros Nex+1 times
Msig = zeros(Nex);   % Generate a Nex x Nex array of zeros 


% Call free precess function for TR and TE 	
[Ate,Bte] = freeprecess(TE,T1,T2,0);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,0);

M = [zeros(2,Nf);ones(1,Nf)]; % Create 3 x Nf matrix with first two rows filled with zeros, and last row with 1
on = ones(1,Nf); % Create a 1 x Nf row vector filled with ones 
	
Rfph = 0;  % Initial excitation pulse 
Rfinc = inc; % Redefine the increment 

for n=1:Nex 

	A = Ate * throt(flip,Rfph); % Set A to matrix A of freeprecess for t = TE applied to throt
	B = Bte; % Restate Bte 
	M = A*M+B*on; % Express in matrix form 

	Msig(n) = mean( squeeze(M(1,:)+1i*M(2,:)) ) * exp(-1i*Rfph);
	
	M=Atr*M+Btr*on;

	for k=1:Nf
		M(:,k) = zrot(phi(k))*M(:,k);
    end

	Rfph = Rfph+Rfinc;
	Rfinc = Rfinc+inc;
end


time = (0:Nex-1)*TR+TE;
subplot(2,1,1);
plot(time,abs(Msig));
xlabel('Time (ms)');
ylabel('Magnitude');
grid on;


subplot(2,1,2);
plot(time,angle(Msig));
xlabel('Time (ms)');
ylabel('Phase');
grid on;

%% B-5b) Redoing part B-5a but defining a function 
df = 0;		% off-resonance frequency (Hz)
T1 = 600;	% T1 relaxation time (ms)
T2 = 100;	% T2 relaxation time (ms) 
TE = 2;		% echo time (ms) 
TR = 10;    % repitition time (ms)
flip = (0:0.01:0.5)*pi;	% flip angle (radians) --> listed as alpha on pulse sequence diagram
inc = 117/180*pi; % Increment of excitations (switched to radians by dividing ny 180degrees*pi)
Nex = 100;

for k = 1:length(flip)
    sig1(k) = spgrsignal(flip(k),T1,T2,TE,TR,df,Nex,inc);
    [Mss] = srsignal(flip(k),T1,T2,TE,TR,df);
    sig2(k) = M(1)+1i*M(2);
end 

figure ('Name','RF spoiled GRE and SR','NumberTitle', 'off')
plot (flip, abs(sig1), 'r-', flip, abs(sig2), 'b--');
xlabel('Flip angle (radians)');
ylabel('Signal Magntitude');
grid on;
title ('Signal Amplitude vs Flip Angle for RF-spoiled GRE and SR');
legend('RF-spoiled GRE', 'Saturation-Recovery Approximation');

%% Part C-1) Signal vs Sequence Parameters 
% ---------------------------------------- %
%% C-1a) Plotting spin-echo signal level as a function of TR

figure ('Name','Spin-Echo Signal Level (TR)','NumberTitle', 'off')

% TE = 50ms, TR = 100-4000ms
fplot(@(x)abs(sesignal(600,100,50,x,0))',[100,4000]); 

xlabel('Repitition Time (ms)')
ylabel('Spin-Echo Signal Level')
title ('Spin Echo Signal level vs Repitition time')
grid on;

%% C-1b) Plotting spin-echo signal level as a function of TE

figure ('Name','Spin-Echo Signal Level (TE)','NumberTitle', 'off')

% TE = 0-500ms, TR = 1000ms
fplot(@(x)abs(sesignal(600,100,x,1000,0))',[0,500]); 

xlabel('Echo Time (ms)')
ylabel('Spin-Echo Signal Level')
title ('Spin Echo Signal level vs Repitition time')
grid on;

%% C-1c) Plotting Ratio of signal level to sqaure-root of TR

figure ('Name','Ratio of Spin-Echo Signal Level','NumberTitle', 'off')

% Dividing function by sqrt(TR) for best SNR efficiency (fplot from C-1a)
% SNR is proportional to sqrt of readout time 
% SNR = signal-to-noise ratio 

fplot(@(x)abs(sesignal(600,100,50,x,0))/sqrt(x)',[100,4000]); 

xlabel('Repitition Time (ms)')
ylabel('Spin-Echo Signal Level')
title ('Spin Echo Signal level vs Repitition time')
grid on;

%% Part C-2) Contrast 
% -------------------- % 
%% C-2a) Plotting Signal and Signal Difference for Tissue A and B (TR)

figure ('Name','Signal and Signal Difference at TR for Tissue A and B','NumberTitle', 'off')

% Tissue A (T1 = 600ms, T2 = 100ms, TE = 50ms, TR = x)
fplot(@(x)abs(sesignal(600,100,50,x,0))',[50,4000]);
hold on;

% Tissue B (T1 = 1000ms, T2 = 150ms, TE = 50ms, TR = x)
fplot(@(x)abs(sesignal(1000,150,50,x,0))',[50,4000]);

% Tissue A - Tissue B 
fplot(@(x)abs(sesignal(600,100,50,x,0))-abs(sesignal(1000,150,50,x,0))',[50,4000]);

xlabel('Repitition Time (ms)');
ylabel('Spin-Echo Signal Level');
title ('Spin Echo Signal level vs Repitition time');
legend ('Tissue A', 'Tissue B', 'TA-TB');
grid on;

%% C-2b) Plotting Signal and Signal Difference for Tissue A and B (TE) 

figure ('Name','Signal and Signal Difference at TE for Tissue A and B','NumberTitle', 'off')

% Tissue A (T1 = 600ms, T2 = 100ms, TE = x, TR = 1000ms)
fplot(@(x)abs(sesignal(600,100,x,1000,0))',[0,500]); 
hold on; 

%Tissue B (T1 = 1000ms, T2 = 150ms, TE = x, TR = 1000ms)
fplot(@(x)abs(sesignal(1000,150,x,1000,0))',[0,500]); 

% Tissue A - Tissue B
fplot(@(x)abs(sesignal(600,100,x,1000,0))-abs(sesignal(1000,150,x,1000,0))', [0,500]);

xlabel('Echo Time (ms)');
ylabel('Spin-Echo Signal Level');
title ('Spin Echo Signal level vs Echo time');
legend ('Tissue A', 'Tissue B', 'TA-TB');
grid on;

%% C-2c) Pure T2 contrast by increasing TR from C-2b

figure ('Name','Pure T2 Contrast','NumberTitle', 'off')

% Tissue A (T1 = 600ms, T2 = 100ms, TE = x, TR = 4000ms)
fplot(@(x)abs(sesignal(600,100,x,4000,0))',[0,500]); 
hold on; 

%Tissue B (T1 = 1000ms, T2 = 150ms, TE = x, TR = 4000ms)
fplot(@(x)abs(sesignal(1000,150,x,4000,0))',[0,500]); 

% Tissue A - Tissue B
fplot(@(x)abs(sesignal(600,100,x,4000,0))-abs(sesignal(1000,150,x,4000,0))', [0,500]);

xlabel('Echo Time (ms)');
ylabel('Spin-Echo Signal Level');
title ('Spin Echo Signal level vs Echo time');
legend ('Tissue A', 'Tissue B', 'TA-TB');
grid on;

%% C-2d) Normalize all three quantities with sqrt(TR) from C-2a (contrast efficiency)

figure ('Name','Contrast Efficiency','NumberTitle', 'off')

% Tissue A (T1 = 600ms, T2 = 100ms, TE = 50ms, TR = x)
fplot(@(x)abs(sesignal(600,100,50,x,0))/sqrt(x)',[50,4000]);
hold on;

% Tissue B (T1 = 1000ms, T2 = 150ms, TE = 50ms, TR = x)
fplot(@(x)abs(sesignal(1000,150,50,x,0))/sqrt(x)',[50,4000]);

% Tissue A - Tissue B 
fplot(@(x)abs(sesignal(600,100,50,x,0))/sqrt(x)-abs(sesignal(1000,150,50,x,0))/sqrt(x)',[50,4000]);

xlabel('Repitition Time (ms)');
ylabel('Spin-Echo Signal Level');
title ('Spin Echo Signal level vs Repitition time');
legend ('Tissue A', 'Tissue B', 'TA-TB');
grid on;

%% Part C-3) Multiple Spin Echo Sequences 
% ---------------------------------------- %
%% C-3a) Contrast-to-noise ratio (CNR) efficiency 

df = 0; % Off-resonance frequency (Hz) 

T1A = 600;  % Longitudinal relaxation time for tissue A (ms) 
T2A = 100;  % Transverse relaxation time for tissue A (ms) 

T1B = 1000; % Longitudinal relaxation time for tissue B (ms) 
T2B = 150;  % Transverse relaxation time for tissue B (ms)

TR = 10:100:5000; % Example TR values 
TE = 0:20:1000;   % Example TE values 

% Create matrix of zeros TE x TR
CNR_efficiency = zeros(length(TE), length(TR));

% Generate for loop to initialize matrix and iterate over TE and TR 
for k = 1:length(TE)

    for i = 1:length(TR)

        if (TE(k)>TR(i)) % Ensure CNR effieiciency is 0 if TE > TR
            CNR_efficiency(k,i) = 0;
        else % Calculate CNR efficiency 
           CNR_efficiency(k,i)= abs(abs(sesignal(T1A,T2A,TE(k),TR(i),df))-abs(sesignal(T1B,T2B,TE(k),TR(i),df)))/sqrt(TR(i));
        end

    end

      tt = sprintf('%d %% complete.', round(100 * k / length(TE))); % Completion percentage during loop
      %disp(tt); % Display percentage if you wish 
end

% Max size of colormap 
Cmx = max(size(colormap));

% Normalize CNR_efficiency to [0,Cmx]
Cp = CNR_efficiency - min(CNR_efficiency(:));
Cp = Cmx*Cp/max(abs(Cp(:)));

% Display colorbar on side of graph
% --> I am also doing this when I plot with 'colorbar' but they are different 
%Cs = (1:length(TE))' * Cmx/length(TE) * ones(1,round(length(TR)/20));
%Cp = [Cp Cs];

% Plot the CNR efficiency 
figure ('Name','Contrast-to-Noise Efficiency','NumberTitle', 'off')
imagesc(TR,TE,Cp);
colorbar;
xlabel('TR(ms)');
ylabel('TE(ms)');
title('CNR efficiency as a function of TR and TE');

%% Max TR and TE for T1-contrast sequence
% max TR = 375ms 
% max TE = 0ms 

%% C-3b) % Max TR and TE for T2-contrast sequence  

% I am a bit confused on how to read the efficiency plot
% Answers written above and below are the answers from bloch sim 

% max TR = 3000ms 
% max TE = 130ms

%% C-3c) Plotting Signal-to-Noise Ratio (SNR) Efficiency 
df = 0; % Off-resonance frequency (Hz) 
 
T1A = 600;  % Longitudinal relaxation time for tissue A (ms) 
T2A = 100;  % Trasverse relaxation time for tissue A (ms) 

T1B = 1000; % Longitudinal relaxation time for tissue B (ms) 
T2B = 150;  % Transverse relaxation time for tissue B (ms)

TR = 10:100:5000; % Example TR values 
TE = 0:20:1000;   % Example TE values 

% Create matrix of zeros TE x TR for Tissue A and B respectively 
SNR_efficiencyA = zeros(length(TE), length(TR));
SNR_efficiencyB = zeros(length(TE), length(TR));
  
% For SNR, Tissue A and B are calculated separately 
for k = 1:length(TE)

    for i = 1:length(TR)

        if (TE(k)>TR(i)) % Ensure CNR effieiciency is 0 if TE > TR for both tissues 
            SNR_efficiencyA(k,i) = 0;
            SNR_efficiencyB(k,i) = 0;
        else % Calculate CNR efficiency for Tissue A and Tissue B 
           SNR_efficiencyA(k,i) = abs(abs(sesignal(T1A,T2A,TE(k),TR(i),df))/sqrt(TR(i)));
           SNR_efficiencyB(k,i) = abs(abs(sesignal(T1B,T2B,TE(k),TR(i),df))/sqrt(TR(i)));
        end

    end
      tt = sprintf('%d %% complete.', round(100 * k / length(TE))); % Completion percentage during loop
      %disp(tt); % Display percentage if you wish     
end

% Max size of colormap
Cmx = max(size(colormap));

% Normalize CNR_efficiencyA to [0,Cmx]
CpA = SNR_efficiencyA - min(SNR_efficiencyA(:));
CpA = Cmx*CpA/max(abs(CpA(:)));

%Normalize CNR_efficiencyB to [0,Cmx]
CpB = SNR_efficiencyB-min(SNR_efficiencyB(:));
CpB = Cmx*CpB/max(abs(CpB(:)));

%Cs = (1:length(TE))' * Cmx/length(TE) * ones(1,round(length(TR)/20));

%Cp = [Cp Cs];

% Plotting SNR Efficiency for Tissue A
figure ('Name','Signal-to-Noise Efficiency Tissue A','NumberTitle', 'off')
imagesc(TR,TE,CpA);
colorbar;
xlabel('TR(ms)');
ylabel('TE(ms)');
title('CNR efficiency as a function of TR and TE for Tissue A');

% Plotting SNR Efficiency for Tissue B 
figure ('Name','Signal-to-Noise Efficiency Tissue B','NumberTitle', 'off')
imagesc(TR,TE,CpB);
colorbar;
xlabel('TR(ms)');
ylabel('TE(ms)');
title('CNR efficiency as a function of TR and TE for Tissue B');

%% C-3d) Plotting SNR Efficiency as a function of Echo Train Length (ETL)

df = 0;     % Off-resonance frequency (Hz) 
ETL = 1:30; % Echo-train length  

T1A = 600;  % Longitudinal relaxation time for tissue A (ms) 
T2A = 100;  % Transverse relaxation time for tissue A (ms) 

TR = 3000;  % Repitition time (ms) 
TE = 15;    % Echo time (ms) 

for k = 1:length(ETL)
    SNR_efficiency(k) = sum(abs(fsesignal(T1A,T2A,TE,TR,df,ETL(k))))/sqrt(ETL(k))/sqrt(TR);
end

% Plot SNR Efficiency vs ETL
figure ('Name','Signal-to-Noise Efficiency as a Function of Echo train length ','NumberTitle', 'off')
plot(ETL, SNR_efficiency);
xlabel('Echo-Train Length');
ylabel('Signal-to-noise ratio efficiency');
title('SNR efficiency as a function og ETL');

%% C-3e) Plotting SNR Efficiency of Tissue A and B and CNR efficiency as a function of ETL
% --> Will need to clear all on command line for this to run
% --> Something in above code is messing with it 

df = 0;       % Off-resonance frequency (Hz) 
ETL = (1:60); % Echo-Train length 

T1A = 600;    % Longitudinal relaxation time for tissue A (ms) 
T2A = 100;    % Transverse relaxation time for tissue A (ms) 

T1B = 1000;   % Longitudinal relaxation time for tissue B (ms) 
T2B = 150;    % Transverse relaxation time for tissue B (ms)

TR = 3000;    % Repitition time (ms)  
TE = 15;      % Echo time (ms) 

% Calculate SNR efficiency for Tissue A and B respectively 
for k = 1:length(ETL)
    SNR_efficiencyA(k) = sum(abs(fsesignal(T1A,T2A,TE,TR,df,ETL(k))))/sqrt(ETL(k))/sqrt(TR);
    SNR_efficiencyB(k) = sum(abs(fsesignal(T1B,T2B,TE,TR,df,ETL(k))))/sqrt(ETL(k))/sqrt(TR);
    %CNR_efficiency(k)= abs(abs(fsesignal(T1A,T2A,TE,TR,df,ETL(k)))-abs(sesignal(T1B,T2B,TE,TR,df,ETL(k))))/sqrt(TR);
end

% Calculate CNR efficiency using SNR efficiency from Tissue A and B 
CNR_efficiency = SNR_efficiencyB - SNR_efficiencyA;

figure ('Name','Signal-to-Noise Efficiency as a Function of Echo train length ','NumberTitle', 'off')

% Plot SNR efficiency of tissue A
plot (ETL, SNR_efficiencyA);
xlabel('Echo train length');
ylabel('SNR Efficiency');

grid on;
title('SNR Efficiency for Tissue A and B vs ETL');

hold on;  % allows all the data for each efficiency to be plotted on same graph 

% Plot SNR efficiency of tissue B and CNR efficiencu 
plot (ETL, SNR_efficiencyB);
plot (ETL, CNR_efficiency);
legend('Tissue A', 'Tissue B', 'CNR Efficiency');

%% C-3f) Reading Exercise 
% --> Looking at FSE knee cartilage mri image 

%% C-3g) Determining T2 in each pixel 
% T2 = TE/ln(0.95(exp(TE/T2*))= 119ms 

%% Part F-1) Basic Selective Excitation 
% -------------------------------------- %
%% F-1a) Plotting transverse magnetization immediately after second delta pulse 

phi = [pi/4 pi/4]; % tip angle (radians)
df = (-500:500);   % off-resonance frequency (Hz) 
sig = 0*df;        % allocate space?? 
T1 = 600;          % Longitudinal relaxation time (ms)
T2 = 100;          % Transverse relaxation time (ms) 


for n = 1:length(df) % Loop will iterate from 1 to entire range of df 
    Mt = [0,0,1]'; % initial magnetization 
    dt = 2.3; % time stamp between delta pulses (ms) 
    Mt = throt(abs(phi(1)),angle(phi(1))) * Mt; % Apply initial magnetization to throt

    for k = 2:length(phi) % Loop will iterate from 2 to entire range of phi 
        [Afp,Bfp]=freeprecess(dt,T1,T2,df(n)); %free precession matrix for plotting with time step dt
        Mt = Afp*Mt+Bfp; % apply free precession at Mt
        Mt = throt(abs(phi(k)),angle(phi(k))) * Mt; % Apply precessed Mt to throt 
    end 
    sig(n) = Mt(1)+1i*Mt(2); % Apply complex number representation to retain both phase and magnitude info 
end 

figure ('Name','Transverse Magnetization Immediately After Second Delta Pulse ','NumberTitle', 'off')

% Plotting magnitude of transverse magentization 
subplot(2,1,1);
plot(df,abs(sig));
xlabel('Frequency (Hz)');
ylabel('Signal (Fraction of M0)');
grid on;

% Plotting phase of transverse magnetization 
subplot(2,1,2);
plot(df,angle(sig));
xlabel('Frequency (Hz)');
ylabel('Signal (Radians)');
grid on;

%% F-1b) Plotting signal amplitude as a function of position, assuming no variation in resonant frequency 
% --> Turning on a gradient 

phi = [pi/4 pi/4]; % tip angle (radians)
x = (-20:20);      % Position range (mm) 
grad = 0.1;        % Gradient strength (G/cm)
sig = 0*x;         % allocate space?? 
T1 = 600;          % Longitudinal relaxation time (ms)
T2 = 100;          % Transverse relaxation time (ms) 
df = 0;            % Off-resonance frequency (Hz) 

for n = 1:length(x)
    Mt = [0,0,1]'; % initial magnetization 
    dt = 2.3; % time stamp between delta pulses (ms) 
    Mt = throt(abs(phi(1)),angle(phi(1))) * Mt; % Apply initial magnetization to throt

    for k = 2:length(phi)
        [Afp,Bfp]=freeprecess(dt,T1,T2,df); %free precession matrix for plotting with time step dt
        Mt = Afp*Mt+Bfp; %apply free precession at Mt 
        grad_rotation = 4258*2*pi*(dt/1000)*grad*(x(n)/10); % Apply rotation from gradient 
        %grad_position = gamma*2*pi*(dt/1000)*grad*(x(n)/10); % Why does it not work if I have it written as gamma instead of 4258??
        %grad_poistion = 4258 * grad * x(n);

        % Apply RF rotation 
        Mt = zrot(grad_rotation)*Mt; 
        Mt = throt(abs(phi(k)),angle(phi(k)))*Mt;
    end 

    sig(n) = Mt(1)+1i*Mt(2); % Apply complex number representation to retain both phase and magnitude info 

end 

figure ('Name','Applying a Gradient to Plot Signal Amplitude vs Position ','NumberTitle', 'off')

% Plot phase of transverse magnetization with a gradient 
subplot(2,1,1);
plot(x,abs(sig));
xlabel('Position range (mm)');
ylabel('Signal (Fraction of M0)');
grid on;

% Plot magnitude of transverse magnetization with a gradient 
subplot(2,1,2);
plot(x,angle(sig));
xlabel('Position range (mm)');
ylabel('Signal (Radians)');
grid on;

%% F-1c) Apply half the negative gradient area after the last RF pulse 
% --> G = Gauss, 1 Tesla = 10^4 Gauss

phi = [pi/4 pi/4]; % tip angle (radians)
x = (-20:20);      % Position range (mm) 
grad = 0.1;        % Gradient strength (G/cm)
sig = 0*x;         % allocate space?? 
T1 = 600;          % T1 relaxation time (ms)
T2 = 100;          % T2 relaxation time (ms) 
df = 0;            % Off-resonance frequency (Hz) 

for n = 1:length(x)
    Mt = [0,0,1]'; % initial magnetization 
    dt = 2.3; % time stamp between delta pulses (ms) 
    Mt = throt(abs(phi(1)),angle(phi(1))) * Mt; % Apply initial magnetization to throt

    for k = 2:length(phi)
        [Afp,Bfp]=freeprecess(dt,T1,T2,df); % free precession matrix for plotting with time step dt
        Mt = Afp*Mt+Bfp; % apply free precession at Mt 
        grad_rotation = 4258*2*pi*(dt/1000)*grad*(x(n)/10); % Apply rotation from gradient 

        % Apply RF rotation 
        Mt = zrot(grad_rotation)*Mt;
        Mt = throt(abs(phi(k)),angle(phi(k)))*Mt;
    end 
     [Afp,Bfp]=freeprecess(dt,T1,T2,df); % free precession matrix for plotting with time step dt
        Mt = Afp*Mt+Bfp; % apply free precession at Mt 
        grad_rotation = 4258*2*pi*(dt/1000)*(-0.5*grad)*(x(n)/10); % applying half the negative gradient rotation after last RF pulse 
        
        Mt = zrot(grad_rotation)*Mt;
       
    sig(n) = Mt(1)+1i*Mt(2); % apply complex number representation to retain both phase and magnitude info 
end 

figure ('Name','Applying a Gradient to Plot Signal Amplitude vs Position ','NumberTitle', 'off')

subplot(2,1,1);
plot(x,abs(sig));
xlabel('Frequency (Hz)');
ylabel('Signal (Fraction of M0)');
grid on;


subplot(2,1,2);
plot(x,angle(sig));
xlabel('Frequency (Hz)');
ylabel('Signal (Radians)');
grid on;

%% Part F-2) The Hard Pulse Approximation
% ---------------------------------------- %
%% F-2a) Plotting discrete rotations vs time (no gradient applied) 

gamma = 4258;   % gyromagnetic ratio (Hz/Gauss) 
% gamma = 4258e4; % gyromagnetic ratio (Hz/Tesla)
t = (0:0.1:6);  % Time B1 field is applied with sampling interval of 0.1 (ms)
dt = t(2)-t(1); % Time stamp

% Define B1(t)  
B1 = 0.06*sinc((t-3)/1); % 00.6 in Gauss, t-3/1 in ms 
% B1 = 0.06e4*sinc((t-3)/1); % 0.06 in Tesla 
% Why does it not work when I convert to Tesla? angle very large 

figure ('Name','Discrete Rotations vs Time for the RF pulse ','NumberTitle', 'off')
rotations = gamma*2*pi*B1*(dt/1000);
 

stem (t,rotations); % Using stem to show the bubbles of each point--> WHY??? 
xlabel('time (ms)');
ylabel('Flip Angle (Radians)')

flip = sum(rotations)*180/pi;
flip = sprintf('flip is %f degrees',round(flip));
disp('Flip angle for on-resonant mag with no applied gradient')
disp(flip)

% I am getting a flip angle of 92 degrees, using 0.06G as stated in quest
% Ans gives flip = 82 but they have constant as 0.05 

%% F-2b) Repeating 2a but with sampling interval of 4 us

gamma = 4258;    % gyromagnetic ratio (Hz/Gauss) 
% gamma = 4258e4; % gyromagnetic ratio (Hz/Tesla)
t = (0:0.004:6); % Time B1 field is applied with sampling interval of 0.1 (ms)
dt = t(2)-t(1);  % Time stamp 

% Define B1(t)  
B1 = 0.06*sinc((t-3)/1); % 00.6 in Gauss, t-3/1 in ms 
% B1 = 0.06e4*sinc((t-3)/1); % 0.06 in Tesla 
% Why does it not work when I convert to Tesla? angle very large 

figure ('Name','Discrete Rotations vs Time for the RF pulse ','NumberTitle', 'off')
rotations = gamma*2*pi*B1*(dt/1000);
 

stem (t,rotations); % Using stem to show the bubbles of each point 
xlabel('time (ms)');
ylabel('Flip Angle (Radians)')

flip = sum(rotations)*180/pi;
flip = sprintf('flip is %f degrees',round(flip));
disp('Flip angle for on-resonant mag with no applied gradient')
disp(flip)

% Angle stays the same as before but much smaller amplitudes (x10^-3)
% Area under curves are also filled in 

%% F-2c) Simulating the hard-pulse approximation 

T1 = 600; 
T2 = 100; 

gamma = 4258;   % gyromagnetic ratio (Hz/Gauss) 
% gamma = 4258e4; % gyromagnetic ratio (Hz/Tesla)
t = (0:0.04:6); % Time B1 field is applied with sampling interval of 0.04 (ms)
dt = t(2)-t(1); % time stamp

B1 = 0.06*sinc((t-3)/1); % Define B1 field 
rotations = gamma*2*pi*B1*(dt/1000);

df = (-1000:1000); % Off-resonance frequency 
sig = 0*df;

for n = 1:length(df)
    M = [0,0,1]'; % Initial magnetization
    M = throt(abs(rotations(1)),angle(rotations(1)))*M; % Apply RF rotations

    for k = 2:length(rotations)
        [Afp,Bfp]=freeprecess(dt,T1,T2,df(n)); %free precession matrix for plotting with time step dt
        M = Afp*M+Bfp; %apply free precession at Mt 
        M = throt(abs(rotations(k)),angle(rotations(k)))*M; % Apply RF rotations
    end 
    sig(n) = M(1)+1i*M(2); % Apply complex number representation to retain both magnitude and phase info 
end

figure ('Name','Simulating Hard Pulse approximation','NumberTitle', 'off')

% Plotting magnitude
subplot(2,1,1);
plot(df,abs(sig));
xlabel('Frequency (Hz)');
ylabel('Signal (Fraction of M0)');
grid on;

% Plotting phase
subplot(2,1,2);
plot(df,angle(sig));
xlabel('Frequency (Hz)');
ylabel('Phase (Radians)');
grid on;

%% Part F-2) Spatially Selective Simulations
% ------------------------------------------- %
%% F-3a) Gradient x0.5 --> RF x1 --> Gradient x0.5

T1 = 600;                     % Longitudinal relaxation time 
T2 = 100;                     % Transverse relaxation time 
t = (0:0.0001:0.006);         % Sample time (s) 
x = (-20:0.1:20);             % positions (mm) 
B1 = 0.06*sinc((1000*t-3)/1); % B1 field 
df = 0;                       % Off-resonance frequency (Hz) 


[msig] = sliceprofile(B1,0.1*ones(size(t)),t,T1,T2,x,df); % Calling the function 

figure ('Name','Simulating Hard Pulse approximation using sliceprofile.m','NumberTitle', 'off')

% Plotting magnitude 
subplot(2,1,1);
plot(x,abs(msig));
xlabel('Position (mm)');
ylabel('Magnitude');
grid on;

% Plotting phase 
subplot(2,1,2);
plot(x,angle(msig));
xlabel('Position (mm)');
ylabel('Phase');
grid on;

%% F-3b) Using F-3a, appens a 'refocusing' gradient to gradient waveform 
% --> Refocusing gradient bring dephased spins into coherence 
% --> Refocus by reversing direction of spins 

T1 = 600;             % Longitudinal relaxation time 
T2 = 100;             % Transverse relaxation time 
df = 0;               % Off-resonance frequency (Hz) 
t = (0:0.0001:0.006); % Sample time (s) 
x = (-20:0.1:20);     % positions (mm) 

B1 = 0.06*sinc((1000*t-3)/1); % B1 field 
grad = 0.1*ones(size(t)); % gradient 
phase_flat = -0.5; % Refocusing phase factor/ratio
grad = [grad grad*phase_flat]; % refocusing gradient 
t = [t t+0.006]; % shift time to account for refocusing gradient 
B1 = [B1 0*B1]; % making RF zero for refocusing gradient 

[m,msig] = sliceprofile(B1,grad,t,T1,T2,x,df); % Call the function 

figure ('Name','Applying Refocusing Gradients (functions of position)','NumberTitle', 'off')

% Plot the magntidue 
subplot(1,3,1);
plot(x,abs(msig));
xlabel('Position (mm)');
ylabel('Magnitude');
title('Magnitude');
grid on;

% Plot the phase 
subplot(1,3,2);
plot(x,angle(msig));
xlabel('Position (mm)');
ylabel('Phase');
title('Phase');
grid on;

% Plot z-component of magnetization 
subplot(1,3,3);
plot(x,m(3,:));
xlabel('Position (mm)');
ylabel('Z-component of Magnetization');
title('Z-component of Magnetization');
grid on;

figure ('Name','Applying Refocusing Gradients (functions of time)','NumberTitle', 'off')

% Plot RF Waveform 
subplot(1,2,1);
plot(t,B1);
xlabel('time (s)');
ylabel('RF Waveform');
title('RF Waveform');
grid on;

% Plot Gradient Waveform 
subplot(1,2,2);
plot(t,grad);
xlabel('time (s)');
ylabel('Gradient Waveform ');
title('Gradient Waveform');
grid on;

%% F-3c) Determining Best refocusing to slice gradient ratio to get the phase flat 


T1 = 600;             % Longitudinal relaxation time 
T2 = 100;             % Transverse relaxation time 
df = 0;               % Off-resonance frequency (Hz) 
t = (0:0.0001:0.006); % Sample time (s) 
x = (-20:0.1:20);     % positions (mm) 

B1 = 0.06*sinc((1000*t-3)/1); % B1 field 
grad = 0.1*ones(size(t)); % gradient 
phase_flat = -0.52; % I am not sure if I am mathematically supposed to find this, but I tests with different values from 0.48-0.53 and ~0.52 was game the flatest option 
grad = [grad grad*phase_flat]; % refocusing gradient 
t = [t t+0.006]; % shift time to account for refocusing gradient 
B1 = [B1 0*B1]; % making RF zero for refocusing gradient 

[m,msig] = sliceprofile(B1,grad,t,T1,T2,x,df); % Call the function 

figure ('Name','Applying Refocusing Gradients (functions of position)','NumberTitle', 'off')

% Plot the magntidue 
subplot(1,3,1);
plot(x,abs(msig));
xlabel('Position (mm)');
ylabel('Magnitude');
title('Magnitude');
grid on;

% Plot the phase 
subplot(1,3,2);
plot(x,angle(msig));
xlabel('Position (mm)');
ylabel('Phase');
title('Phase');
grid on;

% Plot z-component of magnetization 
subplot(1,3,3);
plot(x,m(3,:));
xlabel('Position (mm)');
ylabel('Z-component of Magnetization');
title('Z-component of Magnetization');
grid on;

figure ('Name','Applying Refocusing Gradients (functions of time)','NumberTitle', 'off')

% Plot RF Waveform 
subplot(1,2,1);
plot(t,B1);
xlabel('time (s)');
ylabel('RF Waveform');
title('RF Waveform');
grid on;

% Plot Gradient Waveform 
subplot(1,2,2);
plot(t,grad);
xlabel('time (s)');
ylabel('Gradient Waveform ');
title('Gradient Waveform');
grid on;

%% F-3d) Shift resonance freq in F-3c to 100 Hz

T1 = 600;             % Longitudinal relaxation time 
T2 = 100;             % Transverse relaxation time 
df = 100;             % Off-resonance frequency (Hz) 
t = (0:0.0001:0.006); % Sample time (s) 
x = (-20:0.1:20);     % positions (mm) 

B1 = 0.06*sinc((1000*t-3)/1); % B1 field 
grad = 0.1*ones(size(t)); % gradient 
phase_flat = -0.52; % Refocusing phase factor/ratio
grad = [grad grad*phase_flat]; % refocusing gradient 
t = [t t+0.006]; % shift time to account for refocusing gradient 
B1 = [B1 0*B1]; % making RF zero for refocusing gradient 

[m,msig] = sliceprofile(B1,grad,t,T1,T2,x,df); % Call the function 

figure ('Name','Applying Refocusing Gradients (functions of position)','NumberTitle', 'off')

% Plot the magntidue 
subplot(1,3,1);
plot(x,abs(msig));
xlabel('Position (mm)');
ylabel('Magnitude');
title('Magnitude');
grid on;

% Plot the phase 
subplot(1,3,2);
plot(x,angle(msig));
xlabel('Position (mm)');
ylabel('Phase');
title('Phase');
grid on;

% Plot z-component of magnetization 
subplot(1,3,3);
plot(x,m(3,:));
xlabel('Position (mm)');
ylabel('Z-component of Magnetization');
title('Z-component of Magnetization');
grid on;

figure ('Name','Applying Refocusing Gradients (functions of time)','NumberTitle', 'off')

% Plot RF Waveform 
subplot(1,2,1);
plot(t,B1);
xlabel('time (s)');
ylabel('RF Waveform');
title('RF Waveform');
grid on;

% Plot Gradient Waveform 
subplot(1,2,2);
plot(t,grad);
xlabel('time (s)');
ylabel('Gradient Waveform ');
title('Gradient Waveform');
grid on;

% Slice profile is shifted to the left a few mm 

%% F-3e) Modulate RF pulse by 900Hz exponential 

T1 = 600;             % Longitudinal relaxation time 
T2 = 100;             % Transverse relaxation time 
df = 0;               % Off-resonance frequency (Hz) 
t = (0:0.0001:0.006); % Sample time (s) 
x = (-50:1:50);       % positions (mm) 

RFmodulate = exp(2*pi*900*i*t); % Factor to modulate RF pulse 
B1 = 0.06*sinc((1000*t-3)/1); % B1 field 
B1 = B1.*RFmodulate; % Modulating the RF pulse --> want B1 just as an element ('.') not as matrix 
grad = 0.1*ones(size(t)); % gradient 
phase_flat = -0.52; % Refocusing phase factor/ratio
grad = [grad grad*phase_flat]; % refocusing gradient 
t = [t t+0.006]; % shift time to account for refocusing gradient 
B1 = [B1 0*B1]; % making RF zero for refocusing gradient 

[m,msig] = sliceprofile(B1,grad,t,T1,T2,x,df); % Call the function 

figure ('Name','Applying Refocusing Gradients (functions of position)','NumberTitle', 'off')

% Plot the magntidue 
subplot(1,3,1);
plot(x,abs(msig));
xlabel('Position (mm)');
ylabel('Magnitude');
title('Magnitude');
grid on;

% Plot the phase 
subplot(1,3,2);
plot(x,angle(msig));
xlabel('Position (mm)');
ylabel('Phase');
title('Phase');
grid on;

% Plot z-component of magnetization 
subplot(1,3,3);
plot(x,m(3,:));
xlabel('Position (mm)');
ylabel('Z-component of Magnetization');
title('Z-component of Magnetization');
grid on;

figure ('Name','Applying Refocusing Gradients (functions of time)','NumberTitle', 'off')

% Plot RF Waveform 
subplot(1,2,1);
plot(t,real(B1),'--',t,imag(B1),':',t,abs(B1),'-');
xlabel('time (s)');
ylabel('RF Waveform');
title('RF Waveform');
legend('Real','Imaginary','Abs');
grid on;

% Plot Gradient Waveform 
subplot(1,2,2);
plot(t,grad);
xlabel('time (s)');
ylabel('Gradient Waveform ');
title('Gradient Waveform');
grid on;

% Plots shifted to the right by ~ 10mm. 
% 'peaks' of plots are not off-center
% phase is more sparatic 

%% Replace RF with RF of 3e and RF modulated by -900Hz pure exponential

T1 = 600;             % Longitudinal relaxation time 
T2 = 100;             % Transverse relaxation time 
df = 0;               % Off-resonance frequency (Hz) 
t = (0:0.0001:0.006); % Sample time (s) 
x = (-50:1:50);       % positions (mm) 

RFmodulate = exp(2*pi*900*i*t); % Factor to modulate RF pulse 
B1 = 0.06*sinc((1000*t-3)/1); % B1 field 
grad = 0.1*ones(size(t)); % gradient 
phase_flat = -0.52; % Refocusing phase factor/ratio
grad = [grad grad*phase_flat]; % refocusing gradient 
t = [t t+0.006]; % shift time to account for refocusing gradient 
B1 = (B1.*RFmodulate + B1.*conj(RFmodulate)); % making RF zero for refocusing gradient 
B1 = [B1 0*B1];
[m,msig] = sliceprofile(B1,grad,t,T1,T2,x,df); % Call the function 

figure ('Name','Applying Refocusing Gradients (functions of position)','NumberTitle', 'off')

% Plot the magntidue 
subplot(1,3,1);
plot(x,abs(msig));
xlabel('Position (mm)');
ylabel('Magnitude');
title('Magnitude');
grid on;

% Plot the phase 
subplot(1,3,2);
plot(x,angle(msig));
xlabel('Position (mm)');
ylabel('Phase');
title('Phase');
grid on;

% Plot z-component of magnetization 
subplot(1,3,3);
plot(x,m(3,:));
xlabel('Position (mm)');
ylabel('Z-component of Magnetization');
title('Z-component of Magnetization');
grid on;

figure ('Name','Applying Refocusing Gradients (functions of time)','NumberTitle', 'off')

% Plot RF Waveform 
subplot(1,2,1);
plot(t,real(B1),'--',t,imag(B1),':',t,abs(B1),'-');
xlabel('time (s)');
ylabel('RF Waveform');
title('RF Waveform');
legend('Real','Imaginary','Abs');
grid on;

% Plot Gradient Waveform 
subplot(1,2,2);
plot(t,grad);
xlabel('time (s)');
ylabel('Gradient Waveform ');
title('Gradient Waveform');
grid on;

 
