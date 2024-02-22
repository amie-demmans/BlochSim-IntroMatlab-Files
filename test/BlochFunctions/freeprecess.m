function [Afp,Bfp]=freeprecess(t,T1,T2,df)
phi = 2*pi*df*t/1000; % Resonant precession
E1 = exp(-t/T1); % z-component magnetization vector (Mz)
E2 = exp(-t/T2); % x/y magnetization vectore (Mx/My)

% Combine relaxation and precession effects 
Afp = [E2,0,0;0,E2,0;0,0,E1]*zrot(phi); 
Bfp = [0,0,1-E1]'; 