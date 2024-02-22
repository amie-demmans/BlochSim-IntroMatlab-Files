function [Msig,Mss] = spgrsignal(flip,T1,T2,TE,TR,df,Nex,inc)

if (nargin < 8)
	inc = 117/180*pi;
end

if (nargin < 7)
	Nex = 100;
end

if (nargin < 6)
	df = 0;
end

Nf = 100;	% Simulate 50 different gradient-spoiled spins.
phi = (1:Nf)/Nf*2*pi;  % phase of excitation

M=zeros(3,Nf,Nex+1); % Generate a 3 x Nf array of zeros Nex+1 times
%Msig = zeros(Nex);   % Generate a Nex x Nex array of zeros 


% Call free precess function for TR and TE 	
[Ate,Bte] = freeprecess(TE,T1,T2,df);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,df);

M = [zeros(2,Nf);ones(1,Nf)]; % Create 3 x Nf matrix with first two rows filled with zeros, and last row with 1
on = ones(1,Nf); % Create a 1 x Nf row vector filled with ones 
	
Rfph = 0;  % Initial excitation pulse 
Rfinc = inc; % Redefine the increment 

for n=1:Nex % Generate loop for n = 1 to number of excitations 

	A = Ate * throt(flip,Rfph);
	B = Bte;
	M = A*M+B*on;

	Msig = mean( squeeze(M(1,:)+1i*M(2,:)) ) * exp(-1i*Rfph);
    Mss = M;
	
	M=Atr*M+Btr*on;

	for k=1:Nf
		M(:,k) = zrot(phi(k))*M(:,k);
    end

	Rfph = Rfph+Rfinc;
	Rfinc = Rfinc+inc;
end