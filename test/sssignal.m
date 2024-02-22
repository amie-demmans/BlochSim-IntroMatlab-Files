function [Mss] = sssignal(flip,T1,T2,TE,TR,df)

Rflip = yrot(flip); %flip matrix w rotation along yaxis

[Atr,Btr] = freeprecess(TR,T1,T2,df); %free precession matrix for time = TR
[Ate,Bte] = freeprecess(TE,T1,T2,df); %free precession matric for time = TE


Mss = (inv(eye(3)-Atr*Rflip)*Btr); % apply steady state magnetization at point prior to excitation
Mss = Rflip*Mss; % apply flip angle to Mss
Mss = Ate*Mss+Bte; % apply free precession at TE

%% Other ways I was trying that didn't work, kept for reference can delete later

%Msig = Mss(1)+1i*Mss(2);

%Mss = Rflip * M;
%	Mss = Ate * Mss + Bte;
%	Mss = Atr * Mss + Btr;
%Mss = Ate*Rflip*Atr*Mss + (Ate*Rflip*Btr+Bte);

%Mss = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte);
%Msig = Mss(1)+1i*Mss(2);


