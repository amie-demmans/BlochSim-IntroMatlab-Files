function [Mss] = sesignal(T1,T2,TE,TR,df)

Rflip = yrot(pi/2); %flip matrix w rotation along yaxis
Refocus = xrot(pi); %refocusing pulse with rotation along x
[Atr,Btr] = freeprecess(TR-TE,T1,T2,df); %free precession matrix for time = TR
%[Ate,Bte] = freeprecess(TR-(TE/2),T1,T2,df);
[Ate2,Bte2] = freeprecess(TE/2,T1,T2,df); %free precession matric for time = TE/2
%[Ate3,Bte3] = freeprecess(TE,T1,T2,df);
Atr = [0,0,0;0,0,0;0,0,1]*Atr; % define new matrix Atr for steady state magentization neglecting transverse magnetization

%%  My way that makes sense to me in terms of the pulse sequence diagram 

Mss = (inv(eye(3)-Atr*Rflip)*Btr); % apply steady state magnetization to neglect transverse magnetization
Mss = Rflip*Mss;                   % apply 90 degree flip to Mss 
Mss = Ate2*Mss+Bte2;               % apply free precession at TE/2
Mss = Refocus*Mss;                 % apply 180 degree flip to Mss 
Mss = Ate2*Mss+Bte2;               % apply free precession at TE/2 
Mss = Mss(1)+1i*Mss(2);            % complex number representation to display both magnitude and phase info

%% Bloch sim answer

%Mss = inv(eye(3)-Ate2*Refocus*Ate2*Rflip*Atr) * (Bte2+Ate2*Refocus*(Bte2+Ate2*Rflip*Btr));
%Msig = Mss(1)+1i*Mss(2);

%% Other way I was trying that didn't work out but keeping for reference can delete later

%Mss = (inv(eye(3)-Atr*Rflip)*Btr); % apply steady state magnetization at point prior to excitation to neglect transverse magnetization 
%Mss = Rflip*Mss; % apply 90 degree flip angle to Mss
%Mss = Ate*Mss+Bte; % apply free precession at TE/2
%Mss = (inv(eye(3)-Ate*Refocus)*Bte); 
%Mss = Refocus*Mss; %apply refocusing pulse (180 degree) to Mss
%Mss = Ate2*Mss+Bte2; % apply free precession at TR-TE/2
%Mss = Rflip*Mss;
%Mss = Ate3*Mss+Bte3;

%Mss = (inv(eye(3)-Ate*Rflip)*Bte);
%Mss = Rflip*Mss;
%Mss = Ate2*Mss+Bte2;
%Mss = Refocus*Mss;
%Mss = Ate*Mss+Bte;
%Mss = Rflip*Mss;