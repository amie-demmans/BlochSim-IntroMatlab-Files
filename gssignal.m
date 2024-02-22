function Mss = gssignal(flip,T1,T2,TE,TR,df,phi)

Rflip = yrot(flip); %flip matrix w rotation along yaxis
[Atr,Btr] = freeprecess(TR,T1,T2,df); %free precession matrix for time = TR
[Ate,Bte] = freeprecess(TE,T1,T2,df); %free precession matric for time = TE

%M = [0,0,1]';
Atr = zrot(phi)*Atr;
Mss = (inv(eye(3)-Atr*Rflip)*Btr); % apply steady state magnetization at point prior to excitation
Mss = Rflip*Mss; % apply flip angle to Mss 
Mss = Ate*Mss+Bte; % apply free precession at TE