function Mss = srsignal(flip,T1,T2,TR,df)

Rflip = yrot(flip); %flip matrix w rotation along yaxis
[Atr,Btr] = freeprecess(TR,T1,T2,df); %free precession matrix for time = TR

Atr = [0,0,0;0,0,0;0,0,1]*Atr; % define new matrix Atr for steady state magentization neglecting transverse magnetization

Mss = (inv(eye(3)-Atr*Rflip)*Btr); % apply steady state magnetization with new Atr to neglect transverse magnetization
