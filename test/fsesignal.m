function [Msig,Mss] = fsesignal(T1,T2,TE,TR,df,ETL)

Rflip = yrot(pi/2); %flip matrix w rotation along yaxis
Refocus = xrot(pi); %refocusing pulse with rotation along x
[Atr,Btr] = freeprecess(TR-ETL*TE,T1,T2,df); %free precession matrix for time = TR
[Ate,Bte] = freeprecess(TE/2,T1,T2,df);

Atr = [0,0,0;0,0,0;0,0,1]*Atr;

A = eye(3);
B = [0,0,0]';

% For each echo, 'propogate' A and B
for k = 1:ETL
    A = Ate*Refocus*Ate*A;
    B = Bte+Ate*Refocus*(Bte+Ate*B);
end

% Propogate A and B through to just after flip, and calculate steady state 
A = Rflip*Atr*A;
B = Rflip*(Btr+Atr*B);

Mss = inv(eye(3)-A)*B;
M = Mss;

 %Calculate signal on each echo 
for k = 1:ETL
    M = Ate*Refocus*Ate*M + Bte+Ate*Refocus*Bte;
    Mss(:,k) = M;
    Msig(k) = M(1)+1i*M(2);
 end


%%  Another way similar to sesignal  
%for k = 1:ETL 
%Mss(:,ETL) = (inv(eye(3)-Atr*Rflip)*Btr);
%Mss(:,ETL) = Rflip*Mss(:,ETL);
%Mss(:,ETL) = Ate*Mss(:,ETL)+Bte;
%Mss(:,ETL) = Refocus*Mss(:,ETL); 
%Mss(:,ETL) = Ate*Mss(:,ETL)+Bte;
%Mss(:,ETL) = Mss(1)+1i*Mss(2);
%end 

