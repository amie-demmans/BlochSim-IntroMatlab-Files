function [m,msig] = sliceprofile(rotations,grad,t,T1,T2,pos,df)

gamma = 4258; % gyromagnetic ratio (Hz/G)
dt = t(2)-t(1);
RFrotations = 2*pi*gamma*rotations*dt;

pos = pos(:).';
m = zeros(3,length(pos)); % 3xN array 
msig = zeros(1,length(pos)); % 1xN array 
%msig = 0*pos;
%m = [msig,msig,msig]';

for n = 1:length(pos)
    
   M = [0,0,1]';
   [Afp,Bfp] = freeprecess(1000*dt/2, T1, T2, df);
   
   for k = 1:length(RFrotations)
       % Apply gradient-induced precessions for half the sample time 
       M = Afp*M+Bfp;
       grad_rotations = gamma*2*pi*(dt/2)*grad(k)*(pos(n)/10);
       M =  zrot(grad_rotations)*M;

       % Apply RF rotation over full sample time 
       M = throt(abs(RFrotations(k)),angle(RFrotations(k))) *M;

       % Apply gradient-induced precessions for half the sample time 
       M = Afp*M+Bfp;
       grad_rotations = gamma*2*pi*(dt/2)*grad(k)*(pos(n)/10);
       M =  zrot(grad_rotations)*M;
   end 
   m(:,n) = M;
   msig(n) = M(1)+1i*M(2); % Apply complex number representation to retain both magnitude and phase info 
end 





       


