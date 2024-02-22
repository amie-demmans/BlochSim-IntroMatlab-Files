function [Mss,Msig] = gresignal(flip,T1,T2,TE,TR,df)

N = 100; % 100 points the average magnetization is calculated over 
M = zeros(3,N); % create a 3xN matrix of zeros 
%phi = linspace(-2*pi,2*pi,N); %generates row vector of 100 points with phase varying from -2pi to 2pi

phi = ((1:N)/N-0.5 ) * 4*pi;
for k = 1:N % From k = 1 to N 
    [Mss] = gssignal(flip,T1,T2,TE,TR,df,phi(k)); % Call the function gssignal with phi being from k = 1 to N
    M(:,k) = Mss; % 
end

%spoiler_angle = 4*pi;
%spoiler_matrix = xrot(spoiler_angle);
%M = spoiler_matrix*M;

Mss = mean(M')';
Msig = Mss(1)+1i*Mss(2);

