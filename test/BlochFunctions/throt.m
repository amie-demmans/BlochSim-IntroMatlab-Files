function Rth = throt(phi,theta)
Rz = zrot(-theta); % apply zrot to -theta
Rx = xrot(phi);    % apply zrot to phi 
Rth = (Rz')*Rx*Rz; % Rth = Rz(theta)*Rx(phi)*Rz(-theta)


% Returning rotation matrix of phi about axis defined by y = x*tan(theta)
