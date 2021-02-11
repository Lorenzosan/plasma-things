function pol_mat = cart_to_pol(cart_mat,N,M)

% This program takes as input a density matrix on a cartesian grid
% and returns the equivalent on a polar grid.

% cart_mat is the input matrix of plasma electric density on a cartesian grid
% Input matrix has dimension 512x672
% N ad M are the divisions you want in radial and azimuthal direction

pol_mat = zeros(N,M);

nrow = size(cart_mat,1);
ncol = size(cart_mat,2);

m_to_pix = 1.875e-04;

% Trap center co-ordinates, calculated separately
T = [365,269]; 

for i = 1:N
    for j = 1:M
        xp = T(2) + 240*(i-1)*cos(2*pi*(j-1)/M)/N;
        yp = T(1) + 240*(i-1)*sin(2*pi*(j-1)/M)/N;
        fxp = floor(xp);
        fyp = floor(yp);
        %[fxp, fyp];
        delta_x = xp - fxp;
        delta_y = yp - fyp;
        pol_mat(i,j) = cart_mat(fxp,fyp+1)*delta_y*(1-delta_x) + ...
            cart_mat(fxp+1,fyp+1)*delta_x*delta_y + ...
            cart_mat(fxp,fyp)*(1-delta_y)*(1-delta_x) + ...
            cart_mat(fxp+1,fyp)*delta_x*(1-delta_y);
    end

end
