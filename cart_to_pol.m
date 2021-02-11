function pol_mat = cart_to_pol(cart_mat,N,M)

% La matrice di input e' 512x672
% N ed M sono le suddivisioni che si vogliono radialmente e azimutalmente

pol_mat = zeros(N,M);

nrow = size(cart_mat,1);
ncol = size(cart_mat,2);

m_to_pix = 1.875e-04;
T = [365,269]; % centro trappola

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
