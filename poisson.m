clear all

%% 0. Constants definitions
epsilon0 = 8.85*10^-12; % Dielectric vacuum constant
N = 2048;               % Radial grid division
M = 2048;               % Azimuthal grid division
Rw = 0.045;             % Penning-Malmberg trap's radius [m]
r = linspace(0,Rw,N);
r = r';

% Minimum separation between points (radial and azimuthal)
dr  = Rw/N;
dth = 2*pi/M;

%% 1. Definition of boundary conditions on scalar potential phi(R,th)
phiBC = zeros(1,M); 

%% 2. Sectoring of trap boundaries and re-definition of boundary conditions on sectors
sectno = 2;
phi_sect = [0 0];% 80 -80];

phiBC(1,1:M/2) = phi_sect(1);
phiBC(1,M/2:end) = phi_sect(2);

% for i = 1:sectno
%     phiBC((i-1)*(M/sectno)+1:i*(M/sectno)) = phi_sect(i);
% end

%% 3. Definition of electrical density matrix and conversion to polar grid
dm = load('rho1.txt');
density_matrix = cart_to_pol(dm,N,M);

%% 4. Discrete Fourier transform of polar matrix, potential and boundary conditions 
density_matrix_ft = fft(density_matrix,M,2);
phiBC_ft  = fft(phiBC,M);
phi_matrix_ft = zeros(N,M);

%% 5. Definition of matrices for m=0 and m>0
%  Diagonals are defined as in wikipedia's
%  page on tridiagonal systems
%
%  b1 c1 0  0            d1
%  a2 b2 c2 0            d2
%  0  a3 b3 c3 ....  x = d3
%     :
%     :

A = zeros(1,N);
B = zeros(1,N);
C = zeros(1,N);
rhotilde = zeros(1,N);

for m = 0:(M-1)  
    A(1:N-1) = ( 1/(dr^2) - 1./(2*dr*r(1:N-1)) );
    A(N) = 0;
    B(N) = 1;
    C(2:N) = ( 1/(dr^2) + 1./(2*r(2:N)*dr) );
    rhotilde(2:N-1) = -density_matrix_ft(2:N-1,m+1)/epsilon0;
    rhotilde(N) = phiBC_ft(m+1);
    
if m == 0
    B(1) = -4/(dr^2);
    B(2:N-1) = -2/(dr^2);
    C(1) = 4/(dr^2);
    rhotilde(1) = -density_matrix_ft(1,m+1)/epsilon0;
else
    B(1) = 1;
    B(2:N-1) = ( -2/(dr^2) - m^2./(r(2:N-1).^2) );
    C(1) = 0;
    rhotilde(1) = 0;
end

    phi_matrix_ft(:,m+1) = ThomasAlgorithm(A,B,C,rhotilde);
end

%% 6. Output is inverse transformed, in order to retrieve the potential
phi_matrix = ifft(phi_matrix_ft, M, 2);
phi = real(phi_matrix);


%% 7. Outcomes plot
% NB: A slice of potential is added in order to close the circle
phi(:,M+1) = phi(:,1); 

% Generazione vettore r e theta e creazione della matrice (r,theta)
rr = linspace(0,Rw,N);
thh=linspace(0,2*pi,M+1);
[myr,theta] = meshgrid(rr,thh);

% pol2cart converts a pair of cylindrical coordinates in cartesian form
[x,y] = pol2cart(theta,myr);

figure
colormap('jet')
% I need 'r' on columns and th on rows -> transposition
surf(x,y,phi','FaceAlpha',1,'LineStyle','none','FaceColor','flat')
% view(0,90) 
hold on
aspz = max(max(abs(phi)));
if aspz == 0
aspz = 1;
end
% daspect([1 1 aspz])
daspect auto
title("Electric potential inside ELTRAP");
colorbar

figure
colormap('jet')
contour(x,y,phi')
hold on
aspz = max(max(abs(phi)));
if aspz == 0
aspz = 1;
end
daspect([1 1 aspz])
title("Equipotential lines");

figure
colormap('jet')
density_matrix(:,M+1) = density_matrix(:,1); % se omesso perdo un settore
surf(x,y,density_matrix','FaceAlpha',1,'LineStyle','none','FaceColor','flat')
hold on
aspz = max(max(abs(density_matrix)));
if aspz == 0
aspz = 1;
end
daspect auto
title("Plasma density");
colorbar


