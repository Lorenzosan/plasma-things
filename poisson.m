clear all

%% 0. Definizione costanti utili
epsilon0 = 8.85*10^-12; % F/m
N = 2048;           % suddivisione radiale
M = 2048;           % suddivisione azimutale 
Rw = 0.045;        % raggio trappola P-M [m]
r = linspace(0,Rw,N);
r = r';

% Spaziatura dei punti in r e azimutalmente
dr  = Rw/N;
dth = 2*pi/M;

%% 1. Dichiarazione del vettore con le boundary conditions
%  Le condizioni al contorno sul potenziale scalare phi(R,th)
phiBC = zeros(1,M); %zero per comodità

%% 2. Definizione del settore interessato e ridefinizione delle BC
sectno = 2;
phi_sect = [0 0];% 80 -80];

phiBC(1,1:M/2) = phi_sect(1);
phiBC(1,M/2:end) = phi_sect(2);

% for i = 1:sectno
%     phiBC((i-1)*(M/sectno)+1:i*(M/sectno)) = phi_sect(i);
% end

%% 3. Defininzione della matrice polare e matrice del potenziale
% density_matrix = zeros(N,M);
% density_matrix(580:640,1:20) = 0.0001;
% density_matrix(580:640,1:200) = 0.000001;
% density_matrix(580:640,600:825) = -0.000001;

dm = load('rho1.txt');
density_matrix = cart_to_pol(dm,N,M);


%% 4. Calcolo delle DFT della matrice polare, potenziale e BC
density_matrix_ft = fft(density_matrix,M,2);
phiBC_ft  = fft(phiBC,M);
phi_matrix_ft = zeros(N,M);

%% 5. Impostazione delle matrici del sistema per m=0 e m>0
% Definizione delle diagonali come da wikipedia
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
%     beta_m = sin(dth*m/2)/(dth/2);

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

%% 6. Antitrasformo il risultato 
phi_matrix = ifft(phi_matrix_ft, M, 2);
phi = real(phi_matrix);


%% 7. Visualizzazione risultati
% NB: si aggiunge una fetta perchè altrimenti non chiude il cerchio
phi(:,M+1) = phi(:,1); 

% Generazione vettore r e theta e creazione della matrice (r,theta)
rr = linspace(0,Rw,N);
thh=linspace(0,2*pi,M+1);
[myr,theta] = meshgrid(rr,thh);

% pol2cart converte una coppia di coord. cilindriche in cartesiane
[x,y] = pol2cart(theta,myr);

figure
colormap('jet')
% Devo avere r sulle colonne e th sulle righe -> traspongo
% surface(x,y,phi','FaceAlpha',1,'LineStyle','none','FaceColor','flat')
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
% Devo avere r sulle colonne e th sulle righe -> traspongo
contour(x,y,phi')
% view(0,90) 
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
% daspect([1 1 aspz])
daspect auto
title("Plasma density");
colorbar


