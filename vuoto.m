%potenziale di vuoto

%dimensioni della trappola
Rw=0.09/2;
l=0.15;
ls=0.15; %lunghezza singolo elettrodo S
lc=0.09; %lunghezza singolo elettrodo C
dz=0.001; %separazione elettrodi
Nel=12; %numero elettrodi
zn=zeros(1,2*Nel);
zn(1) = 0;
zn(2) = 0.15; %GND
zn(3) = 0.151; 
zn(4) = 0.241; %C1
zn(5) = 0.242;
zn(6) = 0.332; %C2
zn(7) = 0.333;
zn(8) = 0.483; %S8
zn(9) = 0.484;
zn(10)= 0.574; %C4
zn(11)= 0.575;
zn(12)= 0.725; %S2
zn(13)= 0.726;
zn(14)= 0.816; %C5
zn(15)= 0.817;
zn(16)= 0.966; %S4
zn(17)= 0.967;
zn(18)= 1.057; %C6
zn(19)= 1.058;
zn(20)= 1.148; %C7
zn(21)= 1.149;
zn(22)= 1.239; %C8
zn(23)= 1.240;
zn(24)= 1.390; %SH

L = 1.39;

% for j = 2:2*Nel
%     if mod(j,2)==0
%     zn(j)=zn(j-1)+l;
%     else
%     zn(j)=zn(j-1)+dz;    
%     end
% end
% L=l*Nel+dz*(Nel-1);


%assegno i potenziali
prompt={'GND','C1','C2','S8','C4','S2','C5','S4','C6','C7','C8','SH'};
default={'0','0','0','0','0','0','0','0','0','0','0','0'};
dlgtitle='Electrode voltages';
answer=inputdlg(prompt,dlgtitle,[1 50],default);

v(1)=str2double(answer{1});
v(2)=str2double(answer{2});
v(3)=str2double(answer{3});
v(4)=str2double(answer{4});
v(5)=str2double(answer{5});
v(6)=str2double(answer{6});
v(7)=str2double(answer{7});
v(8)=str2double(answer{8});
v(9)=str2double(answer{9});
v(10)=str2double(answer{10});
v(11)=str2double(answer{11});
v(12)=str2double(answer{12});

k=pi/L;
%r=linspace(-Rw,Rw,10000);
z=linspace(0,L,5000);

prompt={'D'};
default={'0'};
dlgtitle='offset';
answ=inputdlg(prompt,dlgtitle,[1 50],default);

D=str2double(answ{1});
a=0;
b=0;
phi_vuoto_z=zeros(1,length(z));

for zind = 1:length(z)
    for n = 1:1000
        k_n = n*pi/L;
        ybess = besseli(0,D*k_n)/( besseli(0,Rw*k_n) * L/2) *sin(k_n*z(zind));
        phi_vuoto_z(zind) = phi_vuoto_z(zind) + ...
            ( v(1)*cos(k_n*zn(1)) - v(12)*(k_n*zn(24)) )*ybess / k_n;
        
        for j = 1:length(v)-1
            phi_vuoto_z(zind) = phi_vuoto_z(zind) + (v(j+1)-v(j))/(zn(2*j+1)-zn(2*j)) * ...
                (sin(k_n*zn(2*j+1))-sin(k_n*zn(2*j))) * ybess / k_n^2;
        end

    end
end

% for zind = 1:length(z)
%     for n = 1:200
%         a=0;
%         for p = 1:length(v)
%             a=a+v(p)*(cos(k*n*((l*p)-l))-cos(k*n*l*p))/(n*pi);
%         end
%         ybess=besseli(0,D*k*n)/besseli(0,Rw*k*n)*sin(k*n*z(zind));
%         phi_vuoto_z(zind)=phi_vuoto_z(zind)+a*ybess;
%     end
%     phi_vuoto_z(zind)=2*phi_vuoto_z(zind);
% end
plot(z, phi_vuoto_z, '-'); hold on; grid on;