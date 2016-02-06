%SGE potential in 2d: sigma.B;
%evolution in time
function SGE2d
global X Z Rho

fig=figure;
set(fig,'DoubleBuffer','on');
set(gca,'xlim',[-80 80],'ylim',[-80 80],...
    'nextplot','replace','Visible','off');

aviobj = avifile('5.avi','fps',1);

dt = 0.1;       Lt = 2.0;        Nt = Lt/dt;  
dx = 1/16;      Lx = 32;        Nx = Lx/dx;
dz = dx;        Lz = 32;        Nz = Lz/dx;
B = 1.0;        B0 = 0.0;
Amplup =1.0/sqrt(2);           Ampldn=1.0/sqrt(2);
wx=1;wz=1;
%
%sigma = 1.0;           %w packet's width
x = (-Nx/2:(Nx/2-1))*dx;
x0 =0.0;
z = (-Nz/2:(Nz/2-1))*dz;
z0 =1.0;
[X, Z] = meshgrid(x, z); 
%Magnetic field and plot:
Bx = -B*X;
By =0;
Bz = B0 + B*Z;
figure(1);
plotB(Nx, Nz, Bx, Bz);
%pause
Babs = sqrt(Bx.^2 + By.^2 + Bz.^2)+10^(-5);
%gaussian w pkt
Wpktup = 1/(pi*wx*wz)^0.5*Amplup*exp(-0.5*((X - x0).^2/wx^2 + (Z - z0).^2/wz^2));
Wpktdn = 1/(pi*wx*wz)^0.5*Amplup*exp(-0.5*((X - x0).^2/wx^2 + (Z - z0).^2/wz^2));
Rho = (abs(Wpktup)).^2 - (abs(Wpktdn)).^2;
figure(2);
plotrho
pause(0.05)
%exponentiation w/ Pauli matrices
%Cs = cos(0.5*dt*Babs);  % Trotter
%Sn = sin(0.5*dt*Babs);
Cs = cos(0.5*dt*Babs);
Sn = sin(0.5*dt*Babs);
ExPot11 = (Cs + 1i*Sn.*Bz./Babs);
ExPot12 = 1i*Sn.*(Bx - 1i*By)./Babs;
ExPot21 = 1i*Sn.*(Bx + 1i*By)./Babs;
ExPot22 = (Cs - 1i*Sn.*Bz./Babs);
%kin en
Psiup = Wpktup;
Psidn = Wpktdn;
ExKinFT = ExpKE(dt, Lx, Nx, Lz, Nz);
entropy=[0];
for k = 1:Nt
%    t = k*dt;
%analytic    
%     den = 1 + i*t;
%     Exn1 = 0.5*i*(Zsq + i*t*zpq)/t;
%     Exn2 = -0.5*i*(z + i*t*z0 - 0.5*B*t*t).^2/(t*den);
%     Exn3 = 0.5*i*B*t*z - i*B*B*t^3/24.0;
%     PsiA = Ampl*exp(Exn1 + Exn2 + Exn3)/den;
%     RhoA = (abs(PsiA)).^2;
%numeric    
    WpupFT = fft2(ExPot11.*Psiup+ExPot12.*Psidn); %Trotter
    WpdnFT = fft2(ExPot21.*Psiup+ExPot22.*Psidn);
    %WpupFT = fft2(Psiup);
    %WpdnFT = fft2(Psidn);
    Wpup = ifft2(ExKinFT.*WpupFT);
    Wpdn = ifft2(ExKinFT.*WpdnFT);
    Psiup = ExPot11.*Wpup+ExPot12.*Wpdn;
    Psidn = ExPot21.*Wpup+ExPot22.*Wpdn;
    %Rho  = 1i*(conj(Psidn).*Psiup-conj(Psiup).*Psidn);
    Rho =conj(Psiup).*Psiup-conj(Psidn).*Psidn;
    %Rho = conj(Psiup).*Psidn+conj(Psidn).*Psiup ;
    %Rho = (abs(Psidn)+abs(Psiup)).^2;
    %Rho=conj(Psiup).*Psiup+conj(Psidn).*Psidn;
    
    
    Psiupup=conj(Psiup).*Psiup;
    Psiupdn=conj(Psiup).*Psidn;
    Psidnup=conj(Psidn).*Psiup;
    Psidndn=conj(Psidn).*Psidn;
    
    
    a=(sum(sum(Psiupup))).*dx^2;
    b=(sum(sum(Psidndn))).*dx^2;
    % c and d is similar to a and b but change of integration range, the
    % overlapping is not seen in the mixed state.
   % c=sum(sum(Psidndn(97:192,:))).*dx^2+sum(sum(Psidndn1(97:192,:))).*dx^2;
   % d=sum(sum(Psiupup1(1:96,:))).*dx^2+sum(sum(Psiupup(1:96,:))).*dx^2;
    
    c=(sum(sum(Psiupdn))).*dx^2;
    d=(sum(sum(Psidnup))).*dx^2;
    
    ev=eig([a d ; c b]);
    
    a1=ev(1);
    a2=ev(2);
     
    entropy=[entropy,-a1*log2(a1)-a2*log2(a2)];
    
    
    
    
    plotrho
    pause(0.01)
%A = moviein(Nx/4);

frame = getframe(gca);
  
aviobj = addframe(aviobj,frame);

end
aviobj = close(aviobj);


fid = fopen('ip2.txt','w');
fprintf(fid,'%12.8f\n',entropy);
status = fclose(fid);
time = 0:dt:Lt; 
plot(time,entropy,'ro'); 
axis([0 Lt 0.0 1.02]) 
xlabel('time') 
ylabel('entropy') 
title('entropy versus Time')






%legend('Evolved W Pkt', 2); %'Location', Northwest);
function plotrho
global X Z Rho
%contour(X, Z, Rho, 10);
surf(X, Z, Rho);
shading interp;
colormap jet;
set(gca,'FontSize',28)
colorbar('FontSize',28)

%axis ([-5 5 -5 5 0 10]);
axis ([-8 8 -8 8]);
%axis square;
xlabel('');  ylabel('');
topline = sprintf('');
title(topline);
%
function plotB(Nx, Nz, Bx, Bz)
global X Z
indx = 1:30:Nx;
indz = 1:30:Nz;
xx = X(indx, indx);
zz = Z(indz, indz);
Bxx = Bx(indx, indx);
Bzz = Bz(indz, indz);
quiver(xx, zz, Bxx, Bzz, 5.0,'LineWidth',5.0);
axis ([-5 5 -5 5]);
%axis square;
xlabel('x');  ylabel('z');
topline = sprintf('(a)');
set(gca,'FontSize',40)
title(topline);
%
function Res = ExpKE(dt, Lx, Nx, Lz, Nz)
Kx = (2*pi/Lx)*[0:Nx/2-1 -Nx/2:-1];     %kx grid
Kz = (2*pi/Lz)*[0:Nz/2-1 -Nz/2:-1];     %kx grid
[KX, KZ] = meshgrid(Kx, Kz);
Kin = 0.5*(KX.*KX + KZ.*KZ);            %kinetic energy 
Res = exp(-i*dt*Kin);
%