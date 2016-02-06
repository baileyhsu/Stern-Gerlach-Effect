%SGE potential in 2d: sigma.B;
%evolution in time
function SGE2d
global X Z Rho

fig=figure;
set(fig,'DoubleBuffer','on');
set(gca,'xlim',[-80 80],'ylim',[-80 80],...
    'nextplot','replace','Visible','off');

aviobj = avifile('5.avi','fps',1);

dt = 0.2;      Lt = 2.0;       Nt = Lt/dt;  
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
z0 =4.0;
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
Wpktup = (Amplup/sqrt(pi*wx*wz))*exp(-1*((X - x0).^2/wx^2 + (Z - z0).^2/wz^2));   
Wpktdn = (Ampldn/sqrt(pi*wx*wz))*exp(-1*((X - x0).^2/wx^2 + (Z - z0).^2/wz^2));
Rho = (abs(Wpktup)).^2 - (abs(Wpktdn)).^2;
figure(2);
plotrho
pause(0.05)
%exponentiation w/ Pauli matrices
%Cs = cos(0.5*dt*Babs);  % Trotter
%Sn = sin(0.5*dt*Babs);
Cs = cos(dt*Babs);
Sn = sin(dt*Babs);
ExPot11 = (Cs + 1i*Sn.*Bz./Babs);
ExPot12 = 1i*Sn.*(Bx - 1i*By)./Babs;
ExPot21 = 1i*Sn.*(Bx + 1i*By)./Babs;
ExPot22 = (Cs - 1i*Sn.*Bz./Babs);
%kin en
Psiup = Wpktup;
Psidn = Wpktdn;
ExKinFT = ExpKE(dt, Lx, Nx, Lz, Nz);
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
    %WpupFT = fft2(ExPot11.*Psiup+ExPot12.*Psidn); %Trotter
    %WpdnFT = fft2(ExPot21.*Psiup+ExPot22.*Psidn);
    WpupFT = fft2(Psiup);
    WpdnFT = fft2(Psidn);
    Wpup = ifft2(ExKinFT.*WpupFT);
    Wpdn = ifft2(ExKinFT.*WpdnFT);
    PsiupFT = fft2(ExPot11.*Wpup+ExPot12.*Wpdn);
    PsidnFT = fft2(ExPot21.*Wpup+ExPot22.*Wpdn);
    
    Psiup=ifft2(ExKinFT.*PsiupFT);
    Psidn=ifft2(ExKinFT.*PsidnFT);
    
    %Rho  = 1i*(conj(Psidn).*Psiup-conj(Psiup).*Psidn);
    Rho =conj(Psiup).*Psiup-conj(Psidn).*Psidn;
    %Rho = conj(Psiup).*Psidn+conj(Psidn).*Psiup ;
    %Rho = (abs(Psidn)+abs(Psiup)).^2;
    %Rho=conj(Psiup).*Psiup+conj(Psidn).*Psidn;
    plotrho
    pause(0.01)
%A = moviein(Nx/4);

frame = getframe(gca);
  
aviobj = addframe(aviobj,frame);

end
aviobj = close(aviobj);
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
axis ([-8 8 -8 8]);
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
Kin = 0.5*0.5*(KX.*KX + KZ.*KZ);            %kinetic energy 
Res = exp(-i*dt*Kin);
%