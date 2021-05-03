function FieldPlots(TP,FSP,CLP,RMP,VXP,VYP,PP)
%Plot all field variables: T, FL, CL, rho, VX, VY, P, dFS 
%suffix 'P'-->Plot
%Created 2019-10-17

global NIX
global NIY
global NMX
global NMY
global dxI
global dxM
global dyI
global dyM
global x
global y

W=120.0/1000.0;%ingot width [m]
H=100.0/1000.0;%ingot height [m]

% X1=[0.5*dxI*1000.0:dxI*1000.0:0.5*W-0.5*dxI*1000.0];
% X2=[0.5*W+0.5*dxM*1000.0:dxM:];
% Y1=[0.5*dyI*1000.0:1000.0*dyI:H-0.5*dyI*1000.0];
% [XI,YI]=meshgrid(X1,Y1);
% 
% [XM,YM]=meshgrid();

%NOTE: meshgrid with different intervals
% x = [0,10,20,30:40,50,60];
% y = 0:10:80;
% [x,y] = meshgrid(x,y);
% plot(x,y,x',y')
% axis([0,60,0,80])

%All T field
TPP=zeros(NIY+NMY+3,NIX+NMX+3);
TPP(1:NIY+1,1:NIX+1)=TP(1:NIY+1,1:NIX+1);
TPP(1:NIY+1,NIX+2)=TPP(1:NIY+1,NIX+1);
TPP(NIY+2,1:NIX+2)=TPP(NIY+1,1:NIX+2);
TPP(1:NIY+1,NIX+3:NIX+3+NMX)=TP(1:NIY+1,NIX+2:NIX+NMX+2);
TPP(NIY+2,NIX+3:NIX+3+NMX)=0.5*(TP(NIY+1,NIX+2:NIX+NMX+2)+TP(NIY+2,NIX+2:NIX+NMX+2));
TPP(NIY+3:NIY+3+NMY,1:NIX+1)=TP(NIY+2:NIY+2+NMY,1:NIX+1);
TPP(NIY+3:NIY+3+NMY,NIX+2)=0.5*(TP(NIY+2:NIY+2+NMY,NIX+1)+TP(NIY+2:NIY+2+NMY,NIX+2));
TPP(NIY+3:NIY+3+NMY,NIX+3:NIX+3+NMX)=TP(NIY+2:NIY+2+NMY,NIX+2:NIX+2+NMX);
subplot(2,3,1)
pcolor(x*1000.0,y*1000.0,TPP);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Temperature [K]');
axis([0.0 x(NIX+NMX+3)*1000.0 0.0 y(NIY+NMY+3)*1000.0]);
shading interp
colorbar;

%FS distribution
subplot(2,3,2)
pcolor(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,FSP*100.0);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Solid Volume Fraction [%]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
shading interp
colorbar;

%Cu concentration distribution
subplot(2,3,3)
pcolor(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,CLP*100.0);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Cu Concentration [wt%]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
shading interp
colorbar;

%density distribution
subplot(2,3,4)
pcolor(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,RMP);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Mean Density [kg/m^3]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
shading interp
colorbar;

%velocity field
VXPP=zeros(NIY,NIX);
VYPP=zeros(NIY,NIX);
for i=1:NIX
    for j=1:NIY
        VXPP(j,i)=0.5*(VXP(j+1,i)+VXP(j+1,i+1));
    end
end
for i=1:NIX
    for j=1:NIY
        VYPP(j,i)=0.5*(VYP(j,i+1)+VYP(j+1,i+1));
    end
end

subplot(2,3,5)
vxy=quiver([0.5*dxI*1000.0:dxI*1000.0:0.5*W-0.5*dxI*1000.0],[0.5*dyI*1000.0:dyI*1000.0:H-0.5*dyI*1000.0],VXPP,VYPP);
%vxy.AutoScaleFactor=0.1;
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Velocity [m/sec]');
axis([0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
colorbar;

%pressure field
subplot(2,3,6)
pcolor(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,PP);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Pressure [Pa]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
shading interp
colorbar;

drawnow
end

