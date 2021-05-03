function FieldsPlot(TP,FSP,CLP,RMP,VXP,VYP,PP,t,n)
%Plot all field variables: T, FL, CL, rho, VX, VY, P, dFS
%suffix 'P'-->Plot
%Created 2019-10-17

%Modified 2020-1-1

global NIX
global NIY
global dxI
global dyI
global x
global y

W=120.0;
H=100.0;

%All T field
subplot(2,3,1)
pcolor(x*1000.0,y*1000.0,TP);%TP(1.5,:)=0.5(TP(1,:)+TP(2,:))=TP(1,:)=TP(2,:); hereinafter
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
title(['Temperature [K]   Time=',num2str(t),' sec','   Step= ',num2str(n)]);
shading interp
colorbar;

%FS distribution
subplot(2,3,2)
pcolor(x*1000.0,y*1000.0,FSP*100.0);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Solid Volume Fraction [%]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
shading interp
colorbar;

%Cu concentration distribution
subplot(2,3,3)
pcolor(x*1000.0,y*1000.0,CLP*100.0);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Cu Concentration [wt%]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
shading interp
colorbar;

%density distribution
subplot(2,3,4)
pcolor(x*1000.0,y*1000.0,RMP);
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
vxmax=max(max(abs(VXPP)));
vxmin=min(min(abs(VXPP)));
x_scale=0.5*(vxmax+vxmin);

for i=1:NIX
    for j=1:NIY
        VYPP(j,i)=0.5*(VYP(j,i+1)+VYP(j+1,i+1));
    end
end
vymax=max(max(abs(VYPP)));
vymin=min(min(abs(VYPP)));
y_scale=0.5*(vymax+vymin);

scale=sqrt(x_scale^2+y_scale^2);
refvx=zeros(1,NIX);
refvx(3)=scale;
refvy=zeros(1,NIX);
VX=zeros(NIY+1,NIX);
VX(1,:)=refvx;
VX(2:NIY+1,:)=VXPP;
VY=zeros(NIY+1,NIX);
VY(1,:)=refvy;
VY(2:NIY+1,:)=VYPP;

subplot(2,3,5)
quiver([0.5*dxI*1000.0:dxI*1000.0:W-0.5*dxI*1000.0],[-0.5*dyI*1000.0:dyI*1000.0:H-0.5*dyI*1000.0],VX,VY);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Velocity [m/sec]');
axis([0 x(NIX+2)*1000.0 -2.0*dyI*1000.0 y(NIY+2)*1000.0]);
scaletext=['Velocity scale: ',num2str(scale),' m/s'];
text(0.5*dxI*1000.0,-1.4*dyI*1000.0,scaletext);

%pressure field
subplot(2,3,6)
pcolor(x*1000.0,y*1000.0,PP);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Pressure [Pa]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
shading interp
colorbar;

drawnow

end

