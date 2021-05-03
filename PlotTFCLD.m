function PlotTFCLD(VXP,VYP,PP)
%Plot all field variables: T, FL, CL, rho, VX, VY, P, dFS 
%suffix 'P'-->Plot
%Created 2019-10-17

%Modified for heat and species diffusion verification! 2019-11-8

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
global gap

W=1.0;%ingot width [mm]
H=1.0;%ingot height [mm]


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

subplot(2,2,1)
x1=[0.5*dxI:dxI:W-0.5*dxI];
y1=[0.5*dyI:dyI:H-0.5*dyI]';
quiver(x1,y1,VXPP,VYPP);
%vxy.AutoScaleFactor=0.1;
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Velocity [m/sec]');
axis([0.0 1.0 0.0 1.0]);
colorbar;

%pressure field
subplot(2,2,2)
pcolor(x(1:NIX+2),y(1:NIY+2),PP);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Pressure [Pa]');
axis([0.0 1.0 0.0 1.0]);
shading interp
colorbar;

subplot(2,2,3)
[x2,y2]=meshgrid(x(2:NIX+1),y(2:NIY+1));
starty=0.0:0.05:1.0;
startx=ones(size(starty))*0.85;
streamline(x2,y2,VXPP,VYPP,startx,starty);

xlswrite('VY.xlsx',VYPP);
xlswrite('VX.xlsx',VXPP);
xlswrite('PP.xlsx',PP);
drawnow

end

