VXP=xlsread('VX.xlsx');
VYP=xlsread('VY.xlsx');

W=1.0;%ingot width [m]
H=1.0;%ingot height [m]
[rows,NIX]=size(VXP);
[NIY,cols]=size(VYP);%number of unit in ingot, y-axis
dxI=W/NIX;%%x-axis interval in ingot [m]
dyI=H/NIY;%y-axis interval in ingot [m]

x=zeros(NIX+2,1);
x(1)=0.0;
x(2)=0.5*dxI;
x(3:NIX+1)=x(2)+[1:NIX-1]*dxI;
x(NIX+2)=x(NIX+1)+0.5*dxI;

y=zeros(NIY+2,1);
y(1)=0.0;
y(2)=0.5*dyI;
y(3:NIY+1)=y(2)+[1:NIY-1]*dyI;
y(NIY+2)=y(NIY+1)+0.5*dyI;

figure
[x2,y2]=meshgrid(x(2:NIX+1),y(2:NIY+1));
starty=0.0:0.05:1.0;
startx=ones(size(starty))*0.85;
streamline(x2,y2,VXP,VYP,startx,starty);


figure
contour(x2,y2,VXP,100);
axis([0.0 1.0 0.0 1.0]);
axis ij;

figure
quiver([0.5*dxI:dxI:W-0.5*dxI],[0.5*dyI:dyI:H-0.5*dyI],VXP,VYP);
%vxy.AutoScaleFactor=0.1;
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Velocity [m/sec]');
axis([0 x(NIX+2) 0.0 y(NIY+2)]);
colorbar;

