function PlotTC(TP,FSP,CLP,RMP,Time,Step)
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

%All T field
subplot(2,3,1)
xx=zeros(NIX+NMX+2,1);
yy=zeros(NIY+NMY+2,1);
xx(1)=0.0;
xx(2)=0.5*dxI;
xx(3:NIX+1)=xx(2)+[1:NIX-1]*dxI;
xx(NIX+2)=xx(NIX+1)+0.5*dxI+0.5*dxM+gap;
xx(NIX+3:NIX+1+NMX)=xx(NIX+2)+[1:NMX-1]*dxM;
xx(NIX+2+NMX)=xx(NIX+1+NMX)+0.5*dxM;

yy(1)=0.0;
yy(2)=0.5*dyI;
yy(3:NIY+1)=yy(2)+[1:NIY-1]*dyI;
yy(NIY+2)=yy(NIY+1)+0.5*dyI+0.5*dyM+gap;
yy(NIY+3:NIY+1+NMY)=yy(NIY+2)+[1:NMY-1]*dyM;
yy(NIY+1+NMY+1)=yy(NIY+1+NMY)+0.5*dyM;

contourf(xx*1000.0,yy*1000.0,TP);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
titles=['Temperature [K]      Time=',num2str(Time),' sec       Step=',num2str(Step)];
title(titles);
axis([0.0 x(NIX+NMX+2)*1000.0 0.0 (y(NIY+NMY+2)-gap)*1000.0]);
%shading interp
colorbar;

%FS distribution
subplot(2,3,2)
contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,FSP*100.0);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Solid Volume Fraction [%]');
axis([0.0 x(NIX+2)*1000.0 0.0 (y(NIY+2)-gap)*1000.0]);
%shading interp
colorbar;

%Cu concentration distribution
subplot(2,3,3)
contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,CLP*100.0);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Cu Concentration [wt%]');
axis([0.0 x(NIX+2)*1000.0 0.0 (y(NIY+2)-gap)*1000.0]);
%shading interp
colorbar;

%density distribution
subplot(2,3,4)
contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,RMP);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Mean Density [kg/m^3]');
axis([0.0 x(NIX+2)*1000.0 0.0 (y(NIY+2)-gap)*1000.0]);
%shading interp
colorbar;

drawnow

end

