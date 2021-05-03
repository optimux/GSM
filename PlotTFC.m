function PlotTFC(TP,FSP,CLP,CSMP,RMP,VXP,VYP,PRP,Time,Step)
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
global dx
global dy
global x
global y
global gap

%% ------------------------ Scalar Field ----------------------------------
figure(1)
% %All T field
% subplot(3,3,1)
% xx=zeros(NIX+NMX+2,1);
% yy=zeros(NIY+NMY+2,1);
% xx(1)=0.0;
% xx(2)=0.5*dxI;
% xx(3:NIX+1)=xx(2)+[1:NIX-1]*dxI;
% xx(NIX+2)=xx(NIX+1)+0.5*dxI+0.5*dxM+gap;
% xx(NIX+3:NIX+1+NMX)=xx(NIX+2)+[1:NMX-1]*dxM;
% xx(NIX+2+NMX)=xx(NIX+1+NMX)+0.5*dxM;
%
% yy(1)=0.0;
% yy(2)=0.5*dyI;
% yy(3:NIY+1)=yy(2)+[1:NIY-1]*dyI;
% yy(NIY+2)=yy(NIY+1)+0.5*dyI+0.5*dyM+gap;
% yy(NIY+3:NIY+1+NMY)=yy(NIY+2)+[1:NMY-1]*dyM;
% yy(NIY+1+NMY+1)=yy(NIY+1+NMY)+0.5*dyM;
%
% contourf(xx*1000.0,yy*1000.0,TP);
% axis ij
% xlabel('X [mm]');
% ylabel('Y [mm]');
% titles=['Temperature [K]      Time=',num2str(Time),' sec       Step=',num2str(Step)];
% title(titles);
% axis([0.0 x(NIX+NMX+2)*1000.0 0.0 (y(NIY+NMY+2)-gap)*1000.0]);
% %shading interp
% colorbar;
%
% %FS distribution
% subplot(3,3,2)
% contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,FSP*100.0);
% axis ij
% xlabel('X [mm]');
% ylabel('Y [mm]');
% title('Solid Volume Fraction [%]');
% axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
% %shading interp
% colorbar;
%
% %Cu concentration distribution
% subplot(3,3,3)
% contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,CLP*100.0);
% axis ij
% xlabel('X [mm]');
% ylabel('Y [mm]');
% title('Cu Concentration (Liq) [wt%]');
% axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
% %shading interp
% colorbar;
%
% %density distribution
% subplot(3,3,4)
% contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,RMP);
% axis ij
% xlabel('X [mm]');
% ylabel('Y [mm]');
% title('Mean Density [kg/m^3]');
% axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
% %shading interp
% colorbar;
%
% %old mean solid concentration
% subplot(3,3,5)
% contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,CSMP);
% axis ij
% xlabel('X [mm]');
% ylabel('Y [mm]');
% title('Cu Concentration (S) [kg/m^3]');
% axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
% %shading interp
% colorbar;
%
% %pressure field
% subplot(3,3,6)
% contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,PRP);
% axis ij
% xlabel('X [mm]');
% ylabel('Y [mm]');
% title('P-P_0 [Pa]');
% axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
% %shading interp
% colorbar;
% %
% %T in ingot
% subplot(3,3,7)
% TPI=TP(1:NIY+2,1:NIX+2);
% TPI(1:NIY+1,NIX+2)=TPI(1:NIY+1,NIX+1);
% TPI(NIY+2,1:NIX+2)=TPI(NIY+1,1:NIX+2);
% contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,TPI);
% axis ij
% xlabel('X [mm]');
% ylabel('Y [mm]');
% title('T [K]');
% axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
% %shading interp
% colorbar;
%
% %% ----------------- Pressure & Gradient Smoothness Check -----------------
% PRGX=zeros(NIY,NIX-1);%x-axis relative pressure gradient [Pa/m]
% PRGY=zeros(NIY-1,NIX);%x-axis relative pressure gradient [Pa/m]
%
% for i=1:NIX-1
%     for j=1:NIY
%       PRGX(j,i)=2.0*(PRP(j+1,i+2)-PRP(j+1,i+1))/(dx(i)+dx(i+1));
%     end
% end
%
% for i=1:NIX
%     for j=1:NIY-1
%       PRGY(j,i)=2.0*(PRP(j+2,i+1)-PRP(j+1,i+1))/(dy(j)+dy(j+1));
%     end
% end
%
% subplot(3,3,8)
% pcolor([dxI:dxI:NIX*dxI-dxI],[0.5*dyI:dyI:NIY*dyI-0.5*dyI],PRGX);
% axis ij
% xlabel('X [mm]');
% ylabel('Y [mm]');
% title('(P-P_0)_X  [Pa/m]');
% shading interp
% colorbar;
%
% subplot(3,3,9)
% pcolor([0.5*dxI:dxI:NIX*dxI-0.5*dxI],[dyI:dyI:NIY*dyI-dyI],PRGY);
% axis ij
% xlabel('X [mm]');
% ylabel('Y [mm]');
% title('(P-P_0)_Y  [Pa/m]');
% shading interp
% colorbar;

%==========================================================================

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
FSlevel=[0.0,0.003,0.01,0.1,0.35,0.65,0.8];
contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,FSP*100.0,FSlevel);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Solid Volume Fraction [%]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
%shading interp
colorbar;

FSXu=xlsread('G:\PhD\Img\Xu1991TFS.xlsx','Xu1991FS','A1:B376');%in Celsus
hold on
plot(FSXu(:,1),200.0-FSXu(:,2),'k-.');%200 is height of Xu1991
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
hold off

%Cu concentration distribution
subplot(2,3,3)
contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,CLP*100.0);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Cu Concentration (Liq) [wt%]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
%shading interp
colorbar;

%old mean solid concentration
subplot(2,3,4)
contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,CSMP);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Cu Concentration (S) [wt%]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
%shading interp
colorbar;

%pressure field
subplot(2,3,5)
contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,PRP);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('P-P_0 [Pa]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
%shading interp
colorbar;

%T in ingot
subplot(2,3,6)
Tlevel=[590.0,605.0,635.0,647.0,647.6,648.0,655.0];
TPI=TP(1:NIY+2,1:NIX+2);
TPI(1:NIY+1,NIX+2)=TPI(1:NIY+1,NIX+1);
TPI(NIY+2,1:NIX+2)=TPI(NIY+1,1:NIX+2);
contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,TPI-273.15,Tlevel);
axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('T [K]');
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
%shading interp
colorbar;

TXu=xlsread('G:\PhD\Img\Xu1991TFS.xlsx','Xu1991T','A1:B376');%in Celsus
hold on
plot(TXu(:,1),200.0-TXu(:,2),'k-.');%200 is height of Xu1991
axis([0.0 x(NIX+2)*1000.0 0.0 y(NIY+2)*1000.0]);
hold off

% figure(2)
% %% ------------------------- Velocity Field -------------------------------
% %velocity field
% VXPP=zeros(NIY,NIX);
% VYPP=zeros(NIY,NIX);
% for i=1:NIX
%     for j=1:NIY
%         VXPP(j,i)=0.5*(VXP(j+1,i)+VXP(j+1,i+1));
%     end
% end
% vxmax=max(max(abs(VXPP)));
% vxmin=min(min(abs(VXPP)));
% x_scale=0.5*(vxmax+vxmin);
%
% for i=1:NIX
%     for j=1:NIY
%         VYPP(j,i)=0.5*(VYP(j,i+1)+VYP(j+1,i+1));
%     end
% end
% vymax=max(max(abs(VYPP)));
% vymin=min(min(abs(VYPP)));
% y_scale=0.5*(vymax+vymin);
%
% scale=sqrt(x_scale^2+y_scale^2);
% v_max=sqrt(vxmax^2+vymax^2);
%
% refvx=zeros(1,NIX);
% refvx(3)=scale;
% refvy=zeros(1,NIX);
% VX=zeros(NIY+1,NIX);
% VX(1,:)=refvx;
% VX(2:NIY+1,:)=VXPP;
% VY=zeros(NIY+1,NIX);
% VY(1,:)=refvy;
% VY(2:NIY+1,:)=VYPP;
%
% ya=y;
% ya(1)=0.0-0.5*dyI;
%
% quiver(x(2:NIX+1)*1000.0,ya(1:NIY+1)*1000.0,VX,VY);
% axis ij
% xlabel('X [mm]');
% ylabel('Y [mm]');
% title('Velocity [m/sec]');
% axis([0 x(NIX+2)*1000.0 -2.0*dyI*1000.0 y(NIY+2)*1000.0]);
% %scaletext=['Scale: ',num2str(scale),' m/s'];
% scaletext=['V_m_a_x: ',num2str(v_max,'%E'),' m/s'];
% text(2.5*dxI*1000.0,-1.2*dyI*1000.0,scaletext);

%% -------------------------- Velocity + FS -------------------------------
figure(3)
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
v_max=sqrt(vxmax^2+vymax^2);

refvx=zeros(1,NIX);
refvx(3)=scale;
refvy=zeros(1,NIX);
VX=zeros(NIY+1,NIX);
VX(1,:)=refvx;
VX(2:NIY+1,:)=VXPP;
VY=zeros(NIY+1,NIX);
VY(1,:)=refvy;
VY(2:NIY+1,:)=VYPP;

ya=y;
ya(1)=0.0-0.5*dyI;

contourf(x(1:NIX+2)*1000.0,y(1:NIY+2)*1000.0,FSP*100.0);
colorbar;
hold on
quiver(x(2:NIX+1)*1000.0,ya(1:NIY+1)*1000.0,VX,VY,'r');

axis ij
xlabel('X [mm]');
ylabel('Y [mm]');
title('Velocity [m/sec]');
axis([0 x(NIX+2)*1000.0 -2.0*dyI*1000.0 y(NIY+2)*1000.0]);
%scaletext=['Scale: ',num2str(scale),' m/s'];
scaletext=['V_m_a_x: ',num2str(v_max,'%E'),' m/s'];
text(2.5*dxI*1000.0,-1.2*dyI*1000.0,scaletext);

hold off

drawnow

%print VFS.eps -depsc2 -r600

%% -------------------------- GNUPLOT Velocity -------------------------------
% %velocity field
% VXPP=zeros(NIY,NIX);
% VYPP=zeros(NIY,NIX);
% for i=1:NIX
%     for j=1:NIY
%         VXPP(j,i)=0.5*(VXP(j+1,i)+VXP(j+1,i+1));
%     end
% end
% vxmax=max(max(abs(VXPP)));
% vxmin=min(min(abs(VXPP)));
% x_scale=0.5*(vxmax+vxmin);
% 
% for i=1:NIX
%     for j=1:NIY
%         VYPP(j,i)=0.5*(VYP(j,i+1)+VYP(j+1,i+1));
%     end
% end
% vymax=max(max(abs(VYPP)));
% vymin=min(min(abs(VYPP)));
% y_scale=0.5*(vymax+vymin);
% 
% scale=sqrt(x_scale^2+y_scale^2);
% v_max=sqrt(vxmax^2+vymax^2);
% 
% 
% TPI=TP(1:NIY+2,1:NIX+2);
% TPI(1:NIY+1,NIX+2)=TPI(1:NIY+1,NIX+1);
% TPI(NIY+2,1:NIX+2)=TPI(NIY+1,1:NIX+2);
% 
% fid=fopen('G:\PhD\Img\CuAlVT.txt','w');
% for i=1:NIX+2
%     for j=1:NIY+2
%         if(i==1||i==NIX+2)%first and last column, only T in celsus
%             fprintf(fid,'%f %f %f\n',x(i)*1000.0,y(j)*1000.0,TPI(j,i)-273.15);%x,y in meter, TPI in celsus degree
%         elseif(j==1||j==NIY+2)
%             fprintf(fid,'%f %f %f\n',x(i)*1000.0,y(j)*1000.0,TPI(j,i)-273.15);%x,y in meter, TPI in celsus degree
%         else
%             fprintf(fid,'%f %f %f %f %f\n',x(i)*1000.0,y(j)*1000.0,TPI(j,i)-273.15,4.0*VXPP(j-1,i-1)/scale,4.0*VYPP(j-1,i-1)/scale);%x,y in meter
%         end
%     end
%     fprintf(fid,'\n');
% end
% fclose(fid);
% 
% 
% fid=fopen('G:\PhD\Img\CuAlVP.txt','w');
% for i=1:NIX+2
%     for j=1:NIY+2
%         if(i==1||i==NIX+2)%first and last column, only T in celsus
%             fprintf(fid,'%f %f %f\n',x(i)*1000.0,y(j)*1000.0,PRP(j,i));%x,y in meter, TPI in celsus degree
%         elseif(j==1||j==NIY+2)
%             fprintf(fid,'%f %f %f\n',x(i)*1000.0,y(j)*1000.0,PRP(j,i));%x,y in meter, TPI in celsus degree
%         else
%             fprintf(fid,'%f %f %f %f %f\n',x(i)*1000.0,y(j)*1000.0,PRP(j,i),4.0*VXPP(j-1,i-1)/scale,4.0*VYPP(j-1,i-1)/scale);%x,y in meter
%         end
%     end
%         fprintf(fid,'\n');
% end
% fclose(fid);

end

