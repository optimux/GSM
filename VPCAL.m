function [VX,VY,P,iters] = VPCAL(VXOTM,VYOTM,POTM,FLOTM,FLNTM,RLOTM,RLNTM,FROTM,FRNTM,RSOTM,RSNTM,dFSTM,dtb)
%To solve V-P momentum conservation
%Created 2019-10-20

%O==Old
%N==New

global NIX
global NIY
global mu
global dx
global dy
global g

%======================= VERY IMPORTANT NOTE ==============================
%It is highly recommended the usage of packages for complex formulae, i.e.,
%temporary variables. These packages can be utilized repetedly, easily
%modified, and produce a concise numerical equation.
%======================= VERY IMPORTANT NOTE ==============================


%% ................======== BASIC VELOCITY FORMULAE ========...............

%----------------------------- 1. RFVX ------------------------------------
%main grid (j,k)
RFVXI=zeros(NIY+2,NIX+2);%I==integer grid(main grid) 
for i=2:NIX+1
    %VX(1:NIY+2,1)=0.0 && VX(1:NIY+1,image_1)=0.0 --> RFVXI(1:NIY+2,1)=0.0 [image_1==0]
    %VX(1:NIY+2,NIX+1)=0.0 && VX(1:NIY+2,image_2)=0.0 --> RFVXI(1:NIY+2,NIX+2)=0.0 [image_2==NIX+2]
    %VX(1,image_1:image_2)=0.0 --> RFVXI(1,1:NIX+2)=0.0 [No slip boundary]
    %VX(NIY+2,image_1:image_2)=0.0 --> RFVXI(NIY+2,1:NIX+2)=0.0 [No slip boundary]
    for j=2:NIY+1
        RFVXI(j,i)=RLOTM(j,i)*FLOTM(j,i)*0.5*(VXOTM(j,i-1)+VXOTM(j,i));%[kg/m^2/sec]
    end
end

%staggered grid (j+0.5,k+0.5), (j-0.5,k+0.5)
RFVXH=zeros(NIY+1,NIX+1);%H==half grid(staggered grid) 
for i=1:NIX+1
    for j=1:NIY+1
        RFVXH(j,i)=0.25*(RLOTM(j,i)+RLOTM(j+1,i)+RLOTM(j,i+1)+RLOTM(j+1,i+1))*0.25*(FLOTM(j,i)+FLOTM(j+1,i)+FLOTM(j,i+1)+FLOTM(j+1,i+1))*0.5*(VXOTM(j+1,i)+VXOTM(j,i));%[kg/m^2/sec]
    end
end
%----------------------------- 1. RFVX ------------------------------------

%----------------------------- 2. RFVY ------------------------------------
%main grid (j,k), (j,k+1)
RFVYI=zeros(NIY+2,NIX+2);%I==integer grid(main grid)
for i=2:NIX+1
    for j=2:NIY+1
        %VY(1,1:NIX+2)=0.0 && VY(image_3,1:NIX+2)=0.0 --> RFVYI(1,1:NIX+2)=0.0 [image_3==0]
        %VY(NIY+1,1:NIX+2)=0.0 && VY(image_4,1:NIX+2)=0.0 --> RFVYI(NIY+1,1:NIX+2)=0.0 [image_4==NIY+2]
        %VY(1:NIY+2,1)=Free --> RFVYI(1:NIY+2,1)!=0.0 [Slip boundary]
        %VY(1:NIY+2,NIX+2)=0.0 --> RFVYI(1:NIY+2,NIX+2)=0.0 [No slip boundary]
        RFVYI(j,i)=RLOTM(j,i)*FLOTM(j,i)*0.5*(VYOTM(j-1,i)+VYOTM(j,i));%[kg/m^2/sec]
    end
end

%VY(1:NIY+2,1)=Free --> RFVYI(1:NIY+1,1)!=0.0 [Slip boundary]
%This sentence can be added to last FOR loop
for i=2:NIY+1
    RFVYI(i,1)=RLOTM(i,1)*FLOTM(i,1)*0.5*(VYOTM(j-1,1)+VYOTM(j,1));%[kg/m^2/sec]
end

%staggered grid (j+0.5,k+0.5), (j+0.5,k-0.5)
RFVYH=zeros(NIY+1,NIX+1);%H==half grid(staggered grid)
for i=1:NIX+1
    for j=1:NIY+1
        RFVYH(j,i)=0.25*(RLOTM(j,i)+RLOTM(j+1,i)+RLOTM(j,i+1)+RLOTM(j+1,i+1))*0.25*(FLOTM(j,i)+FLOTM(j+1,i)+FLOTM(j,i+1)+FLOTM(j+1,i+1))*0.5*(VYOTM(j,i+1)+VYOTM(j,i));%[kg/m^2/sec]
    end
end
%----------------------------- 2. RFVY ------------------------------------

%---------------------- 3. x-flow momentum --------------------------------
%Upwind scheme for x-axis mass flow
VVXI=zeros(NIY+2,NIX);%I==integer grid(main grid)
for i=1:NIX
    for j=1:NIY+2
        VVXI(j,i)=VXOTM(j,i)*max(RFVXI(j,i+1),0.0)+VXOTM(j,i+1)*min(RFVXI(j,i+1),0.0);%including VVX1, VVX2 [kg/m/sec^2]
    end
end

%Upwind scheme for y-axis mass flow
VVXH=zeros(NIY+1,NIX+1);%H==half grid(staggered grid)
for i=1:NIX+1
    for j=1:NIY+1
        VVXH(j,i)=VXOTM(j,i)*max(RFVYH(j,i),0.0)+VXOTM(j+1,i)*min(RFVYH(j,i),0.0);%[kg/m/sec^2]
    end
end
%---------------------- 3. x-flow momentum --------------------------------

%---------------------- 4. y-flow momentum --------------------------------
%Upwind scheme for x-axis mass flow
VVYI=zeros(NIY,NIX+2);%I==integer grid(main grid)
for i=1:NIX+2
    for j=1:NIY
        VVYI(j,i)=VYOTM(j,i)*max(RFVYI(j+1,i),0.0)+VYOTM(j+1,i)*min(RFVYI(j+1,i),0.0);%[kg/m/sec^2]
    end
end

%Upwind scheme for y-axis mass flow
VVYH=zeros(NIY+1,NIX+1);%H==half grid(staggered grid)
for i=1:NIX+1
    for j=1:NIY+1
        VVYH(j,i)=VYOTM(j,i)*max(RFVXH(j,i),0.0)+VYOTM(j,i+1)*min(RFVXH(j,i),0.0);%[kg/m/sec^2]
    end
end
%---------------------- 4. y-flow momentum --------------------------------

%-------------------- 5. x-viscous momentum -------------------------------
VDXI=zeros(NIY+2,NIX);
for i=1:NIX
    for j=1:NIY+2
        VDXI(j,i)=mu*(VXOTM(j,i+1)*0.5*(FLOTM(j,i+1)+FLOTM(j,i+2))-VXOTM(j,i)*0.5*(FLOTM(j,i)+FLOTM(j,i+1)))/dx(i);%[Pa.sec times 1/sec == Pa]
    end
end
%VDXI(1:NIY+2,image_1)=0.0 at image_1 and image_2
%VDXI(1:NIY+2,image_2)=0.0 at image_1 and image_2

VDXH=zeros(NIY+1,NIX+1);
for i=1:NIX+1
    for j=2:NIY
        VDXH(j,i)=2.0*mu*(0.5*(FLOTM(j+1,i)+FLOTM(j+1,i+1))*VXOTM(j+1,i)-0.5*(FLOTM(j,i)+FLOTM(j,i+1))*VXOTM(j,i))/(dy(j-1)+dy(j));%[Pa.sec times 1/sec == Pa]
    end
end

for i=1:NIX+1
    VDXH(1,i)=2.0*mu*(0.5*(FLOTM(2,i)+FLOTM(2,i+1))*VXOTM(2,i)-0.5*(FLOTM(1,i)+FLOTM(1,i+1))*VXOTM(1,i))/dy(1);%[Pa.sec times 1/sec == Pa]
    VDXH(NIY+1,i)=2.0*mu*(0.5*(FLOTM(NIY+2,i)+FLOTM(NIY+2,i+1))*VXOTM(NIY+2,i)-0.5*(FLOTM(NIY+1,i)+FLOTM(NIY+1,i+1))*VXOTM(NIY+1,i))/dy(NIY);%[Pa.sec times 1/sec == Pa]
end

%-------------------- 5. x-viscous momentum -------------------------------

%-------------------- 6. y-viscous momentum -------------------------------
VDYI=zeros(NIY,NIX+2);
for i=1:NIX+2
    for j=1:NIY
        VDYI(j,i)=mu*(VYOTM(j+1,i)*0.5*(FLOTM(j+2,i)+FLOTM(j+1,i))-VYOTM(j,i)*0.5*(FLOTM(j+1,i)+FLOTM(j,i)))/dy(j);%[Pa.sec times 1/sec == Pa]
    end
end

VDYH=zeros(NIY+1,NIX+1);
for i=2:NIX
    for j=1:NIY+1
        VDYH(j,i)=2.0*mu*(0.5*(FLOTM(j+1,i+1)+FLOTM(j,i+1))*VYOTM(j,i+1)-0.5*(FLOTM(j+1,i)+FLOTM(j,i))*VYOTM(j,i))/(dx(i-1)+dx(i));%[Pa.sec times 1/sec == Pa]
    end
end

for i=1:NIY+1
    VDYH(i,1)=2.0*mu*(0.5*(FLOTM(i+1,2)+FLOTM(i,2))*VYOTM(i,2)-0.5*(FLOTM(i+1,1)+FLOTM(i,1))*VYOTM(i,1))/dx(1);%[Pa.sec times 1/sec == Pa]
    VDYH(i,NIX+1)=2.0*mu*(0.5*(FLOTM(i+1,NIX+2)+FLOTM(i,NIX+2))*VYOTM(i,NIX+2)-0.5*(FLOTM(i+1,NIX+1)+FLOTM(i,NIX+1))*VYOTM(i,NIX+1))/dx(NIX);%[Pa.sec times 1/sec == Pa]
end

%-------------------- 6. y-viscous momentum -------------------------------

%--------------------- 7. x-convection part -------------------------------
VXCONV=zeros(NIY,NIX+1);
for i=2:NIX
    for j=1:NIY
        VXCONV(j,i)=-2.0*dtb*(VVXI(j+1,i)-VVXI(j+1,i-1))/(dx(i)+dx(i-1))-dtb*(VVXH(j+1,i)-VVXH(j,i))/dy(j);%[kg/m^2/sec]
    end
end

for i=1:NIY
    VXCONV(i,1)=0.0-dtb*(VVXH(i+1,1)-VVXH(i,1))/dy(i);%VVXI have no value at image_1 sites while RFVXI=0.0 at image_1 sites [kg/m^2/sec]
    VXCONV(i,NIX+1)=0.0-dtb*(VVXH(i+1,NIX+1)-VVXH(i,NIX+1))/dy(i);%%VVXI have no value at image_2 sites while RFVXI=0.0 at image_2 sites [kg/m^2/sec]
end
%--------------------- 7. x-convection part -------------------------------

%---------------------- 8. x-diffusion part -------------------------------
VXDIF=zeros(NIY,NIX+1);
for i=2:NIX
    for j=1:NIY
        VXDIF(j,i)=2.0*dtb*(VDXI(j+1,i)-VDXI(j+1,i-1))/(dx(i)+dx(i-1))+dtb*(VDXH(j+1,i)-VDXH(j,i))/dy(j);%[Pa.sec/m]
    end
end

for i=1:NIY
    VXDIF(i,1)=2.0*dtb*(VDXI(i+1,1)-0.0)/dx(1)+dtb*(VDXH(i+1,1)-VDXH(i,1))/dy(i);%[Pa.sec/m]
    VXDIF(i,NIX+1)=2.0*dtb*(0.0-VDXI(i+1,NIX))/dx(NIX)+dtb*(VDXH(j+1,NIX+1)-VDXH(j,NIX+1))/dy(i);%[Pa.sec/m]
end
%---------------------- 8. x-diffusion part -------------------------------

%---------------------- 9. y-convection part ------------------------------
VYCONV=zeros(NIY+1,NIX);
for i=1:NIX
    for j=2:NIY
        VYCONV(j,i)=-dtb*(VVYH(j,i+1)-VVYH(j,i))/dx(i)-2.0*dtb*(VVYI(j,i+1)-VVYI(j-1,i+1))/(dy(j)+dy(j-1));%[Pa.sec/m]
    end
end

for i=1:NIX
    VYCONV(1,i)=-dtb*(VVYH(1,i+1)-VVYH(1,i))/dx(i)-2.0*dtb*(VVYI(1,i+1)-0.0)/dy(1);%VVYI(imag_3,1:NIX+2)=0.0 have no image_3 sites [Pa.sec/m]
    VYCONV(NIY+1,i)=-dtb*(VVYH(NIY+1,i+1)-VVYH(NIY+1,i))/dx(i)-2.0*dtb*(0.0-VVYI(NIY,i+1))/dy(NIY);%VVYI(imag_4,1:NIX+2)=0.0 have no image_4 sites [Pa.sec/m]
end

%---------------------- 9. y-convection part ------------------------------

%---------------------- 10. y-diffusion part ------------------------------
VYDIF=zeros(NIY+1,NIX);
for i=1:NIX
    for j=2:NIY
        VYDIF(j,i)=dtb*(VDYH(j,i+1)-VDYH(j,i))/dx(i)+2.0*dtb*(VDYI(j,i+1)-VDYI(j-1,i+1))/(dy(j)+dy(j-1));%[Pa.sec/m]
    end
end

for i=1:NIX
    VYDIF(1,i)=dtb*(VDYH(1,i+1)-VDYH(1,i))/dx(i)+2.0*dtb*(VDYI(1,i+1)-0.0)/dy(1);%[Pa.sec/m]
    VYDIF(NIY+1,i)=dtb*(VDYH(NIY+1,i+1)-VDYH(NIY+1,i))/dx(i)+2.0*dtb*(0.0-VDYI(NIY,i))/dy(NIY);%[Pa.sec/m]
end
%---------------------- 10. y-diffusion part ------------------------------

%-------------------- 11. x- pressure gradient ----------------------------
FPDX=zeros(NIY+2,NIX+1);%D==difference
for i=2:NIX
    for j=1:NIY+2
        FPDX(j,i)=-2.0*dtb*(FLOTM(j,i+1)*POTM(j,i+1)-FLOTM(j,i)*POTM(j,i))/(dx(i)+dx(i-1));%[Pa.sec/m == kg/m^2/sec]
    end
end

for i=1:NIY+2
    FPDX(i,1)=-2.0*dtb*(FLOTM(i,2)*POTM(i,2)-FLOTM(i,1)*POTM(i,1))/dx(1);%[Pa.sec/m == kg/m^2/sec]
    FPDX(i,NIX+1)=-2.0*dtb*(FLOTM(i,NIX+2)*POTM(i,NIX+2)-FLOTM(i,NIX+1)*POTM(i,NIX+1))/dx(NIX);%[Pa.sec/m == kg/m^2/sec]
end

%-------------------- 11. x- pressure gradient ----------------------------

%---------------------- 12. y- pressure part ------------------------------
FPDY=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=2:NIY
        FPDY(j,i)=-2.0*dtb*(FLOTM(j+1,i)*POTM(j+1,i)-FLOTM(j,i)*POTM(j,i))/(dy(j)+dy(j-1));%[Pa.sec/m == kg/m^2/sec]
    end
end

for i=1:NIX+2
    FPDY(1,i)=-2.0*dtb*(FLOTM(2,i)*POTM(2,i)-FLOTM(1,i)*POTM(1,i))/dy(1);%[Pa.sec/m == kg/m^2/sec]
    FPDY(NIY+1,i)=-2.0*dtb*(FLOTM(NIY+2,i)*POTM(NIY+2,i)-FLOTM(NIY+1,i)*POTM(NIY+1,i))/dy(NIY);%[Pa.sec/m == kg/m^2/sec]
end

%---------------------- 12. y- pressure part ------------------------------

%----------------------- 13. y- gravity part ------------------------------
FRG=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=1:NIY+1
        FRG(j,i)=0.5*g*dtb*(0.5*(FROTM(j,i)+FROTM(j+1,i))+0.5*(FRNTM(j,i)+FRNTM(j+1,i)));%[kg/m^2/sec]
    end
end

%----------------------- 13. y- gravity part ------------------------------

%--------------------------- 14. Old RFVX ---------------------------------
RFVX=zeros(NIY+2,NIX+1);
for i=1:NIX+1
    for j=1:NIY+2
        RFVX(j,i)=0.5*(RLOTM(j,i)+RLOTM(j,i+1))*0.5*(FLOTM(j,i)+FLOTM(j,i+1))*VXOTM(j,i);%[kg/m^2/sec]
    end
end

%--------------------------- 14. Old RFVX ---------------------------------

%--------------------------- 15. Old RFVY ---------------------------------
RFVY=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=1:NIY+1
        RFVY(j,i)=0.5*(RLOTM(j,i)+RLOTM(j+1,i))*0.5*(FLOTM(j,i)+FLOTM(j+1,i))*VYOTM(j,i);%[kg/m^2/sec]
    end
end

%--------------------------- 15. Old RFVY ---------------------------------


%% ...............======== DFLP COEFFICIENT MATRIX ========................

%--------------------------- Permeability K -------------------------------
FSNTM=1.0-FLNTM;
KFLTM=zeros(NIY+2,NIX+2);%permeability from D-17
for i=1:NIX+2
    for j=1:NIY+2
        if(FLNTM(j,i)>=1.0/3.0)
            KFLTM(j,i)=2.6e-5*(1.923e-2*FLNTM(j,i)^2+(4.0+3.0*FSNTM(j,i)-3.0*sqrt(FSNTM(j,i)*(8.0-3.0*FSNTM(j,i))))/FSNTM(j,i));%[mm^2]
        else
            KFLTM(j,i)=5.0e-7*FLNTM(j,i)^2;%[mm^2]
        end
    end
end
KFLTM=KFLTM*10^-6;%[m^2]

%--------------------------- Permeability K -------------------------------

%------------------- mu*FL/(RL*K) in denomenator --------------------------
UFRKY=zeros(NIY+1,NIX+2);%Y==y-axis
for i=1:NIX+2
    for j=1:NIY+1
        UFRKY(j,i)=mu*0.5*(FLNTM(j,i)+FLNTM(j+1,i))/(0.5*(RLNTM(j,i)+RLNTM(j+1,i))*0.5*(KFLTM(j,i)+KFLTM(j+1,i)));%[1/sec]
    end
end

UFRKX=zeros(NIY+2,NIX+1);%X==x-axis
for i=1:NIX+1
    for j=1:NIY+2
        UFRKX(j,i)=mu*0.5*(FLNTM(j,i)+FLNTM(j,i+1))/(0.5*(RLNTM(j,i)+RLNTM(j,i+1))*0.5*(KFLTM(j,i)+KFLTM(j,i+1)));%[1/sec]
    end
end

%------------------- mu*FL/(RL*K) in denomenator --------------------------

%x-axis parameter in matrix
% AXB=zeros(NIY,NIX+1);%B==Backward
% for i=2:NIX
%     for j=1:NIY
%         AXB(j,i)=2.0*dy(j)/((dx(i-1)+dx(i))*(1.0+dtb*UFRKX(j+1,i)));
%     end
% end
%
% for i=1:NIY
%     AXB(i,1)=2.0*dy(i)/(dx(1)*(1.0+dtb*UFRKX(i+1,1)));
%     AXB(i,NIX+1)=2.0*dy(i)/(dx(NIX)*(1.0+dtb*UFRKX(i+1,NIX+1)));
% end

AX=zeros(NIY,NIX+1);
for i=2:NIX
    for j=1:NIY
        AX(j,i)=2.0*dy(j)/((dx(i-1)+dx(i))*(1.0+dtb*UFRKX(j+1,i)));%[1]
    end
end

for i=1:NIY
    AX(i,1)=2.0*dy(i)/(dx(1)*(1.0+dtb*UFRKX(i+1,1)));%[1]
    AX(i,NIX+1)=2.0*dy(i)/(dx(NIX)*(1.0+dtb*UFRKX(i+1,NIX+1)));%[1]
end

%y-axis parameter in matrix
AY=zeros(NIY+1,NIX);
for i=1:NIX
    for j=2:NIY
        AY(j,i)=2.0*dx(i)/((dy(j)+dy(j-1))*(1.0+dtb*UFRKY(j,i+1)));%[1]
    end
end

for i=1:NIX
    AY(1,i)=2.0*dx(i)/(dy(1)*(1.0+dtb*UFRKY(1,i+1)));%[1]
    AY(NIY+1,i)=2.0*dx(i)/(dy(NIY)*(1.0+dtb*UFRKY(NIY+1,i+1)));%[1]
end

%% ...............=========== THE BOMB CYCLE ============..................

err0DRFV=0.001;
errDRFV=1.0;

%NOTE: DFLP at 4 physical walls or boundaries are useless since coefficients are set zero there, thus DFLP has onlt NIX*NIY elements.
DFLP=zeros(NIY*NIX,1);%consider delta(FL*P) as one variable which has (NIX*NIY) elements in one column
DFLPT=zeros(NIY*NIX,1);%T=total, total pressure correction

DRFVX=zeros(NIY+2,NIX+1);
DRFVY=zeros(NIY+1,NIX+2);
VX=zeros(NIY+2,NIX+1);
VY=zeros(NIY+1,NIX+2);
iters=0;

while(errDRFV>err0DRFV)
    
    %===================== APPROXIMATION OF RFVXY =========================
    %--------------------------- RFVXTM -----------------------------------
    RFVXTM=zeros(NIY+2,NIX+1);
    for i=1:NIX+1
        for j=2:NIY+1
            RFVXTM(j,i)=(RFVX(j,i)+VXCONV(j-1,i)+VXDIF(j-1,i)+FPDX(j,i))/(1.0+dtb*UFRKX(j,i));%[kg/m^2/sec]
        end
    end
    %RFVXTM(1,1:NIX+1)=0.0 are available but have no effect on VX
    %RFVXTM(NIY+2,1:NIX+1)=0.0 are available but have no effect on VX
    
    %--------------------------- RFVXTM -----------------------------------
    
    %--------------------------- RFVYTM -----------------------------------
    RFVYTM=zeros(NIY+1,NIX+2);
    for i=2:NIX+1
        for j=1:NIY+1
            RFVYTM(j,i)=(RFVY(j,i)+VYCONV(j,i-1)+VYDIF(j,i-1)+FPDY(j,i)+FRG(j,i))/(1.0+dtb*UFRKY(j,i));%[kg/m^2/sec]
        end
    end
    
    for i=2:NIY
        %1-st column is in slip boundary condition
        RFVYTM(i,1)=(RFVY(i,1)-0.0-2.0*dtb*(VVYI(i,1)-VVYI(i-1,1))/(dy(i)+dy(i-1))+0.0+2.0*dtb*(VDYI(i,1)-VDYI(i-1,1))/(dy(i)+dy(i-1))+FPDY(i,1)+FRG(i,1))/(1.0+dtb*UFRKY(i,1));%[kg/m^2/sec]
        %NOTE: VVYH(i,image_1)=VVYH(i,1)=0.0 since RFVX(j+0.5,k+0.5)==0.0 --> first '0.0' in above formula; 
        %VDYH(i,1)=VDYH(i,image_1)=0.0 since VY(i,1)==VY(i,image_1)==VY(i,image_image) --> second '0.0' in above formula
        
        %last column is in no slip boundary condition
        RFVYTM(i,NIX+2)=(RFVY(i,NIX+2)-0.0-2.0*dtb*(VVYI(i,NIX+2)-VVYI(i-1,NIX+2))/(dy(i)+dy(i-1))+0.0+2.0*dtb*(VDYI(i,NIX+2)-VDYI(i-1,NIX+2))/(dy(i)+dy(i-1))+FPDY(i,NIX+2)+FRG(i,NIX+2))/(1.0+dtb*UFRKY(i,NIX+2));%no x -axis mass flow
    end
    
    %--------------------------- RFVYTM -----------------------------------
    %===================== APPROXIMATION OF RFVXY =========================
    
    %================== MASS CONSERVATION RESIDUAL  =======================
    %----------------------------- RES ------------------------------------
    RES=zeros(NIY,NIX);
    for i=1:NIX
        for j=1:NIY
            RES(j,i)=dy(j)*(RFVXTM(j+1,i+1)-RFVXTM(j+1,i))+dx(i)*(RFVYTM(j+1,i+1)-RFVYTM(j,i+1))+0.5*dFSTM(j+1,i+1)*(RSOTM(j+1,i+1)+RSNTM(j+1,i+1))*dx(i)*dy(j)/dtb+(FRNTM(j+1,i+1)-FROTM(j+1,i+1))*dx(i)*dy(j)/dtb;%[kg/m/sec == Pa.sec]
        end
    end
    
    RESTM=zeros(NIY*NIX,1);%set as [B] which has NIX*NIY elements in one column
    for i=1:NIX
        RESTM(1+(i-1)*NIY:i*NIY)=RES(1:NIY,i)/dtb;%[kg/m/sec == Pa.sec]
    end
    
%     for i=1:NIY*NIX
%     if(RESTM(i)<10^-5)
%         RESTM(i)=0.0;
%     end
%     end
    %----------------------------- RES ------------------------------------
    %=================== MASS CONSERVATION RESIDUAL  ======================
    
    %=========================== [A] MATRIX ===============================
    %NOTE: In [A][X]=[B], [A], [X] and [B] have dimension of [NIY]*[NIX]; to solve this 2D
    %matrix, we reshape [X] and [B] into vector of length [NIX*NIY], this will give a
    %large sparse matrix of [A] of dimension [NIX*NIY]*[NIX*NIY]
    
    %coefficient matrix, a(j,k), a(j+0.5,k), a(j-0.5,k), a(j,k-0.5), a(j,k+0.5)
    A=zeros(NIY*NIX,NIX*NIY);
    for i=2:NIX-1
        for j=2:NIY-1
            k=(i-1)*NIY+j;
            A(k,k-NIY)=-AX(j,i);%[1]
            A(k,k-1)=-AY(j,i);%[1]
            A(k,k)=AX(j,i)+AY(j,i)+AY(j+1,i)+AX(j,i+1);%[1]
            A(k,k+1)=-AY(j+1,i);%[1]
            A(k,k+NIY)=-AX(j,i+1);%[1]
        end
    end
    
    %special boundary condition at [1,1]
    A(1,1)=AX(1,2)+AY(2,1);%[1]
    A(1,2)=-AY(2,1);%[1]
    A(1,NIY+1)=-AX(1,2);%[1]
    
    %special boundary condition at [NIY,1]
    A(NIY,NIY)=AX(NIY,2)+AY(NIY,1);%[1]
    A(NIY,NIY-1)=-AY(NIY,1);%[1]
    A(NIY,2*NIY)=-AX(NIY,2);%[1]
    
    %special boundary condition at [1,NIX]
    A((NIX-1)*NIY+1,(NIX-1)*NIY+1)=AX(1,NIX)+AY(2,NIX);%[1]
    A((NIX-1)*NIY+1,(NIX-1)*NIY+2)=-AY(2,NIX);%[1]
    A((NIX-1)*NIY+1,(NIX-2)*NIY+1)=-AX(1,NIX);%[1]
    
    %special boundary condition at [NIY,NIX]
    A(NIX*NIY,NIX*NIY)=AX(NIY,NIX)+AY(NIY,NIX);%[1]
    A(NIX*NIY,NIX*NIY-1)=-AY(NIY,NIX);%[1]
    A(NIX*NIY,(NIX-1)*NIY)=-AX(NIY,NIX);%[1]
    
    %Left boundary
    for j=2:NIY-1
        A(j,j)=AY(j,1)+AY(j+1,1)+AX(j,2);%[1]
        A(j,j-1)=-AY(j,1);%[1]
        A(j,j+1)=-AY(j+1,1);%[1]
        A(j,j+NIY)=-AX(j,2);%[1]
    end
    
    %Right boundary
    for i=2:NIY-1
        j=(NIX-1)*NIY+i;
        A(j,j)=AY(i,NIX)+AY(i+1,NIX)+AX(i,NIX);%[1]
        A(j,j-1)=-AY(i,NIX);%[1]
        A(j,j+1)=-AY(i+1,NIX);%[1]
        A(j,j-NIY)=-AX(i,NIX);%[1]
    end
    
    %Top boundary
    for i=2:NIX-1
        j=(i-1)*NIY+1;
        A(j,j)=AX(1,i)+AY(2,i)+AX(1,i+1);%[1]
        A(j,j-NIY)=-AX(1,i);%[1]
        A(j,j+1)=-AY(2,i);%[1]
        A(j,j+NIY)=-AX(1,i+1);%[1]
    end
    
    %Bottom boundary
    for i=2:NIX-1
        j=i*NIY;
        A(j,j)=AX(NIY,i)+AY(NIY,i)+AX(NIY,i+1);%[1]
        A(j,j-NIY)=-AX(NIY,i);%[1]
        A(j,j-1)=-AY(NIY,i);%[1]
        A(j,j+NIY)=-AX(NIY,i+1);%[1]
    end
    
    %=========================== [A] MATRIX ===============================
    
    %======================= PRESSURE CORRECTION ==========================
    %Solve linear equations [A][X]=[B]
    DFLP=-A\RESTM;%[kg/m/sec^2 == Pa]
    DFLPT=DFLPT+DFLP;
    %======================= PRESSURE CORRECTION ==========================
    
    %====================== MASS FLOW CORRECTION ==========================
    %----------------------------- DRFVX ----------------------------------
    for i=2:NIX
        %P(1,1:NIX+2)=0.0(always) --> DRFVX(1,1:NIX+1)=0.0
        %VX(NIY+2,1:NIX+1)=0.0(always) --> DRFVX(NIY+2,1:NIX+1)?=?0.0
        for j=2:NIY+1
            %VX(1:NIY+2,1)=0.0(always) --> DRFVX(1:NIY+2,1)=0.0
            %VX(1:NIY+2,NIX+1)=0.0(always) --> DRFVX(1:NIY+2,NIX+1)=0.0
            k=(i-2)*NIY+j-1;
            DRFVX(j,i)=-2.0*dtb*(DFLP(k+NIY)-DFLP(k))/((dx(i)+dx(i-1))*(1.0+dtb*UFRKX(j,i)));%[Pa.sec/m == kg/m^2/sec]
        end
    end
    %DFLP=[NIY*NIX,1]
    
    %----------------------------- DRFVX ----------------------------------
    
    %----------------------------- DRFVY ----------------------------------
    for i=2:NIX+1
        %DRFVY(1:NIY+1,1)=FREE, !!!!!!!!!!!!!!!!!!!!@#$%^&*()_+  USELESS  +_)(*&^%$#@!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        %VY(1:NIY+1,NIX+2)=0.0(always) --> DRFVY(1:NIY+1,NIX+2)=0.0
        for j=2:NIY
            %VY(1,1:NIX+2)=0.0(always) --> DRFVY(1,1:NIX+2)=0.0
            %VY(NIY+1,1:NIX+2)=0.0(always) --> DRFVY(NIY+1,1:NIX+2)=0.0
            k=(i-2)*NIY+j-1;
            DRFVY(j,i)=-2.0*dtb*(DFLP(k+1)-DFLP(k))/((dy(j)+dy(j-1))*(1.0+dtb*UFRKY(j,i)));%[Pa.sec/m == kg/m^2/sec]
        end
    end
       
    %----------------------------- DRFVY ----------------------------------
    %======================= MASS FLOW CORRECTION =========================
    
    %------------------------- NEW RFV ESTIMATION -------------------------
    RFVX=RFVXTM+DRFVX;
    RFVY=RFVYTM+DRFVY;
    
    %------------------------- NEW RFV ESTIMATION -------------------------
    
    %------------------------- RFVX RFVY error ----------------------------
    errDRFVX=abs(max(max(DRFVX)));
    errDRFVY=abs(max(max(DRFVY)));
    errDRFV=max(errDRFVY,errDRFVX);
    %------------------------- RFVX RFVY error ----------------------------
      
    iters=iters+1;
end

    %-------------------------- New Pressure ------------------------------    
    for i=2:NIX+1
        for j=2:NIY+1
            k=(i-2)*NIY+j-1;
            POTM(j,i)=DFLPT(k);
        end
    end
    
    for i=2:NIY+2
        POTM(i,1)=POTM(i,2);%FPDX(1:NIY+2,1)=0.0
        POTM(i,NIX+2)=POTM(i,NIX+1);%FPDX(1:NIY+2,NIX+1)=0.0
    end
    for i=1:NIX+2
        POTM(NIY+2,i)=POTM(NIY+1,i)+0.5*g*dy(NIY)*RLNTM(NIY+2,i);
    end
    
    P=POTM;
    
    %-------------------------- New Pressure ------------------------------    

%------------------------------ NEW VX ------------------------------------
for i=2:NIX
    %VX(1,1:NIX+1)=0.0 always
    %VX(NIY+2,1:NIX+1)=0.0 always
    for j=2:NIY+1
        %VX(1:NIY+2,1)=0.0 always
        %VX(1:NIY+2,NIX+1)=0.0 always
        VX(j,i)=RFVX(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j,i+1))*0.5*(FLNTM(j,i)+FLNTM(j,i+1)));%[m/sec]
    end
end

%------------------------------ NEW VX ------------------------------------

%------------------------------ NEW VY ------------------------------------
for i=1:NIX+1
    %VY(1,1:NIX+2)=0.0 always
    %VY(NIY+1,1:NIX+2)=0.0 always
    for j=2:NIY
        %VY(2:NIY,1)=Free [Slip boundary condition]
        %VY(1:NIY+1,NIX+2)=0.0 always
        VY(j,i)=RFVY(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j+1,i))*0.5*(FLNTM(j,i)+FLNTM(j+1,i)));%[m/sec]
    end
end
%------------------------------ NEW VY ------------------------------------

end

