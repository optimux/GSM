function [VX,VY,P,RESM] = VPGM(VXOTM,VYOTM,PRTM,MUOTM,MUNTM,FLOTM,FLNTM,RLOTM,RLNTM,RSOTM,RSNTM,dFSTM)
%To solve V-P momentum conservation
%Created 2019-10-20

%Modified for iteration logical 2019-11-10

%Modified for boundary condition settings 2019-11-13

%Modified for old terms corrected from new terms 2020-2-14

%Modified for mass check 2020-2-16

%O==Old
%N==New
%G-->Gravity full consideration

global NIX
global NIY
global dx
global dy
global g
global RL0
global dtb

%======================= VERY IMPORTANT NOTE ==============================
%It is highly recommended the usage of packages for complex formulae, i.e.,
%temporary variables. These packages can be utilized repetedly, easily
%modified, and produce a concise numerical equation.
%======================= VERY IMPORTANT NOTE ==============================


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
        UFRKY(j,i)=0.5*(MUNTM(j,i)+MUNTM(j+1,i))*0.5*(FLNTM(j,i)+FLNTM(j+1,i))/(0.5*(RLNTM(j,i)+RLNTM(j+1,i))*0.5*(KFLTM(j,i)+KFLTM(j+1,i)));%[1/sec]
    end
end

UFRKX=zeros(NIY+2,NIX+1);%X==x-axis
for i=1:NIX+1
    for j=1:NIY+2
        UFRKX(j,i)=0.5*(MUNTM(j,i)+MUNTM(j,i+1))*0.5*(FLNTM(j,i)+FLNTM(j,i+1))/(0.5*(RLNTM(j,i)+RLNTM(j,i+1))*0.5*(KFLTM(j,i)+KFLTM(j,i+1)));%[1/sec]
    end
end

%------------------- mu*FL/(RL*K) in denomenator --------------------------

AX=zeros(NIY,NIX+1);
for i=2:NIX
    for j=1:NIY
        AX(j,i)=2.0*dy(j)/((dx(i-1)+dx(i))*(1.0+dtb*UFRKX(j+1,i)));%main domain, [1]
    end
end

for i=1:NIY
    AX(i,1)=2.0*dy(i)/(2.0*dx(1)*(1.0+dtb*UFRKX(i+1,1)));%1st column, [1]
    AX(i,NIX+1)=2.0*dy(i)/(2.0*dx(NIX)*(1.0+dtb*UFRKX(i+1,NIX+1)));%last column, [1]
end

%y-axis parameter in matrix
AY=zeros(NIY+1,NIX);
for i=1:NIX
    for j=2:NIY
        AY(j,i)=2.0*dx(i)/((dy(j)+dy(j-1))*(1.0+dtb*UFRKY(j,i+1)));%main domain, [1]
    end
end

for i=1:NIX
    AY(1,i)=2.0*dx(i)/(2.0*dy(1)*(1.0+dtb*UFRKY(1,i+1)));%1st row, [1]
    AY(NIY+1,i)=2.0*dx(i)/(2.0*dy(NIY)*(1.0+dtb*UFRKY(NIY+1,i+1)));%last row, [1]
end

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
A(1,2)=-AY(2,1)*0.0;%[1]
A(1,NIY+1)=-AX(1,2)*0.0;%[1]
A(1,1)=1.0;

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


%% ...............=========== THE BOMB CYCLE ============..................

%NOTE: DFLP at 4 physical walls or boundaries are useless since coefficients are set zero there, thus DFLP has onlt NIX*NIY elements.
DFLP=zeros(NIY*NIX,1);%consider delta(FL*P) as one variable which has (NIX*NIY) elements in one column
B=zeros(NIY*NIX,1);%set as [B] which has NIX*NIY elements in one column

DRFVX=zeros(NIY+2,NIX+1);
DRFVY=zeros(NIY+1,NIX+2);
VX=zeros(NIY+2,NIX+1);
VY=zeros(NIY+1,NIX+2);

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



%%................======== BASIC VELOCITY FORMULAE ========...............

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
for i=2:NIY+1
    RFVXI(i,1)=-RFVXI(i,2);%left VX=0.0 boundary
    RFVXI(i,NIX+2)=-RFVXI(i,NIX+1);%rihgt VX=0.0 boundary
end

for i=2:NIX+1
    RFVXI(1,i)=-RFVXI(2,i);%top VX=0.0 boundary
    RFVXI(NIY+2,i)=-RFVXI(NIY+1,i);%bottom VX=0.0 boundary
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
for i=1:NIX+1
    for j=2:NIY+1
        %VY(1,1:NIX+2)=0.0 && VY(image_3,1:NIX+2)=0.0 --> RFVYI(1,1:NIX+2)=0.0 [image_3==0]
        %VY(NIY+1,1:NIX+2)=0.0 && VY(image_4,1:NIX+2)=0.0 --> RFVYI(NIY+2,1:NIX+2)=0.0 [image_4==NIY+2]
        %VY(1:NIY+2,1)=Free --> RFVYI(1:NIY+2,1)!=0.0 [Slip boundary]
        %VY(1:NIY+2,NIX+2)=0.0 --> RFVYI(1:NIY+2,NIX+2)=0.0 [No slip boundary]
        RFVYI(j,i)=RLOTM(j,i)*FLOTM(j,i)*0.5*(VYOTM(j-1,i)+VYOTM(j,i));%[kg/m^2/sec]
    end
end

%VY(1:NIY+2,1)=Free --> RFVYI(1:NIY+1,1)!=0.0 [Slip boundary]
%This sentence can be added to last FOR loop
for i=2:NIX+1
    RFVYI(1,i)=-RFVYI(2,i);%top VY=0.0 boundary
    RFVYI(NIY+2,i)=-RFVYI(NIY+1,i);%bottom VY=0.0 boundary
end

for i=2:NIY+1
    RFVYI(i,NIX+2)=-RFVYI(i,NIX+1);%right VY=0.0 boundary
    %left VY=FREE
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
VVXI=zeros(NIY+2,NIX+2);%I==integer grid(main grid)
for i=2:NIX+1
    for j=1:NIY+2
        VVXI(j,i)=VXOTM(j,i-1)*max(RFVXI(j,i),0.0)+VXOTM(j,i)*min(RFVXI(j,i),0.0);%including VVX1, VVX2 [kg/m/sec^2]
    end
end

for i=2:NIY+1
    VVXI(i,1)=-VXOTM(i,2)*max(-RLOTM(i,1)*FLOTM(i,1)*0.5*VXOTM(i,2),0.0)+0.0;%left boundary
    VVXI(i,NIX+2)=0.0-VXOTM(i,NIX)*min(-RLOTM(i,NIX+2)*FLOTM(i,NIX+2)*0.5*VXOTM(i,NIX),0.0);%right boundary
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
VVYI=zeros(NIY+2,NIX+2);%I==integer grid(main grid)
for i=1:NIX+2
    for j=2:NIY+1
        VVYI(j,i)=VYOTM(j-1,i)*max(RFVYI(j,i),0.0)+VYOTM(j,i)*min(RFVYI(j,i),0.0);%[kg/m/sec^2]
    end
end

for i=2:NIX+1
    VVYI(1,i)=-VYOTM(2,i)*max(-RLOTM(1,i)*FLOTM(1,i)*0.5*VYOTM(2,i),0.0)+0.0;
    VVYI(NIY+2,i)=0.0-VYOTM(NIY,i)*min(-RLOTM(NIY+2,i)*FLOTM(NIY+2,i)*0.5*VYOTM(NIY,i),0.0);
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
VDXI=zeros(NIY+2,NIX+2);
for i=2:NIX+1
    for j=1:NIY+2
        VDXI(j,i)=MUOTM(j,i)*(VXOTM(j,i)*0.5*(FLOTM(j,i)+FLOTM(j,i+1))-VXOTM(j,i-1)*0.5*(FLOTM(j,i-1)+FLOTM(j,i)))/dx(i-1);%[Pa.sec times 1/sec == Pa]
    end
end

for i=2:NIY+1
    VDXI(i,1)=MUOTM(i,1)*(0.0-(-VXOTM(i,2)*0.5*(FLOTM(i,3)+FLOTM(i,2))))/dx(1);%left boundary
    VDXI(i,NIY+2)=MUOTM(i,NIX+2)*(0.5*(FLOTM(i,NIX)+FLOTM(i,NIX+1))*(-VXOTM(i,NIX))-0.0)/dx(NIX);%right boundary
end

%VDXI(1:NIY+2,image_1)=0.0 at image_1 and image_2
%VDXI(1:NIY+2,image_2)=0.0 at image_1 and image_2

VDXH=zeros(NIY+1,NIX+1);
for i=1:NIX+1
    for j=2:NIY
        VDXH(j,i)=2.0*0.25*(MUOTM(j,i)+MUOTM(j,i+1)+MUOTM(j+1,i)+MUOTM(j+1,i+1))*(0.5*(FLOTM(j+1,i)+FLOTM(j+1,i+1))*VXOTM(j+1,i)-0.5*(FLOTM(j,i)+FLOTM(j,i+1))*VXOTM(j,i))/(dy(j-1)+dy(j));%[Pa.sec times 1/sec == Pa]
    end
end

for i=1:NIX+1
    VDXH(1,i)=2.0*0.25*(MUOTM(1,i)+MUOTM(1,i+1)+MUOTM(2,i)+MUOTM(2,i+1))*(0.5*(FLOTM(2,i)+FLOTM(2,i+1))*VXOTM(2,i)-0.5*(FLOTM(1,i)+FLOTM(1,i+1))*VXOTM(1,i))/(2.0*dy(1));%VXOTM(1,i)=-VXOTM(2,i); [Pa.sec times 1/sec == Pa]
    VDXH(NIY+1,i)=2.0*0.25*(MUOTM(NIY+1,i)+MUOTM(NIY+1,i+1)+MUOTM(NIY+2,i)+MUOTM(NIY+2,i+1))*(0.5*(FLOTM(NIY+2,i)+FLOTM(NIY+2,i+1))*VXOTM(NIY+2,i)-0.5*(FLOTM(NIY+1,i)+FLOTM(NIY+1,i+1))*VXOTM(NIY+1,i))/(2.0*dy(NIY));%VXOTM(NIY+1,i)=-VXOTM(NIY+2,i);[Pa.sec times 1/sec == Pa]
end

%-------------------- 5. x-viscous momentum -------------------------------

%-------------------- 6. y-viscous momentum -------------------------------
VDYI=zeros(NIY+2,NIX+2);
for i=1:NIX+2
    for j=2:NIY+1
        VDYI(j,i)=MUOTM(j,i)*(VYOTM(j,i)*0.5*(FLOTM(j+1,i)+FLOTM(j,i))-VYOTM(j-1,i)*0.5*(FLOTM(j,i)+FLOTM(j-1,i)))/dy(j-1);%[Pa.sec times 1/sec == Pa]
    end
end

for i=2:NIX+1
    VDYI(1,i)=MUOTM(1,i)*(0.0-0.5*(FLOTM(2,i)+FLOTM(3,i))*(-VYOTM(2,i)))/dy(1);%top boundary
    VDYI(NIY+2,i)=MUOTM(NIY+2,i)*(0.5*(FLOTM(NIY+1,i)+FLOTM(NIY,i))*(-VYOTM(NIY,i))-0.0)/dy(NIY);%bottom boundary
end

VDYH=zeros(NIY+1,NIX+1);
for i=2:NIX
    for j=1:NIY+1
        VDYH(j,i)=2.0*0.25*(MUOTM(j,i)+MUOTM(j,i+1)+MUOTM(j+1,i)+MUOTM(j+1,i+1))*(0.5*(FLOTM(j+1,i+1)+FLOTM(j,i+1))*VYOTM(j,i+1)-0.5*(FLOTM(j+1,i)+FLOTM(j,i))*VYOTM(j,i))/(dx(i-1)+dx(i));%[Pa.sec times 1/sec == Pa]
    end
end

for i=1:NIY+1
    VDYH(i,1)=2.0*0.25*(MUOTM(i,1)+MUOTM(i,2)+MUOTM(i+1,1)+MUOTM(i+1,2))*(0.5*(FLOTM(i+1,2)+FLOTM(i,2))*VYOTM(i,2)-0.5*(FLOTM(i+1,1)+FLOTM(i,1))*VYOTM(i,1))/(2.0*dx(1));%VYOTM(i,1)=VYOTM(i,2) for slip boundary; [Pa.sec times 1/sec == Pa]
    VDYH(i,NIX+1)=2.0*0.25*(MUOTM(i,NIX+1)+MUOTM(i,NIX+2)+MUOTM(i+1,NIX+1)+MUOTM(i+1,NIX+2))*(0.5*(FLOTM(i+1,NIX+2)+FLOTM(i,NIX+2))*VYOTM(i,NIX+2)-0.5*(FLOTM(i+1,NIX+1)+FLOTM(i,NIX+1))*VYOTM(i,NIX+1))/(2.0*dx(NIX));%VYOTM(i,NIX+2)=-VYOTM(i,NIX+1); [Pa.sec times 1/sec == Pa]
end

%-------------------- 6. y-viscous momentum -------------------------------

%--------------------- 7. x-convection part -------------------------------
VXCONV=zeros(NIY,NIX+1);
for i=2:NIX
    for j=1:NIY
        VXCONV(j,i)=-2.0*dtb*(VVXI(j+1,i+1)-VVXI(j+1,i))/(dx(i)+dx(i-1))-dtb*(VVXH(j+1,i)-VVXH(j,i))/dy(j);%[kg/m^2/sec]
    end
end

for i=1:NIY
    VXCONV(i,1)=-2.0*dtb*(VVXI(i+1,2)-VVXI(i+1,1))/(2.0*dx(1))-dtb*(VVXH(i+1,1)-VVXH(i,1))/dy(i);%left boundary [kg/m^2/sec]
    VXCONV(i,NIX+1)=-2.0*dtb*(VVXI(i+1,NIX+2)-VVXI(i+1,NIX+1))/(2.0*dx(NIX))-dtb*(VVXH(i+1,NIX+1)-VVXH(i,NIX+1))/dy(i);%right boundary [kg/m^2/sec]
end
%--------------------- 7. x-convection part -------------------------------

%---------------------- 8. x-diffusion part -------------------------------
VXDIF=zeros(NIY,NIX+1);
for i=2:NIX
    for j=1:NIY
        VXDIF(j,i)=2.0*dtb*(VDXI(j+1,i+1)-VDXI(j+1,i))/(dx(i)+dx(i-1))+dtb*(VDXH(j+1,i)-VDXH(j,i))/dy(j);%[Pa.sec/m]
    end
end

for i=1:NIY
    VXDIF(i,1)=2.0*dtb*(VDXI(i+1,2)-VDXI(i+1,1))/(2.0*dx(1))+dtb*(VDXH(i+1,1)-VDXH(i,1))/dy(i);%left boundary [Pa.sec/m]
    VXDIF(i,NIX+1)=2.0*dtb*(VDXI(i+1,NIX+2)-VDXI(i+1,NIX+1))/(2.0*dx(NIX))+dtb*(VDXH(i+1,NIX+1)-VDXH(i,NIX+1))/dy(i);%right boundary [Pa.sec/m]
end
%---------------------- 8. x-diffusion part -------------------------------

%---------------------- 9. y-convection part ------------------------------
VYCONV=zeros(NIY+1,NIX);
for i=1:NIX
    for j=2:NIY
        VYCONV(j,i)=-dtb*(VVYH(j,i+1)-VVYH(j,i))/dx(i)-2.0*dtb*(VVYI(j+1,i+1)-VVYI(j,i+1))/(dy(j)+dy(j-1));%[Pa.sec/m]
    end
end

for i=1:NIX
    VYCONV(1,i)=-dtb*(VVYH(1,i+1)-VVYH(1,i))/dx(i)-2.0*dtb*(VVYI(2,i+1)-VVYI(1,i+1))/(2.0*dy(1));%VVYI(imag_3,1:NIX+2)=0.0 have no image_3 sites [Pa.sec/m]
    VYCONV(NIY+1,i)=-dtb*(VVYH(NIY+1,i+1)-VVYH(NIY+1,i))/dx(i)-2.0*dtb*(VVYI(NIY+2,i+1)-VVYI(NIY+1,i+1))/(2.0*dy(NIY));%VVYI(imag_4,1:NIX+2)=0.0 have no image_4 sites [Pa.sec/m]
end

%---------------------- 9. y-convection part ------------------------------

%---------------------- 10. y-diffusion part ------------------------------
VYDIF=zeros(NIY+1,NIX);
for i=1:NIX
    for j=2:NIY
        VYDIF(j,i)=dtb*(VDYH(j,i+1)-VDYH(j,i))/dx(i)+2.0*dtb*(VDYI(j+1,i+1)-VDYI(j,i+1))/(dy(j)+dy(j-1));%[Pa.sec/m]
    end
end

for i=1:NIX
    VYDIF(1,i)=dtb*(VDYH(1,i+1)-VDYH(1,i))/dx(i)+2.0*dtb*(VDYI(2,i+1)-VDYI(1,i+1))/(2.0*dy(1));%top bounadry [Pa.sec/m]
    VYDIF(NIY+1,i)=dtb*(VDYH(NIY+1,i+1)-VDYH(NIY+1,i))/dx(i)+2.0*dtb*(VDYI(NIY+2,i+1)-VDYI(NIY+1,i+1))/(2.0*dy(NIY));%bottom boundary [Pa.sec/m]
end
%---------------------- 10. y-diffusion part ------------------------------

%-------------------- 11. x- pressure gradient ----------------------------
FPDX=zeros(NIY+2,NIX+1);%D==difference
for i=2:NIX
    for j=1:NIY+2
        %FPDX(j,i)=-2.0*dtb*(FLOTM(j,i+1)*vpa(POTM(j,i+1),13)-FLOTM(j,i)*vpa(POTM(j,i),13))/(dx(i)+dx(i-1));%[Pa.sec/m == kg/m^2/sec]
        FPDX(j,i)=-2.0*dtb*(FLOTM(j,i+1)*PRTM(j,i+1)-FLOTM(j,i)*PRTM(j,i))/(dx(i)+dx(i-1));%[Pa.sec/m == kg/m^2/sec]
    end
end
%NOTE: pressure with 13 digits precision; otherwise numerical precision
%will give non-zero results at t=0.0 sec

for i=1:NIY+2
    FPDX(i,1)=-2.0*dtb*(FLOTM(i,2)*PRTM(i,2)-FLOTM(i,1)*PRTM(i,1))/(2.0*dx(1));%[Pa.sec/m == kg/m^2/sec]
    FPDX(i,NIX+1)=-2.0*dtb*(FLOTM(i,NIX+2)*PRTM(i,NIX+2)-FLOTM(i,NIX+1)*PRTM(i,NIX+1))/(2.0*dx(NIX));%[Pa.sec/m == kg/m^2/sec]
end

PSFX=zeros(NIY+2,NIX+1);%static pressure * FL, to extract static pressure from actual pressure
for i=2:NIX
    for j=1:NIY+2
        PSFX(j,i)=-2.0*dtb*(FLNTM(j,i+1)*P0(j,i+1)-FLNTM(j,i)*P0(j,i))/(dx(i)+dx(i-1));%[Pa.sec/m == kg/m^2/sec]
    end
end

for i=1:NIY+2
    PSFX(i,1)=-2.0*dtb*(FLNTM(i,2)*P0(i,2)-FLNTM(i,1)*P0(i,1))/(2.0*dx(1));%left boundary [Pa.sec/m == kg/m^2/sec]
    PSFX(i,NIX+1)=-2.0*dtb*(FLNTM(i,NIX+2)*P0(i,NIX+2)-FLNTM(i,NIX+1)*P0(i,NIX+1))/(2.0*dx(NIX));%right boundary [Pa.sec/m == kg/m^2/sec]
end

%-------------------- 11. x- pressure gradient ----------------------------

%---------------------- 12. y- pressure part ------------------------------

FPDY=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=2:NIY
        FPDY(j,i)=-2.0*dtb*(FLOTM(j+1,i)*PRTM(j+1,i)-FLOTM(j,i)*PRTM(j,i))/(dy(j)+dy(j-1));%[Pa.sec/m == kg/m^2/sec]
    end
end

for i=1:NIX+2
    FPDY(1,i)=-2.0*dtb*(FLOTM(2,i)*PRTM(2,i)-FLOTM(1,i)*PRTM(1,i))/(2.0*dy(1));%[Pa.sec/m == kg/m^2/sec]
    FPDY(NIY+1,i)=-2.0*dtb*(FLOTM(NIY+2,i)*PRTM(NIY+2,i)-FLOTM(NIY+1,i)*PRTM(NIY+1,i))/(2.0*dy(NIY));%[Pa.sec/m == kg/m^2/sec]
end

PSFY=zeros(NIY+1,NIX+2);%static pressure * FL, to extract static pressure from actual pressure
for i=1:NIX+2
    for j=2:NIY
        PSFY(j,i)=-2.0*dtb*(FLNTM(j+1,i)*P0(j+1,i)-FLNTM(j,i)*P0(j,i))/(dy(j)+dy(j-1));%[Pa.sec/m == kg/m^2/sec]
    end
end

for i=1:NIX+2
    PSFY(1,i)=-2.0*dtb*(FLNTM(2,i)*P0(2,i)-FLNTM(1,i)*P0(1,i))/(2.0*dy(1));%upper boundary [Pa.sec/m == kg/m^2/sec]
    PSFY(NIY+1,i)=-2.0*dtb*(FLNTM(NIY+2,i)*P0(NIY+2,i)-FLNTM(NIY+1,i)*P0(NIY+1,i))/(2.0*dy(NIY));%bottom boundary [Pa.sec/m == kg/m^2/sec]
end

%---------------------- 12. y- pressure part ------------------------------

%----------------------- 13. y- gravity part ------------------------------

FRLGOTM=zeros(NIY+1,NIX+2);%G==gravity related terms in RFVYTM
FRLGNTM=zeros(NIY+1,NIX+2);

% %************** Scheme 1: use approximation, see PPT P.45 *****************
for i=1:NIX+2
    for j=1:NIY+1
        FRLGOTM(j,i)=0.5*(FLOTM(j,i)+FLOTM(j+1,i))*(0.5*(RLOTM(j,i)+RLOTM(j+1,i))-0.5*(RL0(j,i)+RL0(j+1,i)));
        FRLGNTM(j,i)=0.5*(FLNTM(j,i)+FLNTM(j+1,i))*(0.5*(RLNTM(j,i)+RLNTM(j+1,i))-0.5*(RL0(j,i)+RL0(j+1,i)));
        
        % FRLGNTM(j,i)=0.5*(FLNTM(j,i)*RLNTM(j,i)+FLNTM(j+1,i)*RLNTM(j+1,i))-0.5*(FLNTM(j,i)*RL0(j,i)+FLNTM(j+1,i)*RL0(j+1,i));
        % FRLGOTM(j,i)=0.5*(FLOTM(j,i)*RLOTM(j,i)+FLOTM(j+1,i)*RLOTM(j+1,i))-0.5*(FLOTM(j,i)*RL0(j,i)+FLOTM(j+1,i)*RL0(j+1,i));
        
        %FRLNTM(j,i)=FRLOTM(j,i);%For constant density. If this is set,
        %then sum(B)=0.0 like lid driven cavity flow
    end
end

FRG=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=1:NIY+1
        FRG(j,i)=0.5*g*dtb*(FRLGOTM(j,i)+FRLGNTM(j,i));%[kg/m^2/sec]
    end
end
% %**************************************************************************


% %******** Scheme 2: Do not use approximation, see PPT P.43 & P.64 *********
% for i=1:NIX+2
%     for j=1:NIY+1
%         FRLGOTM(j,i)=0.5*(FLOTM(j,i)+FLOTM(j+1,i))*0.5*(RLOTM(j,i)+RLOTM(j+1,i));
%         FRLGNTM(j,i)=0.5*(FLNTM(j,i)+FLNTM(j+1,i))*0.5*(RLNTM(j,i)+RLNTM(j+1,i));
%         
%         %FRLGNTM(j,i)=0.5*(FLNTM(j,i)*RLNTM(j,i)+FLNTM(j+1,i)*RLNTM(j+1,i));
%         %FRLGOTM(j,i)=0.5*(FLOTM(j,i)*RLOTM(j,i)+FLOTM(j+1,i)*RLOTM(j+1,i));
%         
%     end
% end
% 
% %FRLNTM(j,i)=FRLOTM(j,i);%For constant density. If this is set,
% %then sum(B)=0.0 like lid driven cavity flow
% 
% 
% FRG=zeros(NIY+1,NIX+2);
% for i=1:NIX+2
%     for j=1:NIY+1
%         FRG(j,i)=g*dtb*(0.5*FRLGOTM(j,i)+0.5*FRLGNTM(j,i))+PSFY(j,i);%[kg/m^2/sec]
%     end
% end
% %************************************************************************

%NOTE: In Boussinesq assumption, only change of density in gravity term
%has effect; in present program, we include all density change effects!
%----------------------- 13. y- gravity part ------------------------------

%===================== APPROXIMATION OF RFVXY =========================
%--------------------------- RFVXTM -----------------------------------
RFVXTM=zeros(NIY+2,NIX+1);
for i=1:NIX+1
    for j=2:NIY+1
        RFVXTM(j,i)=(RFVX(j,i)+VXCONV(j-1,i)+VXDIF(j-1,i)+FPDX(j,i)+PSFX(j,i))/(1.0+dtb*UFRKX(j,i));%[kg/m^2/sec]
    end
end
RFVXTM(:,1)=0.0;% are available but have no effect on VX
RFVXTM(:,NIX+1)=0.0;% are available but have no effect on VX

%--------------------------- RFVXTM -----------------------------------

%--------------------------- RFVYTM -----------------------------------
RFVYTM=zeros(NIY+1,NIX+2);

for i=2:NIX+1
    for j=2:NIY
        RFVYTM(j,i)=(RFVY(j,i)+VYCONV(j,i-1)+VYDIF(j,i-1)+FPDY(j,i)+FRG(j,i))/(1.0+dtb*UFRKY(j,i));%[kg/m^2/sec]
    end
end

for i=2:NIX+1
    RFVYTM(1,i)=(RFVY(1,i)+VYCONV(1,i-1)+VYDIF(1,i-1)+FPDY(1,i)+FRG(1,i))/(1.0+dtb*UFRKY(1,i));%1st row, [kg/m^2/sec]
    RFVYTM(NIY+1,i)=(RFVY(NIY+1,i)+VYCONV(NIY+1,i-1)+VYDIF(NIY+1,i-1)+FPDY(NIY+1,i)+FRG(NIY+1,i))/(1.0+dtb*UFRKY(NIY+1,i));%last row, [kg/m^2/sec]
end
RFVYTM(1,:)=0.0;
RFVYTM(NIY+1,:)=0.0;
%     for i=2:NIY
%         %1-st column is in slip boundary condition
%         RFVYTM(i,1)=(RFVY(i,1)-0.0-2.0*dtb*(VVYI(i,1)-VVYI(i-1,1))/(dy(i)+dy(i-1))+0.0+2.0*dtb*(VDYI(i,1)-VDYI(i-1,1))/(dy(i)+dy(i-1))+FPDY(i,1)+FRG(i,1))/(1.0+dtb*UFRKY(i,1));%[kg/m^2/sec]
%         %NOTE: VVYH(i,image_1)=VVYH(i,1)=0.0 since RFVX(j+0.5,k+0.5)==0.0 --> first '0.0' in above formula;
%         %VDYH(i,1)=VDYH(i,image_1)=0.0 since VY(i,1)==VY(i,image_1)==VY(i,image_image) --> second '0.0' in above formula
%
%         %last column is in no slip boundary condition
%         RFVYTM(i,NIX+2)=(RFVY(i,NIX+2)-0.0-2.0*dtb*(VVYI(i,NIX+2)-VVYI(i-1,NIX+2))/(dy(i)+dy(i-1))+0.0+2.0*dtb*(VDYI(i,NIX+2)-VDYI(i-1,NIX+2))/(dy(i)+dy(i-1))+FPDY(i,NIX+2)+FRG(i,NIX+2))/(1.0+dtb*UFRKY(i,NIX+2));%no x -axis mass flow
%     end

%--------------------------- RFVYTM -----------------------------------
%===================== APPROXIMATION OF RFVXY =========================

%================== MASS CONSERVATION RESIDUAL  =======================
%----------------------------- RES ------------------------------------
RES=zeros(NIY,NIX);
for i=1:NIX
    for j=1:NIY
        RES(j,i)=dy(j)*(RFVXTM(j+1,i+1)-RFVXTM(j+1,i))+dx(i)*(RFVYTM(j+1,i+1)-RFVYTM(j,i+1))+0.5*dFSTM(j+1,i+1)*(RSOTM(j+1,i+1)+RSNTM(j+1,i+1))*dx(i)*dy(j)/dtb+(FLNTM(j+1,i+1)*RLNTM(j+1,i+1)-FLOTM(j+1,i+1)*RLOTM(j+1,i+1))*dx(i)*dy(j)/dtb;%[kg/m/sec == Pa.sec]
        %NOTE: This RES is absolutely universal to all conditions,
        %especially for those related to significant density change which
        %can be seen from RSOTM, RSNTM
        
    end
end

%     for i=2:NIX-1
%         RES(1,i)=dy(1)*(RFVXTM(2,i+1)-RFVXTM(2,i))+dx(i)*(RFVYTM(2,i+1)-0.0)+0.5*dFSTM(2,i+1)*(RSOTM(2,i+1)+RSNTM(2,i+1))*dx(i)*dy(1)/dtb+(FRLNTM(2,i+1)-FRLOTM(2,i+1))*dx(i)*dy(1)/dtb;%top boundary
%         RES(NIY,i)=dy(NIY)*(RFVXTM(NIY+1,i+1)-RFVXTM(NIY+1,i))+dx(i)*(0.0-RFVYTM(NIY,i+1))+0.5*dFSTM(NIY+1,i+1)*(RSOTM(NIY+1,i+1)+RSNTM(NIY+1,i+1))*dx(i)*dy(NIY)/dtb+(FRLNTM(NIY+1,i+1)-FRLOTM(NIY+1,i+1))*dx(i)*dy(NIY)/dtb;%bottom boundary
%     end
%
%     for i=2:NIX-1
%         RES(i,1)=dy(i)*(RFVXTM(i+1,2)-0.0)+dx(1)*(RFVYTM(i+1,2)-RFVYTM(i,2))+0.5*dFSTM(i+1,2)*(RSOTM(i+1,2)+RSNTM(i+1,2))*dx(1)*dy(i)/dtb+(FRLNTM(i+1,2)-FRLOTM(i+1,2))*dx(1)*dy(i)/dtb;%left boundary
%         RES(i,NIX)=dy(i)*(0.0-RFVXTM(i+1,NIX))+dx(NIX)*(RFVYTM(i+1,NIX+1)-RFVYTM(i,NIX+1))+0.5*dFSTM(i+1,NIX+1)*(RSOTM(i+1,NIX+1)+RSNTM(i+1,NIX+1))*dx(NIX)*dy(i)/dtb+(FRLNTM(i+1,NIX+1)-FRLOTM(i+1,NIX+1))*dx(NIX)*dy(i)/dtb;%right boundary
%     end
%
%     RES(1,1)=dy(1)*(RFVXTM(2,2)-0.0)+dx(1)*(RFVYTM(2,2)-0.0)+0.5*dFSTM(2,2)*(RSOTM(2,2)+RSNTM(2,2))*dx(1)*dy(1)/dtb+(FRLNTM(2,2)-FRLOTM(2,2))*dx(1)*dy(1)/dtb;
%     RES(1,NIX)=dy(1)*(0.0-RFVXTM(2,NIX))+dx(NIX)*(RFVYTM(2,NIX+1)-0.0)+0.5*dFSTM(2,NIX+1)*(RSOTM(2,NIX+1)+RSNTM(2,NIX+1))*dx(NIX)*dy(1)/dtb+(FRLNTM(2,NIX+1)-FRLOTM(2,NIX+1))*dx(NIX)*dy(1)/dtb;
%     RES(NIY,1)=dy(NIY)*(RFVXTM(NIY+1,2)-0.0)+dx(1)*(0.0-RFVYTM(NIY,2))+0.5*dFSTM(NIY+1,2)*(RSOTM(NIY+1,2)+RSNTM(NIY+1,2))*dx(1)*dy(NIY)/dtb+(FRLNTM(NIY+1,2)-FRLOTM(NIY+1,2))*dx(1)*dy(NIY)/dtb;
%     RES(NIY,NIX)=dy(NIY)*(0.0-RFVXTM(NIY+1,NIX))+dx(NIX)*(0.0-RFVYTM(NIY,NIX+1))+0.5*dFSTM(NIY+1,NIX+1)*(RSOTM(NIY+1,NIX+1)+RSNTM(NIY+1,NIX+1))*dx(NIX)*dy(NIY)/dtb+(FRLNTM(NIY+1,NIX+1)-FRLOTM(NIY+1,NIX+1))*dx(NIX)*dy(NIY)/dtb;
%
for i=1:NIX
    B(1+(i-1)*NIY:i*NIY)=-RES(1:NIY,i)/dtb;%[kg/m/sec == Pa.sec]
end
B(1)=0.0;

%     RHS=sum(B);
%     [u,s]=eig(A);%s is eigenvalue which is used to determine wether matrix A is positive-definite or not
%     kk=rank(A);
%in this case, s has both positive and negative values, though negative
%value is rather small which is can be taken as 0. Thus, matrix A here
%is positive-semidefinite.

%----------------------------- RES ------------------------------------
%=================== MASS CONSERVATION RESIDUAL  ======================

%======================= PRESSURE CORRECTION ==========================
%Solve linear equations [A][X]=[B]
%DFLP=A\B;%[kg/m/sec^2 == Pa]
%DFLP=GaussSeidel(A,B);

AS=sparse(A);%convert sparse A into sparse storage format
DFLP=PardisoSolver(AS,B);

%     DFLP2=pentaDiag_solve(A,B);
%     DFLP=DFLP2-DFLP1;
%     [DFLP,errr,iterr,flag]=sor(A,B,B,0.5,2000,1e-4);
%     DFLP=DFLP1;
%     DFLP=conjgrad(A,B);
%DFLP=ConjuGrad(A,B,POTM-P0);
% [DFLP,flags]=PCG(A,B);
%---------------------- New  Relative Pressure ------------------------
for i=2:NIX+1
    for j=2:NIY+1
        k=(i-2)*NIY+j-1;
        PRTM(j,i)=(PRTM(j,i)*FLOTM(j,i)+DFLP(k))/FLNTM(j,i);
    end
end

for i=2:NIY+1
    %left boundary: VX(1:NIY+2,1)=0, RFVXTM(1:NIY+2,1)=0,so DFLP(1:NIY)=DFLP_left_ghost_cell
    PRTM(i,1)=(PRTM(i,1)*FLOTM(i,1)+DFLP(i-1))/FLNTM(i,1);
    
    %right boundary: VX(1:NIY+2,NIX+1)=0, RFVXTM(1:NIY+2,NIX+1)=0,so DFLP((NIX-1)*NIY+1:NIX*NIY)=DFLP_righ_ghost_cell
    PRTM(i,NIX+2)=(PRTM(i,NIX+2)*FLOTM(i,NIX+2)+DFLP((NIX-1)*NIY+i-1))/FLNTM(i,NIX+2);
end
for i=1:NIX
    %top boundary: VY(1,1:NIX+2)=0, RFVYTM(1,1:NIX+2)=0,so DFLP(1:NIY:(NIX-1)*NIY+1)=DFLP_upper_ghost_cell
    PRTM(1,i+1)=(PRTM(1,i+1)*FLOTM(1,i+1)+DFLP(1+(i-1)*NIY))/FLNTM(1,i+1);
    
    %bottom boundary: VY(NIY+1,1:NIX+2)=0, RFVYTM(NIY+1,1:NIX+2)=0,so DFLP(NIY:NIY:NIX*NIY)=DFLP_bottom_ghost_cell
    PRTM(NIY+2,i+1)=(PRTM(NIY+2,i+1)*FLOTM(NIY+2,i+1)+DFLP(i*NIY))/FLNTM(NIY+2,i+1);
end
PRTM(1,1)=PRTM(1,2);
PRTM(NIY+2,1)=PRTM(NIY+2,2);
PRTM(1,NIX+2)=PRTM(1,NIX+1);
PRTM(NIY+2,NIX+2)=PRTM(NIY+2,NIX+1);

%---------------------- New Relative Pressure -------------------------

%======================= PRESSURE CORRECTION ==========================

%====================== MASS FLOW CORRECTION ==========================
%-------------------------- DRFVX RFVXN--------------------------------
RFVXN=zeros(NIY+2,NIX+1);%new RFVX
for i=2:NIX
    %P(1,1:NIX+2)=0.0(always) --> DRFVX(1,1:NIX+1)=0.0
    %VX(NIY+2,1:NIX+1)=0.0(always) --> DRFVX(NIY+2,1:NIX+1)?=?0.0
    for j=2:NIY+1
        %VX(1:NIY+2,1)=0.0(always) --> DRFVX(1:NIY+2,1)=0.0
        %VX(1:NIY+2,NIX+1)=0.0(always) --> DRFVX(1:NIY+2,NIX+1)=0.0
        k=(i-2)*NIY+j-1;
        DRFVX(j,i)=-2.0*dtb*(DFLP(k+NIY)-DFLP(k))/((dx(i)+dx(i-1))*(1.0+dtb*UFRKX(j,i)));%[Pa.sec/m == kg/m^2/sec]
        RFVXN(j,i)=RFVXTM(j,i)+DRFVX(j,i);
    end
end
RFVXN(1,:)=-RFVXN(2,:);
RFVXN(NIY+2,:)=-RFVXN(NIY+1,:);
%here, RFVXN(1:NIY+2,1)=0.0 & RFVXN(1:NIY+2,NIX+1)=0.0
%----------------------------- DRFVX ----------------------------------

%-------------------------- DRFVY RFVYN -------------------------------
RFVYN=zeros(NIY+1,NIX+2);%new RFVY
for i=2:NIX+1
    %DRFVY(1:NIY+1,1)=FREE
    %VY(1:NIY+1,NIX+2)=0.0(always) --> DRFVY(1:NIY+1,NIX+2)=0.0
    for j=2:NIY
        %VY(1,1:NIX+2)=0.0(always) --> DRFVY(1,1:NIX+2)=0.0
        %VY(NIY+1,1:NIX+2)=0.0(always) --> DRFVY(NIY+1,1:NIX+2)=0.0
        k=(i-2)*NIY+j-1;
        DRFVY(j,i)=-2.0*dtb*(DFLP(k+1)-DFLP(k))/((dy(j)+dy(j-1))*(1.0+dtb*UFRKY(j,i)));%[Pa.sec/m == kg/m^2/sec]
        RFVYN(j,i)=RFVYTM(j,i)+DRFVY(j,i);
    end
end
RFVYN(:,1)=RFVYN(:,2);
RFVYN(:,NIX+2)=-RFVYN(:,NIX+1);
%here, RFVYN(1,1:NIX+2)=0.0 & RFVYN(NIX+1,1:NIX+2)=0.0
%----------------------------- DRFVY ----------------------------------
%======================= MASS FLOW CORRECTION =========================

%------------------------- NEW RFV ESTIMATION -------------------------
%     RFVX=RFVXTM+DRFVX;
%     RFVY=RFVYTM+DRFVY;

%------------------------- NEW RFV ESTIMATION -------------------------

%------------------------------ NEW VX ------------------------------------
for i=2:NIX
    %VX(1,1:NIX+1)=0.0 always
    %VX(NIY+2,1:NIX+1)=0.0 always
    for j=2:NIY+1
        %VX(1:NIY+2,1)=0.0 always
        %VX(1:NIY+2,NIX+1)=0.0 always
        %             DVX(j,i)=DRFVX(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j,i+1))*0.5*(FLNTM(j,i)+FLNTM(j,i+1)));%[m/sec]
        %             VXOTM(j,i)=VXOTM(j,i)+DVX(j,i);
        VXOTM(j,i)=RFVXN(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j,i+1))*0.5*(FLNTM(j,i)+FLNTM(j,i+1)));%[m/sec]
    end
end

%NOTE:
%     VXOTM(2:NIY+1,1)=0.0;
%     VXOTM(2:NIY+1,NIX+1)=0.0;

for i=1:NIX+1
    VXOTM(1,i)=-VXOTM(2,i);
    VXOTM(NIY+2,i)=-VXOTM(NIY+1,i);
end

%------------------------------ NEW VX ------------------------------------

%------------------------------ NEW VY ------------------------------------
for i=2:NIX+1
    %VY(1,1:NIX+2)=0.0 always
    %VY(NIY+1,1:NIX+2)=0.0 always
    for j=2:NIY
        %VY(2:NIY,1)=Free [Slip boundary condition]
        %VY(1:NIY+1,NIX+2)=0.0 always
        %             DVY(j,i)=DRFVY(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j+1,i))*0.5*(FLNTM(j,i)+FLNTM(j+1,i)));%[m/sec]
        %             VYOTM(j,i)=VYOTM(j,i)+DVY(j,i);
        VYOTM(j,i)=RFVYN(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j+1,i))*0.5*(FLNTM(j,i)+FLNTM(j+1,i)));
    end
end

%NOTE:
%     VYOTM(1,2:NIX+1)=0.0;
%     VYOTM(NIY+1,2:NIX+1)=0.0;

for i=1:NIY+1
    VYOTM(i,1)=VYOTM(i,2);%left slip boundary
    VYOTM(i,NIX+2)=-VYOTM(i,NIX+1);%right boundary
end

%------------------------------ NEW VY ------------------------------------


%------------------------- MASS BALANCE CHECK -----------------------------
RESM=zeros(NIY+2,NIX+2);%absolute residual of mass
RESMR=zeros(NIY+2,NIX+2);%relative residual of mass
mass_check=1.0;%mass check marker
for i=2:NIX+1
    for j=2:NIY+1
        RESM(j,i)=dtb*dy(j-1)*(RFVXN(j,i)-RFVXN(j,i-1))+dtb*dx(i-1)*(RFVYN(j,i)-RFVYN(j-1,i))+0.5*dFSTM(j,i)*(RSOTM(j,i)+RSNTM(j,i))*dx(i-1)*dy(j-1)+(FLNTM(j,i)*RLNTM(j,i)-FLOTM(j,i)*RLOTM(j,i))*dx(i-1)*dy(j-1);%absolute mass residual [kg/sec]
        %NOTE: This RESM is absolutely universal to all conditions,
        %especially for those related to significant density change which
        %can be seen from RSOTM, RSNTM
        
        RESMR(j,i)=RESM(j,i)/(dx(i-1)*dy(j-1)*RL0(j,i)/dtb);%relative mass residual with respect to initial mass flow [(kg/sec)/(kg/sec)==1]
        if(abs(RESMR(j,i))>1.0e-9)
            fprintf('Mass check: (%2d,%2d) not balanced!\n',j-1,i-1);
            mass_check=-1.0;
        end
        
    end
end

RESM(1,2:NIX+1)=RESM(2,2:NIX+1);
RESM(NIY+2,2:NIX+1)=RESM(NIY+1,2:NIX+1);
RESM(1:NIY+2,1)=RESM(1:NIY+2,2);
RESM(1:NIY+2,NIX+2)=RESM(1:NIY+2,NIX+1);

if(mass_check>0.0)
    fprintf('Mass Balanced -- Max error %E\n',max(max(abs(RESM))));
end
%------------------------- MASS BALANCE CHECK -----------------------------

%------------------------ RETURN NEW VARIABLES ----------------------------
P=PRTM;
VX=VXOTM;
VY=VYOTM;
%------------------------ RETURN NEW VARIABLES ----------------------------

end
