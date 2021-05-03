function [VX,VY,P,iters] = VPC1(VXOTM,VYOTM,POTM,FLOTM,FLNTM,RLOTM,RLNTM,RSOTM,RSNTM,dFSTM)
%To solve V-P momentum conservation
%Created 2019-10-20

%Modified for iteration logical 2019-11-10

%Modified for boundary condition settings 2019-11-13

%Carefully rechecked all 2019-11-20

%O==Old
%N==New

global NIX
global NIY
global mu
global dx
global dy
global g
global P0
global RL0
global dtb
global bsq

%======================= VERY IMPORTANT NOTE ==============================
%It is highly recommended the usage of packages for complex formulae, i.e.,
%temporary variables. These packages can be utilized repetedly, easily
%modified, and produce a concise numerical equation.
%======================= VERY IMPORTANT NOTE ==============================

% VXOTM=0.01+zeros(NIY+2,NIX+1);%x-axis velocity [m/sec]
% VXOTM(2:NIY+1,1)=0.0;
% VXOTM(2:NIY+1,NIX+1)=0.0;
% VXOTM(1,:)=-VXOTM(2,:);
% VXOTM(NIY+2,:)=-VXOTM(NIY+1,:);
%
% VYOTM=0.02+zeros(NIY+1,NIX+2);%y-axis velocity [m/sec]
% VYOTM(1,2:NIX+1)=0.0;
% VYOTM(NIY+1,2:NIX+1)=0.0;
% VYOTM(:,1)=VYOTM(:,2);
% VYOTM(:,NIX+2)=-VYOTM(:,NIX+1);

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
A(1,1)=AX(1,2)+AY(2,1)+AX(1,1)+AY(1,1);%[1]
A(1,2)=-AY(2,1);%[1]
A(1,NIY+1)=-AX(1,2);%[1]
% A(1,1)=1.0;

%special boundary condition at [NIY,1]
A(NIY,NIY)=AX(NIY,2)+AY(NIY,1)+AX(NIY,1)+AY(NIY+1,1);%[1]
A(NIY,NIY-1)=-AY(NIY,1);%[1]
A(NIY,2*NIY)=-AX(NIY,2);%[1]

%special boundary condition at [1,NIX]
A((NIX-1)*NIY+1,(NIX-1)*NIY+1)=AX(1,NIX)+AY(2,NIX)+AY(1,NIX)+AX(1,NIX+1);%[1]
A((NIX-1)*NIY+1,(NIX-1)*NIY+2)=-AY(2,NIX);%[1]
A((NIX-1)*NIY+1,(NIX-2)*NIY+1)=-AX(1,NIX);%[1]

%special boundary condition at [NIY,NIX]
A(NIX*NIY,NIX*NIY)=AX(NIY,NIX)+AY(NIY,NIX)+AX(NIY,NIX+1)+AY(NIY+1,NIX);%[1]
A(NIX*NIY,NIX*NIY-1)=-AY(NIY,NIX);%[1]
A(NIX*NIY,(NIX-1)*NIY)=-AX(NIY,NIX);%[1]

%Left boundary
for j=2:NIY-1
    A(j,j)=AY(j,1)+AY(j+1,1)+AX(j,2)+AX(j,1);%[1]
    A(j,j-1)=-AY(j,1);%[1]
    A(j,j+1)=-AY(j+1,1);%[1]
    A(j,j+NIY)=-AX(j,2);%[1]
end

%Right boundary
for i=2:NIY-1
    j=(NIX-1)*NIY+i;
    A(j,j)=AY(i,NIX)+AY(i+1,NIX)+AX(i,NIX)+AX(i,NIX+1);%[1]
    A(j,j-1)=-AY(i,NIX);%[1]
    A(j,j+1)=-AY(i+1,NIX);%[1]
    A(j,j-NIY)=-AX(i,NIX);%[1]
end

%Top boundary
for i=2:NIX-1
    j=(i-1)*NIY+1;
    A(j,j)=AX(1,i)+AY(2,i)+AX(1,i+1)+AY(1,i);%[1]
    A(j,j-NIY)=-AX(1,i);%[1]
    A(j,j+1)=-AY(2,i);%[1]
    A(j,j+NIY)=-AX(1,i+1);%[1]
end

%Bottom boundary
for i=2:NIX-1
    j=i*NIY;
    A(j,j)=AX(NIY,i)+AY(NIY,i)+AX(NIY,i+1)+AY(NIY+1,i);%[1]
    A(j,j-NIY)=-AX(NIY,i);%[1]
    A(j,j-1)=-AY(NIY,i);%[1]
    A(j,j+NIY)=-AX(NIY,i+1);%[1]
end

condition=cond(A)
ddd=zeros(NIX*NIY,1);
for i=1:NIX*NIY
    ddd(i)=sum(A(:,i));
end

% A=zeros(5,NIX*NIY);
% for i=1:NIX*NIY
%
% end

%=========================== [A] MATRIX ===============================


%% ...............=========== THE BOMB CYCLE ============..................

err0RES=0.002;
errRES=1.0;

%NOTE: DFLP at 4 physical walls or boundaries are useless since coefficients are set zero there, thus DFLP has onlt NIX*NIY elements.
DFLP=zeros(NIY*NIX,1);%consider delta(FL*P) as one variable which has (NIX*NIY) elements in one column
B=zeros(NIY*NIX,1);%set as [B] which has NIX*NIY elements in one column

DRFVX=zeros(NIY+2,NIX+1);
DRFVY=zeros(NIY+1,NIX+2);
DVX=zeros(NIY+2,NIX+1);
DVY=zeros(NIY+1,NIX+2);
iters=0;

%while(errRES>err0RES)
    
    %--------------------------- 14. Old RFVX ---------------------------------
    RFVX=zeros(NIY+2,NIX+1);
    for i=1:NIX+1
        for j=1:NIY+2
            RFVX(j,i)=0.5*(RLNTM(j,i)+RLNTM(j,i+1))*0.5*(FLNTM(j,i)+FLNTM(j,i+1))*VXOTM(j,i);%[kg/m^2/sec]
        end
    end
    
    %--------------------------- 14. Old RFVX ---------------------------------
    
    %--------------------------- 15. Old RFVY ---------------------------------
    RFVY=zeros(NIY+1,NIX+2);
    for i=1:NIX+2
        for j=1:NIY+1
            RFVY(j,i)=0.5*(RLNTM(j,i)+RLNTM(j+1,i))*0.5*(FLNTM(j,i)+FLNTM(j+1,i))*VYOTM(j,i);%[kg/m^2/sec]
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
            RFVXI(j,i)=RLNTM(j,i)*FLNTM(j,i)*0.5*(VXOTM(j,i-1)+VXOTM(j,i));%[kg/m^2/sec]
        end
    end
    for i=2:NIY+1
        RFVXI(i,1)=-RFVXI(i,2);%left VX=0.0 boundary
        RFVXI(i,NIX+2)=-RFVXI(i,NIX+1);%rihgt VX=0.0 boundary
    end
    
    for i=1:NIX+2
        RFVXI(1,i)=-RFVXI(2,i);%top VX=0.0 boundary
        RFVXI(NIY+2,i)=-RFVXI(NIY+1,i);%bottom VX=0.0 boundary
    end
    
    %staggered grid (j+0.5,k+0.5), (j-0.5,k+0.5)
    RFVXH=zeros(NIY+1,NIX+1);%H==half grid(staggered grid)
    for i=1:NIX+1
        for j=1:NIY+1
            RFVXH(j,i)=0.25*(RLNTM(j,i)+RLNTM(j+1,i)+RLNTM(j,i+1)+RLNTM(j+1,i+1))*0.25*(FLNTM(j,i)+FLNTM(j+1,i)+FLNTM(j,i+1)+FLNTM(j+1,i+1))*0.5*(VXOTM(j+1,i)+VXOTM(j,i));%[kg/m^2/sec]
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
            RFVYI(j,i)=RLNTM(j,i)*FLNTM(j,i)*0.5*(VYOTM(j-1,i)+VYOTM(j,i));%[kg/m^2/sec]
        end
    end
    
    %VY(1:NIY+2,1)=Free --> RFVYI(1:NIY+1,1)!=0.0 [Slip boundary]
    %This sentence can be added to last FOR loop
    for i=1:NIX+1
        RFVYI(1,i)=-RFVYI(2,i);%top VY=0.0 boundary
        RFVYI(NIY+2,i)=-RFVYI(NIY+1,i);%bottom VY=0.0 boundary
    end
    
    for i=1:NIY+2
        RFVYI(i,NIX+2)=-RFVYI(i,NIX+1);%right VY=0.0 boundary
        %left VY=FREE
    end
    
    %staggered grid (j+0.5,k+0.5), (j+0.5,k-0.5)
    RFVYH=zeros(NIY+1,NIX+1);%H==half grid(staggered grid)
    for i=1:NIX+1
        for j=1:NIY+1
            RFVYH(j,i)=0.25*(RLNTM(j,i)+RLNTM(j+1,i)+RLNTM(j,i+1)+RLNTM(j+1,i+1))*0.25*(FLNTM(j,i)+FLNTM(j+1,i)+FLNTM(j,i+1)+FLNTM(j+1,i+1))*0.5*(VYOTM(j,i+1)+VYOTM(j,i));%[kg/m^2/sec]
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
    
    for i=1:NIY+2
        VVXI(i,1)=-VXOTM(i,2)*max(RFVXI(i,1),0.0)+VXOTM(i,1)*min(RFVXI(i,1),0.0);%VVXI(i,2)-VVXI(i,1)==2.0*VXOTM(i,1)*max(RFVXI(i,2),0.0)==0    left boundary
        VVXI(i,NIX+2)=VXOTM(i,NIX+1)*max(RFVXI(i,NIX+2),0.0)+(-VXOTM(i,NIX))*min(RFVXI(i,NIX+2),0.0);%VVXI(i,NIX+2)-VVXI(i,NIX+1)==-2.0*VXOTM(i,NIX+1)*min(RFVXI(i,NIX+1),0.0)=0    right boundary
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
    
    for i=1:NIX+2
        VVYI(1,i)=-VYOTM(2,i)*max(RFVYI(1,i),0.0)+VYOTM(1,i)*min(RFVYI(1,i),0.0);%VVYI(2,i)-VVYI(1,i)=2.0*VYOTM(1,i)*max(RFVYI(2,i),0.0)=0   top boundary
        VVYI(NIY+2,i)=VYOTM(NIY+1,i)*max(RFVYI(NIY+2,i),0.0)+(-VYOTM(NIY,i)*min(RFVYI(NIY+2,i),0.0));%VVYI(NIY+2,i)-VVYI(NIY+1,i)=-2.0*VYOTM(NIY+1,i)*min(RFVYI(NIY+1,i),0.0)=0   bottom boundary
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
            VDXI(j,i)=mu*(VXOTM(j,i)*0.5*(FLNTM(j,i)+FLNTM(j,i+1))-VXOTM(j,i-1)*0.5*(FLNTM(j,i-1)+FLNTM(j,i)))/dx(i-1);%[Pa.sec times 1/sec == Pa]
        end
    end

    for i=1:NIY+2
        VDXI(i,1)=mu*(0.5*(FLNTM(i,1)+FLNTM(i,2))*VXOTM(i,1)-(-VXOTM(i,2)*0.5*(FLNTM(i,3)+FLNTM(i,2))))/dx(1);%VDXI(i,2)-VDXI(i,1)=-VXOTM(i,1)*(FLNTM(i,2)+FLNTM(i,1))/dx(1)=0   left boundary
        VDXI(i,NIX+2)=mu*(0.5*(FLNTM(i,NIX+1)+FLNTM(i,NIX))*(-VXOTM(i,NIX))-0.5*(FLNTM(i,NIX+1)+FLNTM(i,NIX+2))*VXOTM(i,NIX+1))/dx(NIX);%VDXI(i,NIX+2)-VDXI(i,NIX+1)=-(FLNTM(i,NIX+1)+FLNTM(i,NIX+2))*VXOTM(i,NIX+1)/dx(NIX)=0   right boundary
    end
    %NOTE: you have to set tau_xx|(j+1,k)==tau_xx|(j,k) since partial^2(VX)/partial(x)^2==0.0 at (j+0.5,k)
    
    %VDXI(1:NIY+2,image_1)=0.0 at image_1 and image_2
    %VDXI(1:NIY+2,image_2)=0.0 at image_1 and image_2
    
    VDXH=zeros(NIY+1,NIX+1);
    for i=1:NIX+1
        for j=2:NIY
            VDXH(j,i)=2.0*mu*(0.5*(FLNTM(j+1,i)+FLNTM(j+1,i+1))*VXOTM(j+1,i)-0.5*(FLNTM(j,i)+FLNTM(j,i+1))*VXOTM(j,i))/(dy(j-1)+dy(j));%[Pa.sec times 1/sec == Pa]
        end
    end
    
    for i=1:NIX+1
        VDXH(1,i)=2.0*mu*(0.5*(FLNTM(2,i)+FLNTM(2,i+1))*VXOTM(2,i)-0.5*(FLNTM(1,i)+FLNTM(1,i+1))*VXOTM(1,i))/(2.0*dy(1));%top boundary [Pa.sec times 1/sec == Pa]
        VDXH(NIY+1,i)=2.0*mu*(0.5*(FLNTM(NIY+2,i)+FLNTM(NIY+2,i+1))*VXOTM(NIY+2,i)-0.5*(FLNTM(NIY+1,i)+FLNTM(NIY+1,i+1))*VXOTM(NIY+1,i))/(2.0*dy(NIY));%bottom boundary [Pa.sec times 1/sec == Pa]
    end
    
    %-------------------- 5. x-viscous momentum -------------------------------
    
    %-------------------- 6. y-viscous momentum -------------------------------
    VDYI=zeros(NIY+2,NIX+2);
    for i=1:NIX+2
        for j=2:NIY+1
            VDYI(j,i)=mu*(VYOTM(j,i)*0.5*(FLNTM(j+1,i)+FLNTM(j,i))-VYOTM(j-1,i)*0.5*(FLNTM(j,i)+FLNTM(j-1,i)))/dy(j-1);%[Pa.sec times 1/sec == Pa]
        end
    end
    
    for i=1:NIX+2
        VDYI(1,i)=mu*(VYOTM(1,i)*0.5*(FLNTM(1,i)+FLNTM(2,i))-0.5*(FLNTM(2,i)+FLNTM(3,i))*(-VYOTM(2,i)))/dy(1);%VDYI(2,i)-VDYI(1,i)=-(FLNTM(2,i)+FLNTM(1,i))*VYOTM(1,i)/dy(1)=0   top boundary
        VDYI(NIY+2,i)=mu*(0.5*(FLNTM(NIY+1,i)+FLNTM(NIY,i))*(-VYOTM(NIY,i))-VYOTM(NIY+1,i)*0.5*(FLNTM(NIY+1,i)+FLNTM(NIY+2,i)))/dy(NIY);%VDYI(NIY+2,i)-VDYI(NIY+1,i)=-VYOTM(NIY+1,i)*(FLNTM(NIY+1,i)+FLNTM(NIY+2,i))/dy(NIY)=0   bottom boundary
    end
    %NOTE: you have to set tau_yy|(j,k+1)==tau_yy|(j,k) since partial^2(VY)/partial(y)^2==0.0 at (j,k+0.5)
    
    VDYH=zeros(NIY+1,NIX+1);
    for i=2:NIX
        for j=1:NIY+1
            VDYH(j,i)=2.0*mu*(0.5*(FLNTM(j+1,i+1)+FLNTM(j,i+1))*VYOTM(j,i+1)-0.5*(FLNTM(j+1,i)+FLNTM(j,i))*VYOTM(j,i))/(dx(i-1)+dx(i));%[Pa.sec times 1/sec == Pa]
        end
    end
    
    for i=1:NIY+1
        VDYH(i,1)=2.0*mu*(0.5*(FLNTM(i+1,2)+FLNTM(i,2))*VYOTM(i,2)-0.5*(FLNTM(i+1,1)+FLNTM(i,1))*VYOTM(i,1))/(2.0*dx(1));%left slip boundary  [Pa.sec times 1/sec == Pa]
        VDYH(i,NIX+1)=2.0*mu*(0.5*(FLNTM(i+1,NIX+2)+FLNTM(i,NIX+2))*VYOTM(i,NIX+2)-0.5*(FLNTM(i+1,NIX+1)+FLNTM(i,NIX+1))*VYOTM(i,NIX+1))/(2.0*dx(NIX));%right boundary  [Pa.sec times 1/sec == Pa]
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
        VYCONV(1,i)=-dtb*(VVYH(1,i+1)-VVYH(1,i))/dx(i)-2.0*dtb*(VVYI(2,i+1)-VVYI(1,i+1))/(2.0*dy(1));%top boundary [Pa.sec/m]
        VYCONV(NIY+1,i)=-dtb*(VVYH(NIY+1,i+1)-VVYH(NIY+1,i))/dx(i)-2.0*dtb*(VVYI(NIY+2,i+1)-VVYI(NIY+1,i+1))/(2.0*dy(NIY));%bottom boundary [Pa.sec/m]
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
            FPDX(j,i)=-2.0*dtb*(FLOTM(j,i+1)*(POTM(j,i+1)-P0(j,i+1))-FLOTM(j,i)*(POTM(j,i)-P0(j,i)))/(dx(i)+dx(i-1));%[Pa.sec/m == kg/m^2/sec]
        end
    end
    %NOTE: pressure with 13 digits precision; otherwise numerical precision
    %will give non-zero results at t=0.0 sec
    
    for i=1:NIY+2
        FPDX(i,1)=-2.0*dtb*(FLOTM(i,2)*(POTM(i,2)-P0(i,2))-FLOTM(i,1)*(POTM(i,1)-P0(i,1)))/(2.0*dx(1));%left boundary [Pa.sec/m == kg/m^2/sec]
        FPDX(i,NIX+1)=-2.0*dtb*(FLOTM(i,NIX+2)*(POTM(i,NIX+2)-P0(i,NIX+2))-FLOTM(i,NIX+1)*(POTM(i,NIX+1)-P0(i,NIX+1)))/(2.0*dx(NIX));%right boundary [Pa.sec/m == kg/m^2/sec]
    end
    
    %-------------------- 11. x- pressure gradient ----------------------------
    
    %---------------------- 12. y- pressure part ------------------------------
    
    FPDY=zeros(NIY+1,NIX+2);
    for i=1:NIX+2
        for j=2:NIY
            FPDY(j,i)=-2.0*dtb*(FLOTM(j+1,i)*(POTM(j+1,i)-P0(j+1,i))-FLOTM(j,i)*(POTM(j,i)-P0(j,i)))/(dy(j)+dy(j-1));%[Pa.sec/m == kg/m^2/sec]
        end
    end
    
    for i=1:NIX+2
        FPDY(1,i)=-2.0*dtb*(FLOTM(2,i)*(POTM(2,i)-P0(2,i))-FLOTM(1,i)*(POTM(1,i)-P0(1,i)))/(2.0*dy(1));%top boundary  [Pa.sec/m == kg/m^2/sec]
        FPDY(NIY+1,i)=-2.0*dtb*(FLOTM(NIY+2,i)*(POTM(NIY+2,i)-P0(NIY+2,i))-FLOTM(NIY+1,i)*(POTM(NIY+1,i)-P0(NIY+1,i)))/(2.0*dy(NIY));%bottom boundary  [Pa.sec/m == kg/m^2/sec]
    end
    
    %---------------------- 12. y- pressure part ------------------------------
    
    %----------------------- 13. y- gravity part ------------------------------
    
    FRLOTM=zeros(NIY+1,NIX+2);
    FRLNTM=zeros(NIY+1,NIX+2);
    for i=1:NIX+2
        for j=1:NIY+1
            FRLOTM(j,i)=0.5*(FLOTM(j,i)+FLOTM(j+1,i))*(0.5*(RLOTM(j,i)+RLOTM(j+1,i))-0.5*(RL0(j,i)+RL0(j+1,i)));
            FRLNTM(j,i)=0.5*(FLNTM(j,i)+FLNTM(j+1,i))*(0.5*(RLNTM(j,i)+RLNTM(j+1,i))-0.5*(RL0(j,i)+RL0(j+1,i)));
            %FRLNTM(j,i)=FRLOTM(j,i);%For constant density. If this is set,
            %then sum(B)=0.0 like lid driven cavity flow
        end
    end
    
    FRG=zeros(NIY+1,NIX+2);
    for i=1:NIX+2
        for j=2:NIY
            FRG(j,i)=0.5*g*dtb*(FRLOTM(j,i)+FRLNTM(j,i));%[kg/m^2/sec]
        end
    end
    
    for i=1:NIX+2
        FRG(1,i)=0.0*0.5*(0.5*g*dtb*(FRLOTM(1,i)+FRLNTM(1,i)));%only half cell at top boundary
        FRG(NIY+1,i)=0.0*0.5*(0.5*g*dtb*(FRLOTM(NIY+1,i)+FRLNTM(NIY+1,i)));%only half cell at bottom boundary
    end
    %NOTE: In Boussinesq assumption, only change of density in gravity term
    %has effect; in present program, we include all density change effects!
    %----------------------- 13. y- gravity part ------------------------------
    
    %===================== APPROXIMATION OF RFVXY =========================
    %--------------------------- RFVXTM -----------------------------------
    RFVXTM=zeros(NIY+2,NIX+1);
    for i=1:NIX+1
        for j=2:NIY+1
            RFVXTM(j,i)=(RFVX(j,i)+VXCONV(j-1,i)+VXDIF(j-1,i)+FPDX(j,i))/(1.0+dtb*UFRKX(j,i));%[kg/m^2/sec]
        end
    end
    RFVXTM(1:NIY+2,1)=0.0;%left boundary, to make sure
    RFVXTM(1:NIY+2,NIX+1)=0.0;%right boundary, to make sure
    
    %--------------------------- RFVXTM -----------------------------------
    
    %--------------------------- RFVYTM -----------------------------------
    RFVYTM=zeros(NIY+1,NIX+2);
    
    for i=2:NIX+1
        for j=2:NIY
            RFVYTM(j,i)=(RFVY(j,i)+VYCONV(j,i-1)+VYDIF(j,i-1)+FPDY(j,i)+FRG(j,i))/(1.0+dtb*UFRKY(j,i));%[kg/m^2/sec]
        end
    end
    
    for i=1:NIX+2
        RFVYTM(1,i)=0.0;%(RFVY(1,i)+VYCONV(1,i-1)+VYDIF(1,i-1)+FPDY(1,i)+FRG(1,i))/(1.0+dtb*UFRKY(1,i));%top boundary, [kg/m^2/sec]
        RFVYTM(NIY+1,i)=0.0;%(RFVY(NIY+1,i)+VYCONV(NIY+1,i-1)+VYDIF(NIY+1,i-1)+FPDY(NIY+1,i)+FRG(NIY+1,i))/(1.0+dtb*UFRKY(NIY+1,i));%bottom boundary, [kg/m^2/sec]
    end
    
%     for i=2:NIY
%         RFVYTM(i,1)=RFVYTM(i,2);%left boundary
%         RFVYTM(i,NIX+2)=-RFVYTM(i,NIX+1);%right boundary
%     end
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
            RES(j,i)=dy(j)*(RFVXTM(j+1,i+1)-RFVXTM(j+1,i))+dx(i)*(RFVYTM(j+1,i+1)-RFVYTM(j,i+1))+0.5*dFSTM(j+1,i+1)*(RSOTM(j+1,i+1)+RSNTM(j+1,i+1))*dx(i)*dy(j)/dtb+bsq*(FRLNTM(j+1,i+1)-FRLOTM(j+1,i+1))*dx(i)*dy(j)/dtb;%[kg/m/sec == Pa.sec]
        end
    end
    
    for i=2:NIX-1
        RES(1,i)=dy(1)*(RFVXTM(2,i+1)-RFVXTM(2,i))+dx(i)*(RFVYTM(2,i+1)-RFVYTM(1,i+1))+0.5*dFSTM(2,i+1)*(RSOTM(2,i+1)+RSNTM(2,i+1))*dx(i)*dy(1)/dtb+bsq*(FRLNTM(2,i+1)-FRLOTM(2,i+1))*dx(i)*dy(1)/dtb;%top boundary
        RES(NIY,i)=dy(NIY)*(RFVXTM(NIY+1,i+1)-RFVXTM(NIY+1,i))+dx(i)*(RFVYTM(NIY+1,i+1)-RFVYTM(NIY,i+1))+0.5*dFSTM(NIY+1,i+1)*(RSOTM(NIY+1,i+1)+RSNTM(NIY+1,i+1))*dx(i)*dy(NIY)/dtb+bsq*(FRLNTM(NIY+1,i+1)-FRLOTM(NIY+1,i+1))*dx(i)*dy(NIY)/dtb;%bottom boundary
    end
    
    for i=2:NIY-1
        RES(i,1)=dy(i)*(RFVXTM(i+1,2)-RFVXTM(i+1,1))+dx(1)*(RFVYTM(i+1,2)-RFVYTM(i,2))+0.5*dFSTM(i+1,2)*(RSOTM(i+1,2)+RSNTM(i+1,2))*dx(1)*dy(i)/dtb+bsq*(FRLNTM(i+1,2)-FRLOTM(i+1,2))*dx(1)*dy(i)/dtb;%left boundary
        RES(i,NIX)=dy(i)*(RFVXTM(i+1,NIX+1)-RFVXTM(i+1,NIX))+dx(NIX)*(RFVYTM(i+1,NIX+1)-RFVYTM(i,NIX+1))+0.5*dFSTM(i+1,NIX+1)*(RSOTM(i+1,NIX+1)+RSNTM(i+1,NIX+1))*dx(NIX)*dy(i)/dtb+bsq*(FRLNTM(i+1,NIX+1)-FRLOTM(i+1,NIX+1))*dx(NIX)*dy(i)/dtb;%right boundary
    end
    
    RES(1,1)=dy(1)*(RFVXTM(2,2)-RFVXTM(2,1))+dx(1)*(RFVYTM(2,2)-RFVYTM(1,2))+0.5*dFSTM(2,2)*(RSOTM(2,2)+RSNTM(2,2))*dx(1)*dy(1)/dtb+bsq*(FRLNTM(2,2)-FRLOTM(2,2))*dx(1)*dy(1)/dtb;
    RES(1,NIX)=dy(1)*(RFVXTM(2,NIX+1)-RFVXTM(2,NIX))+dx(NIX)*(RFVYTM(2,NIX+1)-RFVYTM(1,NIX+1))+0.5*dFSTM(2,NIX+1)*(RSOTM(2,NIX+1)+RSNTM(2,NIX+1))*dx(NIX)*dy(1)/dtb+bsq*(FRLNTM(2,NIX+1)-FRLOTM(2,NIX+1))*dx(NIX)*dy(1)/dtb;
    RES(NIY,1)=dy(NIY)*(RFVXTM(NIY+1,2)-RFVXTM(NIY+1,1))+dx(1)*(RFVYTM(NIY+1,2)-RFVYTM(NIY,2))+0.5*dFSTM(NIY+1,2)*(RSOTM(NIY+1,2)+RSNTM(NIY+1,2))*dx(1)*dy(NIY)/dtb+bsq*(FRLNTM(NIY+1,2)-FRLOTM(NIY+1,2))*dx(1)*dy(NIY)/dtb;
    RES(NIY,NIX)=dy(NIY)*(RFVXTM(NIY+1,NIX+1)-RFVXTM(NIY+1,NIX))+dx(NIX)*(RFVYTM(NIY+1,NIX+1)-RFVYTM(NIY,NIX+1))+0.5*dFSTM(NIY+1,NIX+1)*(RSOTM(NIY+1,NIX+1)+RSNTM(NIY+1,NIX+1))*dx(NIX)*dy(NIY)/dtb+bsq*(FRLNTM(NIY+1,NIX+1)-FRLOTM(NIY+1,NIX+1))*dx(NIX)*dy(NIY)/dtb;
    %RES(1,1)=2.0*RES(1,1);
    
    for i=1:NIX
        B(1+(i-1)*NIY:i*NIY)=-RES(1:NIY,i)/dtb;%+A(1+(i-1)*NIY:i*NIY,1+(i-1)*NIY:i*NIY)*DFLP(1+(i-1)*NIY:i*NIY);%[kg/m/sec == Pa.sec]
    end
    %B(1)=0.0;
    
    RHS=sum(B)
    %----------------------------- RES ------------------------------------
    %=================== MASS CONSERVATION RESIDUAL  ======================
    
    %======================= PRESSURE CORRECTION ==========================
    %Solve linear equations [A][X]=[B]
    %DFLP=A\B;%[kg/m/sec^2 == Pa]
    DFLP=GaussSeidel(A,B);
    
    %-------------------------- New Pressure ------------------------------
    for i=2:NIX+1
        for j=2:NIY+1
            k=(i-2)*NIY+j-1;
            POTM(j,i)=POTM(j,i)+DFLP(k)/FLNTM(j,i);
        end
    end
    
    %Boundary pressure extrapolation
    for i=2:NIY+2
        POTM(i,1)=POTM(i,2)-(POTM(i,3)-POTM(i,2))*dx(1)/((dx(1)+dx(2))/2.0);%left boundary
        POTM(i,NIX+2)=POTM(i,NIX+1)-(POTM(i,NIX)-POTM(i,NIX+1))*dx(NIX)/((dx(NIX)+dx(NIX-1))/2.0);%right boundary
    end
    for i=1:NIX+2
        POTM(1,i)=POTM(2,i)-dy(1)*(POTM(3,i)-POTM(2,i))/((dy(1)+dy(2))/2.0);%top boundary
        POTM(NIY+2,i)=POTM(NIY+1,i)-dy(NIY)*(POTM(NIY,i)-POTM(NIY+1,i))/((dy(NIY)+dy(NIY-1))/2.0);%bottom boundary
    end
    %-------------------------- New Pressure ------------------------------
    
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
    
    %----------------------------- DRFVX ----------------------------------
    
    %----------------------------- DRFVY ----------------------------------
    for i=2:NIX+1
        %DRFVY(1:NIY+1,1)=FREE
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
%         RFVX=RFVXTM+DRFVX;
%         RFVY=RFVYTM+DRFVY;

LL=zeros(NIY,NIX);%to test whether mass conservation is met
for i=1:NIX
    for j=1:NIY
        LL(j,i)=dy(j)*(DRFVX(j+1,i+1)-DRFVX(j+1,i))+dx(i)*(DRFVY(j+1,i+1)-DRFVY(j,i+1));
    end
end
    
    %------------------------- NEW RFV ESTIMATION -------------------------
    
    %------------------------------ NEW VX ------------------------------------
    for i=2:NIX
        %VX(1,1:NIX+1)=0.0 always
        %VX(NIY+2,1:NIX+1)=0.0 always
        for j=2:NIY+1
            %VX(1:NIY+2,1)=0.0 always
            %VX(1:NIY+2,NIX+1)=0.0 always
            DVX(j,i)=DRFVX(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j,i+1))*0.5*(FLNTM(j,i)+FLNTM(j,i+1)));%[m/sec]
            VXOTM(j,i)=VXOTM(j,i)+DVX(j,i);
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
    for i=1:NIX+1
        %VY(1,1:NIX+2)=0.0 always
        %VY(NIY+1,1:NIX+2)=0.0 always
        for j=2:NIY
            %VY(2:NIY,1)=Free [Slip boundary condition]
            %VY(1:NIY+1,NIX+2)=0.0 always
            DVY(j,i)=DRFVY(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j+1,i))*0.5*(FLNTM(j,i)+FLNTM(j+1,i)));%[m/sec]
            VYOTM(j,i)=VYOTM(j,i)+DVY(j,i);
        end
    end
    
    %NOTE:
    VYOTM(1,1:NIX+2)=0.0;
    VYOTM(NIY+1,1:NIX+2)=0.0;
    
    for i=1:NIY+1
        VYOTM(i,1)=VYOTM(i,2);%slip boundary
        VYOTM(i,NIX+2)=-VYOTM(i,NIX+1);
    end
    
    %------------------------------ NEW VY ------------------------------------
%     [T,FL,CL,dFS]=TFCPV(T,FL,FS,KL,KS,CL,CS,VXOTM,VYOTM,FRCPS,FRCS,FK,RSCPS,RLCPL,RS,RL,CpL);
    
    
    %------------------------- RFVX RFVY error ----------------------------
    %     errDVX=abs(max(max(DVX)));
    %     errDVY=abs(max(max(DVY)));
    %     errDV=max(errDVY,errDVX);
    errRES=max(abs(B));
    %------------------------- RFVX RFVY error ----------------------------
    
    iters=iters+1;
%end

P=POTM;
VX=VXOTM;
VY=VYOTM;

% %------------------------------ RETURN NEW VX ------------------------------------
% for i=2:NIX
%     %VX(1,1:NIX+1)=0.0 always
%     %VX(NIY+2,1:NIX+1)=0.0 always
%     for j=2:NIY+1
%         %VX(1:NIY+2,1)=0.0 always
%         %VX(1:NIY+2,NIX+1)=0.0 always
%         VX(j,i)=VXOTM(j,i);%RFVX(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j,i+1))*0.5*(FLNTM(j,i)+FLNTM(j,i+1)));%[m/sec]
%     end
% end
%
% %------------------------------ RETURN NEW VX ------------------------------------
%
% %------------------------------ RETURN NEW VY ------------------------------------
% for i=1:NIX+1
%     %VY(1,1:NIX+2)=0.0 always
%     %VY(NIY+1,1:NIX+2)=0.0 always
%     for j=2:NIY
%         %VY(2:NIY,1)=Free [Slip boundary condition]
%         %VY(1:NIY+1,NIX+2)=0.0 always
%         VY(j,i)=VYOTM(j,i);%RFVY(j,i)/(0.5*(RLNTM(j,i)+RLNTM(j+1,i))*0.5*(FLNTM(j,i)+FLNTM(j+1,i)));%[m/sec]
%     end
% end
% %------------------------------ RETURN NEW VY ------------------------------------


end
%
