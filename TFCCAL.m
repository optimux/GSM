function [T,FL,CL,dFS,iters]=TFCCAL(TTM,FLTM,FSTM,KLTM,KSTM,CLTM,CSTM,VXTM,VYTM,FRCPSTM,FRCSTM,FKTM,RSCPSTM,RLCPLTM,RSTM,RLTM,CPLTM,dtb)%FRCSTM,FRHTM,
%To iteratively solve energy and species conservation equations with finite
%volume method. Model, discretization and iteration can be found "A
%Continuum Model for Computer Simulation of Macrosegregations in Ingots
%during solidification" by Daming Xu 1989 and his companion papers in 1991

%Created 2019-10-17

%TTM: [NIY+NMY+2,NIX+NMX+2]
%FLTM: [NIY+2,NIX+2]
%FSTM: [NIY+2,NIX+2]
%KLTM: [NIY+2,NIX+2]
%KSTM: [NIY+2,NIX+2]

%FKTM: [NIY+2,NIX+2]

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
global dx
global dy
global gap
global DL
global CE

err0dFS=0.0001;%controlled accuracy for dFS, T, CL
err0TTM=0.0001;
err0CLTM=0.0001;

TVB=zeros(NIY+NMY+2,1);%TVB: T of vertical boundary [K]
THB=zeros(NIX+NMX+2,1);%THB: T of horizontal boundary [K]

%======================= VERY IMPORTANT NOTE ==============================
%It is highly recommended the usage of packages for complex formulae, i.e.,
%temporary variables. These packages can be utilized repetedly, easily
%modified, and produce a concise numerical equation.
%======================= VERY IMPORTANT NOTE ==============================


%################################################################### ENERGY BALANCE #########################################################################
%% .............========= DIFFUSION HEAT FLUX ==========...................
%[FRCPSTM,FRCSTM,FKTM,FR] = INTSP(RSX);%integrate solid properties: intg(FS*RS*CpS), intg(FS*RS*CS), intg(FS*KS), intg(FS*RS); [NIY+2,NIX+2]

%predefine mould thermal conductivity from D-9 [W/mm/K]
KM=zeros(NIY+1+NMY+1,NIX+1+NMX+1);
for i=1:NIX+NMX+2
    for j=1:NIY+NMY+2
        if(TTM(j,i)<601.61)%601.61 [K]=328.46+273.15
            KM(j,i)=7.9718e-2-1.0387e-4*(TTM(j,i)-273.15);
        else
            KM(j,i)=5.43102e-2-2.6516e-5*(TTM(j,i)-273.15);
        end
    end
end
KM=KM*1000.0;%[W/m/K]
KA=KM;%thermal conductivity of system: mould part [W/m/K]
KA(1:NIY+1,1:NIX+1)=FSTM(1:NIY+1,1:NIX+1).*KSTM(1:NIY+1,1:NIX+1)+FLTM(1:NIY+1,1:NIX+1).*KLTM(1:NIY+1,1:NIX+1);%thermal conductivity of system: Ingot part [W/m/K]; it should be FS.*KS+FL.*KL, but CS is NaN
%NOTE: mean thermal conductivity is based on volume fraction, not mass
%fraction. This is also described in Table 1 of Bennon and Incropera 1988!

%------------- thermal conductivity on finite volume faces ----------------
%....................... all faces, x-axis [W/m/K] ........................
KX=zeros(NIY+NMY+2,NIX+NMX+1);
%=== (1)Ingot ===
for i=1:NIX
    for j=1:NIY+1
        KX(j,i)=0.5*(FKTM(j,i)+FKTM(j,i+1)+FLTM(j,i)*KLTM(j,i)+FLTM(j,i+1)*KLTM(j,i+1));%[(FS*KS)+(FL*KL)]|(j+-0.5,k)
    end
end

%=== (2)Vertical Gas gap ===
%KX(1:NIY+1,NIX+1) is vertical gas gap, heat flux is directly calculated, no KX needed!
%KX(NIY+2:NIY+NMY+2,NIX+1)=0.5*(KM(NIY+2:NIY+NMY+2,NIX+1)+KM(NIY+2:NIY+NMY+2,NIX+2));

%=== (3)Mould ===
for i=NIX+2:NIX+NMX+1
    for j=1:NIY+1
        KX(j,i)=0.5*(KA(j,i)+KA(j,i+1));
    end
end

for i=1:NIX+NMX+1
    for j=NIY+2:NIY+NMY+2
        KX(j,i)=0.5*(KA(j,i)+KA(j,i+1));
    end
end

%=== (4)Vertical boundary ===
%KX(1:NIY+NMY+2,NIX+NMX+1) is also the mean of TTM(:,NIX+NMX+1) and TTM(:,NIX+NMX+2); used for radiation-convection boundary condition
%........................all faces, x-axis [W/m/K] ........................

%........................all faces, y-axis [W/m/K] ........................
KY=zeros(NIY+NMY+1,NIX+NMX+2);
%=== (1)Ingot ===
for i=1:NIX+1
    for j=1:NIY
        KY(j,i)=0.5*(FKTM(j,i)+FKTM(j+1,i)+FLTM(j,i)*KLTM(j,i)+FLTM(j+1,i)*KLTM(j+1,i));%[(FS*KS)+(FL*KL)]|(j,k+-0.5)
    end
end

%=== (2)Horizontal gas gap ===
%KY(NIY+1,1:NIX+1) is herizontal gas gap, heat flux is directly calculated, no KY needed!

%=== (3)Mould ===
for i=NIX+2:NIX+NMX+2
    for j=1:NIY+1
        KY(j,i)=0.5*(KA(j,i)+KA(j+1,i));
    end
end

for i=1:NIX+NMX+2
    for j=NIY+2:NIY+NMY+1
        KY(j,i)=0.5*(KA(j,i)+KA(j+1,i));
    end
end

%=== (4)Herizontal boundary ===
%KY(NIY+NMY+1,:) is also the mean of TTM(NIY+NMY+1,:) and TTM(NIY+NMY+2,:); used for radiation-convection boundary condition
%........................all faces, y-axis [W/m/K] ........................

%------------- thermal conductivity on finite volume faces ----------------

%----------------------------- T gradient ---------------------------------
[TVB,THB]= BoundTemp(TTM(:,NIX+NMX+1),TTM(NIY+NMY+1,:),KA(:,NIX+NMX+1),KA(NIY+NMY+1,:));%calculate boundary temperature via conduction==radiation + convection
TTM(:,NIX+NMX+2)=TVB;

DTX=zeros(NIY+NMY+2,NIX+NMX+1);%T gradient between finite volumes in x-axis
for i=2:NIX+NMX
    for j=1:NIY+NMY+2
        DTX(j,i)=(TTM(j,i+1)-TTM(j,i))/(0.5*dx(i-1)+0.5*dx(i));%[K/m]
    end
end
%gas gap is also calculated as this firstly, it will be modified later

for i=1:NIY+NMY+2
    DTX(i,1)=(TTM(i,2)-TTM(i,1))/(0.5*dx(1));%1st column == 0.0
    DTX(i,NIX+NMX+1)=(TTM(i,NIX+NMX+2)-TTM(i,NIX+NMX+1))/(0.5*dx(NIX+NMX));%outmost boundary of mould
end

TTM(NIY+NMY+2,:)=THB;

DTY=zeros(NIY+NMY+1,NIX+NMX+2);%T gradient between finite volumes in y-axis
for i=1:NIX+NMX+2
    for j=2:NIY+NMY
        DTY(j,i)=(TTM(j+1,i)-TTM(j,i))/(0.5*dy(j-1)+0.5*dy(j));%[K/m]
    end
end
%gas gap is also calculated as this firstly, it will be modified later

for i=1:NIX+NMX+2
    DTY(1,i)=(TTM(2,i)-TTM(1,i))/(0.5*dy(1));%1st row == 0.0
    DTY(NIY+NMY+1,i)=(TTM(NIY+NMY+2,i)-TTM(NIY+NMY+1,i))/(0.5*dy(NIY+NMY));%heat loss to the ambient through bottom surface
end

%----------------------------- T gradient ---------------------------------

%----------------------------- Heat flux ----------------------------------
[FVB,FHB] = BoundTFlux(TTM(:,NIX+NMX+1),TTM(NIY+NMY+1,:),KA(:,NIX+NMX+1),KA(NIY+NMY+1,:));

%heat flux on all x faces
TQDX=zeros(NIY+NMY+2,NIX+NMX+1);%T=T, Q=flux, D==Diffusion, X=x-axis
for i=1:NIX+NMX+1
    for j=1:NIY+NMY+2
        TQDX(j,i)=KX(j,i)*DTX(j,i);%[W/m^2]
    end
end

for i=2:NIY+NMY+1
    TQDX(i,NIX+NMX+1)=-FVB(i-1);%right most boundary of mould
end

%vertical gas gap
for i=1:NIY+1
    TIMX=TTM(i,NIX+1)-0.5*dxI*(TTM(i,NIX+1)-TTM(i,NIX+2))/(0.5*dxI+0.5*dxM+gap);%T at the face between NIX+1 and NIX+2 column by linear interpolation
    Katm=2.72835e-5+4.898321e-8*(TIMX-273.15);%atmosphere thermal conductivity from D-11 [W/mm/K]
    Katm=Katm*1000.0;%[W/m/K]
    TQDX(i,NIX+1)=(TTM(i,NIX+2)-TTM(i,NIX+1))/(0.5*dx(NIX)/KA(i,NIX+1)+gap/Katm+0.5*dx(NIX+1)/KA(i,NIX+2));%[W/m^2]
end

%heat flux on all y faces
TQDY=zeros(NIY+NMY+1,NIX+NMX+2);%T=T, Q=flux, D==Diffusion, Y=y-axis
for i=1:NIX+NMX+2
    for j=1:NIY+NMY+1
        TQDY(j,i)=KY(j,i)*DTY(j,i);%[W/m^2]
    end
end

for i=2:NIX+NMX+1
    TQDY(NIY+NMY+1,i)=-FHB(i-1);%bottom boundary of mould
end

%horizontal gas gap
for i=1:NIX+1
    TIMY=TTM(NIY+1,i)-0.5*dyI*(TTM(NIY+1,i)-TTM(NIY+2,i))/(0.5*dyI+0.5*dyM+gap);%T at the face between NIY+1 and NIY+2 column by linear interpolation
    Katm=2.72835e-5+4.898321e-8*(TIMY-273.15);%atmosphere thermal conductivity from D-11c[W/mm/K]
    Katm=Katm*1000.0;%[W/m/K]
    TQDY(NIY+1,i)=(TTM(NIY+2,i)-TTM(NIY+1,i))/(0.5*dy(NIY)/KA(NIY+1,i)+gap/Katm+0.5*dy(NIY+1)/KA(NIY+2,i));%[W/m^2]
end

%total heat flux of each finite volume
TQDT=zeros(NIY+NMY,NIX+NMX);%T=T, Q=flux, D==Diffusion, T=total
for i=1:NIX+NMX
    for j=1:NIY+NMY
        TQDT(j,i)=-2.0*(TQDX(j+1,i)-TQDX(j+1,i+1))*dtb/dx(i)-2.0*(TQDY(j,i+1)-TQDY(j+1,i+1))*dtb/dy(j);%[J/m^3]
    end
end

%Modify vertical gas gap
for i=1:NIY
    TQDT(i,NIX)=-2.0*(TQDX(i+1,NIX)-TQDX(i+1,NIX+1))*dtb/(dx(NIX)+0.5*gap)-2.0*(TQDY(i,NIX+1)-TQDY(i+1,NIX+1))*dtb/dy(i);%[J/m^3]
    TQDT(i,NIX+1)=-2.0*(TQDX(i+1,NIX+1)-TQDX(i+1,NIX+2))*dtb/(dx(NIX+1)+0.5*gap)-2.0*(TQDY(i,NIX+2)-TQDY(i+1,NIX+2))*dtb/dy(i);%[J/m^3]
end

%Modify horizontal gas gap
for i=1:NIX
    TQDT(NIY,i)=-2.0*(TQDX(NIY+1,i)-TQDX(NIY+1,i+1))*dtb/dx(i)-2.0*(TQDY(NIY,i+1)-TQDY(NIY+1,i+1))*dtb/(dy(NIY)+0.5*gap);%[J/m^3]
    TQDT(NIY+1,i)=-2.0*(TQDX(NIY+2,i)-TQDX(NIY+2,i+1))*dtb/dx(i)-2.0*(TQDY(NIY+1,i+1)-TQDY(NIY+2,i+1))*dtb/(dy(NIY+1)+0.5*gap);%[J/m^3]
end

%----------------------------- Heat flux ----------------------------------

%% ..............========= INTERNAL HEAT STORAGE ========..................
IHS=zeros(NIY+NMY+2,NIX+NMX+2);%old T related Internal Heat Storage

%=== (1)ingot ===
for i=1:NIX+1
    for j=1:NIY+1
        IHS(j,i)=FRCPSTM(j,i)*TTM(j,i)+FLTM(j,i)*RLTM(j,i)*CPLTM(j,i)*TTM(j,i);%[J/m^3]
    end
end

%=== (2)mould ===
RCPM=zeros(NIY+NMY+2,NIX+NMX+2);%mould product of density*heat capacity [J/mm^3/K]; ingot part is useless
for i=1:NIX+NMX+2
    for j=1:NIY+NMY+2
        if(TTM(j,i)<=935.15)%935.15 [K]=662+273.15
            RCPM(j,i)=4.00115e-3+1.000725e-6*(TTM(j,i)-273.15);
        else
            RCPM(j,i)=1.733164e-3+4.425e-6*(TTM(j,i)-273.15);
        end
    end
end
RCPM=RCPM*10^9;%[J/m^3/K]

for i=NIX+2:NIX+NMX+2
    for j=1:NIY+NMY+2
        IHS(j,i)=RCPM(j,i)*TTM(j,i);%[J/m^3]
    end
end

for i=1:NIX+1
    for j=NIY+2:NIY+NMY+2
        IHS(j,i)=RCPM(j,i)*TTM(j,i);%[J/m^3]
    end
end


%% ..............========= CONVECTION HEAT FLUX =========..................

%------------------------ Part One: RFVX RFVY -----------------------------
%                          mass flux density

%x-axis velocity
%         -------
%        |       |
%    --> |   *   |  -->
%        |       |
%         -------
%NOTE: in fact, vx at loft and bottom faces should be estimated,
%but they can be interpolated; so only vx normal to vertical faces are
%calculated!

RFVX=zeros(NIY+NMY+2,NIX+NMX+1);%RL*FL*VX, all VX faces
for i=2:NIX
    %VX(1:NIY+1,1)=0.0 --> RFVX(1:NIY+1,1)=0.0 (left boundary condition in ingot)
    %VX(1:NIY+1,NIX+1)=0.0 --> RFVX(1:NIY+1,NIX+1)=0.0 (right boundary condition in ingot)
    %VX(1:NIY+NMY+2,NIX+2:NIX+NMX+1)=0.0 --> RFVX(1:NIY+NMY+2,NIX+2:NIX+NMX+1)=0.0 (in right mould)
    %VX(NIY+2:NIY+NMY+2,1:NIX+NIX+1)=0.0 --> RFVX(NIY+2:NIY+NMY+2,1:NIX+NIX+1)=0.0 (in bottom mould)
    for j=2:NIY+1
        %VX(1,1:NIX+1)=0.0 --> RFVX(1,1:NIX+1)=0.0 (top boundary condition in ingot)
        RFVX(j,i)=VXTM(j,i)*0.5*(RLTM(j,i)+RLTM(j,i+1))*0.5*(FLTM(j,i)+FLTM(j,i+1));%[kg/m^2/sec]
    end
end

%y-axis velocity
%            ^
%            |
%         -------
%        |       |
%        |   *   |
%        |       |
%         -------
%            ^
%            |
%NOTE: in fact, vy at left and right faces should be estimated,
%but they can be interpolated; so only vy normal to horizontal faces are
%calculated!

RFVY=zeros(NIY+NMY+1,NIX+NMX+2);%RL*FL*VY, all VY faces
for i=2:NIX+1
    for j=2:NIY
        %VY(1,1:NIX+2)=0.0 --> RFVY(1,1:NIX+2)=0.0 (top boundary condition in ingot)
        %VY(NIY+1,1:NIX+2)=0.0 --> RFVY(NIY+1,1:NIX+2)=0.0 (bottom boundary condition in ingot)
        %VY(NIY+2:NIY+NMY+1,1:NIX+2)=0.0 --> RFVY(NIY+2:NIY+NMY+1,1:NIX+2)=0.0 (in bottom  mould)
        %VY(1:NIY+NMY+1,NIX+2:NIX+NMX+2)=0.0 --> RFVY(1:NIY+NMY+1,NIX+2:NIX+NMX+2)=0.0 (in right mould)
        RFVY(j,i)=VYTM(j,i)*0.5*(RLTM(j,i)+RLTM(j+1,i))*0.5*(FLTM(j,i)+FLTM(j+1,i));%[kg/m^2/sec]
    end
end

for i=1:NIY+1
    RFVY(i,1)=RFVY(i,2);%FREE (left boundary in ingot)
end
%------------------------ Part One: RFVX RFVY -----------------------------

%----------------------- Part Two: TRFVX TRFVY ----------------------------
TRFVX=zeros(NIY+NMY+2,NIX+NMX+1);
for i=2:NIX
    %VX(1:NIY+1,1)=0.0 --> TRFVX(1:NIY+1,1)=0.0 (left boundary condition in ingot)
    %VX(1:NIY+1,NIX+1)=0.0 --> TRFVX(1:NIY+1,NIX+1)=0.0 (right boundary condition in ingot)
    %VX(1:NIY+NMY+2,NIX+2:NIX+NMX+1)=0.0 --> TRFVX(1:NIY+NMY+2,NIX+2:NIX+NMX+1)=0.0 (in right mould)
    %VX(NIY+2:NIY+NMY+2,1:NIX+NIX+1)=0.0 --> TRFVX(NIY+2:NIY+NMY+2,1:NIX+NIX+1)=0.0 (in bottom mould)
    for j=2:NIY+1
        %VX(1,1:NIX+1)=0.0 --> TRFVX(1,1:NIX+1)=0.0 (top boundary condition in ingot)
        TRFVX(j,i)=TTM(j,i)*max(RFVX(j,i),0.0)+TTM(j,i+1)*min(RFVX(j,i),0.0);%[kg.K/m^2/sec]
    end
end

TRFVY=zeros(NIY+NMY+1,NIX+NMX+2);
for i=2:NIX+1
    for j=2:NIY
        %VY(1,1:NIX+2)=0.0 --> TRFVY(1,1:NIX+2)=0.0 (top boundary condition in ingot)
        %VY(NIY+1,1:NIX+2)=0.0 --> TRFVY(NIY+1,1:NIX+2)=0.0 (bottom boundary condition in ingot)
        %VY(NIY+2:NIY+NMY+1,1:NIX+2)=0.0 --> TRFVY(NIY+2:NIY+NMY+1,1:NIX+2)=0.0 (in bootom  mould)
        %VY(1:NIY+NMY+1,NIX+2:NIX+NMX+2)=0.0 --> TRFVY(1:NIY+NMY+1,NIX+2:NIX+NMX+2)=0.0 (in right mould)
        TRFVY(j,i)=TTM(j,i)*max(RFVY(j,i),0.0)+TTM(j+1,i)*min(RFVY(j,i),0.0);%[kg.K/m^2/sec]
    end
end
for i=1:NIY+1
    TRFVY(i,1)=TRFVY(i,2);%FREE (left boundary in ingot)
end

%----------------------- Part Two: TRFVX TRFVY ----------------------------

%----------------------- Part Three: Heat flux ----------------------------
TQVX=zeros(NIY+NMY+2,NIX+NMX);%T=T, Q=flux, VX=VX
for i=1:NIX
    %TRFVX(1:NIY+NMY+2,NIX+1:NIX+NMX+1)=0.0 --> TQVX(1:NIY+NMY+2,NIX+1:NIX+NMX)=0.0 (right part in mould)
    %TRFVX(NIY+2:NIY+NMY+2,1:NIX+NMX+1)=0.0 --> TQVX(NIY+2:NIY+NMY+2,1:NIX+NMX)=0.0 (bottom part in mould)
    for j=2:NIY+1
        %TRFVX(1,1:NIX+1)=0.0 --> TQVX(1,1:NIX+NMX)=0.0 (top boundary in ingot)
        %TQVX(j,i)=-dtb*0.5*(CPLTM(j+1,i+1)+CPLTM(j+1,i+2))*TRFVX(j,i+1)/dx(i)+dtb*0.5*(CPLTM(j+1,i+1)+CPLTM(j+1,i))*TRFVX(j,i)/dx(i);%[J/m^3] --> This CPLTM average is not necessary anymore!
        TQVX(j,i)=-dtb*CPLTM(j,i+1)*(TRFVX(j,i+1)-TRFVX(j,i))/dx(i);%[J/m^3]
    end
end

TQVY=zeros(NIY+NMY,NIX+NMX+2);%T=T, Q=flux, VY=VY
for i=1:NIX+1
    for j=1:NIY
        %TRFVY(NIY+1:NIY+1+NMY,1:NIX+1)=0.0 --> TQVY(NIY+1:NIY+NMY,1:NIX+1)=0.0 (bottom part in mould)
        %TRFVY(1:NIY+NMY+1,NIX+2:NIX+NMX+2)=0.0 --> TQVY(1:NIY+NMY,NIX+2:NIX+NMX+2)=0.0 (right part in mould)
        %TQVY(j,i)=-dtb*0.5*(CPLTM(j+2,i+1)+CPLTM(j+1,i+1))*TRFVY(j+1,i)/dy(j)+dtb*0.5*(CPLTM(j+1,i+1)+CPLTM(j,i+1))*TRFVY(j,i)/dy(j);%[J/m^3] --> This CPLTM average is not necessary anymore!
        TQVY(j,i)=-dtb*CPLTM(j+1,i)*(TRFVY(j+1,i)-TRFVY(j,i))/dy(j);%[J/m^3]
        
    end
end

TQVT=zeros(NIY+NMY,NIX+NMX);%T=T, Q=flux, T=total
TQVT(1:NIY,1:NIX)=TQVX(2:NIY+1,1:NIX)+TQVY(1:NIY,2:NIX+1);%only in ingot [J/m^3]

%------------------------ Part Three: Heat flux ---------------------------

%################################################################### ENERGY BALANCE #########################################################################


%################################################################## SPECIES BALANCE #########################################################################

%% .............========= DIFFUSION SPECIES FLUX =========.................

%x-axis FL*RL*CL
FRCLX=zeros(NIY+2,NIX+2);
for i=1:NIX+2
    for j=1:NIY+2
        FRCLX(j,i)=FLTM(j,i)*RLTM(j,i)*CLTM(j,i);%[kg/m^3]
    end
end

%y-axis FL*RL*CL
FRCLY=zeros(NIY+2,NIX+2);
for i=1:NIX+2
    for j=1:NIY+2
        FRCLY(j,i)=FLTM(j,i)*RLTM(j,i)*CLTM(j,i);%[kg/m^3]
    end
end

%x-axis diffusion
CLQDX=zeros(NIY+2,NIX+1);%CL=CL, Q=flux, D=diffusion, X=x-axis
for i=2:NIX
    %CL(1:NIY+2,1)=CL(1:NIY+2,2) --> CLQDX(1:NIY+2,1)=0.0
    %CL(1:NIY+2,NIX+1)=CL(1:NIY+2,NIX+2) --> CLQDX(1:NIY+2,NIX+1)=0.0
    for j=1:NIY+2
        CLQDX(j,i)=2.0*dtb*DL*(FRCLX(j,i+1)-FRCLX(j,i))/((dx(i-1)+dx(i))*dx(i-1));%[kg/m^3]
    end
end

%y-axis diffusion
CLQDY=zeros(NIY+1,NIX+2);%CL=CL, Q=flux, D=diffusion, Y=y-axis
for i=1:NIX+2
    for j=2:NIY
        %CL(1,1:NIX+2)=CL(2,1:NIX+2) --> CLQDY(1,1:NIX+2)=0.0
        %CL(NIY+2,1:NIX+2)=CL(NIY+1,1:NIX+2) --> CLQDY(NIY+1,1:NIX+2)=0.0
        CLQDY(j,i)=2.0*dtb*DL*(FRCLY(j+1,i)-FRCLY(j,i))/((dy(j-1)+dy(j))*dy(j-1));%[kg/m^3]
    end
end

%total species flux of each finite volume
CLQDT=zeros(NIY,NIX);%CL=CL, Q=flux, D=diffusion, T=total
for i=1:NIX
    for j=1:NIY
        CLQDT(j,i)=CLQDX(j+1,i+1)-CLQDX(j+1,i)+CLQDY(j+1,i+1)-CLQDY(j,i+1);%[kg/m^3]
    end
end

%% .............========= INTERNAL SPECIES STORAGE =========...............

KP=zeros(NIY+2,NIX+2);%solid/liquid partition coefficient
for i=1:NIX+2
    for j=1:NIY+2
        if(CLTM(j,i)~=CE)
            KP(j,i)=0.12824+5.699124e-5*CLTM(j,i)*100.0+3.728777e-5*(CLTM(j,i)*100.0)^2;%solid/liquid partition coefficient from D-15 [1]
        else
            KP(j,i)=1.0;
        end
    end
end

ISS=zeros(NIY+2,NIX+2);%internal species storage
FRLTM=zeros(NIY+2,NIX+2);%old time FL*RL
RSKPTM=zeros(NIY+2,NIX+2);%old time RS*KP
dFSTM=zeros(NIY+2,NIX+2);%temporary dFS
phi=0.5;%integration coefficient (0,1)

for i=1:NIX+2
    for j=1:NIY+2
        FRLTM(j,i)=FLTM(j,i)*RLTM(j,i)*CLTM(j,i);%[kg/m^3]
        RSKPTM(j,i)=RSTM(j,i)*KP(j,i)*CLTM(j,i);%[kg/m^3]
        ISS(j,i)=FRLTM(j,i)-(1.0-phi)*dFSTM(j,i)*RSKPTM(j,i);%[kg/m^3]
    end
end

%% ...........========== CONVECTION SPECIES FLUX ===========...............

%---------------------- Part One: CLRFVX CLRFVY ---------------------------

%x-axis species convection flux
CLRFVX=zeros(NIY+2,NIX+1);
for i=2:NIX
    %VX(1:NIY+2,1)=0.0 --> CLRFVX(1:NIY+2,1)=0.0
    %VX(1:NIY+2,NIX+1)=0.0 --> CLRFVX(1:NIY+2,NIX+1)=0.0
    for j=1:NIY+2
        %VX(1,1:NIX+1)=0.0 --> CLRFVX(1,1:NIX+1)=0.0
        %VX(NIY+2,1:NIX+1)=0.0 --> CLRFVX(NIY+2,1:NIX+1)=0.0
        CLRFVX(j,i)=CLTM(j,i)*max(RFVX(j,i),0.0)+CLTM(j,i+1)*min(RFVX(j,i),0.0);%[kg/m^2/sec]
    end
end

%y-axis species convection flux
CLRFVY=zeros(NIY+1,NIX+2);
for i=2:NIX+1
    %VY(1:NIY+1,1)=FREE --> CLRFVY(1:NIY+1,1)==CLRFVY(1:NIY+1,2)
    %VY(1:NIY+1,NIX+2)=0.0 --> CLRFVY(1:NIY+1,NIX+2)=0.0
    for j=2:NIY
        %VY(1,1:NIX+2)=0.0 --> CLRFVY(1,1:NIX+2)=0.0
        %VY(NIY+1,1:NIX+2)=0.0 --> CLRFVY(NIY+1,1:NIX+2)=0.0
        CLRFVY(j,i)=CLTM(j,i)*max(RFVY(j,i),0.0)+CLTM(j+1,i)*min(RFVY(j,i),0.0);%[kg/m^2/sec]
    end
end

for i=1:NIY+1
    CLRFVY(i,1)=CLRFVY(i,2);%left slip boundary
end

%---------------------- Part One: CLRFVX CLRFVY ---------------------------

%------------------------ Part Two: Species flux --------------------------

%x-axis species convection flux
CLQVX=zeros(NIY+2,NIX);%CL=CL, Q=flux, VX=VX
for i=1:NIX
    for j=1:NIY+2
        CLQVX(j,i)=-dtb*(CLRFVX(j,i+1)-CLRFVX(j,i))/dx(i);%[kg/m^3]
    end
end

%y-axis species convection flux
CLQVY=zeros(NIY,NIX+2);%CL=CL, Q=flux, VY=VY
for i=1:NIX+2
    for j=1:NIY
        CLQVY(j,i)=-dtb*(CLRFVY(j+1,i)-CLRFVY(j,i))/dy(j);%[kg/m^3]
    end
end

%total species convection flux
CLQVT=zeros(NIY,NIX);%CL=CL, Q=flux, V=velocity, T=total
CLQVT=CLQVX(2:NIY+1,1:NIX)+CLQVY(1:NIY,2:NIX+1);%[kg/m^3]

%------------------------ Part Two: Species flux --------------------------

%################################################################## SPECIES BALANCE #########################################################################

%################################################################ T-FS-CL ITERATION #########################################################################
errdFS=ones(NIY+1,NIX+1);%initial error of dFS, to triger while loop
errTTM=ones(NIY+NMY+2,NIX+NMX+2);
errCLTM=ones(NIY,NIX);

errdFSmax=max(max(errdFS));
errTTMmax=max(max(errTTM));
errCLTMmax=max(max(errCLTM));

TFN=zeros(NIY+NMY+2,NIX+NMX+2);%T Factor of Next step
TN=zeros(NIY+NMY+2,NIX+NMX+2);%T of Next step
CLFN=zeros(NIY,NIX);%CL Factor of Next step
CLN=zeros(NIY+2,NIX+2);%CL of Next step
CSN=zeros(NIY+2,NIX+2);%CS of Next step
FLN=zeros(NIY+2,NIX+2);%FL of Next step
RLN=zeros(NIY+2,NIX+2);%RL of Next step
RSN=zeros(NIY+2,NIX+2);%RS of Next step
KPTM=zeros(NIY+2,NIX+2);%KP of Next step
KPN=zeros(NIY+2,NIX+2);%KP of Next step
TLN=zeros(NIY+2,NIX+2);%Liquidus of Next step
HLTM=zeros(NIY+2,NIX+2);%latent heat of Next step
HLN=zeros(NIY+2,NIX+2);%latent heat of Next step
dFSN=zeros(NIY+2,NIX+2);%dFS of Next step
CPSN=zeros(NIY+2,NIX+2);%CPS of Next step
CPLN=zeros(NIY+2,NIX+2);%CPL of Next step

%TFN: T Factor of Next step
%=== (1)ingot ===
for i=1:NIX+1
    for j=1:NIY+1
        TFN(j,i)=FRCPSTM(j,i)+FLTM(j,i)*RLTM(j,i)*CPLTM(j,i)+(RSCPSTM(j,i)-RLCPLTM(j,i))*dFSTM(j,i);%[J/m^3/K]
    end
end

%=== (2)mould ===
for i=NIX+2:NIX+NMX+2
    for j=1:NIY+NMY+2
        TFN(j,i)=RCPM(j,i);%[J/m^3/K]
    end
end

for i=1:NIX+1
    for j=NIY+2:NIY+NMY+2
        TFN(j,i)=RCPM(j,i);%[J/m^3/K]
    end
end

%----------------------------- STEP ONE -----------------------------------
%Let the initial value of dFSTM equal zero and calculate the approximation of TN using Eq. (8).
for j=2:NIX+NMX+1
    for k=2:NIY+NMY+1
        TN(k,j)=(IHS(k,j)+TQVT(k-1,j-1)+TQDT(k-1,j-1))/TFN(k,j);%[J/m^3 divided by J/m^3/K == K]
    end
end
%RSCPSTM: RS*CPS at n step, as approxiamation of RS*CPS at n+1 step
%RLCPLTM: RS*CPL at n step, as approxiamation of RS*CPL at n+1 step
TN(1,1:NIX+NMX+2)=TN(2,1:NIX+NMX+2);%abdiabatic top boundary
TN(1:NIY+NMY+2,1)=TN(1:NIY+NMY+2,2);%abdiabatic left boundary
TN(1:NIY+NMY+2,NIX+NMX+2)=TTM(1:NIY+NMY+2,NIX+NMX+2);%right mould boundary, set from BoundTemp.m
TN(NIY+NMY+2,1:NIX+NMX+2)=TTM(NIY+NMY+2,1:NIX+NMX+2);%bottom mould boundary, set from BoumdTemp.m
%TTM=TN;

%----------------------------- STEP ONE -----------------------------------

%----------------------------- STEP TWO -----------------------------------
%Take the approximations TN and CLTM(n) as the initial values of temperature and
%liquid composition at time n+1, and calculate the approximation of CLTM(n+1) using Eq. (9) and the approximate liquidus
%temperature at n+1 from Eq. (3).

%FL of Next step
%FLN=FLTM-dFSTM;

%Partition coefficient of old step
for i=1:NIX+2
    for j=1:NIY+2
        if(CLTM(j,i)~=CE)
            KPTM(j,i)=0.12824+5.699124e-5*CLTM(j,i)*100.0+3.728777e-5*(CLTM(j,i)*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
        else
            KPTM(j,i)=1.0;
        end
    end
end

%CL of Next step
for j=1:NIX
    for k=1:NIY
        %CLFN--Factor ahead CL(n+1): FL, RL, dFS, RS, KP should use NEW values,
        %but for the first iteration, set old values as new values
        CLFN(k,j)=FLTM(k+1,j+1)*RLTM(k+1,j+1)+phi*dFSTM(k+1,j+1)*RSTM(k+1,j+1)*KPTM(k+1,j+1);%[kg/m^3]
        CLN(k+1,j+1)=(ISS(k+1,j+1)+CLQVT(k,j)+CLQDT(k,j))/CLFN(k,j);%[kg/m^3 divided by kg/m^3 ==1]; ISS should update with NEW dFS ,but no NEW dFS available now, set old dFSTM==0.0 as new one
        errCLTM(k,j)=abs((CLN(k+1,j+1)-CLTM(k+1,j+1))/CLN(k+1,j+1));
    end
end
CLN(1,2:NIX+1)=CLN(2,2:NIX+1);%top boundary
CLN(NIY+2,2:NIX+1)=CLN(NIY+1,2:NIX+1);%bottom boundary
CLN(2:NIY+1,1)=CLN(2:NIY+1,2);%left boundary
CLN(2:NIY+1,NIX+2)=CLN(2:NIY+1,NIX+1);%right boundary
CLN(1,1)=CLN(1,2);
CLN(1,NIX+2)=CLN(1,NIX+1);
CLN(NIY+2,1)=CLN(NIY+2,2);
CLN(NIY+2,NIX+2)=CLN(NIY+2,NIX+1);
%CLTM=CLN;

%Liquidus of Next step
TLN=660.37-2.34581*CLN*100.0-3.129e-2*(CLN*100.0).^2+273.15;%ingot liquidus [K]; CLN here not in wt% while in wt% from D-16

% %Partition coefficient of Next step
% for i=1:NIX+2
%     for j=1:NIY+2
%         if(CLTM(j,i)~=CE)
%             KPN(j,i)=0.12824+5.699124e-5*CLTM(j,i)*100.0+3.728777e-5*(CLTM(j,i)*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
%         else
%             KPN(j,i)=1.0;
%         end
%     end
% end

%----------------------------- STEP TWO -----------------------------------

%logical variable for while loop
err=logical((err0dFS<=errdFSmax)||(err0TTM<=errTTMmax)||(err0CLTM<=errCLTMmax));

%main loop for T-FS-CL
iters=0;
while(err)
    %-------------------------- STEP THREE --------------------------------
    %Calculate the approximation of dFSN by Eq. (10) and take the average:
    %dFSN=0.5*(dFSN+dFSTM)
    
    for i=1:NIX+2
        for j=1:NIY+2
            if(CSTM(j,i)<0.1)
                HLTM(j,i)=397.67-2.3288*CSTM(j,i)*100.0;%latent heat from D-7 [J/g]
            else
                HLTM(j,i)=333.59;%[J/g]
            end
        end
    end
    HLTM=HLTM*1000.0;%[J/kg]
    
    %RL of Next step
    % RLN=2.5222e-3+2.703e-5*CLTM*100.0-3.16e-7*(TTM(1:NIY+2,1:NIX+2)-273.15);%liquid density of Next step from D-6 [g/mm^3]
    % RLN=RLN*10^6;%[kg/m^3]
    %
    % %RS of Next step
    % for i=1:NIX+2
    %     for j=1:NIY+2
    %         if(CSN(j,i)<CE)
    %             RSN(j,i)=2.58e-3;%[g/mm^3]
    %         else
    %             RSN(j,i)=3.4e-3;%[g/mm^3]
    %         end
    %     end
    % end
    % RSN=RSN*10^6;%[kg/m^3]
    
    
    for j=1:NIX+1
        for k=1:NIY+1
            dFSN(k,j)=(FLTM(k,j)*RLTM(k,j)*CLTM(k,j)+FRCSTM(k,j))*(TLN(k,j)-TN(k,j))/(RSTM(k,j)*HLTM(k,j));%[1]
        end
    end
    %NOTE: when TL>T, ingot is super liquidus and no solid forms and thus
    %dFSN<0 (unphysical), so we correct it as following:
    for i=1:NIX+2
        for j=1:NIY+2
            if(dFSN(j,i)<0.0)
                dFSN(j,i)=0.0;
            end
        end
    end
    
    dFSN=0.5*(dFSN+dFSTM);
    
    for j=1:NIX+1
        for k=1:NIY+1
            errdFS(k,j)=abs((dFSN(k,j)-dFSTM(k,j))/dFSN(k,j));
            if(isnan(errdFS(k,j)))
                %errdFS==NaN --> dFSN==0.0 --> no solid forms
                errdFS(k,j)=0.0;
            end
        end
    end
    
    %-------------------------- STEP THREE --------------------------------
    
    %--------------------------- STEP FOUR --------------------------------
    %Again using Eq. (8), calculate the new approximation of T(n+1) with new initial
    %values, dFSN and CLN, and take TN=0.5*(TN+TTM).
    
    % %Partition coefficient of Next step
    for i=1:NIX+2
        for j=1:NIY+2
            if(CLN(j,i)~=CE)
                KPN(j,i)=0.12824+5.699124e-5*CLN(j,i)*100.0+3.728777e-5*(CLN(j,i)*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
            else
                KPN(j,i)=1.0;
            end
        end
    end
    
    %RS of Next step
    for i=1:NIX+2
        for j=1:NIY+2
            CSN(j,i)=CLN(j,i)*KPN(j,i);
            if(CSN(j,i)<CE)
                RSN(j,i)=2.58e-3;%[g/mm^3]
            else
                RSN(j,i)=3.4e-3;%[g/mm^3]
            end
        end
    end
    RSN=RSN*10^6;%[kg/m^3]
    
    %CpS of Next step
    for i=1:NIX+2
        for j=1:NIY+2
            CPSN(j,i)=0.88+4.446e-4*(TN(j,i)-273.15)-2.274e-3*CSN(j,i)*100.0;%solid specific heat capacity from D-3 [J/g/K]
        end
    end
    CPSN=CPSN*1000.0;%[J/kg/K]
    
    %RL of Next step
    RLN=2.5222e-3+2.703e-5*CLN*100.0-3.16e-7*(TN(1:NIY+2,1:NIX+2)-273.15);%liquid density of Next step from D-6 [g/mm^3]
    RLN=RLN*10^6;%[kg/m^3]
    
    %CpL of Next step
    for i=1:NIX+2
        for j=1:NIY+2
            CPLN(j,i)=1.086-5.928e-3*CLN(j,i)*100.0;%liquid specific heat capacity from D-4 [J/g/K]
        end
    end
    CPLN=CPLN*1000.0;%[J/kg/K]
    
    for i=1:NIX+2
        for j=1:NIY+2
            RSCPSTM(j,i)=RSN(j,i)*CPSN(j,i);%[J/m^3/K]
            RLCPLTM(j,i)=RLN(j,i)*CPLN(j,i);%[J/m^3/K]
        end
    end
    
    for j=1:NIX+1
        for k=1:NIY+1
            TFN(k,j)=FRCPSTM(k,j)+FLTM(k,j)*RLTM(k,j)*CPLTM(k,j)+(RSCPSTM(k,j)-RLCPLTM(k,j))*dFSN(k,j);%[J/m^3/K]
        end
    end
    
    
    for j=2:NIX+NMX+1
        %TN(1:NIY+NMY+2,NIX+NMX+2)=determined right mould suface
        %TN(1:NIY+NMY+2,1)=adiabatic left ingot
        for k=2:NIY+NMY+1
            %TN(NIY+NMY+2,1:NIX+NMX+2)=determined bottom mould surface
            %TN(1,1:NIX+NMX+2)=adiabatic top ingot
            TN(k,j)=(IHS(k,j)+TQVT(k-1,j-1)+TQDT(k-1,j-1))/TFN(k,j);%[J/m^3 divided by J/m^3/K == K]
        end
    end
    TN=0.5*(TN+TTM);
    TN(1,1:NIX+NMX+2)=TN(2,1:NIX+NMX+2);%abdiabatic top boundary
    TN(1:NIY+NMY+2,1)=TN(1:NIY+NMY+2,2);%abdiabatic left boundary

    for j=2:NIX+NMX+1
        for k=2:NIY+NMY+1
            errTTM(k,j)=abs((TN(k,j)-TTM(k,j))/TN(k,j));
        end
    end
    
    %--------------------------- STEP FOUR --------------------------------
    
    %--------------------------- STEP FIVE --------------------------------
    %With the new initial value TN (as well as dFSN and CLN), calculate the new approximation
    %of CLN(n+1) using Eq. (9), and take CLN=0.5*(CLN+CLTM).
    
    %FL of of Next step
    FLN=FLTM-dFSN;
    
    %RL of Next step
    RLN=2.5222e-3+2.703e-5*CLN*100.0-3.16e-7*(TN(1:NIY+2,1:NIX+2)-273.15);%liquid density of Next step from D-6 [g/mm^3]
    RLN=RLN*10^6;%[kg/m^3]
    
    %RS: CLN is updated once, KPN is also updated from CLN, thus CSN keeps
    %constant --> RSN does not change
    %         for i=1:NIX
    %             for j=1:NIY
    %                 if(CSN(j,i)<CE)
    %                     RSN(j,i)=2.58e-3;%[g/mm^3]
    %                 else
    %                     RSN(j,i)=3.4e-3;%[g/mm^3]
    %                 end
    %             end
    %         end
    %         RSN=RSN*10^6;%[kg/m^3]
    
    %KP of Next step
    %     for i=1:NIX+2
    %         for j=1:NIY+2
    %             if(CLTM(j,i)~=CE)
    %                 KPN(j,i)=0.12824+5.699124e-5*CLTM(j,i)*100.0+3.728777e-5*(CLTM(j,i)*100.0)^2;%solid/liquid partition coefficient from D-15 [1]
    %             else
    %                 KPN(j,i)=1.0;
    %             end
    %         end
    %     end
    
    %New ISS due to updated dFSN
    for i=1:NIX+2
        for j=1:NIY+2
            ISS(j,i)=FRLTM(j,i)-(1.0-phi)*dFSN(j,i)*RSKPTM(j,i);%[kg/m^3]
        end
    end
    
    
    %CL of Next step
    for j=1:NIX
        for k=1:NIY
            %Factor ahead CL(n+1)
            CLFN(k,j)=FLN(k+1,j+1)*RLN(k+1,j+1)+phi*dFSN(k+1,j+1)*RSN(k+1,j+1)*KPN(k+1,j+1);%[kg/m^3]
            CLN(k+1,j+1)=(ISS(k+1,j+1)+CLQVT(k,j)+CLQDT(k,j))/CLFN(k,j);%[kg/m^3 divided by kg/m^3 ==1]
        end
    end
    
    CLN=0.5*(CLN+CLTM);
    
    for j=1:NIX
        for k=1:NIY
            errCLTM(k,j)=abs((CLN(k+1,j+1)-CLTM(k+1,j+1))/CLN(k+1,j+1));
        end
    end
    %--------------------------- STEP FIVE --------------------------------
    
    %--------------------------- STEP SIX ---------------------------------
    %Judge if all of the following conditions are true:
    errdFSmax=max(max(errdFS(2:NIY+1,2:NIX+1)));
    errTTMmax=max(max(errTTM(2:NIY+NMY+1,2:NIX+NMX+1)));
    errCLTMmax=max(max(errCLTM));
    %NOTE: errdFS=ones(NIY+1,NIX+1) errTTM=ones(NIY+NMY+2,NIX+NMX+2) errCLTM=ones(NIY,NIX);
    
    
    err=logical((err0dFS<=errdFSmax)||(err0TTM<=errTTMmax)||(err0CLTM<=errCLTMmax));%logical variable for while loop
    
    %If the conditions are not true, take TLN(n+1)=TLN(CLN(n+1)) and repeat steps 3-6 until
    %all the above conditions are satisfied
    %Liquidus of Next step
    
    %--------------------------- STEP SIX ---------------------------------
    
    %--------------------------- STEP SEVEN -------------------------------
    %If jugments are not true, then repeat step 3-6 with updated TN, CLN,
    %dFSN
    
    %Partition coefficient to update CSTM with CLN
    for i=1:NIX+2
        for j=1:NIY+2
            if(CLN(j,i)~=CE)
                KPN(j,i)=0.12824+5.699124e-5*CLN(j,i)*100.0+3.728777e-5*(CLN(j,i)*100.0)^2;%solid/liquid partition coefficient from D-15 [1]
            else
                KPN(j,i)=1.0;
            end
            CSTM(j,i)=CLN(j,i)*KPN(j,i);
        end
    end
    
    %Update FL
    FLTM=FLN;
    
    %Update RL
    RLTM=2.5222e-3+2.703e-5*CLN*100.0-3.16e-7*(TN(1:NIY+2,1:NIX+2)-273.15);%liquid density of Next step from D-6 [g/mm^3]
    RLTM=RLTM*10^6;%[kg/m^3]
    
    %Update TL
    TLN=660.37-2.34581*CLN*100.0-3.129e-2*(CLN*100.0).^2+273.15;%ingot liquidus from D-16 [K]; CLN here not in wt% while in wt% from D-16
    
    %Update RSTM from CSTM
    for i=1:NIX+2
        for j=1:NIY+2
            if(CSTM(j,i)<CE)
                RSTM(j,i)=2.58e-3;%[g/mm^3]
            else
                RSTM(j,i)=3.4e-3;%[g/mm^3]
            end
        end
    end
    RSTM=RSTM*10^6;%[kg/m^3]
    
    %Update dFSTM
    dFSTM=dFSN;
    
    %Update CPLTM from CLN
    for i=1:NIX+2
        for j=1:NIY+2
            CPLTM(j,i)=1.086-5.928e-3*CLN(j,i)*100.0;%liquid specific heat capacity from D-4 [J/g/K]
        end
    end
    CPLTM=CPLTM*1000.0;%[J/kg/K]
    
    %Update TTM
    TTM=TN;
    
    %Update CL
    CLTM=CLN;
    
    %-------------------------- STEP SEVEN --------------------------------
    
    iters=iters+1;
    fprintf('Iteration:%5d\n',iters);
    
end
%################################################################ T-FS-CL ITERATION #########################################################################

T=TN;
FL=FLN;
CL=CLN;
dFS=dFSN;
end

