function [T,FL,CL,dFS]=TFCPV(TTM,FLTM,FSTM,KLTM,KSTM,CLTM,CSTM,VXTM,VYTM,FRCPSTM,FRCSTM,FKTM,RSCPSTM,RLCPLTM,RSTM,RLTM,CPLTM)%FRCSTM,FRHTM,
%To iteratively solve energy and species conservation equations with finite
%volume method. Model, discretization and iteration can be found "A
%Continuum Model for Computer Simulation of Macrosegregations in Ingots
%during solidification" by Daming Xu 1989 and his companion papers in 1991

%Created 2019-10-17

%Modified for T-FS-CL iteration verification! 2019-11-10

%Modified in P-V iteration 2019-11-24; only ingot are useful

%TTM: [NIY+NMY+2,NIX+NMX+2]
%FLTM: [NIY+2,NIX+2]
%FSTM: [NIY+2,NIX+2]
%KLTM: [NIY+2,NIX+2]
%KSTM: [NIY+2,NIX+2]

%FRCPSTM: temporary FRCPS=FS*RS*CPS, of [NIY+2,NIX+2]
%FRCSTM: temporary FRCS=FS*RS*CS, of [NIY+2,NIX+2]
%FKTM: temporary FK=FS*KS, of [NIY+2,NIX+2]
%RSCPSTM: temporary RSCPS=RS*CPS, of [NIY+2,NIX+2]
%RLCPLTM: temporary RLCPL=RL*CPL, of [NIY+2,NIX+2]

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
global dtb

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
[TVB,THB]= BoundT(TTM(:,NIX+NMX+1),TTM(NIY+NMY+1,:),KA(:,NIX+NMX+1),KA(NIY+NMY+1,:));%calculate boundary temperature via conduction==radiation + convection
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

%right most boundary of mould
for i=2:NIY+NMY+1
    %Situation 1: Eq. (3-31), (D-12) with DT=T(i,NIX+NMX+1)-TR, but DT should be TS-TR (TS=surface T)
    TQDX(i,NIX+NMX+1)=-FVB(i-1);
    
    %Situation 2: solve TS first then calculate flux
    %     if(TTM(i,NIX+NMX+2)<601.61)
    %         KBM=7.9718e-2-1.0387e-4*(TTM(i,NIX+NMX+2)-273.15);%temporary boundary thermal conductivity [W/mm/K]
    %     else
    %         KBM=5.43102e-2-2.6516e-5*(TTM(i,NIX+NMX+2)-273.15);%temporary boundary thermal conductivity [W/mm/K]
    %     end
    %
    %     KBM=KBM*1000.0;%[W/m/K]
    %     TQDX(i,NIX+NMX+1)=0.5*(KX(i,NIX+NMX+1)+KBM)*DTX(i,NIX+NMX+1);
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

%bottom boundary of mould
for i=2:NIX+NMX+1
    %Situation 1: Eq. (3-31), (D-14) with DT=T(NIY+NMY+1,i)-TR, but DT should be TS-TR (TS=surface T)
    TQDY(NIY+NMY+1,i)=-FHB(i-1);
    
    %Situation 2: solve TS first then calculate flux
    %     if(TTM(NIY+NMY+2,i)<601.61)
    %         KBM=7.9718e-2-1.0387e-4*(TTM(NIY+NMY+2,i)-273.15);%temporary boundary thermal conductivity [W/mm/K]
    %     else
    %         KBM=5.43102e-2-2.6516e-5*(TTM(NIY+NMY+2,i)-273.15);%temporary boundary thermal conductivity [W/mm/K]
    %     end
    %
    %     KBM=KBM*1000.0;%[W/m/K]
    %     TQDY(NIY+NMY+1,i)=0.5*(KY(NIY+NMY+1,i)+KBM)*DTY(NIY+NMY+1,i);
    
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
        %TRFVX(1,1:NIX+1)=0.0 --> TQVX(1,1:NIX)=0.0 (top boundary in ingot)
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
errdFS=ones(NIY+2,NIX+2);%initial error of dFS, to triger while loop
errTTM=ones(NIY+2,NIX+2);
errCLTM=ones(NIY+2,NIX+2);
err=zeros(NIY+2,NIX+2);
% errdFSmax=max(max(errdFS));
% errTTMmax=max(max(errTTM));
% errCLTMmax=max(max(errCLTM));

TFN=zeros(NIY+NMY+2,NIX+NMX+2);%T Factor of Next step
TN=zeros(NIY+NMY+2,NIX+NMX+2);%T of Next step
TLTM=zeros(NIY+2,NIX+2);%Liquidus of old step
HLTM=zeros(NIY+2,NIX+2);%latent heat of Next step

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
TN(1:NIY+NMY+2,NIX+NMX+2)=TN(1:NIY+NMY+2,NIX+NMX+1);%right mould boundary
TN(NIY+NMY+2,1:NIX+NMX+2)=TN(NIY+NMY+1,1:NIX+NMX+2);%bottom mould boundary
TTM=TN;

%----------------------------- STEP ONE -----------------------------------
TLTM=660.37-2.34581*CLTM*100.0-3.129e-2*(CLTM*100.0).^2+273.15;%old liquidus to determine superliquidus or subliquidus [K]
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

modifymarker=0;

%logical variable for while loop
for i=1:NIX+2
    for j=1:NIY+2
        err(j,i)=logical((err0dFS<=errdFS(j,i))||(err0TTM<=errTTM(j,i))||(err0CLTM<=errCLTM(j,i)));
    end
end

%main loop for T-FS-CL
            CLFN=zeros(NIY,NIX);%CL Factor of Next step
            CLN=zeros(NIY+2,NIX+2);%CL of Next step
            CSN=zeros(NIY+2,NIX+2);%CS of Next step
            FLN=zeros(NIY+2,NIX+2);%FL of Next step
            RLN=zeros(NIY+2,NIX+2);%RL of Next step
            RSN=zeros(NIY+2,NIX+2);%RS of Next step
            KPN=zeros(NIY+2,NIX+2);%KP of Next step
            TLN=zeros(NIY+2,NIX+2);%Liquidus of Next step
            dFSN=zeros(NIY+2,NIX+2);%dFS of Next step
            CPSN=zeros(NIY+2,NIX+2);%CPS of Next step
            CPLN=zeros(NIY+2,NIX+2);%CPL of Next step

for i=1:NIX
    for j=1:NIY
        %SITUATION 1: SUPERLIQUIDUS
        if(TN(j+1,i+1)>=TLTM(j+1,i+1))
            dFSN(j+1,i+1)=dFSTM(j+1,i+1);%For superliquidus, we do nothing!
            FLN(j+1,i+1)=FLTM(j+1,i+1);
            CLN(j+1,i+1)=CLTM(j+1,i+1);
            %SITUATION 2: CRYSTALLIZATION
        else
            %Once some cell are subliquidus, we mark ingot T has been
            %solved by T-FS-CL iteration. Even one cell!
            modifymarker=modifymarker+1;
            
            %----------------------------- STEP TWO -----------------------------------
            %Take the approximations TN and CLTM(n) as the initial values of temperature and
            %liquid composition at time n+1, and calculate the approximation of CLTM(n+1) using Eq. (9) and the approximate liquidus
            %temperature at n+1 from Eq. (3).
            
            %FL of Next step
            FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSTM(j+1,i+1);
            
            %NOTE: new T has effect on RL,RS,KP
            %New RL uodated by new T
            RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLTM(j+1,i+1)*100.0-3.16e-7*(TN(j+1,i+1)-273.15);%liquid density of Next step from D-6 [g/mm^3]
            RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
            
            %New KP updated by new T
            
            if(CLTM(j+1,i+1)~=CE)
                KPN(j+1,i+1)=0.12824+5.699124e-5*CLTM(j+1,i+1)*100.0+3.728777e-5*(CLTM(j+1,i+1)*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
            else
                KPN(j+1,i+1)=1.0;
            end
            
            %New RS updated by new T
            
            CSN(j+1,i+1)=CLTM(j+1,i+1)*KPN(j+1,i+1);
            if(CSN(j+1,i+1)<CE)
                RSN(j+1,i+1)=2.58e-3;%[g/mm^3]
            else
                RSN(j+1,i+1)=3.4e-3;%[g/mm^3]
            end
            
            RSN(j+1,i+1)=RSN(j+1,i+1)*10^6;%[kg/m^3]
            
            %CL of Next step
            
            %CLFN--Factor ahead CL(n+1): FL, RL, dFS, RS, KP should use NEW
            %values FLN, RLN, RSN, KPN; dFSTM==dFSN==0.0 for the first step
            CLFN(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1)+phi*dFSTM(j+1,i+1)*RSN(j+1,i+1)*KPN(j+1,i+1);%[kg/m^3]
            CLN(j+1,i+1)=(ISS(j+1,i+1)+CLQVT(j,i)+CLQDT(j,i))/CLFN(j,i);%[kg/m^3 divided by kg/m^3 ==1]; ISS should update with NEW dFS ,but no NEW dFS available now, set old dFSTM==0.0 as new one
            errCLTM(j+1,i+1)=abs((CLN(j+1,i+1)-CLTM(j+1,i+1))/CLN(j+1,i+1));
            
            CLTM(j+1,i+1)=CLN(j+1,i+1);
            
            TLN(j+1,i+1)=660.37-2.34581*CLN(j+1,i+1)*100.0-3.129e-2*(CLN(j+1,i+1)*100.0).^2+273.15;%new liquidus [K]
            
            %----------------------------- STEP TWO -----------------------------------
            while(err(j+1,i+1))
                
                %-------------------------- STEP THREE --------------------------------
                %Calculate the approximation of dFSN by Eq. (10) and take the average:
                %dFSN=0.5*(dFSN+dFSTM)
         
                dFSN(j+1,i+1)=(FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+FRCPSTM(j+1,i+1))*(TLN(j+1,i+1)-TN(j+1,i+1))/(RSTM(j+1,i+1)*HLTM(j+1,i+1));%[1]
                
                dFSN(j+1,i+1)=0.5*(dFSN(j+1,i+1)+dFSTM(j+1,i+1));
                
                errdFS(j+1,i+1)=abs((dFSN(j+1,i+1)-dFSTM(j+1,i+1))/dFSN(j+1,i+1));

                dFSTM(j+1,i+1)=dFSN(j+1,i+1);
                %-------------------------- STEP THREE --------------------------------
                
                %--------------------------- STEP FOUR --------------------------------
                %Again using Eq. (8), calculate the new approximation of T(n+1) with new initial
                %values, dFSN and CLN, and take TN=0.5*(TN+TTM).
                
                %NOTE: new T, CL have effect on RL,RS,KP,CpS,CpL
                %Partition coefficient of Next step
                
                if(CLN(j+1,i+1)~=CE)
                    KPN(j+1,i+1)=0.12824+5.699124e-5*CLN(j+1,i+1)*100.0+3.728777e-5*(CLN(j+1,i+1)*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
                else
                    KPN(j+1,i+1)=1.0;
                end
                
                %RS of Next step
                CSN(j+1,i+1)=CLN(j+1,i+1)*KPN(j+1,i+1);
                if(CSN(j+1,i+1)<CE)
                    RSN(j+1,i+1)=2.58e-3;%[g/mm^3]
                else
                    RSN(j+1,i+1)=3.4e-3;%[g/mm^3]
                end
                
                RSN(j+1,i+1)=RSN(j+1,i+1)*10^6;%[kg/m^3]
                
                %CpS of Next step
                
                CPSN(j+1,i+1)=0.88+4.446e-4*(TN(j+1,i+1)-273.15)-2.274e-3*CSN(j+1,i+1)*100.0;%solid specific heat capacity from D-3 [J/g/K]
                
                CPSN(j+1,i+1)=CPSN(j+1,i+1)*1000.0;%[J/kg/K]
                
                %RL of Next step
                RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLN(j+1,i+1)*100.0-3.16e-7*(TN(j+1,i+1)-273.15);%liquid density of Next step from D-6 [g/mm^3]
                RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
                
                %CpL of Next step
                
                CPLN(j+1,i+1)=1.086-5.928e-3*CLN(j+1,i+1)*100.0;%liquid specific heat capacity from D-4 [J/g/K]
                
                CPLN(j+1,i+1)=CPLN(j+1,i+1)*1000.0;%[J/kg/K]
                
                
                RSCPSTM(j+1,i+1)=RSN(j+1,i+1)*CPSN(j+1,i+1);%[J/m^3/K]
                RLCPLTM(j+1,i+1)=RLN(j+1,i+1)*CPLN(j+1,i+1);%[J/m^3/K]
                
                
                TFN(j+1,i+1)=FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*dFSN(j+1,i+1);%[J/m^3/K]
                
                
                %     for j=2:NIX+NMX+1
                %         %TN(1:NIY+NMY+2,NIX+NMX+2)=determined right mould suface
                %         %TN(1:NIY+NMY+2,1)=adiabatic left ingot
                %         for k=2:NIY+NMY+1
                %             %TN(NIY+NMY+2,1:NIX+NMX+2)=determined bottom mould surface
                %             %TN(1,1:NIX+NMX+2)=adiabatic top ingot
                TN(j+1,i+1)=(IHS(j+1,i+1)+TQVT(j,i)+TQDT(j,i))/TFN(j+1,i+1);%[J/m^3 divided by J/m^3/K == K]
                
                
                TN(j+1,i+1)=0.5*(TN(j+1,i+1)+TTM(j+1,i+1));
                %     TN(1,1:NIX+NMX+2)=TN(2,1:NIX+NMX+2);%abdiabatic top boundary
                %     TN(1:NIY+NMY+2,1)=TN(1:NIY+NMY+2,2);%abdiabatic left boundary
                
                errTTM(j+1,i+1)=abs((TN(j+1,i+1)-TTM(j+1,i+1))/TN(j+1,i+1));
                
                TTM(j+1,i+1)=TN(j+1,i+1);
                %--------------------------- STEP FOUR --------------------------------
                
                %--------------------------- STEP FIVE --------------------------------
                %With the new initial value TN (as well as dFSN and CLN), calculate the new approximation
                %of CLN(n+1) using Eq. (9), and take CLN=0.5*(CLN+CLTM).
                
                %FL of of Next step
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSN(j+1,i+1);
                
                %NOTE: new T have effect on RL,CpS,CpL; RS,KP are functions of CLN
                %which is updated before so remain unchanged!
                
                %RL of Next step
                RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLN(j+1,i+1)*100.0-3.16e-7*(TN(j+1,i+1)-273.15);%liquid density of Next step from D-6 [g/mm^3]
                RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
                
                
                %New ISS due to updated dFSN
                
                ISS(j+1,i+1)=FRLTM(j+1,i+1)-(1.0-phi)*dFSN(j+1,i+1)*RSKPTM(j+1,i+1);%[kg/m^3]
                
                %Factor ahead CL(n+1)
                CLFN(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1)+phi*dFSN(j+1,i+1)*RSN(j+1,i+1)*KPN(j+1,i+1);%[kg/m^3]
                %CL of Next step
                CLN(j+1,i+1)=(ISS(j+1,i+1)+CLQVT(j,i)+CLQDT(j,i))/CLFN(j,i);%[kg/m^3 divided by kg/m^3 ==1]
                
                CLN(j+1,i+1)=0.5*(CLN(j+1,i+1)+CLTM(j+1,i+1));
                
                errCLTM(j+1,i+1)=abs((CLN(j+1,i+1)-CLTM(j+1,i+1))/CLN(j+1,i+1));
                CLTM(j+1,i+1)=CLN(j+1,i+1);
                %--------------------------- STEP FIVE --------------------------------
                
                %--------------------------- STEP SIX ---------------------------------
                
                err(j+1,i+1)=logical((err0dFS<=errdFS(j+1,i+1))||(err0TTM<=errTTM(j+1,i+1))||(err0CLTM<=errCLTM(j+1,i+1)));%logical variable for while loop
                
                %If the conditions are not true, take TLN(n+1)=TLN(CLN(n+1)) and repeat steps 3-6 until
                %all the above conditions are satisfied
                %Liquidus of Next step
                
                %--------------------------- STEP SIX ---------------------------------
                
                %--------------------------- STEP SEVEN -------------------------------
                %If jugments are not true, then repeat step 3-6 with updated TN, CLN,
                %dFSN
                
                TLN(j+1,i+1)=660.37-2.34581*CLN(j+1,i+1)*100.0-3.129e-2*(CLN(j+1,i+1)*100.0).^2+273.15;%ingot liquidus from D-16 [K]; CLN here not in wt% while in wt% from D-16
                
                %-------------------------- STEP SEVEN --------------------------------
                                
            end
        end
    end
end
    %################################################################ T-FS-CL ITERATION #########################################################################
    
    %Ingot T-FS-CL has been iteratively sovled?
    if(modifymarker>0.0)%Yes, ingot changed
        %NOTE: T in mould is the same as T without T-FS-CL iteration, since
        %IHS, TQVT, TQDT and TFN stay as old values
        
        %Second: returen new field variables
        TN(1,1:NIX+NMX+2)=TN(2,1:NIX+NMX+2);%abdiabatic top boundary
        TN(1:NIY+NMY+2,1)=TN(1:NIY+NMY+2,2);%abdiabatic left boundary
        TN(1:NIY+NMY+2,NIX+NMX+2)=TN(1:NIY+NMY+2,NIX+NMX+1);%right mould boundary
        TN(NIY+NMY+2,1:NIX+NMX+2)=TN(NIY+NMY+1,1:NIX+NMX+2);%bottom mould boundary
        T=TN;
        
        FLN(1,2:NIX+1)=FLN(2,2:NIX+1);%top boundary
        FLN(NIY+2,2:NIX+1)=FLN(NIY+1,2:NIX+1);%bottom boundary
        FLN(2:NIY+1,1)=FLN(2:NIY+1,2);%left boundary
        FLN(2:NIY+1,NIX+2)=FLN(2:NIY+1,NIX+1);%right boundary
        FLN(1,1)=FLN(1,2);
        FLN(1,NIX+2)=FLN(1,NIX+1);
        FLN(NIY+2,1)=FLN(NIY+2,2);
        FLN(NIY+2,NIX+2)=FLN(NIY+2,NIX+1);
        FL=FLN;
        
        CLN(1,2:NIX+1)=CLN(2,2:NIX+1);%top boundary
        CLN(NIY+2,2:NIX+1)=CLN(NIY+1,2:NIX+1);%bottom boundary
        CLN(2:NIY+1,1)=CLN(2:NIY+1,2);%left boundary
        CLN(2:NIY+1,NIX+2)=CLN(2:NIY+1,NIX+1);%right boundary
        CLN(1,1)=CLN(1,2);
        CLN(1,NIX+2)=CLN(1,NIX+1);
        CLN(NIY+2,1)=CLN(NIY+2,2);
        CLN(NIY+2,NIX+2)=CLN(NIY+2,NIX+1);
        CL=CLN;
        
        dFSN(1,2:NIX+1)=dFSN(2,2:NIX+1);%top boundary
        dFSN(NIY+2,2:NIX+1)=dFSN(NIY+1,2:NIX+1);%bottom boundary
        dFSN(2:NIY+1,1)=dFSN(2:NIY+1,2);%left boundary
        dFSN(2:NIY+1,NIX+2)=dFSN(2:NIY+1,NIX+1);%right boundary
        dFSN(1,1)=dFSN(1,2);
        dFSN(1,NIX+2)=dFSN(1,NIX+1);
        dFSN(NIY+2,1)=dFSN(NIY+2,2);
        dFSN(NIY+2,NIX+2)=dFSN(NIY+2,NIX+1);
        dFS=dFSN;
        
    else%No, only ingot T changed; FS, CL unchange! return updated T and constant FL,CL
        
        TN(1,1:NIX+NMX+2)=TN(2,1:NIX+NMX+2);%abdiabatic top boundary
        TN(1:NIY+NMY+2,1)=TN(1:NIY+NMY+2,2);%abdiabatic left boundary
        TN(1:NIY+NMY+2,NIX+NMX+2)=TN(1:NIY+NMY+2,NIX+NMX+1);%right mould boundary
        TN(NIY+NMY+2,1:NIX+NMX+2)=TN(NIY+NMY+1,1:NIX+NMX+2);%bottom mould boundary
        T=TN;
        
        FL=FLTM;
        CL=CLTM;
        dFS=dFSTM;
        
    end
    
    
end

