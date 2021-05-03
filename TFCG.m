function [T,FL,CL,dFS,iters,RESE]=TFCG(TTM,FLTM,FSTM,KLTM,KSTM,CLTM,CSTM,VXTM,VYTM,FRCPSTM,FRCSTM,FKTM,RSCPSTM,RLCPLTM,RSTM,RLTM,CPLTM,step)
%To iteratively solve energy and species conservation equations with finite
%volume method. Model, discretization and iteration can be found "A
%Continuum Model for Computer Simulation of Macrosegregations in Ingots
%during solidification" by Daming Xu 1989 and his companion papers in 1991

%Created 2019-10-17

%Modified for T-FS-CL iteration verification! 2019-11-10

%Modified for second paper 2020-1-1

%Modified for a complete solidification process 2020-1-7

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
global dx
global dy
global gap
global DL
global CE
global dtb
global T0
global TR
global TE
global CL0
global FSE
global FSELOG

err0dFS=10^-8;%controlled accuracy for dFS, T, CL
err0TTM=10^-8;
err0CLTM=10^-8;
phi=0.5;%integration coefficient (0,1)

%======================= VERY IMPORTANT NOTE ==============================
%It is highly recommended the usage of packages for complex formulae, i.e.,
%temporary variables. These packages can be utilized repetedly, easily
%modified, and produce a concise numerical equation.
%======================= VERY IMPORTANT NOTE ==============================


%################################################################### ENERGY BALANCE #########################################################################
%% .............========= DIFFUSION HEAT FLUX ==========...................

KA=FSTM.*KSTM+FLTM.*KLTM;%thermal conductivity of system: only ingot part [W/m/K]; it should be FS.*KS+FL.*KL, but CS is NaN
%NOTE: mean thermal conductivity is based on volume fraction, not mass
%fraction. This is also described in Table 1 of Bennon and Incropera 1988!

%------------- thermal conductivity on finite volume faces ----------------
%....................... all faces, x-axis [W/m/K] ........................
KX=zeros(NIY+2,NIX+1);
%=== (1)Ingot ===
for i=1:NIX+1
    for j=1:NIY+2
        KX(j,i)=0.5*(FKTM(j,i)+FKTM(j,i+1)+FLTM(j,i)*KLTM(j,i)+FLTM(j,i+1)*KLTM(j,i+1));%[(FS*KS)+(FL*KL)]|(j+-0.5,k)
    end
end

%=== (2)Vertical Gas gap ===
%KX(1:NIY+2,NIX+1) is vertical gas gap, heat flux is directly calculated, KX is not necessary!

%........................all faces, y-axis [W/m/K] ........................
KY=zeros(NIY+1,NIX+2);
%=== (1)Ingot ===
for i=1:NIX+2
    for j=1:NIY+1
        KY(j,i)=0.5*(FKTM(j,i)+FKTM(j+1,i)+FLTM(j,i)*KLTM(j,i)+FLTM(j+1,i)*KLTM(j+1,i));%[(FS*KS)+(FL*KL)]|(j,k+-0.5)
    end
end

%------------- thermal conductivity on finite volume faces ----------------

%----------------------------- T gradient ---------------------------------

DTX=zeros(NIY+2,NIX+1);%T gradient between finite volumes in x-axis
for i=2:NIX
    for j=1:NIY+2
        DTX(j,i)=(TTM(j,i+1)-TTM(j,i))/(0.5*dx(i-1)+0.5*dx(i));%[K/m]
    end
end

for i=1:NIY+2
    DTX(i,1)=(TTM(i,2)-TTM(i,1))/(2.0*0.5*dx(1));%1st column == 0.0
    DTX(i,NIX+1)=(TTM(i,NIX+2)-TTM(i,NIX+1))/(2.0*0.5*dx(NIX));%last column; this will be modified later
end

DTY=zeros(NIY+1,NIX+2);%T gradient between finite volumes in y-axis
for i=1:NIX+2
    for j=2:NIY
        DTY(j,i)=(TTM(j+1,i)-TTM(j,i))/(0.5*dy(j-1)+0.5*dy(j));%[K/m]
    end
end

for i=1:NIX+2
    DTY(1,i)=(TTM(2,i)-TTM(1,i))/(2.0*0.5*dy(1));%1st row == 0.0 except 2nd, 3rd
    DTY(NIY+1,i)=(TTM(NIY+2,i)-TTM(NIY+1,i))/(2.0*0.5*dy(NIY));%last row
end

%----------------------------- T gradient ---------------------------------

%----------------------------- Heat flux ----------------------------------

%heat flux on all x faces
TQDX=zeros(NIY+2,NIX+1);%T=T, Q=flux, D==Diffusion, X=x-axis
for i=1:NIX+1
    for j=1:NIY+2
        TQDX(j,i)=KX(j,i)*DTX(j,i);%[W/m^2]
    end
end

%vertical gas gap == right boundary
for i=2:NIY+1
    TIMX=TTM(i,NIX+1)+dx(NIX)*(TTM(i,NIX+1)-TTM(i,NIX))/(dx(NIX)+dx(NIX-1));%T at the left face of air gap by linear extrapolation
    Katm=2.72835e-5+4.898321e-8*(0.5*(TIMX-273.15)+0.5*(TR-273.15));%atmosphere thermal conductivity from D-11 [W/mm/K]
    Katm=Katm*1000.0;%[W/m/K]
    TQDX(i,NIX+1)=(TR-TTM(i,NIX+1))/(0.5*dx(NIX)/KA(i,NIX+1)+gap/Katm);%[W/m^2]
end

%heat flux on all y faces
TQDY=zeros(NIY+1,NIX+2);%T=T, Q=flux, D==Diffusion, Y=y-axis
for i=1:NIX+2
    for j=1:NIY+1
        TQDY(j,i)=KY(j,i)*DTY(j,i);%[W/m^2]
    end
end

%riser boundary
for i=2:3
    TQDY(1,i)=KY(1,i)*(TTM(2,i)-T0);
end

%total heat flux of each finite volume
TQDT=zeros(NIY,NIX);%T=T, Q=flux, D==Diffusion, T=total
for i=1:NIX
    for j=1:NIY
        TQDT(j,i)=-2.0*(TQDX(j+1,i)-TQDX(j+1,i+1))*dtb/dx(i)-2.0*(TQDY(j,i+1)-TQDY(j+1,i+1))*dtb/dy(j);%[J/m^3]
    end
end

%----------------------------- Heat flux ----------------------------------

%% ..............========= INTERNAL HEAT STORAGE ========..................
IHS=zeros(NIY+2,NIX+2);%old T related Internal Heat Storage

%=== (1)ingot ===
for i=1:NIX+2
    for j=1:NIY+2
        IHS(j,i)=FRCPSTM(j,i)*TTM(j,i)+FLTM(j,i)*RLTM(j,i)*CPLTM(j,i)*TTM(j,i);%[J/m^3]
    end
end

%% ..............========= CONVECTION HEAT FLUX =========..................

%------------------------ Part One: RFVX RFVY -----------------------------
%                          mass flux density

%x-axis velocity
%        +-------+
%        |       |
%    --> |   *   |  -->
%        |       |
%        +-------+
%NOTE: in fact, vx at loft and bottom faces should be estimated,
%but they can be interpolated; so only vx normal to vertical faces are
%calculated!

RFVX=zeros(NIY+2,NIX+1);%RL*FL*VX, all VX faces
for i=1:NIX+1
    %VX(1:NIY+2,1)=0.0 --> RFVX(1:NIY+2,1)=0.0 (left boundary condition)
    %VX(1:NIY+2,NIX+1)=0.0 --> RFVX(1:NIY+2,NIX+1)=0.0 (right boundary condition)
    for j=1:NIY+2
        %VX(1,1:NIX+1)=-VX(2,1:NIX+1) --> RFVX(1,1:NIX+1)=-RFVX(2,1:NIX+1) (top boundary condition)
        %VX(NIY+2,1:NIX+1)=-VX(NIY+1,1:NIX+1) --> RFVX(NIY+2,1:NIX+1)=-RFVX(NIY+1,1:NIX+1) (bottom boundary condition)
        RFVX(j,i)=VXTM(j,i)*0.5*(RLTM(j,i)+RLTM(j,i+1))*0.5*(FLTM(j,i)+FLTM(j,i+1));%[kg/m^2/sec]
    end
end

%y-axis velocity
%            ^
%            |
%        +-------+
%        |       |
%        |   *   |
%        |       |
%        +-------+
%            ^
%            |
%NOTE: in fact, vy at left and right faces should be estimated,
%but they can be interpolated; so only vy normal to horizontal faces are
%calculated!

RFVY=zeros(NIY+1,NIX+2);%RL*FL*VY, all VY faces
for i=1:NIX+2
    %VY(1:NIY+1,1)=-VY(1:NIY,2) --> RFVY(1:NIY+1,1)=-RFVY(1:NIY+1,2) (left boundary)
    %VY(1:NIY+1,NIX+2)=-VY(1:NIY,NIX+1) --> RFVY(1:NIY+1,NIX+2)=-RFVY(1:NIY+1,NIX+1) (right boundary)
    for j=1:NIY+1
        %VY(1,4:NIX+2)=0.0 --> RFVY(1,4:NIX+2)=0.0 (top boundary condition except 1st, 2rd, 3nd)
        %VY(NIY+1,1:NIX+2)=0.0 --> RFVY(NIY+1,1:NIX+2)=0.0 (bottom boundary condition)
        RFVY(j,i)=VYTM(j,i)*0.5*(RLTM(j,i)+RLTM(j+1,i))*0.5*(FLTM(j,i)+FLTM(j+1,i));%[kg/m^2/sec]
    end
end

%NOTE: VY(1,2), VY(1,3) is riser,recheck it!!!
if(RFVY(1,2)>=0.0&&RFVY(1,3)>=0.0)
    fprintf('Riser flows downward. Recheck: OK!\n');
end
%------------------------ Part One: RFVX RFVY -----------------------------

%----------------------- Part Two: TRFVX TRFVY ----------------------------
TRFVX=zeros(NIY+2,NIX+1);
for i=1:NIX+1
    %VX(1:NIY+2,1)=0.0 --> TRFVX(1:NIY+2,1)=0.0 (left boundary condition)
    %VX(1:NIY+2,NIX+1)=0.0 --> TRFVX(1:NIY+2,NIX+1)=0.0 (right boundary condition)
    for j=1:NIY+2
        %VX(1,1:NIX+1)=-VX(2,1:NIX+1)  (top boundary condition)
        %VX(NIY+2,1:NIX+1)=-VX(NIY+1,1:NIX+1) (bottom boundary condition)
        TRFVX(j,i)=TTM(j,i)*max(RFVX(j,i),0.0)+TTM(j,i+1)*min(RFVX(j,i),0.0);%[kg.K/m^2/sec]
    end
end

TRFVY=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=1:NIY+1
        %VY(1,4:NIX+2)=0.0 --> TRFVY(1,4:NIX+2)=0.0 (top boundary condition except 1st, 2nd, 3rd)
        %VY(NIY+1,1:NIX+2)=0.0 --> TRFVY(NIY+1,1:NIX+2)=0.0 (bottom boundary condition)
        TRFVY(j,i)=TTM(j,i)*max(RFVY(j,i),0.0)+TTM(j+1,i)*min(RFVY(j,i),0.0);%[kg.K/m^2/sec]
    end
end

%----------------------- Part Two: TRFVX TRFVY ----------------------------

%----------------------- Part Three: Heat flux ----------------------------
TQVX=zeros(NIY+2,NIX);%T=T, Q=flux, VX=VX
for i=1:NIX
    for j=1:NIY+2
        TQVX(j,i)=-dtb*CPLTM(j,i+1)*(TRFVX(j,i+1)-TRFVX(j,i))/dx(i);%[J/m^3]
    end
end

TQVY=zeros(NIY,NIX+2);%T=T, Q=flux, VY=VY
for i=1:NIX+2
    for j=1:NIY
        TQVY(j,i)=-dtb*CPLTM(j+1,i)*(TRFVY(j+1,i)-TRFVY(j,i))/dy(j);%[J/m^3]
    end
end

TQVT=zeros(NIY,NIX);%T=T, Q=flux, T=total
TQVT(1:NIY,1:NIX)=TQVX(2:NIY+1,1:NIX)+TQVY(1:NIY,2:NIX+1);%only in ingot [J/m^3]

%------------------------ Part Three: Heat flux ---------------------------

%------------------------ Part Four: Latent heat --------------------------
dFSTM=zeros(NIY+2,NIX+2);%temporary dFS

HSTM=zeros(NIY+2,NIX+2);%latent heat from D-7
for i=1:NIX+2
    for j=1:NIY+2
        if(CSTM(j,i)<0.1)
            HSTM(j,i)=397.67-2.3288*CSTM(j,i)*100.0;%[J/g]
        else
            HSTM(j,i)=333.59;%[J/g]
        end
    end
end
HSTM=HSTM*1000.0;%[J/kg]

LH=zeros(NIY+NMY+2,NIX+NMX+2);%latent heat [J/m^3]
for i=1:NIX+2
    for j=1:NIY+2
        LH(j,i)=dFSTM(j,i)*RSTM(j,i)*HSTM(j,i);
    end
end
%------------------------ Part Four: Latent heat --------------------------


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
    %CL(1:NIY+2,1)=CL(1:NIY+2,2) --> CLQDX(1:NIY+2,1)=0.0 (left boundary)
    %CL(1:NIY+2,NIX+1)=CL(1:NIY+2,NIX+2) --> CLQDX(1:NIY+2,NIX+1)=0.0 (right boundary)
    for j=1:NIY+2
        CLQDX(j,i)=2.0*dtb*DL*(FRCLX(j,i+1)-FRCLX(j,i))/(dx(i-1)+dx(i));%[kg/m^2]
    end
end

for j=1:NIY+2
    CLQDX(j,1)=2.0*dtb*DL*(FRCLX(j,2)-FRCLX(j,1))/(dx(1)+dx(1));%set dx(0)==dx(1); left boundary
    CLQDX(i,NIX+1)=2.0*dtb*DL*(FRCLX(j,NIX+2)-FRCLX(j,NIX+1))/(dx(NIX)+dx(NIX));%set dx(NIX+1)==dx(NIX);right boundary
end

%y-axis diffusion
CLQDY=zeros(NIY+1,NIX+2);%CL=CL, Q=flux, D=diffusion, Y=y-axis
for i=1:NIX+2
    for j=2:NIY
        %CL(1,1:NIX+2)=CL(2,1:NIX+2) --> CLQDY(1,1:NIX+2)=0.0
        %CL(NIY+2,1:NIX+2)=CL(NIY+1,1:NIX+2) --> CLQDY(NIY+1,1:NIX+2)=0.0
        CLQDY(j,i)=2.0*dtb*DL*(FRCLY(j+1,i)-FRCLY(j,i))/(dy(j-1)+dy(j));%[kg/m^2]
    end
end

for i=1:NIX+2
    CLQDY(1,i)=2.0*dtb*DL*(FRCLY(2,i)-FRCLY(1,i))/(dy(1)+dy(1));%set dy(0)==dy(1); top boundary
    CLQDY(NIY+1,i)=2.0*dtb*DL*(FRCLY(NIY+2,i)-FRCLY(NIY+1,i))/(dy(NIY)+dy(NIY));%set dy(NIY)==dy(NIY+1); bottom boundary
end

%total species flux of each finite volume
CLQDT=zeros(NIY,NIX);%CL=CL, Q=flux, D=diffusion, T=total
for i=1:NIX
    for j=1:NIY
        CLQDT(j,i)=(CLQDX(j+1,i+1)-CLQDX(j+1,i))/dx(i)+(CLQDY(j+1,i+1)-CLQDY(j,i+1))/(dy(j));%[kg/m^3]
    end
end

%% .............========= INTERNAL SPECIES STORAGE =========...............

KP=zeros(NIY+2,NIX+2);%solid/liquid partition coefficient
for i=1:NIX+2
    for j=1:NIY+2
        if(CLTM(j,i)<=CE)
            KP(j,i)=0.12824+5.699124e-5*CLTM(j,i)*100.0+3.728777e-5*(CLTM(j,i)*100.0)^2;%solid/liquid partition coefficient from D-15 [1]
        else
            KP(j,i)=1.0;
        end
    end
end

ISS=zeros(NIY+2,NIX+2);%internal species storage
FRLTM=zeros(NIY+2,NIX+2);%old time FL*RL
RSKPTM=zeros(NIY+2,NIX+2);%old time RS*KP

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
for i=1:NIX+1
    %VX(1:NIY+2,1)=0.0 --> CLRFVX(1:NIY+2,1)=0.0
    %VX(1:NIY+2,NIX+1)=0.0 --> CLRFVX(1:NIY+2,NIX+1)=0.0
    for j=1:NIY+2
        CLRFVX(j,i)=CLTM(j,i)*max(RFVX(j,i),0.0)+CLTM(j,i+1)*min(RFVX(j,i),0.0);%[kg/m^2/sec]
    end
end

%y-axis species convection flux
CLRFVY=zeros(NIY+1,NIX+2);
for i=1:NIX+2
    for j=1:NIY+1
        %VY(1,4:NIX+2)=0.0 --> CLRFVY(1,4:NIX+2)=0.0
        %VY(NIY+1,1:NIX+2)=0.0 --> CLRFVY(NIY+1,1:NIX+2)=0.0
        CLRFVY(j,i)=CLTM(j,i)*max(RFVY(j,i),0.0)+CLTM(j+1,i)*min(RFVY(j,i),0.0);%[kg/m^2/sec]
    end
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

TFN=zeros(NIY+2,NIX+2);%T Factor of Next step
TN=zeros(NIY+2,NIX+2);%T of Next step
TLTM=zeros(NIY+2,NIX+2);%Liquidus of old step
%HSTM=zeros(NIY+2,NIX+2);%latent heat of Next step
iters=zeros(NIY+2,NIX+2);
TTME=TTM;%for eutectic T calculation

%TFN: T Factor of Next step
%=== (1)ingot ===
for i=1:NIX+2
    for j=1:NIY+2
        TFN(j,i)=FRCPSTM(j,i)+FLTM(j,i)*RLTM(j,i)*CPLTM(j,i)+(RSCPSTM(j,i)-RLCPLTM(j,i))*dFSTM(j,i);%[J/m^3/K]
    end
end

%----------------------------- STEP ONE -----------------------------------

%Let the initial value of dFSTM equal zero and calculate the approximation of TN using Eq. (8).
for j=2:NIX+1
    for k=2:NIY+1
        TN(k,j)=(IHS(k,j)+TQVT(k-1,j-1)+TQDT(k-1,j-1)+LH(k,j))/TFN(k,j);%[J/m^3 divided by J/m^3/K == K]
    end
end
%RSCPSTM: RS*CPS at n step, as approxiamation of RS*CPS at n+1 step
%RLCPLTM: RS*CPL at n step, as approxiamation of RS*CPL at n+1 step
TN(1,2:NIX+1)=TN(2,2:NIX+1);%top boundary
TN(1,1:3)=T0;
TN(NIY+2,2:NIX+1)=TN(NIY+1,2:NIX+1);%bottom boundary
TN(2:NIY+2,1)=TN(2:NIY+2,2);

%TN(1:NIY+2,1)=TN(1:NIY+2,2);%left boundary
TN(2:NIY+1,NIX+2)=TN(2:NIY+1,NIX+1)+dx(NIX)*(TN(2:NIY+1,NIX+1)-TN(2:NIY+1,NIX))/(dx(NIX)+dx(NIX-1));%right boundary by linear extrapolation, for graphics only
TN(1,NIX+2)=TN(2,NIX+2);%only for graphics
TN(NIY+2,NIX+2)=TN(NIY+1,NIX+2);%only for graphics
TTM=TN;

%----------------------------- STEP ONE -----------------------------------
TLTM=660.37-2.34581*CLTM*100.0-3.129e-2*(CLTM*100.0).^2+273.15;%old liquidus to determine superliquidus or subliquidus [K]

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
TEC=zeros(NIY+2,NIX+2);%New T for Energy Check
dFSN=zeros(NIY+2,NIX+2);%dFS of Next step

LAT=zeros(NIY,NIX);%latent heat for energy check
caseID=3*ones(NIY,NIX);%case ID
dFSBE=zeros(NIY+2,NIX+2);%solid fraction increment before eutectic
dFSAE=zeros(NIY+2,NIX+2);%solid fraction increment after eutectic
dFSSEN=zeros(NIY+2,NIX+2);%used in energy check for sensible heat

RSE=3400.0;%eutectic solid density [kg/m^3]
HSE=333.59*1000.0;%eutectic latent heat [J/kg]
RLE=2.5222e-3+2.703e-5*CE*100.0-3.16e-7*(TE-273.15);%liquid density from D-6 [g/mm^3];  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
RLE=RLE*10^6;%[kg/m^3]
CPLE=1.086-5.928e-3*CE*100.0;%liquid specific heat capacity from D-4 [J/g/K]
CPLE=CPLE*1000.0;%[J/kg/K]

for i=1:NIX
    for j=1:NIY
        
        dFSCAL=(FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1))*(TE-TN(j+1,i+1))/(RSE*HSE-(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*TE);%Eq.(14)
        %NOTE: RSE*HSE is heat source while RSCPSTM-RLCPLTM is heat sink
        %since in this system the heat storage capacity of solid is better than that of liquid
        
        caseID(j,i)=3;%flag of condition
        if(TN(j+1,i+1)>TLTM(j+1,i+1))
            caseID(j,i)=1;
        end
        if((CLTM(j+1,i+1)<CE)&&(TN(j+1,i+1)>=TE)&&(TN(j+1,i+1)<=TLTM(j+1,i+1)))
            caseID(j,i)=2;
        end
        if((CLTM(j+1,i+1)>=CE)&&(FSTM(j+1,i+1)+dFSCAL<=1.0)&&(TN(j+1,i+1)<TE))
            caseID(j,i)=4;
        end
        if((FSTM(j+1,i+1)<1.0)&&(TN(j+1,i+1)<TE)&&(dFSCAL+FSTM(j+1,i+1)>1.0))
            caseID(j,i)=5;
        end
        if(FSTM(j+1,i+1)>=1.0)
            caseID(j,i)=6;
        end
        
        %         if(cases==3)
        %             fprintf('dddd\n');
        %         end
        
        switch caseID(j,i)
            
            case 1%SITUATION 1: SUPERLIQUIDUS
                dFSN(j+1,i+1)=dFSTM(j+1,i+1);%For superliquidus, we do nothing!
                FLN(j+1,i+1)=FLTM(j+1,i+1);
                CLN(j+1,i+1)=CLTM(j+1,i+1);
                iters(j+1,i+1)=1;
                
                LAT(j,i)=dFSN(j+1,i+1)*RSTM(j+1,i+1)*HSTM(j+1,i+1);%latent heat for energy check
                TEC(j+1,i+1)=TN(j+1,i+1);%T for energy check
                
                %SITUATION 2: LIQUID --> MUSHY REGION
                %0.0<=FSTM<FSE, TE<=TN<=TLTM --> CLTM<=CE, TE<=TN<=TLTM
            case 2%&&(TN(j+1,i+1)<=TLTM(j+1,i+1))
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
                %New RL updated by new T
                RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLTM(j+1,i+1)*100.0-3.16e-7*(TN(j+1,i+1)-273.15);%liquid density of Next step from D-6 [g/mm^3]
                RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
                
                %New KP updated by new T
                %NOTE: since CLTM<CE, thus KPN will follow the first
                %statement, and there's no chance that the second statement
                %will be executed!!!
                if(CLTM(j+1,i+1)<CE)
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
                
                TLN(j+1,i+1)=660.37-2.34581*CLN(j+1,i+1)*100.0-3.129e-2*(CLN(j+1,i+1)*100.0)^2+273.15;%new liquidus [K]
                
                iters(j+1,i+1)=0;
                m=1;
                %----------------------------- STEP TWO -----------------------------------
                while(err(j+1,i+1))
                    
                    %-------------------------- STEP THREE --------------------------------
                    %Calculate the approximation of dFSN by Eq. (10) and take the average:
                    %dFSN=0.5*(dFSN+dFSTM)
                    
                    dFSN(j+1,i+1)=(FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+FRCPSTM(j+1,i+1))*(TLN(j+1,i+1)-TN(j+1,i+1))/(RSTM(j+1,i+1)*HSTM(j+1,i+1));%[1]
                    
                    dFSN(j+1,i+1)=0.5*(dFSN(j+1,i+1)+dFSTM(j+1,i+1));
                    
                    errdFS(j+1,i+1)=abs((dFSN(j+1,i+1)-dFSTM(j+1,i+1))/dFSN(j+1,i+1));
                    
                    dFSTM(j+1,i+1)=dFSN(j+1,i+1);
                    %-------------------------- STEP THREE --------------------------------
                    %                                         figure(2)
                    %                     plot(m,TN(j+1,i+1),'g*');
                    %                     hold on
                    %
                    %                     plot(m,TLN(j+1,i+1),'r*');
                    %                     hold on
                    %                     figure(3)
                    %                     plot(m,dFSN(j+1,i+1),'k*');
                    %                     hold on
                    %                     m=m+1;
                    
                    %--------------------------- STEP FOUR --------------------------------
                    %Again using Eq. (8), calculate the new approximation of T(n+1) with new initial
                    %values, dFSN and CLN, and take TN=0.5*(TN+TTM).
                    
                    %NOTE: new T, CL have effect on RL,RS,KP,CpS,CpL
                    %Partition coefficient of Next step
                    
                    %NOTE: since CLTM<CE, thus KPN will follow the first
                    %statement, and there's no chance that the second statement
                    %will be executed!!!
                    %                     if(CLN(j+1,i+1)<CE)
                    %                         KPN(j+1,i+1)=0.12824+5.699124e-5*CLN(j+1,i+1)*100.0+3.728777e-5*(CLN(j+1,i+1)*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
                    %                     else
                    %                         KPN(j+1,i+1)=1.0;
                    %                     end
                    %
                    %                     %RS of Next step
                    %                     CSN(j+1,i+1)=CLN(j+1,i+1)*KPN(j+1,i+1);
                    %                     if(CSN(j+1,i+1)<CE)
                    %                         RSN(j+1,i+1)=2.58e-3;%[g/mm^3]
                    %                     else
                    %                         RSN(j+1,i+1)=3.4e-3;%[g/mm^3]
                    %                     end
                    %                     RSN(j+1,i+1)=RSN(j+1,i+1)*10^6;%[kg/m^3]
                    %
                    %                     if(CSN(j+1,i+1)<0.1)
                    %                         HSN(j+1,i+1)=397.67-2.3288*CSN(j+1,i+1)*100.0;%latent heat from D-7 [J/g]
                    %                     else
                    %                         HSN(j+1,i+1)=333.59;%[J/g]
                    %                     end
                    %                     HSN(j+1,i+1)=HSN(j+1,i+1)*1000.0;%[J/kg]
                    %
                    %                     %CpS of Next step
                    %                     CPSN(j+1,i+1)=0.88+4.446e-4*(TN(j+1,i+1)-273.15)-2.274e-3*CSN(j+1,i+1)*100.0;%solid specific heat capacity from D-3 [J/g/K]
                    %                     CPSN(j+1,i+1)=CPSN(j+1,i+1)*1000.0;%[J/kg/K]
                    %
                    %                     %RL of Next step
                    %                     RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLN(j+1,i+1)*100.0-3.16e-7*(TN(j+1,i+1)-273.15);%liquid density of Next step from D-6 [g/mm^3]
                    %                     RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
                    %
                    %                     %CpL of Next step
                    %
                    %                     CPLN(j+1,i+1)=1.086-5.928e-3*CLN(j+1,i+1)*100.0;%liquid specific heat capacity from D-4 [J/g/K]
                    %                     CPLN(j+1,i+1)=CPLN(j+1,i+1)*1000.0;%[J/kg/K]
                    %
                    %
                    %                     RSCPSTM(j+1,i+1)=RSN(j+1,i+1)*CPSN(j+1,i+1);%[J/m^3/K]
                    %                     RLCPLTM(j+1,i+1)=RLN(j+1,i+1)*CPLN(j+1,i+1);%[J/m^3/K]
                    
                    LH(j+1,i+1)=dFSN(j+1,i+1)*RSTM(j+1,i+1)*HSTM(j+1,i+1);%latent heat [J/m^3]
                    
                    TFN(j+1,i+1)=FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*dFSN(j+1,i+1);%[J/m^3/K]
                    
                    TN(j+1,i+1)=(IHS(j+1,i+1)+TQVT(j,i)+TQDT(j,i)+LH(j+1,i+1))/TFN(j+1,i+1);%[J/m^3 divided by J/m^3/K == K]
                    TN(j+1,i+1)=0.5*(TN(j+1,i+1)+TTM(j+1,i+1));
                    
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
                    
                    TLN(j+1,i+1)=660.37-2.34581*CLN(j+1,i+1)*100.0-3.129e-2*(CLN(j+1,i+1)*100.0)^2+273.15;%ingot liquidus from D-16 [K]; CLN here not in wt% while in wt% from D-16
                    
                    %-------------------------- STEP SEVEN --------------------------------
                    iters(j+1,i+1)=iters(j+1,i+1)+1;
                end
                
                LAT(j,i)=dFSN(j+1,i+1)*RSTM(j+1,i+1)*HSTM(j+1,i+1);%latent heat for energy check
                TEC(j+1,i+1)=TN(j+1,i+1);%T for energy check
                
                %THERMAL EQUILIBRIUM --> T==Tliq
                TN(j+1,i+1)=TLN(j+1,i+1);
                
            case 4%SITUATION 4: EUTECTIC
                modifymarker=modifymarker+1;
                
                dFSN(j+1,i+1)=dFSCAL;
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSN(j+1,i+1);
                TN(j+1,i+1)=TE;
                CLN(j+1,i+1)=CE;
                iters(j+1,i+1)=1;
                
                TEC(j+1,i+1)=TN(j+1,i+1);%T for energy check
                LAT(j,i)=dFSN(j+1,i+1)*RSE*HSE;%%latent heat for energy check
                %NOTE: dFSN(j+1,i+1)*(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*TE is sensible heat
                
            case 5%SITUATION 5: EUTECTIC --> COMPLETE SOLID
                modifymarker=modifymarker+1;
                
                dFSN(j+1,i+1)=1.0-FSTM(j+1,i+1);
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSN(j+1,i+1);
                
                TN(j+1,i+1)=TN(j+1,i+1)+(RSE*HSE-(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*TE)*(1.0-FSTM(j+1,i+1))/(FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1));
                CLN(j+1,i+1)=CE;
                iters(j+1,i+1)=1;
                TEC(j+1,i+1)=TN(j+1,i+1);%T for energy check
                
            case 6%COMPLETE SOLID COOLING
                
                FLN(j+1,i+1)=0.0;
                dFSN(j+1,i+1)=0.0;
                CLN(j+1,i+1)=CE;
                iters(j+1,i+1)=1;
                
                LAT(j,i)=dFSN(j+1,i+1)*RSTM(j+1,i+1)*HSTM(j+1,i+1);%latent heat for energy check
                TEC(j+1,i+1)=TN(j+1,i+1);%T for energy check
                
            otherwise%SITUATION 3: MUSHY REGION --> EUTECTIC
                
                TNCAL=TN(j+1,i+1);
                dFSTMCAL=dFSTM(j+1,i+1);
                TTMCAL=TTM(j+1,i+1);
                %FL of Next step
                FLNCAL=FLTM(j+1,i+1)-dFSTMCAL;
                
                %   if(FLNCAL>0.0)%correspoding to statement in the following: fprintf('Complete solid BUT should be eutectic\n');
                %NOTE: new T has effect on RL,RS,KP
                %New RL updated by new T
                RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLTM(j+1,i+1)*100.0-3.16e-7*(TNCAL-273.15);%liquid density of Next step from D-6 [g/mm^3]
                RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
                
                %New KP updated by new T
                if(CLTM(j+1,i+1)<CE)
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
                CLFNCAL=FLNCAL*RLN(j+1,i+1)+phi*dFSTM(j+1,i+1)*RSN(j+1,i+1)*KPN(j+1,i+1);%[kg/m^3]
                CLNCAL=(ISS(j+1,i+1)+CLQVT(j,i)+CLQDT(j,i))/CLFNCAL;%[kg/m^3 divided by kg/m^3 ==1]; ISS should update with NEW dFS ,but no NEW dFS available now, set old dFSTM==0.0 as new one
                errCLTM(j+1,i+1)=abs((CLNCAL-CLTM(j+1,i+1))/CLNCAL);
                
                CLTMCAL=CLNCAL;
                
                TLNCAL=660.37-2.34581*CLNCAL*100.0-3.129e-2*(CLNCAL*100.0)^2+273.15;%new liquidus [K]
                m=0;
                %----------------------------- STEP TWO -----------------------------------
                while(err(j+1,i+1))
                    m=m+1;
                    
                    %-------------------------- STEP THREE --------------------------------
                    %Calculate the approximation of dFSN by Eq. (10) and take the average:
                    %dFSN=0.5*(dFSN+dFSTM)
                    
                    dFSNCAL=(FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+FRCPSTM(j+1,i+1))*(TLNCAL-TNCAL)/(RSTM(j+1,i+1)*HSTM(j+1,i+1));%[1]
                    
                    dFSNCAL=0.5*(dFSNCAL+dFSTMCAL);
                    
                    errdFS(j+1,i+1)=abs((dFSNCAL-dFSTMCAL)/dFSNCAL);
                    
                    dFSTMCAL=dFSNCAL;
                    %-------------------------- STEP THREE --------------------------------
                    
                    %--------------------------- STEP FOUR --------------------------------
                    %Again using Eq. (8), calculate the new approximation of T(n+1) with new initial
                    %values, dFSN and CLN, and take TN=0.5*(TN+TTM).
                    
                    %NOTE: new T, CL have effect on RL,RS,KP,CpS,CpL
                    %Partition coefficient of Next step
                    %                         KPN(j+1,i+1)=0.12824+5.699124e-5*CLNCAL*100.0+3.728777e-5*(CLNCAL*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
                    %
                    %                         %RS of Next step
                    %                         CSN(j+1,i+1)=CLNCAL*KPN(j+1,i+1);
                    %                         if(CSN(j+1,i+1)<CE)
                    %                             RSN(j+1,i+1)=2.58e-3;%[g/mm^3]
                    %                         else
                    %                             RSN(j+1,i+1)=3.4e-3;%[g/mm^3]
                    %                         end
                    %                         RSN(j+1,i+1)=RSN(j+1,i+1)*10^6;%[kg/m^3]
                    %
                    %                         if(CSN(j+1,i+1)<0.1)
                    %                             HSN(j+1,i+1)=397.67-2.3288*CSN(j+1,i+1)*100.0;%latent heat from D-7 [J/g]
                    %                         else
                    %                             HSN(j+1,i+1)=333.59;%[J/g]
                    %                         end
                    %                         HSN(j+1,i+1)=HSN(j+1,i+1)*1000.0;%[J/kg]
                    %
                    %                         %CpS of Next step
                    %                         CPSN(j+1,i+1)=0.88+4.446e-4*(TNCAL-273.15)-2.274e-3*CSN(j+1,i+1)*100.0;%solid specific heat capacity from D-3 [J/g/K]
                    %                         CPSN(j+1,i+1)=CPSN(j+1,i+1)*1000.0;%[J/kg/K]
                    %
                    %                         %RL of Next step
                    %                         RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLNCAL*100.0-3.16e-7*(TNCAL-273.15);%liquid density of Next step from D-6 [g/mm^3]
                    %                         RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
                    %
                    %                         %CpL of Next step
                    %
                    %                         CPLN(j+1,i+1)=1.086-5.928e-3*CLNCAL*100.0;%liquid specific heat capacity from D-4 [J/g/K]
                    %                         CPLN(j+1,i+1)=CPLN(j+1,i+1)*1000.0;%[J/kg/K]
                    %
                    %
                    %                         RSCPSTM(j+1,i+1)=RSN(j+1,i+1)*CPSN(j+1,i+1);%[J/m^3/K]
                    %                         RLCPLTM(j+1,i+1)=RLN(j+1,i+1)*CPLN(j+1,i+1);%[J/m^3/K]
                    
                    LH(j+1,i+1)=dFSNCAL(j+1,i+1)*RSTM(j+1,i+1)*HSTM(j+1,i+1);%latent heat [J/m^3]
                    
                    TFNCAL=FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*dFSNCAL;%[J/m^3/K]
                    
                    TNCAL=(IHS(j+1,i+1)+TQVT(j,i)+TQDT(j,i)+LH(j+1,i+1))/TFNCAL;%[J/m^3 divided by J/m^3/K == K]
                    TNCAL=0.5*(TNCAL+TTMCAL);
                    
                    errTTM(j+1,i+1)=abs((TNCAL-TTMCAL)/TNCAL);
                    
                    TTMCAL=TNCAL;
                    %--------------------------- STEP FOUR --------------------------------
                    
                    %--------------------------- STEP FIVE --------------------------------
                    %With the new initial value TN (as well as dFSN and CLN), calculate the new approximation
                    %of CLN(n+1) using Eq. (9), and take CLN=0.5*(CLN+CLTM).
                    
                    %FL of of Next step
                    FLNCAL=FLTM(j+1,i+1)-dFSNCAL;
                    
                    %NOTE: new T have effect on RL,CpS,CpL; RS,KP are functions of CLN
                    %which is updated before so remain unchanged!
                    KPN(j+1,i+1)=0.12824+5.699124e-5*CLNCAL*100.0+3.728777e-5*(CLNCAL*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
                    
                    %RS of Next step
                    CSN(j+1,i+1)=CLNCAL*KPN(j+1,i+1);
                    if(CSN(j+1,i+1)<CE)
                        RSN(j+1,i+1)=2.58e-3;%[g/mm^3]
                    else
                        RSN(j+1,i+1)=3.4e-3;%[g/mm^3]
                    end
                    RSN(j+1,i+1)=RSN(j+1,i+1)*10^6;%[kg/m^3]
                    
                    %RL of Next step
                    RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLNCAL*100.0-3.16e-7*(TNCAL-273.15);%liquid density of Next step from D-6 [g/mm^3]
                    RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
                    
                    
                    %New ISS due to updated dFSN
                    
                    ISS(j+1,i+1)=FRLTM(j+1,i+1)-(1.0-phi)*dFSNCAL*RSKPTM(j+1,i+1);%[kg/m^3]
                    
                    %Factor ahead CL(n+1)
                    CLFNCAL=FLNCAL*RLN(j+1,i+1)+phi*dFSNCAL*RSN(j+1,i+1)*KPN(j+1,i+1);%[kg/m^3]
                    
                    %CL of Next step
                    CLNCAL=(ISS(j+1,i+1)+CLQVT(j,i)+CLQDT(j,i))/CLFNCAL;%[kg/m^3 divided by kg/m^3 ==1]
                    CLNCAL=0.5*(CLNCAL+CLTMCAL);
                    
                    %NOTE: it's coupled change of T-FS-CL that makes
                    %dFS very small; see following plots
                    % figure(1)
                    % plot(m,TNCAL,'k*',m,TLNCAL,'r*');
                    % hold on
                    % figure(2)
                    % plot(m,dFSNCAL,'g*');
                    % hold on
                    
                    errCLTM(j+1,i+1)=abs((CLNCAL-CLTMCAL)/CLNCAL);
                    CLTMCAL=CLNCAL;
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
                    
                    TLNCAL=660.37-2.34581*CLNCAL*100.0-3.129e-2*(CLNCAL*100.0)^2+273.15;%ingot liquidus from D-16 [K]; CLN here not in wt% while in wt% from D-16
                    %-------------------------- STEP SEVEN --------------------------------
                    iters(j+1,i+1)=1;
                    
                end
                
                if(TLNCAL>TE)%back to case 2
                    FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSNCAL;
                    dFSN(j+1,i+1)=dFSNCAL;
                    CLN(j+1,i+1)=CLNCAL;
                    
                    LAT(j,i)=dFSNCAL*RSTM(j+1,i+1)*HSTM(j+1,i+1);%latent heat for energy check
                    TEC(j+1,i+1)=TNCAL;%T for energy check
                    
                    %THERMAL EQUILIBRIUM --> T==Tliq
                    TN(j+1,i+1)=TLNCAL;
                    
                    caseID(j,i)=2;%mark as situation 2
                    
                else
                    dFSCAL=dFSNCAL;
                    dFSSEN(j+1,i+1)=dFSCAL;%used in energy check for sensible heat
                    modifymarker=modifymarker+1;
                    
                    dFSBE(j+1,i+1)=dFSCAL*(TTME(j+1,i+1)-TE)/(TTME(j+1,i+1)-TLNCAL);%solid fraction before eutectic, Eq.(11)
                    % dFSAE1=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TTME(j+1,i+1)-TN(j+1,i+1))*(TE-TLNCAL)/((TTME(j+1,i+1)-TLNCAL)*RSE*HSE);%solid fraction after eutectic, Eq.(12)
                    dFSAE(j+1,i+1)=dFSCAL*(TE-TLNCAL)*RSTM(j+1,i+1)*HSTM(j+1,i+1)/((TTME(j+1,i+1)-TLNCAL)*RSE*HSE);
                    dFSN(j+1,i+1)=dFSBE(j+1,i+1)+dFSAE(j+1,i+1);
                    FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSN(j+1,i+1);
                    TN(j+1,i+1)=TE;
                    CLN(j+1,i+1)=CE;
                    
                    LAT(j,i)=dFSBE(j+1,i+1)*RSTM(j+1,i+1)*HSTM(j+1,i+1);%latent heat for energy check
                    LAT(j,i)=LAT(j,i)+dFSAE(j+1,i+1)*RSE*HSE;%latent heat for energy check
                    TEC(j+1,i+1)=TNCAL;%T for energy check
                    
                    if(FSELOG(j+1,i+1))%FSE only modify once
                        FSE(j+1,i+1)=FSTM(j+1,i+1)+dFSBE(j+1,i+1);%real solid fraction when eutectic
                        FSELOG(j+1,i+1)=-1.0;
                    end
                    
                    fprintf('=============== FSE(%2d,%2d)= %7.5f ==============\n',j+1,i+1,FSE(j+1,i+1));
                    
                    %                 else
                    %                     fprintf('Complete solid BUT should be eutectic\n');
                    %                 end
                    
                end
                
        end
    end
    %################################################################ T-FS-CL ITERATION #########################################################################
    
    %Ingot T-FS-CL has been iteratively sovled?
    if(modifymarker>0.0)%Yes, ingot changed (i.e., subliquidus)
        %NOTE: T in mould is the same as T without T-FS-CL iteration, since
        %IHS, TQVT, TQDT and TFN stay as old values
        
        %Second: returen new field variables
        TN(1,2:NIX+1)=TN(2,2:NIX+1);%abdiabatic top boundary
        TN(1,1:3)=T0;
        TN(NIY+2,2:NIX+1)=TN(NIY+1,2:NIX+1);%adiabatic bottom boundary
        TN(2:NIY+2,1)=TN(2:NIY+2,2);%abdiabatic left boundary
        %TN(1:NIY+2,1)=TN(1:NIY+2,2);%abdiabatic left boundary
        TN(2:NIY+1,NIX+2)=TN(2:NIY+1,NIX+1)+dx(NIX)*(TN(2:NIY+1,NIX+1)-TN(2:NIY+1,NIX))/(dx(NIX)+dx(NIX-1));%right boundary by linear extrapolation, for graphics only
        TN(1,NIX+2)=TN(2,NIX+2);%only for graphics
        TN(NIY+2,NIX+2)=TN(NIY+1,NIX+2);%only for graphics
        T=TN;
        
        FLN(1,2:NIX+1)=FLN(2,2:NIX+1);%top boundary
        FLN(1,2:3)=1.0;
        FLN(NIY+2,2:NIX+1)=FLN(NIY+1,2:NIX+1);%bottom boundary
        FLN(2:NIY+1,1)=FLN(2:NIY+1,2);%left boundary
        %FLN(1:NIY+1,1)=FLN(1:NIY+1,2);%left boundary
        FLN(2:NIY+1,NIX+2)=FLN(2:NIY+1,NIX+1);%right boundary
        FLN(1,1)=FLN(1,2);
        FLN(1,NIX+2)=FLN(1,NIX+1);
        FLN(NIY+2,1)=FLN(NIY+2,2);
        FLN(NIY+2,NIX+2)=FLN(NIY+2,NIX+1);
        FL=FLN;
        
        CLN(1,2:NIX+1)=CLN(2,2:NIX+1);%top boundary
        CLN(1,2:3)=CL0;
        CLN(NIY+2,2:NIX+1)=CLN(NIY+1,2:NIX+1);%bottom boundary
        CLN(2:NIY+1,1)=CLN(2:NIY+1,2);%left boundary
        %CLN(1:NIY+1,1)=CLN(1:NIY+1,2);%left boundary
        CLN(2:NIY+1,NIX+2)=CLN(2:NIY+1,NIX+1);%right boundary
        CLN(1,1)=CLN(1,2);
        CLN(1,NIX+2)=CLN(1,NIX+1);
        CLN(NIY+2,1)=CLN(NIY+2,2);
        CLN(NIY+2,NIX+2)=CLN(NIY+2,NIX+1);
        CL=CLN;
        
        dFSN(1,2:NIX+1)=dFSN(2,2:NIX+1);%top boundary
        dFSN(1,2:3)=0.0;
        dFSN(NIY+2,2:NIX+1)=dFSN(NIY+1,2:NIX+1);%bottom boundary
        dFSN(2:NIY+1,1)=dFSN(2:NIY+1,2);%left boundary
        %dFSN(1:NIY+1,1)=dFSN(1:NIY+1,2);%left boundary
        
        dFSN(2:NIY+1,NIX+2)=dFSN(2:NIY+1,NIX+1);%right boundary
        dFSN(1,1)=dFSN(1,2);
        dFSN(1,NIX+2)=dFSN(1,NIX+1);
        dFSN(NIY+2,1)=dFSN(NIY+2,2);
        dFSN(NIY+2,NIX+2)=dFSN(NIY+2,NIX+1);
        dFS=dFSN;
        
    else%No, only ingot T changed (i.e., superliquidus); FS, CL unchange! return updated T and constant FL,CL
        
        TN(1,2:NIX+1)=TN(2,2:NIX+1);%abdiabatic top boundary
        TN(1,1:3)=T0;
        TN(NIY+2,2:NIX+1)=TN(NIY+1,2:NIX+1);%adiabatic bottom boundary
        TN(2:NIY+2,1)=TN(2:NIY+2,2);%abdiabatic left boundary
        %TN(1:NIY+2,1)=TN(1:NIY+2,2);%abdiabatic left boundary
        
        TN(2:NIY+1,NIX+2)=TN(2:NIY+1,NIX+1)+dx(NIX)*(TN(2:NIY+1,NIX+1)-TN(2:NIY+1,NIX))/(dx(NIX)+dx(NIX-1));%right boundary by linear extrapolation, for graphics only
        TN(1,NIX+2)=TN(2,NIX+2);%only for graphics
        TN(NIY+2,NIX+2)=TN(NIY+1,NIX+2);%only for graphics
        T=TN;
        
        FL=FLTM;
        FL(1,2:3)=1.0;
        CL=CLTM;
        CL(1,2:3)=CL0;
        dFS=dFSTM;
        dFS(1,2:3)=0.0;
        
    end
    
    %% -------------------------- ENERGY CHECK --------------------------------
    
    %IMPORTANT NOTE: if convection is not present, then energy check can be applied to
    %every cell and also system collectively; otherwise, only system energy
    %check can be easily made !!!
    
    %Prapare for sensible heat and latent heat
    SEN=zeros(NIY,NIX);%sensible heat
    RESE=zeros(NIY+2,NIX+2);%difference between energy removed and internal energy loss
    energy_check=1.0;%energy check marker
    for i=1:NIX
        for j=1:NIY
            SEN(j,i)=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TEC(j+1,i+1)-TTME(j+1,i+1))+(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*TEC(j+1,i+1)*dFSN(j+1,i+1);%sensible heat
            if(caseID(j,i)==3)%MUSHY --> EUTECTIC
                SEN(j,i)=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TEC(j+1,i+1)-TTME(j+1,i+1))+(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*TEC(j+1,i+1)*dFSSEN(j+1,i+1);%sensible heat
            end
            if(caseID(j,i)==5)%EUTECTIC --> COMPLETE SOLID
                SEN(j,i)=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TEC(j+1,i+1)-TTME(j+1,i+1))+(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*TE*dFSN(j+1,i+1);%sensible heat
            end
            
            RESE(j+1,i+1)=(SEN(j,i)-TQVT(j,i)-LAT(j,i)-TQDT(j,i))/TQDT(j,i);%1, 2, 3 iterms are the energy lost from system, 4th iterm is the energy removed from system; this is relative residual with respect to energy loss
            %NOTE: RESE does not include geometry configuration, i.e., dx, dy
            
            if(abs(RESE(j+1,i+1))>=1.0e-4)
                fprintf('Energy Check: Cell (%2d,%2d) not balanced!\n',j,i);
                energy_check=-1.0;
            end
            if(isnan(RESE(j+1,i+1)))
                RESE(j+1,i+1)=0.0;%RESE=NaN means TQDT is 0.0, i.e., no energy loss
            end
        end
    end
    
    RESE(1,2:NIX+1)=RESE(2,2:NIX+1);
    RESE(NIY+2,2:NIX+1)=RESE(NIY+1,2:NIX+1);
    RESE(1:NIY+2,1)=RESE(1:NIY+2,2);
    RESE(1:NIY+2,NIX+2)=RESE(1:NIY+2,NIX+1);
    
    if(energy_check>0.0)
        fprintf('Energy Balanced!\n');
    end
    %---------------------------- ENERGY CHECK --------------------------------
end

%% =========================== DISCARD CODES ==============================
% %SITUATION 1: SUPERLIQUIDUS
%         if(TN(j+1,i+1)>TLTM(j+1,i+1))
%             dFSN(j+1,i+1)=dFSTM(j+1,i+1);%For superliquidus, we do nothing!
%             FLN(j+1,i+1)=FLTM(j+1,i+1);
%             CLN(j+1,i+1)=CLTM(j+1,i+1);
%             iters(j+1,i+1)=1;
%             continue%calculation below not nacessary!
%         end
%
%         %SITUATION 2: LIQUID --> MUSHY REGION
%         %0.0<=FSTM<FSE, TE<=TN<=TLTM --> CLTM<=CE, TE<=TN<=TLTM
%         if((CLTM(j+1,i+1)<CE)&&(TN(j+1,i+1)>=TE)&&(TN(j+1,i+1)<=TLTM(j+1,i+1)))
%             %Once some cell are subliquidus, we mark ingot T has been
%             %solved by T-FS-CL iteration. Even one cell!
%             modifymarker=modifymarker+1;
%
%             %----------------------------- STEP TWO -----------------------------------
%             %Take the approximations TN and CLTM(n) as the initial values of temperature and
%             %liquid composition at time n+1, and calculate the approximation of CLTM(n+1) using Eq. (9) and the approximate liquidus
%             %temperature at n+1 from Eq. (3).
%
%             %FL of Next step
%             FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSTM(j+1,i+1);
%
%             %NOTE: new T has effect on RL,RS,KP
%             %New RL updated by new T
%             RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLTM(j+1,i+1)*100.0-3.16e-7*(TN(j+1,i+1)-273.15);%liquid density of Next step from D-6 [g/mm^3]
%             RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
%
%             %New KP updated by new T
%             if(CLTM(j+1,i+1)<=CE)
%                 KPN(j+1,i+1)=0.12824+5.699124e-5*CLTM(j+1,i+1)*100.0+3.728777e-5*(CLTM(j+1,i+1)*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
%             else
%                 KPN(j+1,i+1)=1.0;
%             end
%
%             %New RS updated by new T
%             CSN(j+1,i+1)=CLTM(j+1,i+1)*KPN(j+1,i+1);
%             if(CSN(j+1,i+1)<CE)
%                 RSN(j+1,i+1)=2.58e-3;%[g/mm^3]
%             else
%                 RSN(j+1,i+1)=3.4e-3;%[g/mm^3]
%             end
%             RSN(j+1,i+1)=RSN(j+1,i+1)*10^6;%[kg/m^3]
%
%             %CL of Next step
%
%             %CLFN--Factor ahead CL(n+1): FL, RL, dFS, RS, KP should use NEW
%             %values FLN, RLN, RSN, KPN; dFSTM==dFSN==0.0 for the first step
%             CLFN(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1)+phi*dFSTM(j+1,i+1)*RSN(j+1,i+1)*KPN(j+1,i+1);%[kg/m^3]
%             CLN(j+1,i+1)=(ISS(j+1,i+1)+CLQVT(j,i)+CLQDT(j,i))/CLFN(j,i);%[kg/m^3 divided by kg/m^3 ==1]; ISS should update with NEW dFS ,but no NEW dFS available now, set old dFSTM==0.0 as new one
%             errCLTM(j+1,i+1)=abs((CLN(j+1,i+1)-CLTM(j+1,i+1))/CLN(j+1,i+1));
%
%             CLTM(j+1,i+1)=CLN(j+1,i+1);
%
%             TLN(j+1,i+1)=660.37-2.34581*CLN(j+1,i+1)*100.0-3.129e-2*(CLN(j+1,i+1)*100.0).^2+273.15;%new liquidus [K]
%
%             iters(j+1,i+1)=0;
%             %----------------------------- STEP TWO -----------------------------------
%             while(err(j+1,i+1))
%
%                 %-------------------------- STEP THREE --------------------------------
%                 %Calculate the approximation of dFSN by Eq. (10) and take the average:
%                 %dFSN=0.5*(dFSN+dFSTM)
%
%                 dFSN(j+1,i+1)=(FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+FRCPSTM(j+1,i+1))*(TLN(j+1,i+1)-TN(j+1,i+1))/(RSTM(j+1,i+1)*HSTM(j+1,i+1));%[1]
%
%                 dFSN(j+1,i+1)=0.5*(dFSN(j+1,i+1)+dFSTM(j+1,i+1));
%
%                 errdFS(j+1,i+1)=abs((dFSN(j+1,i+1)-dFSTM(j+1,i+1))/dFSN(j+1,i+1));
%
%                 dFSTM(j+1,i+1)=dFSN(j+1,i+1);
%                 %-------------------------- STEP THREE --------------------------------
%
%                 %--------------------------- STEP FOUR --------------------------------
%                 %Again using Eq. (8), calculate the new approximation of T(n+1) with new initial
%                 %values, dFSN and CLN, and take TN=0.5*(TN+TTM).
%
%                 %NOTE: new T, CL have effect on RL,RS,KP,CpS,CpL
%                 %Partition coefficient of Next step
%
%                 if(CLN(j+1,i+1)<=CE)
%                     KPN(j+1,i+1)=0.12824+5.699124e-5*CLN(j+1,i+1)*100.0+3.728777e-5*(CLN(j+1,i+1)*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
%                 else
%                     KPN(j+1,i+1)=1.0;
%                 end
%
%                 %RS of Next step
%                 CSN(j+1,i+1)=CLN(j+1,i+1)*KPN(j+1,i+1);
%                 if(CSN(j+1,i+1)<CE)
%                     RSN(j+1,i+1)=2.58e-3;%[g/mm^3]
%                 else
%                     RSN(j+1,i+1)=3.4e-3;%[g/mm^3]
%                 end
%                 RSN(j+1,i+1)=RSN(j+1,i+1)*10^6;%[kg/m^3]
%
%                 %CpS of Next step
%                 CPSN(j+1,i+1)=0.88+4.446e-4*(TN(j+1,i+1)-273.15)-2.274e-3*CSN(j+1,i+1)*100.0;%solid specific heat capacity from D-3 [J/g/K]
%                 CPSN(j+1,i+1)=CPSN(j+1,i+1)*1000.0;%[J/kg/K]
%
%                 %RL of Next step
%                 RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLN(j+1,i+1)*100.0-3.16e-7*(TN(j+1,i+1)-273.15);%liquid density of Next step from D-6 [g/mm^3]
%                 RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
%
%                 %CpL of Next step
%
%                 CPLN(j+1,i+1)=1.086-5.928e-3*CLN(j+1,i+1)*100.0;%liquid specific heat capacity from D-4 [J/g/K]
%                 CPLN(j+1,i+1)=CPLN(j+1,i+1)*1000.0;%[J/kg/K]
%
%
%                 RSCPSTM(j+1,i+1)=RSN(j+1,i+1)*CPSN(j+1,i+1);%[J/m^3/K]
%                 RLCPLTM(j+1,i+1)=RLN(j+1,i+1)*CPLN(j+1,i+1);%[J/m^3/K]
%
%
%                 TFN(j+1,i+1)=FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*dFSN(j+1,i+1);%[J/m^3/K]
%
%                 TN(j+1,i+1)=(IHS(j+1,i+1)+0.0*TQVT(j,i)+TQDT(j,i))/TFN(j+1,i+1);%[J/m^3 divided by J/m^3/K == K]
%                 TN(j+1,i+1)=0.5*(TN(j+1,i+1)+TTM(j+1,i+1));
%
%                 errTTM(j+1,i+1)=abs((TN(j+1,i+1)-TTM(j+1,i+1))/TN(j+1,i+1));
%
%                 TTM(j+1,i+1)=TN(j+1,i+1);
%                 %--------------------------- STEP FOUR --------------------------------
%
%                 %--------------------------- STEP FIVE --------------------------------
%                 %With the new initial value TN (as well as dFSN and CLN), calculate the new approximation
%                 %of CLN(n+1) using Eq. (9), and take CLN=0.5*(CLN+CLTM).
%
%                 %FL of of Next step
%                 FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSN(j+1,i+1);
%
%                 %NOTE: new T have effect on RL,CpS,CpL; RS,KP are functions of CLN
%                 %which is updated before so remain unchanged!
%
%                 %RL of Next step
%                 RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLN(j+1,i+1)*100.0-3.16e-7*(TN(j+1,i+1)-273.15);%liquid density of Next step from D-6 [g/mm^3]
%                 RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
%
%
%                 %New ISS due to updated dFSN
%
%                 ISS(j+1,i+1)=FRLTM(j+1,i+1)-(1.0-phi)*dFSN(j+1,i+1)*RSKPTM(j+1,i+1);%[kg/m^3]
%
%                 %Factor ahead CL(n+1)
%                 CLFN(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1)+phi*dFSN(j+1,i+1)*RSN(j+1,i+1)*KPN(j+1,i+1);%[kg/m^3]
%
%                 %CL of Next step
%                 CLN(j+1,i+1)=(ISS(j+1,i+1)+CLQVT(j,i)+CLQDT(j,i))/CLFN(j,i);%[kg/m^3 divided by kg/m^3 ==1]
%                 CLN(j+1,i+1)=0.5*(CLN(j+1,i+1)+CLTM(j+1,i+1));
%
%                 errCLTM(j+1,i+1)=abs((CLN(j+1,i+1)-CLTM(j+1,i+1))/CLN(j+1,i+1));
%                 CLTM(j+1,i+1)=CLN(j+1,i+1);
%                 %--------------------------- STEP FIVE --------------------------------
%
%                 %--------------------------- STEP SIX ---------------------------------
%                 err(j+1,i+1)=logical((err0dFS<=errdFS(j+1,i+1))||(err0TTM<=errTTM(j+1,i+1))||(err0CLTM<=errCLTM(j+1,i+1)));%logical variable for while loop
%
%                 %If the conditions are not true, take TLN(n+1)=TLN(CLN(n+1)) and repeat steps 3-6 until
%                 %all the above conditions are satisfied
%                 %Liquidus of Next step
%
%                 %--------------------------- STEP SIX ---------------------------------
%
%                 %--------------------------- STEP SEVEN -------------------------------
%                 %If jugments are not true, then repeat step 3-6 with updated TN, CLN,
%                 %dFSN
%
%                 TLN(j+1,i+1)=660.37-2.34581*CLN(j+1,i+1)*100.0-3.129e-2*(CLN(j+1,i+1)*100.0).^2+273.15;%ingot liquidus from D-16 [K]; CLN here not in wt% while in wt% from D-16
%
%                 %-------------------------- STEP SEVEN --------------------------------
%                 iters(j+1,i+1)=iters(j+1,i+1)+1;
%
%             end
%             continue%calculation below not nacessary!
%         end
%
%         %SITUATION 3: MUSHY REGION --> EUTECTIC
%
%         %FL of Next step
%         FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSTM(j+1,i+1);
%
%         %New RL updated by new T
%         RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLTM(j+1,i+1)*100.0-3.16e-7*(TN(j+1,i+1)-273.15);%liquid density of Next step from D-6 [g/mm^3]
%         RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
%
%         %New KP updated by new T
%         if(CLTM(j+1,i+1)<=CE)
%             KPN(j+1,i+1)=0.12824+5.699124e-5*CLTM(j+1,i+1)*100.0+3.728777e-5*(CLTM(j+1,i+1)*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
%         else
%             KPN(j+1,i+1)=1.0;
%         end
%
%         %New RS updated by new T
%         CSN(j+1,i+1)=CLTM(j+1,i+1)*KPN(j+1,i+1);
%         if(CSN(j+1,i+1)<CE)
%             RSN(j+1,i+1)=2.58e-3;%[g/mm^3]
%         else
%             RSN(j+1,i+1)=3.4e-3;%[g/mm^3]
%         end
%         RSN(j+1,i+1)=RSN(j+1,i+1)*10^6;%[kg/m^3]
%
%         CLFN(j,i)=FLN(j+1,i+1)*RLN(j+1,i+1)+phi*dFSTM(j+1,i+1)*RSN(j+1,i+1)*KPN(j+1,i+1);%[kg/m^3]
%         CLN(j+1,i+1)=(ISS(j+1,i+1)+CLQVT(j,i)+CLQDT(j,i))/CLFN(j,i);%[kg/m^3 divided by kg/m^3 ==1]; ISS should update with NEW dFS ,but no NEW dFS available now, set old dFSTM==0.0 as new one
%         TLN(j+1,i+1)=660.37-2.34581*CLN(j+1,i+1)*100.0-3.129e-2*(CLN(j+1,i+1)*100.0).^2+273.15;%ingot liquidus from D-16 [K]; CLN here not in wt% while in wt% from D-16
%
%         if((CLTM(j+1,i+1)<CE)&&(TLTM(j+1,i+1)>TE)&&(TN(j+1,i+1)<TE)&&(TLN(j+1,i+1)<TE))%T-old>TE
%             modifymarker=modifymarker+1;
%
%             TN(j+1,i+1)=TE;
%             CLN(j+1,i+1)=CE;
%             dFSCAL=(FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+FRCPSTM(j+1,i+1))*(TTME(j+1,i+1)-TN(j+1,i+1))/(RSE*HSE);%total solid corresponding to total heat loss if no solidification is considered
%             dFSBE=dFSCAL*(TTME(j+1,i+1)-TE)/(TTME(j+1,i+1)-TLN(j+1,i+1));%solid fraction before eutectic, Eq.(11)
%             dFSAE=dFSCAL*(TE-TLN(j+1,i+1))/(TTME(j+1,i+1)-TLN(j+1,i+1));%solid fraction after eutectic, Eq.(12)
%             dFSN(j+1,i+1)=dFSBE+dFSAE;
%
%             FSE=FSTM(j+1,i+1)+dFSBE;%real solid fraction when eutectic
%             frpintf('=============== FSE=%4.2f ==============',FSE);
%             iters(j+1,i+1)=1;
%
%             continue%calculation below not nacessary!
%         end
%
%         %SITUATION 4: EUTECTIC
%         if((FSTM(j+1,i+1)>=FSE)&&(FSTM(j+1,i+1)<=1.0)&&(TN(j+1,i+1)<TE))
%             modifymarker=modifymarker+1;
%
%             RLE=2.5222e-3+2.703e-5*CE*100.0-3.16e-7*(TE-273.15);%liquid density from D-6 [g/mm^3];  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
%             RLE=RLE*10^6;%[kg/m^3]
%
%             CPLE=1.086-5.928e-3*CE*100.0;%liquid specific heat capacity from D-4 [J/g/K]
%             CPLE=CPLE*1000.0;%[J/kg/K]
%
%             dFSN(j+1,i+1)=(FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1))*(TE-TN(j+1,i+1))/(RSE*HSE);%Eq.(14)
%             TN(j+1,i+1)=TE;
%             CLN(j+1,i+1)=CE;
%             iters(j+1,i+1)=1;
%
%             continue%calculation below not nacessary!
%         end
%
%         %SITUATION 5: EUTECTIC --> COMPLETE SOLID
%         dFSCAL=(FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+FRCPSTM(j+1,i+1))*(TTME(j+1,i+1)-TN(j+1,i+1))/(RSE*HSE);%total solid corresponding to total heat loss if no solidification is considered
%         if((FSTM(j+1,i+1)<1.0)&&(TN(j+1,i+1)<TE)&&(dFSCAL+FSTM(j+1,i+1)>1.0))
%             modifymarker=modifymarker+1;
%
%             dFSN(j+1,i+1)=1.0-FSTM(j+1,i+1);
%             %CLN(j+1,i+1) has no meaning
%             RLE=2.5222e-3+2.703e-5*CE*100.0-3.16e-7*(TE-273.15);%liquid density from D-6 [g/mm^3];  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
%             RLE=RLE*10^6;%[kg/m^3]
%             CPLE=1.086-5.928e-3*CE*100.0;%liquid specific heat capacity from D-4 [J/g/K]
%             CPLE=CPLE*1000.0;%[J/kg/K]
%
%             TN(j+1,i+1)=TN(j+1,i+1)+RSE*HSE*(1.0-FSTM(j+1,i+1))/(FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1));
%
%             iters(j+1,i+1)=1;
%         end