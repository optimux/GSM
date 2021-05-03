function [T,FL,CL,dFS,iters,RESE]=VollerTFC(TTM,FLTM,FSTM,KLTM,KSTM,CLTM,CSTM,VXTM,VYTM,FRCPSTM,FRCSTM,FKTM,RSCPSTM,RLCPLTM,RSTM,RSOTM,RLTM,CPLTM)%FRCSTM,FRHTM,
%To iteratively solve energy and species conservation equations with finite
%volume method. Model, discretization and iteration can be found "A
%Continuum Model for Computer Simulation of Macrosegregations in Ingots
%during solidification" by Daming Xu 1989 and his companion papers in 1991

%Created 2019-10-17

%Modified for T-FS-CL iteration verification! 2019-11-10

%Modified for second paper 2020-1-1

%Modified for a complete solidification process 2020-1-7

%Modified for Prakash&Voller 1989 2020-2-4

%TTM: [NIY+2,NIX+2]
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
global DL
global CE
global dtb
global T0
global TE
global CL0
global FSE
global FSELOG
global TW
global Tm
global L

err0dFS=0.0001;%controlled accuracy for dFS, T, CL
err0TTM=0.0001;
err0CLTM=0.0001;

%======================= VERY IMPORTANT NOTE ==============================
%It is highly recommended the usage of packages for complex formulae, i.e.,
%temporary variables. These packages can be utilized repetedly, easily
%modified, and produce a concise numerical equation.
%======================= VERY IMPORTANT NOTE ==============================


%################################################################### ENERGY BALANCE #########################################################################
%% .............========= DIFFUSION HEAT FLUX ==========...................

%------------- thermal conductivity on finite volume faces ----------------
%....................... all faces, x-axis [W/m/K] ........................
KX=zeros(NIY+2,NIX+1);
for i=1:NIX+1
    for j=1:NIY+2
        KX(j,i)=0.5*(FKTM(j,i)+FKTM(j,i+1)+FLTM(j,i)*KLTM(j,i)+FLTM(j,i+1)*KLTM(j,i+1));%[(FS*KS)+(FL*KL)]|(j+-0.5,k)
    end
end

%........................all faces, y-axis [W/m/K] ........................
KY=zeros(NIY+1,NIX+2);
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
    DTX(i,1)=(TTM(i,2)-TTM(i,1))/(2.0*0.5*dx(1));%left cold boundary
    DTX(i,NIX+1)=(TTM(i,NIX+2)-TTM(i,NIX+1))/(2.0*0.5*dx(NIX));%right isolated boundary
end

DTY=zeros(NIY+1,NIX+2);%T gradient between finite volumes in y-axis
for i=1:NIX+2
    for j=2:NIY
        DTY(j,i)=(TTM(j+1,i)-TTM(j,i))/(0.5*dy(j-1)+0.5*dy(j));%[K/m]
    end
end

for i=1:NIX+2
    DTY(1,i)=(TTM(2,i)-TTM(1,i))/(2.0*0.5*dy(1));%top boundary
    DTY(NIY+1,i)=(TTM(NIY+2,i)-TTM(NIY+1,i))/(2.0*0.5*dy(NIY));%bottom boundary
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

%heat flux on all y faces
TQDY=zeros(NIY+1,NIX+2);%T=T, Q=flux, D==Diffusion, Y=y-axis
for i=1:NIX+2
    for j=1:NIY+1
        TQDY(j,i)=KY(j,i)*DTY(j,i);%[W/m^2]
    end
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
    %VY(1:NIY+1,1)=-VY(1:NIY+1,2)   (left boundary)
    %VY(1:NIY+1,NIX+2)=-VY(1:NIY+1,NIX+1)    (right boundary)
    for j=1:NIY+1
        %VY(1,1:NIX+2)=0.0 --> TRFVY(1,1:NIX+2)=0.0 (top boundary condition)
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
%################################################################### ENERGY BALANCE #########################################################################


%################################################################## SPECIES BALANCE #########################################################################

%% .............========= DIFFUSION SPECIES FLUX =========.................

%x-axis diffusion
CLQDX=zeros(NIY+2,NIX+1);%CL=CL, Q=flux, D=diffusion, X=x-axis
for i=2:NIX
    %CL(1:NIY+2,1)=CL(1:NIY+2,2) --> CLQDX(1:NIY+2,1)=0.0 (left boundary)
    %CL(1:NIY+2,NIX+1)=CL(1:NIY+2,NIX+2) --> CLQDX(1:NIY+2,NIX+1)=0.0 (right boundary)
    for j=1:NIY+2
        CLQDX(j,i)=2.0*dtb*DL*0.5*(FLTM(j,i)+FLTM(j,i+1))*0.5*(RLTM(j,i)+RLTM(j,i+1))*(CLTM(j,i+1)-CLTM(j,i))/(dx(i-1)+dx(i));%[kg/m^2]
    end
end

for j=1:NIY+2
    CLQDX(j,1)=2.0*dtb*DL*0.5*(FLTM(j,2)+FLTM(j,1))*0.5*(RLTM(j,2)+RLTM(j,1))*(CLTM(j,2)-CLTM(j,1))/(dx(1)+dx(1));%set dx(0)==dx(1); left boundary
    CLQDX(i,NIX+1)=2.0*dtb*DL*0.5*(FLTM(j,NIX+2)+FLTM(j,NIX+1))*0.5*(RLTM(j,NIX+2)+RLTM(j,NIX+1))*(CLTM(j,NIX+2)-CLTM(j,NIX+1))/(dx(NIX)+dx(NIX));%set dx(NIX+1)==dx(NIX); right boundary
end

%y-axis diffusion
CLQDY=zeros(NIY+1,NIX+2);%CL=CL, Q=flux, D=diffusion, Y=y-axis
for i=1:NIX+2
    for j=2:NIY
        %CL(1,1:NIX+2)=CL(2,1:NIX+2) --> CLQDY(1,1:NIX+2)=0.0
        %CL(NIY+2,1:NIX+2)=CL(NIY+1,1:NIX+2) --> CLQDY(NIY+1,1:NIX+2)=0.0
        CLQDY(j,i)=2.0*dtb*DL*0.5*(FLTM(j+1,i)+FLTM(j,i))*0.5*(RLTM(j+1,i)+RLTM(j,i))*(CLTM(j+1,i)-CLTM(j,i))/(dy(j-1)+dy(j));%[kg/m^2]
    end
end

for i=1:NIX+2
    CLQDY(1,i)=2.0*dtb*DL*0.5*(FLTM(2,i)+FLTM(1,i))*0.5*(RLTM(2,i)+RLTM(1,i))*(CLTM(2,i)-CLTM(1,i))/(dy(1)+dy(1));%set dy(0)==dy(1); top boundary
    CLQDY(NIY+1,i)=2.0*dtb*DL*0.5*(FLTM(NIY+2,i)+FLTM(NIY+1,i))*0.5*(RLTM(NIY+2,i)+RLTM(NIY+1,i))*(CLTM(NIY+2,i)-CLTM(NIY+1,i))/(dy(NIY)+dy(NIY));%set dy(NIY)==dy(NIY+1); bottom boundary
end

%total species flux of each finite volume
CLQDT=zeros(NIY,NIX);%CL=CL, Q=flux, D=diffusion, T=total
for i=1:NIX
    for j=1:NIY
        CLQDT(j,i)=(CLQDX(j+1,i+1)-CLQDX(j+1,i))/dx(i)+(CLQDY(j+1,i+1)-CLQDY(j,i+1))/(dy(j));%[kg/m^3]
    end
end

%% .............========= INTERNAL SPECIES STORAGE =========...............

KP=0.3*ones(NIY+2,NIX+2);%solid/liquid partition coefficient

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

iters=zeros(NIY+2,NIX+2);
TTME=TTM;%for eutectic T calculation

%TFN: T Factor of Next step
TFN=zeros(NIY+2,NIX+2);%T Factor of Next step
for i=1:NIX+2
    for j=1:NIY+2
        TFN(j,i)=FRCPSTM(j,i)+FLTM(j,i)*RLTM(j,i)*CPLTM(j,i)+(RSCPSTM(j,i)-RLCPLTM(j,i))*dFSTM(j,i);%[J/m^3/K]
    end
end

%----------------------------- STEP ONE -----------------------------------

%Let the initial value of dFSTM equal zero and calculate the approximation of TN using Eq. (8).
TN=zeros(NIY+2,NIX+2);%T of Next step
for j=2:NIX+1
    for k=2:NIY+1
        TN(k,j)=(IHS(k,j)+TQVT(k-1,j-1)+TQDT(k-1,j-1))/TFN(k,j);%[J/m^3 divided by J/m^3/K == K]
    end
end
%RSCPSTM: RS*CPS at n step, as approxiamation of RS*CPS at n+1 step
%RLCPLTM: RS*CPL at n step, as approxiamation of RS*CPL at n+1 step
TN(1,2:NIX+1)=TN(2,2:NIX+1);%top boundary
TN(NIY+2,2:NIX+1)=TN(NIY+1,2:NIX+1);%bottom boundary
TN(2:NIY+1,1)=2.0*TW-TN(2:NIY+1,2);%left boundary
TN(2:NIY+1,NIX+2)=TN(2:NIY+1,NIX+1);%right boundary

TN(1,1)=TN(2,1);%for better graphics
TN(NIY+2,1)=TN(NIY+1,1);%for better graphics
TN(1,NIX+2)=TN(2,NIX+2);%for better graphics
TN(NIY+2,NIX+2)=TN(NIY+1,NIX+2);%for better graphics

TTM=TN;
%----------------------------- STEP ONE -----------------------------------
RM=zeros(NIY+2,NIX+2);%mean density
MLTM=zeros(NIY+2,NIX+2);%old liquid mass fraction
TLTM=zeros(NIY+2,NIX+2);%Liquidus of old step
CM=zeros(NIY+2,NIX+2);%mean concentration

for i=1:NIX+2
    for j=1:NIY+2
        RM(j,i)=FRTM(j,i)+FLTM(j,i)*RLTM(j,i);%mean density [kg/m^3]
        MLTM(j,i)=FLTM(j,i)*RLTM(j,i)/RM(j,i);%old liquid mass fraction
        CM(j,i)=FRCSTM(j,i)+0.5*(RSTM(j,i)*CLTM(j,i)*KP(j,i)+RSOTM(j,i)*CSTM(j,i))*dFSTM(j,i)+MLTM(j,i)*CLTM(j,i);%mean concentration
        TLTM(j,i)=Tm+(TE-Tm)*CM(j,i)/CE;%old liquidus to determine superliquidus or subliquidus [K]
    end
end

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
HSN=zeros(NIY+2,NIX+2);%H of Next step
HSN=HSTM;

CPM=zeros(NIY+2,NIX+2);%mean specific hea capacity [J/kg/K]
MSTM=zeros(NIY+2,NIX+2);%solid mass fraction
MLTM=zeros(NIY+2,NIX+2);%old liquid mass fraction
MLN=zeros(NIY+2,NIX+2);%new liquid mass fraction
RM=zeros(NIY+2,NIX+2);%mean density [kg/m^3]
HSTM=L*ones(NIY+2,NIX+2);%latent heat of Next step

dFSBE=0.0;%dFS before eutectic
dFSAE=0.0;%dFS after eutectic
RSE=3400.0;%eutectic solid density [kg/m^3]
HSE=333.59*1000.0;%eutectic latent heat [J/kg]
RLE=2.5222e-3+2.703e-5*CE*100.0-3.16e-7*(TE-273.15);%liquid density from D-6 [g/mm^3];  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
RLE=RLE*10^6;%[kg/m^3]
CPLE=1.086-5.928e-3*CE*100.0;%liquid specific heat capacity from D-4 [J/g/K]
CPLE=CPLE*1000.0;%[J/kg/K]

for i=1:NIX
    for j=1:NIY
        
        %CASE 3: EUTECTIC
        cases=3;
        
        %CASE 1: SUPERLIQUIDUS
        if(TN(j+1,i+1)>=TLTM(j+1,i+1))
            cases=1;
        end
        
        %CASE 2: LIQUID --> MUSHY
        if((CLTM(j+1,i+1)<CE)&&(TN(j+1,i+1)>=TE)&&(TN(j+1,i+1)<=TLTM(j+1,i+1)))
            cases=2;
        end
        
        %CASE 4: COMPLETE SOLID
        if(FSTM(j+1,i+1)>=1.0)
            cases=4;
        end
        
        
        switch cases
            %SITUATION 1: SUPERLIQUIDUS
            case 1
                dFSN(j+1,i+1)=dFSTM(j+1,i+1);%For superliquidus, we do nothing!
                FLN(j+1,i+1)=FLTM(j+1,i+1);
                CLN(j+1,i+1)=CLTM(j+1,i+1);
                iters(j+1,i+1)=1;
                
                %SITUATION 2: LIQUID --> MUSHY REGION
            case 2
                %Once some cell are subliquidus, we mark ingot T has been
                %solved by T-FS-CL iteration. Even one cell!
                modifymarker=modifymarker+1;
                
                while(err(j+1,i+1))
                    RM(j+1,i+1)=FRTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1);%mean density [kg/m^3]
                    MLTM(j+1,i+1)=FLTM(j+1,i+1)*RLTM(j+1,i+1)/RM(j+1,i+1);%liquid mass fraction
                    MSTM(j+1,i+1)=1.0-MLTM(j+1,i+1);%solid mass fraction
                    CPM(j+1,i+1)=MLTM(j+1,i+1)*CPLTM(j+1,i+1)+MSTM(j+1,i+1)*
                    
                    a=1.0-KP(j+1,i+1);
                    b=KP(j+1,i+1)-MLTM(j+1,i+1)*(1.0-KP(j+1,i+1))+CPM(j+1,i+1)*(1.0-KP(j+1,i+1))*(Tm-TN(j+1,i+1))/HSTM(j+1,i+1);
                    d=CPM(j+1,i+1)*(TLTM(j+1,i+1)-(1.0-KP(j+1,i+1)*Tm))/HSTM(j+1,i+1)-KP(j+1,i+1)*(MLTM(j+1,i+1)+CPM(j+1,i+1)*TTME(j+1,i+1)/HSTM(j+1,i+1));%TTME should be replaced by TTM
                    
                    MLN(j+1,i+1)=(sqrt(b*b-4.0*a*d)-b)/(2.0*a);
                    if(MLN(j+1,i+1)>1.0)
                        MLN(j+1,i+1)=1.0;
                    end
                    if(MLN(j+1,i+1)<0.0)
                        MLN(j+1,i+1)=0.0;
                    end
                    
                    CLN(j+1,i+1)=
                    
                    errdFS(j+1,i+1)=
                    errTTM(j+1,i+1)=
                    errCLTM(j+1,i+1)=
                    
                    err(j+1,i+1)=logical((err0dFS<=errdFS(j+1,i+1))||(err0TTM<=errTTM(j+1,i+1))||(err0CLTM<=errCLTM(j+1,i+1)));%logical variable for while loop
                    iters(j+1,i+1)=iters(j+1,i+1)+1;
                    
                end
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
                
                TLN(j+1,i+1)=660.37-2.34581*CLN(j+1,i+1)*100.0-3.129e-2*(CLN(j+1,i+1)*100.0).^2+273.15;%new liquidus [K]
                
                iters(j+1,i+1)=0;
                %----------------------------- STEP TWO -----------------------------------
                while(err(j+1,i+1))
                    
                    %-------------------------- STEP THREE --------------------------------
                    %Calculate the approximation of dFSN by Eq. (10) and take the average:
                    %dFSN=0.5*(dFSN+dFSTM)
                    
                    dFSN(j+1,i+1)=(FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+FRCPSTM(j+1,i+1))*(TLN(j+1,i+1)-TN(j+1,i+1))/(RSTM(j+1,i+1)*HSN(j+1,i+1));%[1]
                    
                    dFSN(j+1,i+1)=0.5*(dFSN(j+1,i+1)+dFSTM(j+1,i+1));
                    
                    errdFS(j+1,i+1)=abs((dFSN(j+1,i+1)-dFSTM(j+1,i+1))/dFSN(j+1,i+1));
                    
                    dFSTM(j+1,i+1)=dFSN(j+1,i+1);
                    %-------------------------- STEP THREE --------------------------------
                    
                    %--------------------------- STEP FOUR --------------------------------
                    %Again using Eq. (8), calculate the new approximation of T(n+1) with new initial
                    %values, dFSN and CLN, and take TN=0.5*(TN+TTM).
                    
                    %NOTE: new T, CL have effect on RL,RS,KP,CpS,CpL
                    %Partition coefficient of Next step
                    
                    %NOTE: since CLTM<CE, thus KPN will follow the first
                    %statement, and there's no chance that the second statement
                    %will be executed!!!
                    if(CLN(j+1,i+1)<CE)
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
                    
                    if(CSN(j+1,i+1)<0.1)
                        HSN(j+1,i+1)=397.67-2.3288*CSN(j+1,i+1)*100.0;%latent heat from D-7 [J/g]
                    else
                        HSN(j+1,i+1)=333.59;%[J/g]
                    end
                    HSN(j+1,i+1)=HSN(j+1,i+1)*1000.0;%[J/kg]
                    
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
                    
                    TN(j+1,i+1)=(IHS(j+1,i+1)+0.0*TQVT(j,i)+TQDT(j,i))/TFN(j+1,i+1);%[J/m^3 divided by J/m^3/K == K]
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
                    
                    TLN(j+1,i+1)=660.37-2.34581*CLN(j+1,i+1)*100.0-3.129e-2*(CLN(j+1,i+1)*100.0).^2+273.15;%ingot liquidus from D-16 [K]; CLN here not in wt% while in wt% from D-16
                    
                    %-------------------------- STEP SEVEN --------------------------------
                    iters(j+1,i+1)=iters(j+1,i+1)+1;
                end
                
            case 4%SITUATION 4: EUTECTIC
                modifymarker=modifymarker+1;
                
                dFSN(j+1,i+1)=dFSCAL;
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSN(j+1,i+1);
                TN(j+1,i+1)=TE;
                CLN(j+1,i+1)=CE;
                iters(j+1,i+1)=1;
                
            case 5%SITUATION 5: EUTECTIC --> COMPLETE SOLID
                modifymarker=modifymarker+1;
                
                dFSN(j+1,i+1)=1.0-FSTM(j+1,i+1);
                FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSN(j+1,i+1);
                
                %CLN(j+1,i+1) has no meaning
                RLE=2.5222e-3+2.703e-5*CE*100.0-3.16e-7*(TE-273.15);%liquid density from D-6 [g/mm^3];  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
                RLE=RLE*10^6;%[kg/m^3]
                CPLE=1.086-5.928e-3*CE*100.0;%liquid specific heat capacity from D-4 [J/g/K]
                CPLE=CPLE*1000.0;%[J/kg/K]
                
                TN(j+1,i+1)=TN(j+1,i+1)+RSE*HSE*(1.0-FSTM(j+1,i+1))/(FLTM(j+1,i+1)*RLE*CPLE+FRCPSTM(j+1,i+1));
                CLN(j+1,i+1)=CE;
                iters(j+1,i+1)=1;
                
            case 6%COMPLETE SOLID COOLING
                
                FLN(j+1,i+1)=0.0;
                dFSN(j+1,i+1)=0.0;
                CLN(j+1,i+1)=CE;
                iters(j+1,i+1)=1;
                
            otherwise%SITUATION 3: MUSHY REGION --> EUTECTIC
                
                TNCAL=TN(j+1,i+1);
                dFSTMCAL=dFSTM(j+1,i+1);
                TTMCAL=TTM(j+1,i+1);
                %FL of Next step
                FLNCAL=FLTM(j+1,i+1)-dFSTMCAL;
                
                if(FLNCAL>0.0)%sometimes, FL == 0.0, then no need to calculate CLN, FLN
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
                    
                    TLNCAL=660.37-2.34581*CLNCAL*100.0-3.129e-2*(CLNCAL*100.0).^2+273.15;%new liquidus [K]
                    m=0;
                    %----------------------------- STEP TWO -----------------------------------
                    while(err(j+1,i+1))
                        m=m+1;
                        
                        %-------------------------- STEP THREE --------------------------------
                        %Calculate the approximation of dFSN by Eq. (10) and take the average:
                        %dFSN=0.5*(dFSN+dFSTM)
                        
                        dFSNCAL=(FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+FRCPSTM(j+1,i+1))*(TLNCAL-TNCAL)/(RSTM(j+1,i+1)*HSN(j+1,i+1));%[1]
                        
                        dFSNCAL=0.5*(dFSNCAL+dFSTMCAL);
                        
                        errdFS(j+1,i+1)=abs((dFSNCAL-dFSTMCAL)/dFSNCAL);
                        
                        dFSTMCAL=dFSNCAL;
                        %-------------------------- STEP THREE --------------------------------
                        
                        %--------------------------- STEP FOUR --------------------------------
                        %Again using Eq. (8), calculate the new approximation of T(n+1) with new initial
                        %values, dFSN and CLN, and take TN=0.5*(TN+TTM).
                        
                        %NOTE: new T, CL have effect on RL,RS,KP,CpS,CpL
                        %Partition coefficient of Next step
                        KPN(j+1,i+1)=0.12824+5.699124e-5*CLNCAL*100.0+3.728777e-5*(CLNCAL*100.0)^2;%solid/liquid partition coefficientfrom D-15 [1]
                        
                        %RS of Next step
                        CSN(j+1,i+1)=CLNCAL*KPN(j+1,i+1);
                        if(CSN(j+1,i+1)<CE)
                            RSN(j+1,i+1)=2.58e-3;%[g/mm^3]
                        else
                            RSN(j+1,i+1)=3.4e-3;%[g/mm^3]
                        end
                        RSN(j+1,i+1)=RSN(j+1,i+1)*10^6;%[kg/m^3]
                        
                        if(CSN(j+1,i+1)<0.1)
                            HSN(j+1,i+1)=397.67-2.3288*CSN(j+1,i+1)*100.0;%latent heat from D-7 [J/g]
                        else
                            HSN(j+1,i+1)=333.59;%[J/g]
                        end
                        HSN(j+1,i+1)=HSN(j+1,i+1)*1000.0;%[J/kg]
                        
                        %CpS of Next step
                        CPSN(j+1,i+1)=0.88+4.446e-4*(TNCAL-273.15)-2.274e-3*CSN(j+1,i+1)*100.0;%solid specific heat capacity from D-3 [J/g/K]
                        CPSN(j+1,i+1)=CPSN(j+1,i+1)*1000.0;%[J/kg/K]
                        
                        %RL of Next step
                        RLN(j+1,i+1)=2.5222e-3+2.703e-5*CLNCAL*100.0-3.16e-7*(TNCAL-273.15);%liquid density of Next step from D-6 [g/mm^3]
                        RLN(j+1,i+1)=RLN(j+1,i+1)*10^6;%[kg/m^3]
                        
                        %CpL of Next step
                        
                        CPLN(j+1,i+1)=1.086-5.928e-3*CLNCAL*100.0;%liquid specific heat capacity from D-4 [J/g/K]
                        CPLN(j+1,i+1)=CPLN(j+1,i+1)*1000.0;%[J/kg/K]
                        
                        
                        RSCPSTM(j+1,i+1)=RSN(j+1,i+1)*CPSN(j+1,i+1);%[J/m^3/K]
                        RLCPLTM(j+1,i+1)=RLN(j+1,i+1)*CPLN(j+1,i+1);%[J/m^3/K]
                        
                        
                        TFNCAL=FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1)+(RSCPSTM(j+1,i+1)-RLCPLTM(j+1,i+1))*dFSNCAL;%[J/m^3/K]
                        
                        TNCAL=(IHS(j+1,i+1)+TQVT(j,i)+TQDT(j,i))/TFNCAL;%[J/m^3 divided by J/m^3/K == K]
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
                        
                        TLNCAL=660.37-2.34581*CLNCAL*100.0-3.129e-2*(CLNCAL*100.0).^2+273.15;%ingot liquidus from D-16 [K]; CLN here not in wt% while in wt% from D-16
                        %-------------------------- STEP SEVEN --------------------------------
                    end
                    
                    dFSCAL=dFSNCAL;
                    modifymarker=modifymarker+1;
                    
                    dFSBE=dFSCAL*(TTME(j+1,i+1)-TE)/(TTME(j+1,i+1)-TLNCAL);%solid fraction before eutectic, Eq.(11)
                    % dFSAE1=(FRCPSTM(j+1,i+1)+FLTM(j+1,i+1)*RLTM(j+1,i+1)*CPLTM(j+1,i+1))*(TTME(j+1,i+1)-TN(j+1,i+1))*(TE-TLNCAL)/((TTME(j+1,i+1)-TLNCAL)*RSE*HSE);%solid fraction after eutectic, Eq.(12)
                    dFSAE=dFSCAL*(TE-TLNCAL)*RSTM(j+1,i+1)*HSTM(j+1,i+1)/((TTME(j+1,i+1)-TLNCAL)*RSE*HSE);
                    dFSN(j+1,i+1)=dFSBE+dFSAE;
                    FLN(j+1,i+1)=FLTM(j+1,i+1)-dFSN(j+1,i+1);
                    TN(j+1,i+1)=TE;
                    CLN(j+1,i+1)=CE;
                    
                    if(FSELOG(j+1,i+1))%FSE only modify once
                        FSE(j+1,i+1)=FSTM(j+1,i+1)+dFSBE;%real solid fraction when eutectic
                        FSELOG(j+1,i+1)=0;
                    end
                    
                    fprintf('=============== FSE= %7.5f ==============\n',FSE(j+1,i+1));
                    iters(j+1,i+1)=1;
                    
                else
                    fprintf('Complete solid BUT should be eutectic\n');
                end
                
        end
        
    end
end
%################################################################ T-FS-CL ITERATION #########################################################################


%################################################################ UPDATE VARIABLES ##########################################################################
%T-FS-CL has been iteratively sovled?
if(modifymarker>0.0)%Yes
    
    %Returen new field variables
    TN(1,2:NIX+1)=TN(2,2:NIX+1);%top boundary
    TN(NIY+2,2:NIX+1)=TN(NIY+1,2:NIX+1);%bottom boundary
    TN(2:NIY+1,1)=2.0*TW-TN(2:NIY+1,2);%left boundary
    TN(2:NIY+1,NIX+2)=TN(2:NIY+1,NIX+1);%right boundary
    
    TN(1,1)=TN(2,1);%for better graphics
    TN(NIY+2,1)=TN(NIY+1,1);%for better graphics
    TN(1,NIX+2)=TN(2,NIX+2);%for better graphics
    TN(NIY+2,NIX+2)=TN(NIY+1,NIX+2);%for better graphics
    T=TN;
    
    FLN(1,2:NIX+1)=FLN(2,2:NIX+1);%top boundary
    FLN(NIY+2,2:NIX+1)=FLN(NIY+1,2:NIX+1);%bottom boundary
    FLN(2:NIY+1,1)=FLN(2:NIY+1,2);%left boundary
    FLN(2:NIY+1,NIX+2)=FLN(2:NIY+1,NIX+1);%right boundary
    FLN(1,1)=FLN(2,1);
    FLN(1,NIX+2)=FLN(2,NIX+2);
    FLN(NIY+2,1)=FLN(NIY+1,1);
    FLN(NIY+2,NIX+2)=FLN(NIY+1,NIX+2);
    FL=FLN;
    
    CLN(1,2:NIX+1)=CLN(2,2:NIX+1);%top boundary
    CLN(NIY+2,2:NIX+1)=CLN(NIY+1,2:NIX+1);%bottom boundary
    CLN(2:NIY+1,1)=CLN(2:NIY+1,2);%left boundary
    CLN(2:NIY+1,NIX+2)=CLN(2:NIY+1,NIX+1);%right boundary
    CLN(1,1)=CLN(2,1);
    CLN(1,NIX+2)=CLN(2,NIX+2);
    CLN(NIY+2,1)=CLN(NIY+1,1);
    CLN(NIY+2,NIX+2)=CLN(NIY+1,NIX+2);
    CL=CLN;
    
    dFSN(1,2:NIX+1)=dFSN(2,2:NIX+1);%top boundary
    dFSN(NIY+2,2:NIX+1)=dFSN(NIY+1,2:NIX+1);%bottom boundary
    dFSN(2:NIY+1,1)=dFSN(2:NIY+1,2);%left boundary
    dFSN(2:NIY+1,NIX+2)=dFSN(2:NIY+1,NIX+1);%right boundary
    dFSN(1,1)=dFSN(2,1);
    dFSN(1,NIX+2)=dFSN(2,NIX+2);
    dFSN(NIY+2,1)=dFSN(NIY+1,1);
    dFSN(NIY+2,NIX+2)=dFSN(NIY+1,NIX+2);
    dFS=dFSN;
    
else%No, Return updated T and constant FL,CL
    
    TN(1,2:NIX+1)=TN(2,2:NIX+1);%top boundary
    TN(NIY+2,2:NIX+1)=TN(NIY+1,2:NIX+1);%bottom boundary
    TN(2:NIY+1,1)=2.0*TW-TN(2:NIY+1,2);%left boundary
    TN(2:NIY+1,NIX+2)=TN(2:NIY+1,NIX+1);%right boundary
    
    TN(1,1)=TN(2,1);%for better graphics
    TN(NIY+2,1)=TN(NIY+1,1);%for better graphics
    TN(1,NIX+2)=TN(2,NIX+2);%for better graphics
    TN(NIY+2,NIX+2)=TN(NIY+1,NIX+2);%for better graphics
    T=TN;
    
    FL=FLTM;
    CL=CLTM;
    dFS=dFSTM;
    
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