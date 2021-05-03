%This routine is used to get through some simple issues in Prakash&Voller
%1989

%Created 2020-1-17
%Modified 2020-2-3
clear all

global NIX
global NIY
global x
global y
global g
global CE
global DL
global mu
global dx
global dy
global dtb
global TE

global L
global TW

global W
global H
global Tm

%% ================= CONSTANTS ==================
W=0.025;%width [m]
H=W;%height [m]
NIX=10;%x-axis grids number
NIY=10;%y-axis grids number
dx(1:NIX)=W/NIX;
dy(1:NIY)=H/NIY;
x=zeros(NIX+2,1);
x(1)=0.0;
x(2)=0.5*dx(1);
for i=3:NIX+1
    x(i)=x(i-1)+0.5*(dx(i-2)+dx(i-1));
end
x(NIX+2)=x(NIX+1)+0.5*dx(NIX);

y=zeros(NIY+2,1);
y(1)=0.0;
y(2)=0.5*dy(1);
for i=3:NIY+1
    y(i)=y(i-1)+0.5*(dy(i-2)+dy(i-1));
end
y(NIY+2)=y(NIY+1)+0.5*dy(NIY);

CE=0.8;%eutectic
DL=10^-6/250.0;%liquid diffusion coefficient [m^2/sec]
g=9.8;%gravity [m/sec^2]
L=300.0;%latent heat [J/kg]
HS0=0.0;%solid enthalpy [J/kg]
HL0=L;%liquid enthalpy [J/kg]
K0=5.0e-11;%permeability coefficient [m^2]
Tm=630.0;%melting point of pure substance [K]
TE=250.0;%eutectic temperature [K]
BT=4.0e-5;%thermal expansion [1/K]
BS=0.0;%solutal expansion [1]
T0=600.0;%initial temperature [K]
CL0=0.1;%initial species concentration [1]
TW=400.0;%wall temperature [K]
mu=10^-6;%kinematic viscosity [m^2/sec]
Tref=TE;
Cref=CE;

NS=100;%total steps
t=zeros(NS+1,1);%all time sequence [sec]
dt=zeros(NS,1);%all time intervals between successive step [sec]

%% ================ FIELD VARIABLES =================
KP=0.3*ones(NIY+2,NIX+2);%solid/liquid partition coefficient

FL=ones(NIY+2,NIX+2);%liquid volume fraction
FLRD=zeros(NIY+2,NIX+2,NS+1);%FL records of all steps [1]
FLRD(:,:,1)=FL;%1st record [1]

FS=1.0-FL;%solid volume fraction
FSRD=zeros(NIY+2,NIX+2,NS+1);%FS records of all steps [1]
FSRD(:,:,1)=FS;%1st record [1]

ML=ones(NIY+2,NIX+2);%liquid mass fraction
MLRD=zeros(NIY+2,NIX+2,NS+1);%FL records of all steps [1]
MLRD(:,:,1)=FL;%1st record [1]

MS=1.0-ML;%solid mass fraction
MSRD=zeros(NIY+2,NIX+2,NS+1);%FS records of all steps [1]
MSRD(:,:,1)=FS;%1st record [1]

dFS=zeros(NIY+2,NIX+2);%solid volume fraction increment [1]
dFSRD=zeros(NIY+2,NIX+2,NS+1);%solid volume fraction increment record of all steps
dFSRD(:,:,1)=dFS;%1st record [1]

CpL=3000.0*ones(NIY+2,NIX+2);%liquid specific heat capacity [J/kg/K]
CpLM=zeros(NIY+2,NIX+2);%temperature averaged CPL
CpS=3000.0*ones(NIY+2,NIX+2);%solid specific heat capacity [J/kg/K]
CpSM=zeros(NIY+2,NIX+2);%temperature averaged CPS
CpM=MS.*CpSM+ML.*CpLM;%temperature averaged CP [J/kg/K]

KL=0.4*ones(NIY+2,NIX+2);%liquid thermal conductivity [W/m/K]
KS=0.4*ones(NIY+2,NIX+2);%solid thermal conductivity [W/m/K]
KM=FL.*KL+FS.*KS;%average thermal conductivity [W/m/K]

RL=1000.0*ones(NIY+2,NIX+2);%liquid density [kg/m^3]
RS=1000.0*ones(NIY+2,NIX+2);%solid density [kg/m^3]
RM=FS.*RS+FL.*RL;%mixture density [kg/m^3]

T=T0*ones(NIY+2,NIX+2);
for i=1:NIY+2
    T(i,1)=2.0*TW-T(i,2);%left cold boundary
end
TRD=zeros(NIY+2,NIX+2,NS+1);%T records of all steps [K]
TRD(:,:,1)=T;%1st T record [K]

CL=CL0*ones(NIY+2,NIX+2);
CLRD=zeros(NIY+2,NIX+2,NS+1);%CL records of all steps [1]
CLRD(:,:,1)=CL;%1st record [1]

CS=KP.*CL;
CSRD=zeros(NIY+2,NIX+2,NS+1);%CS records of all steps [1]
CSRD(:,:,1)=CS;%1st record [1]

P=zeros(NIY+2,NIX+2);
P(2,:)=RL(2,:)*g*0.5*dy(1);
P(1,:)=-P(2,:);
for j=3:NIY+1
    P(j,:)=P(j-1,:)+0.5*g*RL(j-1,:)*dy(j-2)+0.5*g*RL(j,:)*dy(j-1);
end
P(NIY+2,:)=P(NIY+1,:)+0.5*g*RL(NIY+1,:)*dy(NIY)+0.5*g*RL(NIY+1,:)*dy(NIY);
PRD=zeros(NIY+2,NIX+2,NS+1);%P records of all steps [Pa]
PRD(:,:,1)=P;%1st record [Pa]

VX=zeros(NIY+2,NIX+1);%liquid x-axis velocity [m/sec]
VX(1,:)=-VX(2,:);%top boundary
VX(NIY+2,:)=-VX(NIY+1,:);%bottom boundary
VXRD=zeros(NIY+2,NIX+1,NS+1);%VX records of all steps [m/sec]
VXRD(:,:,1)=VX;%1st record [m/sec]

VY=zeros(NIY+1,NIX+2);%liquid y-axis velocity [m/sec]
VY(:,1)=-VY(:,2);%left boundary
VY(:,NIX+2)=-VY(:,NIX+1);%right boundary
VYRD=zeros(NIY+1,NIX+2,NS+1);%VX records of all steps [m/sec]
VYRD(:,:,1)=VY;%1st record [m/sec]

VSX=zeros(NIY+2,NIX+1);%solid x-axis velocity [m/sec]; useless
VSY=zeros(NIY+1,NIX+2);%solid y-axis velocity [m/sec]; useless
VXM=zeros(NIY+2,NIX+1);%mean x-axis velocity [m/sec]
VYM=zeros(NIY+1,NIX+2);%mean y-axis velocity [m/sec]

K=zeros(NIY+2,NIX+2);
for i=1:NIX+2
    for j=1:NIY+2
        K(j,i)=K0*FL(j,i)^3/(1.0-FL(j,i))^2;%permeability [m^2]
    end
end

CM=MS.*CS+ML.*CL;%mixture solute concentration [%]

FRCPS=zeros(NIY+2,NIX+2);%FS*RS*CpS
FRCS=zeros(NIY+2,NIX+2);%FS*RS*CS
FK=zeros(NIY+2,NIX+2);%FS*KS
FR=zeros(NIY+2,NIX+2);%FS*RS
RSCPS=zeros(NIY+2,NIX+2);%RS*CPS, current step
RLCPL=zeros(NIY+2,NIX+2);%RL*CPL, current step
iters=zeros(NIY+2,NIX+2);
RESE=zeros(NIY+2,NIX+2);%residual of energy check

for i=1:NIX+2
    for j=1:NIY+2
        FRCPS(j,i)=FS(j,i)*RS(j,i)*CpS(j,i);
        FRCS(j,i)=FS(j,i)*RS(j,i)*CS(j,i);
        FK(j,i)=FS(j,i)*KS(j,i);
        FR(j,i)=FS(j,i)*RS(j,i);
    end
end

fidS=fopen('Solid.txt','w');
fidE=fopen('Energy.txt','w');
aviobj=VideoWriter('MAGTC3.avi');
open(aviobj);

for n=1:NS

    VollerFields(T,FS,CM,ML,VX,VY,P,t(n+1),n);
    fprintf('Calculating Step: %3d   -->  %3d\n',n-1,n);
    
    [FRCPS,FRCS,FK,FR,RSCPS,RLCPL] = INTSP(FRCPS,FRCS,FK,FR,CpS,CpL,CS,KS,RS,RL,dFS);%integrate solid properties: intg(FS*RS*CpS), intg(FS*RS*CS), intg(FS*KS), intg(FS*RS) [NIY+2,NIX+2]
    
    dtb=0.7*TSTEP(FRCPS,FL,RL,CpL,KL,FK,VX,VY,n);%estimate time step [sec]
    dt(n)=dtb;
    t(n+1)=t(n)+dt(n);
    
    %Prepare for VollerVP.m
    FLO=FL;
    RLO=RL; 
    FRO=FR;
    RSO=RS;
    
    %T-F-C iteration
    [T,FL,CL,dFS,iters,RESE]=VollerTFC(T,FL,FS,KL,KS,CL,CS,VX,VY,FRCPS,FRCS,FK,RSCPS,RLCPL,RS,RSO,RL,CpL);
    fprintf('T-FS-CL Loop:%3d iterations!\n',max(max(iters)));
    
    %Update some variables
    FS=1.0-FL;
    RM=FR+FL.*RL;
    ML=FL.*RL./RM;
    CM=FRCS+0.5*(RS.*CL.*KP+RSO.*CS).*dFS+ML.*CL;%mixture solute concentration [%]
    
    %Velocity-Pressure coupling
    %[VX,VY,P,iters]=VollerVP(VX,VY,P,FLO,FL,RLO,RL,RSO,RS,dFS);
    
    %----------------------- RECORD VARIABLES ----------------------------
    %Record main field variables
    TRD(:,:,n+1)=T;
    FLRD(:,:,n+1)=FL;
    FSRD(:,:,n+1)=FS;
    CLRD(:,:,n+1)=CL;
    CSRD(:,:,n+1)=CS;
    dFSRD(:,:,n)=dFS;
    VXRD(:,:,n+1)=VX;
    VYRD(:,:,n+1)=VY;
    PRD(:,:,n+1)=P;
    %----------------------- RECORD VARIABLES ----------------------------
        VollerFields(T,FS,CM,ML,VX,VY,P,t(n+1),n);
        
    %----------------------- WRITE ON DISK --------------------------------
    fprintf(fidS,'Step= %s\n',num2str(n));
    fprintf(fidS,'%6.4f',FS(2:NIY+1,NIX+1));
    fprintf(fidS,'\n');
    
    fprintf(fidE,'Step= %s\n',num2str(n));
    fprintf(fidE,'%6.4f',RESE);
    fprintf(fidE,'\n');
    %----------------------- WRITE ON DISK --------------------------------

    %-------------------------- MAKE MOVIE --------------------------------
    Current_Frame=getframe(gcf);
    writeVideo(aviobj,Current_Frame);
    %-------------------------- MAKE MOVIE --------------------------------

end
fclose(fidS);
fclose(fidE);
close(aviobj);