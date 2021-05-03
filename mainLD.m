%This script is used to reproduce figures in 'A Continuum Model for
%Computer Simulation of Macrosegregations in Ingots during solidification',
%Daming Xu 1989
%Created 2019-10-16

%Modified for lid driven cavity benchmark comparison at 2020-3-22

%======================= VERY IMPORTANT NOTE ==============================
%It is highly recommended the usage of packages for complex formulae, i.e.,
%temporary variables. These packages can be utilized repetedly, easily
%modified, and produce a concise numerical equation.
%======================= VERY IMPORTANT NOTE ==============================

%NOTE: only half volume is computed

clear all

%% ========================= INITIALIZATION ===============================
global NIX
global NIY
global dxI
global dyI
global x
global y
global g
global CE
global DL
global DS
global mu
global dx
global dy
global RL0
global lambda
global dtb
global P0

W=1.0;%ingot width [m]
H=1.0;%ingot height [m]
NIX=40;%number of unit in ingot, x-axis
NIY=40;%number of unit in ingot, y-axis
dxI=W/NIX;%%x-axis interval in ingot [m]
dyI=H/NIY;%y-axis interval in ingot [m]
dx(1:NIX)=dxI;%[m]
dy(1:NIY)=dyI;%[m]

%control point (main grid) x-axis location [m]
x=zeros(NIX+2,1);
x(1)=0.0;
x(2)=0.5*dxI;
x(3:NIX+1)=x(2)+[1:NIX-1]*dxI;
x(NIX+2)=x(NIX+1)+0.5*dxI;

%control point (main grid) y-axis location [m]
y=zeros(NIY+2,1);
y(1)=0.0;
y(2)=0.5*dyI;
y(3:NIY+1)=y(2)+[1:NIY-1]*dyI;
y(NIY+2)=y(NIY+1)+0.5*dyI;

T0=700.0+273.15;%Initial ingot T [K]
TM=26.0+273.15;%Initial T of mould [K]

CL0=0.045;%initial liquid alloy composition as Cu (==4.5 wt%) [1]
CE=0.332;%eutectic composition (33.2 wt%) [1]
DL=3.0e-11;%Cu diffusion coefficient in liquid Aluminum [m^2/sec]
DS=0.0;%Cu diffusion coefficient in solid matter [m^2/sec]

g=0.0;%Earth gravity [m/sec^2]
lambda=0.75;%under-relaxation coefficient: (1-lambda)*old+lambda*new

NS=1000;%total steps
t=zeros(NS+1,1);%all time sequence [sec]
dt=zeros(NS,1);%all time intervals between successive step [sec]

%----------------------- Main field variables -----------------------------
%T=T(t,x,y);
%CL=CL(t,x,y);
%CS=CS(t,x,y);%CS not time dependent in solid
%FL=FL(t,x,y);
%FS=FS(t,x,y);
%dFS=dFS(t,x,y)
%VX=VX(t,x,y);
%VY=VY(t,x,y);
%P=P(t,x,y);

%RD==RecorD

T=TM*ones(NIY+2,NIX+2);%T in ingot and mould [K]
T(1:NIY+1,1:NIX+1)=T0;%initial ingot T [K]
TRD=zeros(NIY+2,NIX+2,NS+1);%T records of all steps [K]
TRD(:,:,1)=T;%1st T record [K]
TI=zeros(NIY+2,NIX+2);%working variable of T [K]
TI(1:NIY+1,1:NIX+1)=T(1:NIY+1,1:NIX+1);
TI(NIY+2,:)=TI(NIY+1,:);
TI(:,NIX+2)=TI(:,NIX+1);

CL=CL0*ones(NIY+2,NIX+2);%Cu concentration in liquid [1]
CLRD=zeros(NIY+2,NIX+2,NS+1);%CL records of all steps [1]
CLRD(:,:,1)=CL;%1st record [1]
%CLI=CL;%working variable of CL [1]

CS=zeros(NIY+2,NIX+2);%Cu concentration in solid [1]
CSRD=zeros(NIY+2,NIX+2,NS+1);%CS records of all steps [1]
CSRD(:,:,1)=CS;%1st record [1]
%CSI=CS;%working variable of CS [1]

FL=ones(NIY+2,NIX+2);%liquid volume fraction [1]; only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
FLRD=zeros(NIY+2,NIX+2,NS+1);%FL records of all steps [1]
FLRD(:,:,1)=FL;%1st record [1]
%FLI=FL;%working variable of FL [1]

FS=zeros(NIY+2,NIX+2);%solid volume fraction [1]; only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
FSRD=zeros(NIY+2,NIX+2,NS+1);%FS records of all steps [1]
FSRD(:,:,1)=FS;%1st record [1]
%FSI=FS;%working variable of FS [1]

dFS=zeros(NIY+2,NIX+2);%solid volume fraction increment [1]
dFSRD=zeros(NIY+2,NIX+2,NS+1);%solid volume fraction increment record of all steps
dFSRD(:,:,1)=dFS;%1st record [1]
%dFSI=dFS;%working variable of dFS [1]

VX=zeros(NIY+2,NIX+1);%x-axis velocity [m/sec]
VX(1,2:NIX)=1.5;
VX(2,2:NIX)=0.5;
VX(NIY+2,:)=-VX(NIY+1,:);
VXRD=zeros(NIY+2,NIX+1,NS+1);%VX records of all steps [m/sec]
VXRD(:,:,1)=VX;%1st record [m/sec]
VXT=1.0;%top lid velocity [m/sec]

VY=zeros(NIY+1,NIX+2);%y-axis velocity [m/sec]
VY(:,1)=-VY(:,2);
VY(:,NIX+2)=-VY(:,NIX+1);
VYRD=zeros(NIY+1,NIX+2,NS+1);%VX records of all steps [m/sec]
VYRD(:,:,1)=VY;%1st record [m/sec]

PA=zeros(NIY+2,NIX+2);%pressure in liquid [Pa]
PR=zeros(NIY+2,NIX+2);%relative pressure in liquid with respect to initial pressure [Pa]
PRRD=zeros(NIY+2,NIX+2,NS+1);%Relative P records of all steps [Pa]
PRRD(:,:,1)=PR;%1st record [Pa]
%NOTE: Initialization after density initialization

% VX=xlsread('F:\HT\Benchmark\Re=1000_120.xlsx','VX','A1:DQ122');
% VY=xlsread('F:\HT\Benchmark\Re=1000_120.xlsx','VY','A1:DR121');
% PR=xlsread('F:\HT\Benchmark\Re=1000_120.xlsx','PR','A1:DR122');

%----------------------- Main field variables -----------------------------

%---------------------- Auxillary field variables -------------------------
%--------- Ingot properties -----------
%KS(T,CS)=KS(t,x,y)
%KL(T,CL)=KL(t,x,y)
%CpS(T,CS)=CpS(t,x,y)
%CpL(T,CL)=CpL(t,x,y)
%RS(T,CS)=RS(t,x,y)
%RL(T,CL)=RL(t,x,y)
%RM=RM(t,x,y)
%CM=CM(t,x,y)
%HS(CS)=HS(t,x,y)
%KP(CL)=KP(t,x,y)
%TL(CL)=TL(t,x,y)
%KFL(FL)=KFL(t,x,y)

%-------- Mould properties ------------
%KM(T)=KM(t,x,y)
%RCPM(T)=RCPM(t,x,y)

%-------- Other variables ------------
%KA system thermal conductivity
%NFSS=NFSS(x,y)=step at which solidification begins
%NFSE=NFSE(x,y)=step at which solidification ends

%Dynamic viscosity of liquid [Pa.sec]
mu=3.0e-3*ones(NIY+2,NIX+2);
% muRD=zeros(NIY+2,NIX+2,NS+1);
% muRD(:,:,1)=mu;

%Partition coefficient
KP=zeros(NIY+2,NIX+2);
for i=1:NIX+2
    for j=1:NIY+2
        if(CL(j,i)<CE)
            KP(j,i)=0.12824+5.699124e-5*CL(j,i)*100.0+3.728777e-5*(CL(j,i)*100.0)^2;%solid/liquid partition coefficient from D-15 [1]
        else
            KP(j,i)=1.0;
        end
    end
end

%CS are valid only if FS>0, otherwise CS should be NaN, however NaN can not
%be used in the following routine, so we use CL*KP
for i=1:NIX+2
    for j=1:NIY+2
        CS(j,i)=CL(j,i)*KP(j,i);
    end
end

KS=0.20983-2.081e-3*CS*100.0;%solid thermal conductivity from D-1 [W/mm/K]; KS=NaN if FS<=0
KS=KS*1000.0;%[W/m/K]; [NIY+2,NIX+2]
% KSRD=zeros(NIY+2,NIX+2,NS+1);%solid thermal conductivity records of all steps [W/m/K]
% KSRD(:,:,1)=KS;%1st record [W/m/K]

KL=6.5923e-2+3.3e-5*(TI-273.15)-6.807e-4*CL*100.0;%liquid thermal conductivity from D-2 [W/mm/K]
KL=1000.0*KL;%[W/m/K]; [NIY+2,NIX+2]
% KLRD=zeros(NIY+2,NIX+2,NS+1);%liquid thermal conductivity records of all steps [W/m/K]
% KLRD(:,:,1)=KL;%1st record [W/m/K]

CpS=0.88+4.446e-4*(TI-273.15)-2.274e-3*CS*100.0;%solid specific heat capacity from D-3 [J/g/K]
CpS=CpS*1000.0;%[J/kg/K]
% CpSRD=zeros(NIY+2,NIX+2,NS+1);%CpS records of all steps [J/kg/K]
% CpSRD(:,:,1)=CpS;%1st record [J/kg/K]

CpL=1.086-5.928e-3*CL*100.0;%liquid specific heat capacity from D-4 [J/g/K]
CpL=CpL*1000.0;%[J/kg/K]
% CpLRD=zeros(NIY+2,NIX+2,NS+1);%CpS records of all steps [J/kg/K]
% CpLRD(:,:,1)=CpL;%1st record [J/kg/K]

RS=zeros(NIY+2,NIX+2);%density of solid from D-5 [kg/m^3]; only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
for i=1:NIX+2
    for j=1:NIY+2
        if(CS(j,i)<CE)
            RS(j,i)=2.58e-3;%[g/mm^3]
        else
            RS(j,i)=3.4e-3;%[g/mm^3]
        end
    end
end
RS=RS*10^6;%[kg/m^3]
% RSRD=zeros(NIY+2,NIX+2,NS+1);%solid density records of all steps [kg/m^3]
% RSRD(:,:,1)=RS;%1st record [kg/m^3]

RL=2.5222e-3+2.703e-5*CL*100.0-3.16e-7*(TI-273.15);%liquid density from D-6 [g/mm^3]
RL=RL*10^6;%[kg/m^3]
% RLRD=zeros(NIY+2,NIX+2,NS+1);%solid density records of all steps [kg/m^3]
% RLRD(:,:,1)=RL;%1st record [kg/m^3]
RL0=RL;

RM=RL.*FL+RS*0.0;%mean density of ingot [kg/m^3];  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
RMRD=zeros(NIY+2,NIX+2,NS+1);%mean density records of all steps [kg/m^3]
RMRD(:,:,1)=RM;%1st record

%---------- Pressure Initialization ----------
PA(2,:)=0.5*RL(2,:).*FL(2,:)*g*dy(1);%1st row in ingot pressure [Pa]; A==absolute pressure
PA(1,:)=-PA(2,:);%this determines liquid surface is exactly 0.0 Pa [Pa]
for i=1:NIX+2
    for j=3:NIY+1
        PA(j,i)=PA(j-1,i)+0.5*g*dy(j-2)*RL(j-1,i)*FL(j-1,i)+0.5*g*dy(j-1)*RL(j,i)*FL(j,i);
    end
end
PA(NIY+2,:)=PA(NIY+1,:)+0.5*g*dy(NIY)*RL(NIY+1,:)+0.5*g*dy(NIY)*RL(NIY+2,:);%set dy(NIY+1)==dy(NIY)
%NOTE: this makes P at real physical boundary equal to
%P(NIY+1,:)+0.5*g*dy(NIY)*RL(NIY+1,:) since cells (NIY+1,:) are identical
%to cells (NIY+2,:)
P0=PA;
PR=PA-PA;%relative pressure in liquid with respect to initial pressure [Pa]
PRRD(:,:,1)=PR;%1st record [Pa]

CM=(FL.*RL.*CL+RS.*FS*0.0)./RM;%mean Cu concentration [1];  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
% CMRD=zeros(NIY+2,NIX+2,NS+1);%mean Cu concentration records of all steps [1]
% CMRD(:,:,1)=CM;%1st record

HS=zeros(NIY+2,NIX+2);%latent heat from D-7;  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
for i=1:NIX+2
    for j=1:NIY+2
        if(CS(j,i)<0.1)
            HS(j,i)=397.67-2.3288*CS(j,i)*100.0;%[J/g]
        else
            HS(j,i)=333.59;%[J/g]
        end
    end
end
HS=HS*1000.0;%[J/kg]
% HSRD=zeros(NIY+2,NIX+2,NS+1);%latent heat records of all steps [wt%]
% HSRD(:,:,1)=HS;%1st record

TL=660.37-2.34581*CL*100.0-3.129*0.01*(CL*100.0).^2+273.15;%ingot liquidus from D-16 [K]

KFL=zeros(NIY+2,NIX+2);%permeability from D-17;  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
for i=1:NIX+2
    for j=1:NIY+2
        if(FL(j,i)>=1.0/3.0)
            KFL(j,i)=2.6e-5*(1.923e-2*FL(j,i)^2+(4.0+3.0*FS(j,i)-3.0*sqrt(FS(j,i)*(8.0-3.0*FS(j,i))))/FS(j,i));%[mm^2]
        else
            KFL(j,i)=5.0e-7*FL(j,i)^2;%[mm^2]
        end
    end
end
KFL=KFL*10^-6;%[m^2]

NFSS=int16(zeros(NIY+2,NIX+2));%record step at which solid begins to form
NFSE=int16(zeros(NIY+2,NIX+2));%record step at which no more solid to form
%NOTE: these two variables above are used to cope with inhomogeneous Cu
%distribution in solid

FRCPS=zeros(NIY+2,NIX+2);%FS*RS*CpS
FRCS=zeros(NIY+2,NIX+2);%FS*RS*CS
FK=zeros(NIY+2,NIX+2);%FS*KS
FR=zeros(NIY+2,NIX+2);%FS*RS
FRL=zeros(NIY+2,NIX+2);%FL*RL
%FRH=zeros(NIY+2,NIX+2);%RS*HS(CS)
RSCPS=zeros(NIY+2,NIX+2);%RS*CPS, current step
RLCPL=zeros(NIY+2,NIX+2);%RL*CPL, current step
iters=zeros(NIY+2,NIX+2);
RESE=zeros(NIY+2,NIX+2);%residual of energy check
RESM=zeros(NIY+2,NIX+2);%residual of mass check

for i=1:NIX+2
    for j=1:NIY+2
        FRCPS(j,i)=FS(j,i)*RS(j,i)*CpS(j,i);
        FRCS(j,i)=FS(j,i)*RS(j,i)*CS(j,i);
        FK(j,i)=FS(j,i)*KS(j,i);
        FR(j,i)=FS(j,i)*RS(j,i);
        FRL(j,i)=FL(j,i)*RL(j,i);
    end
end

Re=10;%Reynolds number [1]
RL=RL./RL;
RL0=RL;
MUO=RL*VXT*W/Re;%dynamic viscosity of liquid [Pa.sec]
%NOTE: mu is [NIY+2,NIX+2]

%---------------------- Auxillary field variables -------------------------

%FieldPlots(T,FS,CL,RM,VX,VY,P);
% fidS=fopen('Solid.txt','w');
% fidE=fopen('Energy.txt','w');
% aviobj=VideoWriter('LidDriven.avi');
% open(aviobj);

%% ============================ Main loop =================================

for n=1:NS
    fprintf('Calculating Step: %5d   -->  %5d\n',n-1,n);
    
    %--------------------------- T-FS-CL ----------------------------------
    %To find step numbers which show solidification starts and ends, repectively
    %This would help if you want to know solidification location and process
%     if(n==1)
%         NFSS=1;
%         NFSE=1;
%     else
%         [NFSS,NFSE]=FindNSE(FSRD(:,:,n-1),FSRD(:,:,n),n);
%     end
    
    [FRCPS,FRCS,FK,FR,RSCPS,RLCPL] = INTSP(FRCPS,FRCS,FK,FR,CpS,CpL,CS,KS,RS,RL,dFS);%integrate solid properties: intg(FS*RS*CpS), intg(FS*RS*CS), intg(FS*KS), intg(FS*RS) [NIY+2,NIX+2]
    dtb=TSTEP(FRCPS,FL,RL,CpL,KL,FK,VX,VY,n);%estimate time step [sec]
    dt(n)=dtb;
    t(n+1)=t(n)+dt(n);
    
    %Prepare for VPM.m
    FLO=FL;
    RLO=RL;
    FRO=FR;
    RSO=RS;
    
    %[T,FL,CL,dFS,iters,RESE]=TFCM(T,FL,FS,KL,KS,CL,CS,VX,VY,FRCPS,FRCS,FK,RSCPS,RLCPL,RS,RL,CpL);%TFCM is original subroutine
    %fprintf('T-FS-CL Loop:%3d iterations!\n',max(max(iters)));
        
    %--------------------------- T-FS-CL ----------------------------------
    
%     %--------------- UPDATE PARAMETERS WITH NO CONVECTION -----------------
%     
%     %New FL, FS, KP (partition coefficient), CS
%     FS=1.0-FL;
%     for i=1:NIX+2
%         for j=1:NIY+2
%             if(CL(j,i)<CE)
%                 KP(j,i)=0.12824+5.699124e-5*CL(j,i)*100.0+3.728777e-5*(CL(j,i)*100.0)^2;%solid/liquid partition coefficient from D-15 [1]
%             else
%                 KP(j,i)=1.0;
%             end
%             CS(j,i)=CL(j,i)*KP(j,i);
%         end
%     end
%     
%     %Update RL,RS
%     for i=1:NIX+2
%         for j=1:NIY+2
%             if(CS(j,i)<CE)
%                 RS(j,i)=2.58e-3;%[g/mm^3]
%             else
%                 RS(j,i)=3.4e-3;%[g/mm^3]
%             end
%         end
%     end
%     RS=RS*10^6;%[kg/m^3]
%     
%     for i=2:NIX+1
%         for j=2:NIY+1
%             RL(j,i)=2.5222e-3+2.703e-5*CL(j,i)*100.0-3.16e-7*(T(j,i)-273.15);%liquid density from D-6 [g/mm^3]
%         end
%     end
%     RL(1,2:NIX+1)=RL(2,2:NIX+1);%top boundary
%     RL(NIY+2,2:NIX+1)=RL(NIY+1,2:NIX+1);%bottom boundary
%     RL(1:NIY+2,1)=RL(1:NIY+2,2);%left boudanry
%     RL(1:NIY+2,NIX+2)=RL(1:NIY+2,NIX+1);%right boundary
%     RL=RL*10^6;%[kg/m^3]
%     
%     RM=RL.*FL+RS.*FS;%mean density of ingot [kg/m^3]
%     RMRD(:,:,n+1)=RM;
%     
%     KS=0.20983-2.081e-3*CS*100.0;%solid thermal conductivity from D-1 [W/mm/K]; KS=NaN if FS<=0
%     KS=KS*1000.0;%[W/m/K]; [NIY+2,NIX+2]
%     
%     KL=6.5923e-2+3.3e-5*(T(1:NIY+2,1:NIX+2)-273.15)-6.807e-4*CL*100.0;%liquid thermal conductivity from D-2 [W/mm/K]
%     KL=1000.0*KL;%[W/m/K]; [NIY+2,NIX+2]
%     
%     CpS=0.88+4.446e-4*(T(1:NIY+2,1:NIX+2)-273.15)-2.274e-3*CS*100.0;%solid specific heat capacity from D-3 [J/g/K]
%     CpS=CpS*1000.0;%[J/kg/K]
%     
%     CpL=1.086-5.928e-3*CL*100.0;%liquid specific heat capacity from D-4 [J/g/K]
%     CpL=CpL*1000.0;%[J/kg/K]
%     
%     CSM=FRCS./FR;%old mean solid concentration [wt %]
%     
      MU=MUO;%Update dynamic viscosity [Pa.sec]
%
%     %--------------- UPDATE PARAMETERS WITH NO CONVECTION -----------------
    
    %---------------------------- VX-VY-P ---------------------------------
    [VX,VY,PR,RESM]=VPLD(VX,VY,PR,MUO,MU,FLO,FL,RLO,RL,RSO,RS,dFS);
    
    MUO=MU;
    %---------------------------- VX-VY-P ---------------------------------
    
    %------------------------ RECORD VARIABLES ----------------------------
    %Record main field variables
%     TRD(:,:,n+1)=T;
%     FLRD(:,:,n+1)=FL;
%     FSRD(:,:,n+1)=FS;
%     CLRD(:,:,n+1)=CL;
%     CSRD(:,:,n+1)=CS;
%     dFSRD(:,:,n)=dFS;
    VXRD(:,:,n+1)=VX;
    VYRD(:,:,n+1)=VY;
    PRRD(:,:,n+1)=PR;
    %----------------------- RECORD VARIABLES -----------------------------
    
    %----------------------- WRITE ON DISK --------------------------------
 %   PlotTFC(T,FS,CL,CSM,RM,VX,VY,PR,t(n+1),n);
%     fprintf(fidS,'Step= %s\n',num2str(n));
%     fprintf(fidS,'%6.4f',FS(2:NIY+1,NIX+1));
%     fprintf(fidS,'\n');
%     
%     fprintf(fidE,'Step= %s\n',num2str(n));
%     fprintf(fidE,'%6.4f',RESE);
%     fprintf(fidE,'\n');
    %----------------------- WRITE ON DISK --------------------------------

    %-------------------------- MAKE MOVIE --------------------------------
%    Current_Frame=getframe(gcf);
%    writeVideo(aviobj,Current_Frame);
    %-------------------------- MAKE MOVIE --------------------------------
end
% PlotTFC(T,FS,CL,CSM,RM,VX,VY,PR,t(n+1),n);
% fclose(fidS);
% fclose(fidE);
%close(aviobj);
xlswrite('F:\HT\Benchmark\Re=10_40.xlsx',VX,'VX');
xlswrite('F:\HT\Benchmark\Re=10_40.xlsx',VY,'VY');
xlswrite('F:\HT\Benchmark\Re=10_40.xlsx',PR,'PR');