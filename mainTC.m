%This script is used to reproduce figures in 'A Continuum Model for
%Computer Simulation of Macrosegregations in Ingots during solidification',
%Daming Xu 1989
%Created 2019-10-16

%Modified from Daming.m created at 2019-9-23

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
global NMX
global NMY
global dxI
global dxM
global dyI
global dyM
global x
global y
global g
global CE
global DL
global DS
global mu
global dx
global dy
global gap
global TR

W=100.0/1000.0;%ingot width [m]
H=120.0/1000.0;%ingot height [m]
WM=240.0/1000.0;%mould width [m]
HM=190.0/1000.0;%mould height [m]
NIX=8;%number of unit in ingot, x-axis
NIY=12;%number of unit in ingot, y-axis
NMX=7;%number of unit in mould, x-axis; in Fig. 3-6, 6 units in mould
NMY=7;%number of unit in mould, y-axis
dxI=0.5*W/NIX;%%x-axis interval in ingot [m]
dyI=H/NIY;%y-axis interval in ingot [m]
dxM=0.5*(WM-W)/NMX;%x-axis interval in mould [m]
dyM=(HM-H)/NMY;%y-axis interval in mould [m]
gap=0.1/1000.0;%gas gap between ingot and mould [m]
dx(1:NIX)=dxI;%[m]
dx(NIX+1:NIX+NMX)=dxM;%[m]
dy(1:NIY)=dyI;%[m]
dy(NIY+1:NIY+NMY)=dyM;%[m]

%control point (main grid) x-axis location [m]
x=zeros(NIX+2+NMX+1,1);
x(1)=0.0;
x(2)=0.5*dxI;
x(3:NIX+1)=x(2)+[1:NIX-1]*dxI;
x(NIX+2)=x(NIX+1)+0.5*dxI;
x(NIX+3)=x(NIX+2)+0.5*dxM+gap;
x(NIX+4:NIX+2+NMX)=x(NIX+3)+[1:NMX-1]*dxM;
x(NIX+3+NMX)=x(NIX+2+NMX)+0.5*dxM;

%control point (main grid) y-axis location [m]
y=zeros(NIY+2+NMY+1,1);
y(1)=0.0;
y(2)=0.5*dyI;
y(3:NIY+1)=y(2)+[1:NIY-1]*dyI;
y(NIY+2)=y(NIY+1)+0.5*dyI;
y(NIY+3)=y(NIY+2)+0.5*dyM+gap;
y(NIY+4:NIY+2+NMY)=y(NIY+3)+[1:NMY-1]*dyM;
y(NIY+2+NMY+1)=y(NIY+2+NMY)+0.5*dyM;

T0=700.0+273.15;%Initial ingot T [K]
TR=25.0+273.15;%T outside of mould [K]
TM=26.0+273.15;%Initial T of mould [K]

CL0=0.045;%initial liquid alloy composition as Cu (==4.5 wt%) [1]
CE=0.332;%eutectic composition (33.2 wt%) [1]
DL=3.0e-11;%Cu diffusion coefficient in liquid Aluminum [m^2/sec]
DS=0.0;%Cu diffusion coefficient in solid matter [m^2/sec]
mu=3.0e-3;%dynamic viscosity of liquid [Pa.sec]

g=9.81;%Earth gravity [m/sec^2]

NS=100;%total steps
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

T=TM*ones(NIY+NMY+2,NIX+NMX+2);%T in ingot and mould [K]
T(1:NIY+1,1:NIX+1)=T0;%initial ingot T [K]
TRD=zeros(NIY+NMY+2,NIX+NMX+2,NS+1);%T records of all steps [K]
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
%VX=-0.001*ones(NIY+1,NIX+1);%x-axis velocity [m/sec]
VXRD=zeros(NIY+2,NIX+1,NS+1);%VX records of all steps [m/sec]
VXRD(:,:,1)=VX;%1st record [m/sec]

VY=zeros(NIY+1,NIX+2);%y-axis velocity [m/sec]
%VY=0.001*ones(NIY+1,NIX+1);%y-axis velocity [m/sec]
VYRD=zeros(NIY+1,NIX+2,NS+1);%VX records of all steps [m/sec]
VYRD(:,:,1)=VY;%1st record [m/sec]

P=zeros(NIY+2,NIX+2);%pressure in liquid [Pa]; only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
PRD=zeros(NIY+2,NIX+2,NS+1);%P records of all steps [Pa]
PRD(:,:,1)=P;%1st record [Pa]
%NOTE: Initialization after density initialization

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
%HL(CS)=HL(t,x,y)
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

%Partition coefficient
KP=zeros(NIY+2,NIX+2);
for i=1:NIX+2
    for j=1:NIY+2
        if(CL(j,i)~=CE)
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

RL=2.5222e-3+2.703e-5*CL*100.0-3.16e-7*(TI-273.15);%liquid density from D-6 [g/mm^3];  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
RL=RL*10^6;%[kg/m^3]
% RLRD=zeros(NIY+2,NIX+2,NS+1);%solid density records of all steps [kg/m^3]
% RLRD(:,:,1)=RL;%1st record [kg/m^3]

RM=RL.*FL+RS*0.0;%mean density of ingot [kg/m^3];  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
RMRD=zeros(NIY+2,NIX+2,NS+1);%mean density records of all steps [kg/m^3]
RMRD(:,:,1)=RM;%1st record

%---------- Pressure: 32 digits precision ----------
% P(1,:)=0.0;%liquid surface [Pa]
% P(2,:)=vpa(0.5,32)*vpa(RL(2,:),32).*vpa(FL(2,:),32)*vpa(g,32)*vpa(dy(1),32);%1st row pressure [Pa]
% dP=vpa(0.5,32)*vpa(g,32)*vpa(dy(1),32)*vpa(RL(1,1),32)*vpa(FL(1,1),32)+vpa(0.5,32)*vpa(g,32)*vpa(dy(1),32)*vpa(RL(1,1),32)*vpa(FL(1,1),32);
% for i=1:NIX+2
%     for j=3:NIY+1
%         P(j,i)=vpa(P(j-1,i),32)+vpa(0.5,32)*vpa(g,32)*vpa(dy(j-2),32)*vpa(RL(j-1,i),32)*vpa(FL(j-1,i),32)+vpa(0.5,32)*vpa(g,32)*vpa(dy(j-1),32)*vpa(RL(j,i),32)*vpa(FL(j,i),32);
%     end
% end
% P(NIY+2,:)=vpa(P(NIY+1,:),32)+vpa(0.5,32)*vpa(g,32)*vpa(dy(NIY),32)*vpa(RL(NIY+1),32);
% PRD(:,:,1)=P;%1st record [Pa]

%---------- Pressure: 13 digits precision is enough ----------
% P(1,:)=0.0;%liquid surface [Pa]
% P(2,:)=vpa(0.5,13)*vpa(RL(2,:),13).*vpa(FL(2,:),13)*vpa(g,13)*vpa(dy(1),13);%1st row pressure [Pa]
% dP=vpa(0.5,13)*vpa(g,13)*vpa(dy(1),13)*vpa(RL(1,1),13)*vpa(FL(1,1),13)+vpa(0.5,13)*vpa(g,13)*vpa(dy(1),13)*vpa(RL(1,1),13)*vpa(FL(1,1),13);
% for i=1:NIX+2
%     for j=3:NIY+1
%         P(j,i)=vpa(P(j-1,i),13)+vpa(0.5,13)*vpa(g,13)*vpa(dy(j-2),13)*vpa(RL(j-1,i),13)*vpa(FL(j-1,i),13)+vpa(0.5,13)*vpa(g,13)*vpa(dy(j-1),13)*vpa(RL(j,i),13)*vpa(FL(j,i),13);
%     end
% end
% P(NIY+2,:)=vpa(P(NIY+1,:),13)+vpa(0.5,13)*vpa(g,13)*vpa(dy(NIY),13)*vpa(RL(NIY+1),13);
% PRD(:,:,1)=P;%1st record [Pa]

%---------- Pressure: with default precision ----------
P(1,:)=0.0;%liquid surface [Pa]
P(2,:)=0.5*RL(2,:).*FL(2,:)*g*dy(1);%1st row pressure [Pa]
for i=1:NIX+2
    for j=3:NIY+1
        P(j,i)=P(j-1,i)+0.5*g*dy(j-2)*RL(j-1,i)*FL(j-1,i)+0.5*g*dy(j-1)*RL(j,i)*FL(j,i);
    end
end
P(NIY+2,:)=P(NIY+1,:)+0.5*g*dy(NIY)*RL(NIY+1);
PRD(:,:,1)=P;%1st record [Pa]

CM=(FL.*RL.*CL+RS.*FS*0.0)./RM;%mean Cu concentration [1];  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
% CMRD=zeros(NIY+2,NIX+2,NS+1);%mean Cu concentration records of all steps [1]
% CMRD(:,:,1)=CM;%1st record

HL=zeros(NIY+2,NIX+2);%latent heat from D-7;  only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
for i=1:NIX+2
    for j=1:NIY+2
        if(CS(j,i)<0.1)
            HL(j,i)=397.67-2.3288*CS(j,i)*100.0;%[J/g]
        else
            HL(j,i)=333.59;%[J/g]
        end
    end
end
HL=HL*1000.0;%[J/kg]
% HLRD=zeros(NIY+2,NIX+2,NS+1);%latent heat records of all steps [wt%]
% HLRD(:,:,1)=HL;%1st record

TL=660.37-2.34581*CL*100.0-3.129e-2*(CL*100.0).^2+273.15;%ingot liquidus from D-16 [K]

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
%FRH=zeros(NIY+2,NIX+2);%RS*HL(CS)
RSCPS=zeros(NIY+2,NIX+2);%RS*CPS, current step
RLCPL=zeros(NIY+2,NIX+2);%RL*CPL, current step
iter=0;

for i=1:NIX+2
    for j=1:NIY+2
        FRCPS(j,i)=FS(j,i)*RS(j,i)*CpS(j,i);
        FRCS(j,i)=FS(j,i)*RS(j,i)*CS(j,i);
        FK(j,i)=FS(j,i)*KS(j,i);
        FR(j,i)=FS(j,i)*RS(j,i);
        FRL(j,i)=FL(j,i)*RL(j,i);
    end
end

%---------------------- Auxillary field variables -------------------------

%FieldPlots(T,FS,CL,RM,VX,VY,P);
%print -depsc2 tt.eps

%% ============================ Main loop =================================
for n=1:NS
    fprintf('Calculating Step: %3d   -->  %3d\n',n-1,n);
    
    %--------------------------- T-FS-CL ----------------------------------
    
    %To find step numbers which show solidification starts and ends, repectively
    %This would help if you want to know solidification location and process
    if(n==1)
        NFSS=1;
        NFSE=1;
    else
        [NFSS,NFSE]=FindNSE(FSRD(:,:,n-1),FSRD(:,:,n),n);%(dFSRD(:,:,n-1),dFSRD(:,:,n),n);
    end
    
    [FRCPS,FRCS,FK,FR,RSCPS,RLCPL] = INTSP(FRCPS,FRCS,FK,FR,CpS,CpL,CS,KS,RS,RL,dFS);%integrate solid properties: intg(FS*RS*CpS), intg(FS*RS*CS), intg(FS*KS), intg(FS*RS) [NIY+2,NIX+2]
    TimeInterval=TSTEP(FRCPS,FL,RL,CpL,KL,FK,VX,VY,n);%estimate time step [sec]
    dt(n)=TimeInterval;
    t(n+1)=t(n)+dt(n);   
    
    [T,FL,CL,dFS,iter]=TFC(T,FL,FS,KL,KS,CL,CS,VX,VY,FRCPS,FRCS,FK,RSCPS,RLCPL,RS,RL,CpL,TimeInterval);
    fprintf('T-FS-CL Loop:%3d iterations!\n',iter);

    %--------------------------- T-FS-CL ----------------------------------
    
    %----------------------- UPDATE PARAMETERS ----------------------------
    
    %New FL, FS, KP (partition coefficient), CS
    FS=1.0-FL;
    for i=1:NIX+2
        for j=1:NIY+2
            if(CL(j,i)~=CE)
                KP(j,i)=0.12824+5.699124e-5*CL(j,i)*100.0+3.728777e-5*(CL(j,i)*100.0)^2;%solid/liquid partition coefficient from D-15 [1]
            else
                KP(j,i)=1.0;
            end
            CS(j,i)=CL(j,i)*KP(j,i);
        end
    end
    
    %Update RL,RS
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
    
    RL=2.5222e-3+2.703e-5*CL*100.0-3.16e-7*(T(1:NIY+2,1:NIX+2)-273.15);%liquid density from D-6 [g/mm^3]
    RL=RL*10^6;%[kg/m^3]
    
    RMRD(:,:,n+1)=RL.*FL+RS.*FS;
    
    KS=0.20983-2.081e-3*CS*100.0;%solid thermal conductivity from D-1 [W/mm/K]; KS=NaN if FS<=0
    KS=KS*1000.0;%[W/m/K]; [NIY+2,NIX+2]
    
    KL=6.5923e-2+3.3e-5*(T(1:NIY+2,1:NIX+2)-273.15)-6.807e-4*CL*100.0;%liquid thermal conductivity from D-2 [W/mm/K]
    KL=1000.0*KL;%[W/m/K]; [NIY+2,NIX+2]
    
    CpS=0.88+4.446e-4*(T(1:NIY+2,1:NIX+2)-273.15)-2.274e-3*CS*100.0;%solid specific heat capacity from D-3 [J/g/K]
    CpS=CpS*1000.0;%[J/kg/K]
    
    CpL=1.086-5.928e-3*CL*100.0;%liquid specific heat capacity from D-4 [J/g/K]
    CpL=CpL*1000.0;%[J/kg/K]
    
    %----------------------- UPDATE PARAMETERS ----------------------------    
    
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
    
    RM=RL.*FL+RS.*FS;%mean density of ingot [kg/m^3]
    PlotTC(T,FS,CS,RM,t(n+1),n);

    
end

