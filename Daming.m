%This script is used to reproduce figures in 'Numerical method for solution
%of strongly coupled binary alloy solidification problems', Daming Xu 1991
%Created 2019-9-23

clear all

W=160.0/1000.0;%ingot width [m]
H=200.0/1000.0;%ingot height [m]
WM=368.0/1000.0;%mold width [m]
HM=320.0/1000.0;%mold height [m]
TR=25.0+273.15;%T outside of mold [K]
TM=25.0+273.15;%Initial T of mold [K]
XS=20;%x-axis intervals in ingot
YS=20;%y-axis intervals in ingot
dx=W/XS;%[m]
dy=H/YS;%[m]
dz=1.0;%assume unit thickness [m]
x=[0.0:dx:W];%x grids [m]
y=[0.0:dy:H];%y grids [m]
xT=[0.0:dx:WM];%coordinate for temperature, mold included
yT=[0.0:dy:HM];%coordinate for temperature, mold included
XST=WM/dx;%x-axis intervals
YST=HM/dy;%y-axis intervals

T0=700.0+273.15;%Initial T [K]
C0=4.5/100.0;%initial liquid alloy composition as Cu [wt %]
CE=0.332;%eutectic composition [wt %]
DL=3.0e-11;%liquid diffusion coefficient [m^2/sec]
DS=0.0;%solid diffusion coefficient [m^2/sec]
g=9.81;%Earth gravity [m/sec^2]
mu=3.0e-9;%dynamic viscosity of liquid, m^2/sec ==> wrong! 2019-10-16
time=300.0;%total time [sec]
dt=0.15;%min time interval [sec]; from Daming Xu 1991
Nt=time/dt+1;%max pages, with initial condition included
phi=1.0;%integration coefficient

TGTM=ones(YS+1,XS+1)*T0;%temporary T grids in ingot [K]; G==ingot
TGTM(:,1)=TM;TGTM(:,XS+1)=TM;TGTM(1,:)=TM;%boundary temperature in ingot [K]
TG=T0*ones(YS+1,XS+1,Nt);%records of T grids in ingot [K]; G==ingot
TG(:,:,1)=TGTM;%initial T distribution in ingot [K]; G==ingot

TTM=TM*ones(YST+1,XST+1);%mold initial temperature [K]
T=TM*ones(YST+1,XST+1,Nt);%records of T [K]
x1=(XST-XS)/2+1;x2=x1+XS;y1=YST-YS+1;y2=YST+1;
TTM(y1:y2,x1:x2)=TGTM;

FLTM=ones(YS+1,XS+1);%temporary liquid volume fraction grids [1]
FL=ones(YS+1,XS+1,Nt);%records of liquid volume fraction grids [1]
FL(:,:,1)=FLTM;%initial liquid volume fraction grids [1]

FSTM=1.0-FLTM;%temporary solid volume fraction grids [1]
FS=1.0-FL;%records of solid volume fraction grids [1]
FS(:,:,1)=FSTM;%initial solid volume fraction grids [1]

CLTM=C0*ones(YS+1,XS+1);%temporary(initial) composition grids [wt %]; liquid composition grids [wt %]
CL=C0*ones(YS+1,XS+1,Nt);%records of solid volume fraction grids [1]
CL(:,:,1)=CLTM;%initial solid volume fraction grids [1]

RhoLTM=2.5838e-3-3.16e-7*TGTM+1.6567e-5*CLTM+2.6891e-7*CLTM.^2;%temporary liquid density from Table 2 [g/mm^3]
RhoLTM=10^6*RhoLTM;%temporary liquid density [kg/m^3]
RhoL=zeros(YS+1,XS+1,Nt);%records of liquid density [kg/m^3]
RhoL(:,:,1)=RhoLTM;%initial liquid density [kg/m^3]

%temporary solid density from Table 2 [g/mm^3]
if(CLTM<CE)
    RhoSTM=2.8e-3+1.279e-6*CLTM-3.89e-7*TGTM;
else
    RhoSTM=3.6302e-3-4.2e7*TGTM;
end
RhoSTM=10^6*RhoSTM;%temporary solid density [kg/m^3]
RhoS=zeros(YS+1,XS+1,Nt);%records of solid density [kg/m^3]
RhoS(:,:,1)=RhoSTM;%initial solid density [kg/m^3]

RhoTM=RhoLTM*FLTM+RhoSTM*FSTM;%density of mixture [kg/m^3]
Rho=zeros(YS+1,XS+1,Nt);%records of mixture density [kg/m^3]
Rho(:,:,1)=RhoTM;%initial mixture density [kg/m^3]

KLTM=6.5923e-2+3.3e-5*TGTM-6.807e-4*CLTM;%temporary liquid thermal conductivity from Table 2 [W/mm/K]
KLTM=KLTM*1000.0;%temporary liquid thermal conductivity [W/m/K]
KL=zeros(YS+1,XS+1,Nt);%records of liquid thermal conductivity [W/m/K]
KL(:,:,1)=KLTM;%initial liquid thermal conductivity [W/m/K]

CpLTM=1.086-5.928e-3*CLTM;%temporary liquid specific heat capacity from Table 2 [J/g/K]
CpLTM=1000.0*CpLTM;%temporary liquid specific heat capacity [J/kg/K]
CpL=zeros(YS+1,XS+1,Nt);%records of liquid specific heat capacity [J/kg/K]
CpL(:,:,1)=CpLTM;%initial liquid specific heat capacity [J/kg/K]

PTM=zeros(YS+1,XS+1);%temporary pressure [Pa]
for i=1:XS+1
    for j=2:YS+1
        PTM(j,i)=PTM(j-1,i)+RhoTM(j-1,i)*g*dy;
    end
end
P=zeros(YS+1,XS+1,Nt);%records of pressure [Pa]
P(:,:,1)=PTM;%initial pressure [Pa]

RPTM=zeros(YS+1,XS+1);%temporary relative pressure [Pa]
RP=zeros(YS+1,XS+1,Nt);%records of relative pressure [Pa]
RP(:,:,1)=RPTM;%initial relative pressure [Pa]

KPTM=0.12824+5.6691e-5*CLTM+3.7288e-5*CLTM.^2;%temporary partition coefficient [1]
KP=zeros(YS+1,XS+1,Nt);%records of partition coefficient [1]
KP(:,:,1)=KPTM;%initial partition coefficient [1]

CSTM=KPTM.*CLTM;%temporary solid composition grids [wt %]
CS=zeros(YS+1,XS+1,Nt);%records of solid composition grids [wt %]
CS(:,:,1)=CSTM;%initial solid composition grids [wt %]

KSTM=0.18983-2.081e-3*CSTM;%temporary solid thermal conductivity from Table 2 [W/mm/K]; should be 0 if no solid
KSTM=1000.0*KSTM;%temporary solid thermal conductivity [W/m/K]
KS=zeros(YS+1,XS+1,Nt);%records of solid thermal conductivity [W/m/K]
KS(:,:,1)=KSTM;%initial solid thermal conductivity [W/m/K]

CpSTM=0.88+4.446e-4*TGTM-2.274e-3*CSTM;%temporary solid specific heat capacity from Table 2 [J/g/K]
CpSTM=CpSTM*1000.0;%temporary solid specific heat capacity [J/kg/K]
CpS=zeros(YS+1,XS+1,Nt);%records of solid specific heat capacity [J/kg/K]
CpS(:,:,1)=CpSTM;%initial solid specific heat capacity [J/kg/K]

RhoSMTM=2.8e-3+1.279e-6*CLTM-3.89e-7*TGTM;%temporary mean solid density [kg/m^3]; for the 1st time, RhoSM==RhoS (mean value==in situ value)
RhoSM=zeros(YS+1,XS+1,Nt);%records of mean solid density [kg/m^3]
RhoSM(:,:,1)=RhoSMTM;%initial mean solid density [kg/m^3]

CpSMTM=0.88+4.66e-4*TGTM-2.274e-3*CSTM;%temporary mean solid specific heat capacity from Table 2 [J/g/K]; for the 1st time, CpSM==CpS (mean value==in situ value)
CpSMTM=1000.0*CpSMTM;%temporary mean solid specific heat capacity [J/kg/K]
CpSM=zeros(YS+1,XS+1,Nt);%records of mean solid specific heat capacity [J/kg/K]
CpSM(:,:,1)=CpSMTM;%initial mean solid specific heat capacity [J/kg/K]

KMTM=FSTM.*KSTM+FLTM.*KLTM;%temporary mixture thermal conductivity [W/m/K]
KM=zeros(YS+1,XS+1,Nt);%records of mixture thermal conductivity [W/m/K]
KM(:,:,1)=KMTM;%initial mixture specific heat capacity [J/kg/K]

DFSTM=zeros(YS+1,XS+1);%temporary solid volume fraction change [1]
DFS=zeros(YS+1,XS+1,Nt);%records of solid volume fraction change [1]
DFS(:,:,2)=DFSTM;%initial solid volume fraction change [1]: let the initial value of DFS(i+1) equal zero (F1=DFS(i+1)=0) and calculate the approximation of T(i+1)' using Eq. (8).

TF0=FS(:,:,1).*RhoSM(:,:,1).*CpSM(:,:,1)+FL(:,:,1).*RhoL(:,:,1).*CpL(:,:,1);
TF1=TF0+(RhoS(:,:,1).*CpS(:,:,1)-RhoL(:,:,1).*CpL(:,:,1)).*DFS(:,:,2);

TDX=zeros(YS+1,XS);%TDY1,2 in Eq.(8), [W/m^2]
for i=1:XS
    for j=1:YS+1
        TDX(j,i)=0.5*(KM(j,i,1)+KM(j,i+1,1))*(TG(j,i+1,1)-TG(j,i,1))/dx;%(x(i+2)-x(i))/2==dx
    end
end

TDY=zeros(YS,XS+1);%TDZ1,2 in Eq.(8), [W/m^2]
for i=1:XS+1
    for j=1:YS
        TDY(j,i)=0.5*(KM(j,i,1)+KM(j+1,i,1))*(TG(j+1,i,1)-TG(j,i,1))/dy;%(y(i+2)-y(i))/2==dy
    end
end

VXTM=zeros(2*YS+1,2*XS+1);%temporary x velocity on all faces [m/sec]; boundary faces included; these meshes are set every half dx and dy, some of them will not be used
%NOTE: all grids have no explicit velocity
VX=zeros(2*YS+1,2*XS+1,Nt);%records of x velocity on all faces [m/sec]
VX(:,:,1)=VXTM;%initial x velocity on all faces [m/sec]

VYTM=zeros(2*YS+1,2*XS+1);%temporary y velocity on all faces [m/sec]; boundary faces included
%NOTE: all grids have no explicit velocity
VY=zeros(2*YS+1,2*XS+1,Nt);%records of y velocity on all faces [m/sec]
VY(:,:,1)=VYTM;%initial y velocity on all faces [m/sec]

RFVXTM=zeros(2*YS+1,2*XS+1);%approximate value of RFVY, x-axis; [kg/m^2/sec]
RFVXC=zeros(2*YS+1,2*XS+1);%corrected value of RFVY, x-axis; C==corrected; [kg/m^2/sec]
RFVX=zeros(2*YS+1,2*XS+1);%records of superficial mass flow, x-axis; [kg/m^2/sec]

RFVYTM=zeros(2*YS+1,2*XS+1);%approximate value of RFVY, y-axis; [kg/m^2/sec]
RFVYC=zeros(2*YS+1,2*XS+1);%corrected value of RFVY, y-axis; C==corrected; [kg/m^2/sec]
RFVY=zeros(2*YS+1,2*XS+1);%records of superficial mass flow, y-axis; [kg/m^2/sec]

%x-axis mid-grid superficial velocity every dy [kg/m^2/sec]
for i=2:2:2*XS
    k1=int32(i/2);
    for j=1:2:2*YS+1
        k2=int32((j+1)/2);
        RFVXTM(j,i)=0.5*(RhoL(k2,k1,1)*FL(k2,k1,1)+RhoL(k2,k1+1,1)*FL(k2,k1+1,1))*VX(j,i,1);
        RFVYTM(j,i)=0.5*(RhoL(k2,k1,1)*FL(k2,k1,1)+RhoL(k2,k1+1,1)*FL(k2,k1+1,1))*VY(j,i,1);
    end
end

%y-axis mid-grid superficial velocity very dx [kg/m^2/sec]
for i=1:2:2*XS+1
    k1=int32((i+1)/2);
    for j=2:2:2*YS
        k2=int32(j/2);
        RFVXTM(j,i)=0.5*(RhoL(k2,k1,1)*FL(k2,k1,1)+RhoL(k2+1,k1,1)*FL(k2+1,k1,1))*VX(j,i,1);
        RFVYTM(j,i)=0.5*(RhoL(k2,k1,1)*FL(k2,k1,1)+RhoL(k2+1,k1,1)*FL(k2+1,k1,1))*VY(j,i,1);
    end
end

%grid superficial velocity very dx [kg/m^2/sec]
for i=3:2:2*XS-1
    k1=int32((i+1)/2);
    for j=3:2:2*YS-1
        k2=int32((j+1)/2);
        RFVXTM(j,i)=RhoL(k2,k1,1)*FL(k2,k1,1)*0.25*(VX(j-1,i,1)+VX(j+1,i,1)+VX(j,i-1,1)+VX(j,i+1,1));
        RFVYTM(j,i)=RhoL(k2,k1,1)*FL(k2,k1,1)*0.25*(VY(j-1,i,1)+VY(j+1,i,1)+VY(j,i-1,1)+VY(j,i+1,1));
    end
end

%NOTE: boudary mid-grid and all grids are static

TVXTM=zeros(2*YS+1,2*XS+1);%part 2 in Eq.(8), x-axis
TVYTM=zeros(2*YS+1,2*XS+1);%part 2 in Eq.(8), y-axis

for i=2:2:2*XS
    k1=int32(i/2);
    for j=3:2:2*YS-1
        k2=int32((j+1)/2);
        TVXTM(j,i)=TG(k2,k1,1)*max([RFVXTM(j,i),0.0])+TG(k2,k1+1,1)*min([RFVXTM(j,i),0.0]);%TVY1,2 in Eq.(8)
    end
end

for i=3:2:2*XS-1
    k1=int32((i+1)/2);
    for j=2:2:2*YS
        k2=int32(j/2);
        TVYTM(j,i)=TG(k2,k1,1)*max([RFVYTM(j,i),0.0])+TG(k2+1,k1,1)*min([RFVYTM(j,i),0.0]);%TVZ1,2 in Eq.(8)
    end
end

for i=1:XST+1
    for j=1:YST+1
        T1=TF0
    end
end


errdFS=0.0001;%iteration convergence limit for dFS
errT=0.0001;%iteration convergence limit for T
errCL=0.0001;%iteration convergence limit for CL
err=max([errdFS,errT,errCL]);%max of iteration error
error=1.0;

iters=0;%iteration counting

while(error>err)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    iters=iters+1;
    
end


