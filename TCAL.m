function [T,FL,CL,dFS,iters]=TCAL(TTM,FLTM,FSTM,KLTM,KSTM,CLTM,CSTM,VXTM,VYTM,FRCPSTM,FRCSTM,FKTM,RSCPSTM,RLCPLTM,RSTM,RLTM,CPLTM,dtb)%FRCSTM,FRHTM,
%Firstly, regard this problem as a 2D pure solid heat conduction problem in
%which no species diffusion, no convection!
%NOTE: we consider fully liquid ingot as pure solid!

%Created 2019-11-8

%TTM: temporary T, of [NIY+NMY+2,NIX+NMX+2]
%FLTM: temporary FL, of [NIY+2,NIX+2]
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
global dx
global dy
global gap

%======================= VERY IMPORTANT NOTE ==============================
%It is highly recommended the usage of packages for complex formulae, i.e.,
%temporary variables. These packages can be utilized repetedly, easily
%modified, and produce a concise numerical equation.
%======================= VERY IMPORTANT NOTE ==============================


%################################################################### ENERGY BALANCE #########################################################################
%% .............========= DIFFUSION HEAT FLUX ==========...................

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
    %TQDT(i,NIX)=-2.0*(TQDX(i+1,NIX)-TQDX(i+1,NIX+1))*dtb/(dx(NIX)+0.5*gap)-2.0*(TQDY(i,NIX+1)-TQDY(i+1,NIX+1))*dtb/dy(i);%[J/m^3]
    TQDT(i,NIX+1)=-2.0*(TQDX(i+1,NIX+1)-TQDX(i+1,NIX+2))*dtb/(dx(NIX+1)+gap)-2.0*(TQDY(i,NIX+2)-TQDY(i+1,NIX+2))*dtb/dy(i);%[J/m^3]
end

%Modify horizontal gas gap
for i=1:NIX
    %TQDT(NIY,i)=-2.0*(TQDX(NIY+1,i)-TQDX(NIY+1,i+1))*dtb/dx(i)-2.0*(TQDY(NIY,i+1)-TQDY(NIY+1,i+1))*dtb/(dy(NIY)+0.5*gap);%[J/m^3]
    TQDT(NIY+1,i)=-2.0*(TQDX(NIY+2,i)-TQDX(NIY+2,i+1))*dtb/dx(i)-2.0*(TQDY(NIY+1,i+1)-TQDY(NIY+2,i+1))*dtb/(dy(NIY+1)+gap);%[J/m^3]
end

%NOTE: boundary heat fluxes of two Situations are almost idenstical within 3% difference
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

%################################################################### ENERGY BALANCE #########################################################################
VXTM=0.1*ones(NIY+2,NIX+1);
VYTM=0.1*ones(NIY+1,NIX+2);
TFN=zeros(NIY+NMY+2,NIX+NMX+2);%T Factor of Next step
TN=zeros(NIY+NMY+2,NIX+NMX+2);%T of Next step
dFSTM=zeros(NIY+2,NIX+2);%always
TQVT=zeros(NIY+NMY,NIX+NMX);%T=T, Q=flux, V=velocity, T=total, always

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

T=TN;
FL=ones(NIY+2,NIX+2);
CL=0.045*ones(NIY+2,NIX+2);
dFS=zeros(NIY+2,NIX+2);
iters=1;
VXTM=0.1-VXTM;
VYTM=0.1-VYTM;

end

