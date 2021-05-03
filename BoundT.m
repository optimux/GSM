function [TVB,THB]= BoundT(TVTM,THTM,KVTM,KHTM)
%To calculate boundary temperature
%Created 2019-10-18

%Modified for 2D pure solid heat conduction verification! 2019-11-8

%TVB: T of vertical boundary
%THB: T of horizontal boundary
%TVTM: T left to vertical boundary
%THTM: T upstair to horizontal boundary
%KVTM: thermal conductivity left to boundary
%KHTM: thermal conductivity upstair to boundary

global TR
global dx
global dy
global NIX
global NIY
global NMX
global NMY

a1=0.10809*1000.0;%thermal conductivity for T<328.46 C from D-9
b1=-1.0387e-4*1000.0;
a2=0.061553*1000.0;%thermal conductivity for T>=328.46 C from D-9
b2=-2.6516e-5*1000.0;
c1=7.9852e-6*10^6;%convective heat transfer coefficient from D-12 [W/m^2/K]
c2=6.1e-5*10^6;%convective heat transfer coefficient from D-14 [W/m^2/K]

n=length(TVTM);
TVB=zeros(n,1);
m=length(THTM);
THB=zeros(m,1);

err=1.0e-6;%T error [K]
sig=5.67e-8;%Boltzmann constant [W/m^2/K^4]
eps=1.0;%emissivity [1]
tri=1.0/3.0;

%% .............============= Vertical Boundary ================...........
for i=2:n-1
    err0=1.0;
    TX=300.0;%initial guess of boundary temperature [K]
    while(err0>err)
        g=0.5*c1*dx(NIX+NMX)*(TVTM(i)-TR)*(TX-TR)^0.25/dy(i-1)^0.25-0.5*KVTM(i)*TVTM(i)-0.5*TVTM(i)*(a1+b1*TX)+0.5*KVTM(i)*TX+0.5*TX*(a1+b1*TX)-0.5*c1*dx(NIX+NMX)*(TVTM(i)-TX)*(TX-TR)^0.25/dy(i-1)^0.25;
        dgdx=0.5*c1*dx(NIX+NMX)*(TVTM(i)-TR)*0.25*(TX-TR)^(-0.75)/dy(i-1)^0.25-0.5*b1*TVTM(i)+0.5*KVTM(i)+0.5*a1+b1*TX-0.5*c1*dx(NIX+NMX)*(0.25*(TX-TR)^(-0.75)*(TVTM(i)-TX)-(TX-TR)^0.25)/dy(i-1)^0.25;
        err0=abs(g/dgdx);
        TX=TX-g/dgdx;
    end
    
    if(TX<328.46+273.15)
        TVB(i)=TX;
    else
        
        err0=1.0;
        TX=300.0;%initial guess of boundary temperature [K]
        while(err0>err)
            g=0.5*c1*dx(NIX+NMX)*(TVTM(i)-TR)*(TX-TR)^0.25/dy(i-1)^0.25-0.5*KVTM(i)*TVTM(i)-0.5*TVTM(i)*(a2+b2*TX)+0.5*KVTM(i)*TX+0.5*TX*(a2+b2*TX)-0.5*c1*dx(NIX+NMX)*(TVTM(i)-TX)*(TX-TR)^0.25/dy(i-1)^0.25;
            dgdx=0.5*c1*dx(NIX+NMX)*(TVTM(i)-TR)*0.25*(TX-TR)^(-0.75)/dy(i-1)^0.25-0.5*b2*TVTM(i)+0.5*KVTM(i)+0.5*a2+b2*TX-0.5*c1*dx(NIX+NMX)*(0.25*(TX-TR)^(-0.75)*(TVTM(i)-TX)-(TX-TR)^0.25)/dy(i-1)^0.25;
            err0=abs(g/dgdx);
            TX=TX-g/dgdx;
        end
        
        TVB(i)=TX;
    end
    
end
TVB(1)=TVB(2);

%% .............============ Horizontal Boundary ===============...........
for i=2:m-1
    err0=1.0;
    TX=300.0;
    while(err0>err)
        g=0.5*c2*dy(NIY+NMY)*(THTM(i)-TR)*(TX-TR)^tri/dx(i-1)^(2.0*tri)-0.5*THTM(i)*KHTM(i)-0.5*THTM(i)*(a1+b1*TX)+0.5*TX*(KHTM(i)+a1+b1*TX)-0.5*dy(NIY+NMY)*c2*(THTM(i)-TX)*(TX-TR)^tri/dx(i-1)^(2.0*tri);
        dgdx=0.5*c2*dy(NIY+NMY)*(THTM(i)-TR)*tri*(TX-TR)^(-2.0*tri)-0.5*b1*THTM(i)+0.5*KHTM(i)+0.5*a1+b1*TX-0.5*c2*dy(NIY+NMY)*(tri*(TX-TR)^(-2.0*tri)*(THTM(i)-TX)-(TX-TR)^tri)/dx(i-1)^(2.0*tri);
        err0=abs(g/dgdx);
        TX=TX-g/dgdx;
    end
    
    if(TX<328.46+273.15)
        THB(i)=TX;
    else
        
        err0=1.0;
        TX=300.0;%initial guess of boundary temperature [K]
        while(err0>err)
            g=0.5*c2*dy(NIY+NMY)*(THTM(i)-TR)*(TX-TR)^tri/dx(i-1)^(2.0*tri)-0.5*THTM(i)*KHTM(i)-0.5*THTM(i)*(a2+b2*TX)+0.5*TX*(KHTM(i)+a2+b2*TX)-0.5*dy(NIY+NMY)*c2*(THTM(i)-TX)*(TX-TR)^tri/dx(i-1)^(2.0*tri);
            dgdx=0.5*c2*dy(NIY+NMY)*(THTM(i)-TR)*tri*(TX-TR)^(-2.0*tri)-0.5*b2*THTM(i)+0.5*KHTM(i)+0.5*a2+b2*TX-0.5*c2*dy(NIY+NMY)*(tri*(TX-TR)^(-2.0*tri)*(THTM(i)-TX)-(TX-TR)^tri)/dx(i-1)^(2.0*tri);
            err0=abs(g/dgdx);
            TX=TX-g/dgdx;
        end
        
        THB(i)=TX;
    end
end
THB(1)=THB(2);


end

