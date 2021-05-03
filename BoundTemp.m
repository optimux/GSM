function [TVB,THB]= BoundTemp(TVTM,TUTM,KLTM,KUTM)
%To calculate boundary temperature
%Created 2019-10-18

%TVB: T of vertical boundary
%THB: T of horizontal boundary
%TVTM: T left to vertical boundary
%TUTM: T upstair to horizontal boundary
%KLTM: thermal conductivity left to boundary
%KUTM: thermal conductivity upstair to boundary

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
m=length(TUTM);
THB=zeros(m,1);

err=0.01;%T error [K]
sig=5.67e-8;%Boltzmann constant [W/m^2/K^4]
eps=1.0;%emissivity [1]

%% .............============= Vertical Boundary ================...........
for i=2:n-1
    err0=1.0;
    TX=300.0;%initial guess of boundary temperature [K]
    while(err0>err)
        g=c1*(TX-TR)^1.25/dy(i-1)^0.25+eps*sig*TX^4+b1*TX^2/dx(NIX+NMX)+(a1-b1*TVTM(i)+KLTM(i))*TX/dx(NIX+NMX)-(TVTM(i)*KLTM(i)+a1*TVTM(i))/dx(NIX+NMX)-eps*sig*TR^4;
        dgdx=1.25*c1*(TX-TR)^0.25/dy(i-1)^0.25+4.0*eps*sig*TX^3+2.0*b1*TX/dx(NIX+NMX)+(a1-b1*TVTM(i)+KLTM(i))/dx(NIX+NMX);
        err0=abs(g/dgdx);
        TX=TX-g/dgdx;
    end
    
    if(TX<328.46+273.15)
        TVB(i)=TX;
    else
        
        err0=1.0;
        TX=300.0;%initial guess of boundary temperature [K]
        while(err0>err)
            g=c1*(TX-TR)^1.25/dy(i-1)^0.25+eps*sig*TX^4+b2*TX^2/dx(NIX+NMX)+(a2-b2*TVTM(i)+KLTM(i))*TX/dx(NIX+NMX)-(TVTM(i)*KLTM(i)+a2*TVTM(i))/dx(NIX+NMX)-eps*sig*TR^4;
            dgdx=1.25*c1*(TX-TR)^0.25/dy(i-1)^0.25+4.0*eps*sig*TX^3+2.0*b2*TX/dx(NIX+NMX)+(a2-b2*TVTM(i)+KLTM(i))/dx(NIX+NMX);
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
        g=eps*sig*TX^4+c2*(TX-TR)^(4.0/3.0)/dx(i-1)^(2.0/3.0)+(a1-b1*TUTM(i)+KUTM(i))*TX/dy(NIY+NMY)-eps*sig*TR^4+b1*TX^2/dy(NIY+NMY)-(TUTM(i)*KUTM(i)+a1*TUTM(i))/dy(NIY+NMY);
        dgdx=4.0*eps*sig*TX^3+4.0*c2*(TX-TR)^(1.0/3.0)/(3.0*dx(i-1)^(2.0/3.0))+(a1-b1*TUTM(i)+KUTM(i))/dy(NIY+NMY)+2.0*b1*TX/dy(NIY+NMY);
        err0=abs(g/dgdx);
        TX=TX-g/dgdx;
    end
    
    if(TX<328.46+273.15)
        THB(i)=TX;
    else
        
        err0=1.0;
        TX=300.0;%initial guess of boundary temperature [K]
        while(err0>err)
            g=eps*sig*TX^4+c2*(TX-TR)^(4.0/3.0)/dx(i-1)^(2.0/3.0)+(a2-b2*TUTM(i)+KUTM(i))*TX/dy(NIY+NMY)-eps*sig*TR^4+b2*TX^2/dy(NIY+NMY)-(TUTM(i)*KUTM(i)+a2*TUTM(i))/dy(NIY+NMY);
            dgdx=4.0*eps*sig*TX^3+4.0*c2*(TX-TR)^(1.0/3.0)/(3.0*dx(i-1)^(2.0/3.0))+(a2-b2*TUTM(i)+KUTM(i))/dy(NIY+NMY)+2.0*b2*TX/dy(NIY+NMY);
            err0=abs(g/dgdx);
            TX=TX-g/dgdx;
        end
        
        THB(i)=TX;
    end
end
THB(1)=THB(2);

%THB(m)=TVB(n)==??
end

