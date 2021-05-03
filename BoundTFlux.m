function [FVB,FHB] = BoundTFlux(TVTM,TUTM,KLTM,KUTM)
%To calculate boundary heat flux
%Created 2019-11-4

%FVB: flux of vertical boundary
%FHB: flux of horizontal boundary

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
FVB=zeros(n-2,1);
HRV=zeros(n-2,1);
m=length(TUTM);
FHB=zeros(m-2,1);
HRB=zeros(m-2,1);

err=0.01;%T error [K]
sig=5.67e-8;%Boltzmann constant [W/m^2/K^4]
eps=1.0;%emissivity [1]

%------------------- Right Boundary ---------------------
for i=1:n-2
    HRV(i)=c1*((TVTM(i+1)-TR)/dy(i))^0.25;
    FVB(i)=(TVTM(i+1)-TR)*HRV(i)/(1.0+HRV(i)*0.5*dx(NIX+NMX)/KLTM(i+1));
end

%------------------- Bottom Boundary ---------------------
for i=1:m-2
    HRB(i)=c2*((TUTM(i+1)-TR)/dx(i)^2)^(1.0/3.0);
    FHB(i)=(TUTM(i+1)-TR)*HRB(i)/(1.0+HRB(i)*0.5*dy(NIY+NMY)/KUTM(i+1));
end

end

