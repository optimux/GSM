%Solidification shrinkage figure 4 in Xu et al., 1991

%Created 2020-3-18
%output name: shrinkage.xlsx

CL=[0.0:0.01:0.35]';
CE=0.332;%eutectic
T=700.0+273.15;%K

RL=2.5838e-3+1.6567e-5*CL*100.0+2.6891e-7*(CL*100).^2-3.16e-7*(T-273.15);%liquid density from D-6 [g/mm^3]
RL=RL*10^6;%[kg/m^3]

RS=zeros(length(CL),1);%density of solid from D-5 [kg/m^3]; only [2:NIY+1,2:NIX+1] valid, outmost layer for plot
for i=1:length(CL)
    if(CL(i)<CE)
        RS(i)=2.8e-3+1.279e-6*CL(i)*100.0-3.89e-7*(T-273.15);
    else
        RS(i)=3.6302e-3-4.2e-7*(T-273.15);
    end
end
RS=RS*10^6;%[kg/m^3]

beta=1.0-RL./RS;

plot(CL*100.0,beta);
hold on
plot(CL*100.0,0.0);
axis([0.0 35.0 -0.3 0.2]);
xlabel('Cu [wt%]');
ylabel('Shrinkage');