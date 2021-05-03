function VRD=RISER(RLOTM,FLOTM,FROTM,RLTM,FLTM,RSTM,dFSTM,VXOTM,VYOTM,RLR)
%This script is used to calculate velocity at riser bottom
%Created 2020-1-2

%RLOTM: old liquid density
%FLOTM: old FL
%FROTM: old FS*RS
%RLTM: new liquid density
%FLTM: new FL
%RSTM: new solid density
%dFSTM: new solid fraction
%VXOTM: old x-axis velocity
%VYOTM: old y-axis velocity
%RFVXOTM: old RFVX
%RFVYOTM: old RFVY
%RLR: liquid density reference at T0

global dtb
global NIY
global NIX
global dx
global dy
global VR

%% ========================  SITUATION 1  =================================
MO=zeros(NIY,NIX);%old mass in control volume
VN=zeros(NIY,NIX);%new volume of cooler mass
FRTM=zeros(NIY+2,NIX+2);%new time-step FS*RS
VTD=0.0;%total volume deficit due to cooling
MOVS=0.0;%mass of old velocity supplied
RFVXOTM=zeros(NIY+2,NIX+1);
RFVYOTM=zeros(NIY+1,NIX+2);

for i=1:NIX+1
    for j=1:NIY+2
        RFVXOTM(j,i)=0.5*(RLOTM(j,i)+RLOTM(i+1))*0.5*(FLOTM(j,i)+FLOTM(j,i+1))*VXOTM(j,i);
    end
end

for i=1:NIX+2
    for j=1:NIY+1
        RFVYOTM(j,i)=0.5*(RLOTM(j,i)+RLOTM(j+1,i))*0.5*(FLOTM(j,i)+FLOTM(j+1,i))*VYOTM(j,i);
    end
end

for i=1:NIX
    for j=1:NIY
        MO(j,i)=dx(i)*dy(j)*(RLOTM(j+1,i+1)*FLOTM(j+1,i+1)+FROTM(j+1,i+1));%old mass in control volume
        FRTM(j+1,i+1)=FROTM(j+1,i+1)+dFSTM(j+1,i+1)*RSTM(j+1,i+1);
        VN(j,i)=MO(j,i)/(RLTM(j+1,i+1)*FLTM(j+1,i+1)+FRTM(j+1,i+1));
        VTD=VTD+(dx(i)*dy(j)-VN(j,i));
        
        MOVS=MOVS+(RFVXOTM(j+1,i)-RFVXOTM(j+1,i+1))*dy(j)*dtb+(RFVYOTM(j,i+1)-RFVYOTM(j+1,i+1))*dx(i)*dtb;
    end
end

VRD=(VTD*RLR-MOVS)/((dx(1)+dx(2))*dtb*RLR);%riser velocity downward
VR=VRD;
% %% ========================  SITUATION 2  =================================
% MO=zeros(NIY,NIX);%old mass in control volume
% VN=zeros(NIY,NIX);%new volume of cooler mass
% FRTM=zeros(NIY+2,NIX+2);%new time-step FS*RS
% VTD=0.0;%total volume deficit due to cooling
% 
% for i=1:NIX
%     for j=1:NIY
%         MO(j,i)=dx(i)*dy(j)*(RLOTM(j+1,i+1)*FLOTM(j+1,i+1)+FROTM(j+1,i+1));%old mass in control volume
%         FRTM(j+1,i+1)=FROTM(j+1,i+1)+dFSTM(j+1,i+1)*RSTM(j+1,i+1);
%         VN(j,i)=MO(j,i)/(RLTM(j+1,i+1)*FLTM(j+1,i+1)+FRTM(j+1,i+1));
%         VTD=VTD+(dx(i)*dy(j)-VN(j,i));
%     end
% end
% 
% VR=VTD/((dx(1)+dx(2))*dtb);%riser velocity downward

end

