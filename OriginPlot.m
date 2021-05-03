function OriginPlot(TP,FSP,VXP,VYP,VSXP,VSYP)
%This function is used to plot in OriginLab format
%Created on 2020-7-3

global NIX
global NIY
global x
global y

VXPP=0.0;%liquid x-axis velocity [m/sec]
VYPP=0.0;%liquid y-axis velocity [m/sec]
VPP=0.0;%liquid absolute magnitute velocity [m/sec]
DIRL=0.0;%liquid velocity direction assuming 0.0 is right horizon arrow [radians]

VSXPP=0.0;%crystal x-axis velocity [m/sec]
VSYPP=0.0;%crystal y-axis velocity [m/sec]
VSPP=0.0;%solid absolute magnitute velocity [m/sec]
DIRS=0.0;%solid velocity direction assuming 0.0 is right horizon arrow [radians]

Frame=zeros(NIY*NIX,24);

for i=1:NIX
    for j=1:NIY
        %--------------------- LIQUID ------------------------
        VXPP=0.5*(VXP(j+1,i)+VXP(j+1,i+1));%x-axis liquid velocity [m/sec]
        VYPP=0.5*(VYP(j,i+1)+VYP(j+1,i+1));%y-axis liquid velocity [m/sec]
        VPP=sqrt(VXPP^2+VYPP^2);%absolute magnitute of liquid velocity [m/sec]
        if((abs(VYPP)<=1.0e-9)&&(abs(VXPP)>1.0e-9))%VYPP=0.0, VXPP!=0.0
            if(VXPP>0.0)
                DIRL=0.0;
            else
                DIRL=pi;
            end
        end
        if((abs(VXPP)<=1.0e-9)&&(abs(VYPP)>1.0e-9))%VXPP=0.0, VYPP!=0.0
            if(VYPP>0.0)
                DIRL=0.5*pi;
            else
                DIRL=1.5*pi;
            end
        end
        if((abs(VXPP)<=1.0e-9)&&(abs(VYPP)<=1.0e-9))
            DIRL=0.0;
        end
        if((VXPP>0.0)&&(VYPP>0.0))
        DIRL=atan(VYPP/VXPP);%direction of liquid velocity assuming 0.0 is right horizon arrow [radians]
        end
        if((VXPP>0.0)&&(VYPP<0.0))
            DIRL=atan(VYPP/VXPP)+2.0*pi;
        end
        if((VXPP<0.0)&&(VYPP>0.0))
            DIRL=atan(VYPP/VXPP)+pi;
        end
        if((VXPP<0.0)&&(VYPP<0.0))
            DIRL=pi+atan(VYPP/VXPP);
        end
        
        Frame(NIY*(i-1)+j,1)=x(i+1);
        Frame(NIY*(i-1)+j,2)=y(j+1);
        Frame(NIY*(i-1)+j,3)=VPP;
        Frame(NIY*(i-1)+j,4)=DIRL;
        
        %--------------------- OLIVINE --------------------------
        VSXPP=0.5*(VSXP.OL(j+1,i)+VSXP.OL(j+1,i+1));%x-axis solid velocity [m/sec]
        VSYPP=0.5*(VSYP.OL(j,i+1)+VSYP.OL(j+1,i+1));%y-axis solid velocity [m/sec]
        VSPP=sqrt(VSXPP^2+VSYPP^2);%absolute magnitute of solid velocity [m/sec]
        if((abs(VSYPP)<=1.0e-9)&&(abs(VSXPP)>1.0e-9))%VSYPP=0.0, VSXPP!=0.0
            if(VSXPP>0.0)
                DIRL=0.0;
            else
                DIRL=pi;
            end
        end
        if((abs(VSXPP)<=1.0e-9)&&(abs(VSYPP)>1.0e-9))%VSXPP=0.0, VSYPP!=0.0
            if(VSYPP>0.0)
                DIRL=0.5*pi;
            else
                DIRL=1.5*pi;
            end
        end
        if((abs(VSXPP)<=1.0e-9)&&(abs(VSYPP)<=1.0e-9))
            DIRL=0.0;
        end
        if((VSXPP>0.0)&&(VSYPP>0.0))
        DIRL=atan(VSYPP/VSXPP);%direction of liquid velocity assuming 0.0 is right horizon arrow [radians]
        end
        if((VSXPP>0.0)&&(VSYPP<0.0))
            DIRL=atan(VSYPP/VXPP)+2.0*pi;
        end
        if((VSXPP<0.0)&&(VSYPP>0.0))
            DIRL=atan(VSYPP/VSXPP)+pi;
        end
        if((VSXPP<0.0)&&(VSYPP<0.0))
            DIRL=pi+atan(VSYPP/VSXPP);
        end
        
        Frame(NIY*(i-1)+j,5)=x(i+1);
        Frame(NIY*(i-1)+j,6)=y(j+1);
        Frame(NIY*(i-1)+j,7)=VSPP;
        Frame(NIY*(i-1)+j,8)=DIRS;
        
        %--------------------- OPX --------------------------
        VSXPP=0.5*(VSXP.OPX(j+1,i)+VSXP.OPX(j+1,i+1));%x-axis solid velocity [m/sec]
        VSYPP=0.5*(VSYP.OPX(j,i+1)+VSYP.OPX(j+1,i+1));%y-axis solid velocity [m/sec]
        VSPP=sqrt(VSXPP^2+VSYPP^2);%absolute magnitute of solid velocity [m/sec]
        if((abs(VSYPP)<=1.0e-9)&&(abs(VSXPP)>1.0e-9))%VSYPP=0.0, VSXPP!=0.0
            if(VSXPP>0.0)
                DIRL=0.0;
            else
                DIRL=pi;
            end
        end
        if((abs(VSXPP)<=1.0e-9)&&(abs(VSYPP)>1.0e-9))%VSXPP=0.0, VSYPP!=0.0
            if(VSYPP>0.0)
                DIRL=0.5*pi;
            else
                DIRL=1.5*pi;
            end
        end
        if((abs(VSXPP)<=1.0e-9)&&(abs(VSYPP)<=1.0e-9))
            DIRL=0.0;
        end
        if((VSXPP>0.0)&&(VSYPP>0.0))
        DIRL=atan(VSYPP/VSXPP);%direction of liquid velocity assuming 0.0 is right horizon arrow [radians]
        end
        if((VSXPP>0.0)&&(VSYPP<0.0))
            DIRL=atan(VSYPP/VXPP)+2.0*pi;
        end
        if((VSXPP<0.0)&&(VSYPP>0.0))
            DIRL=atan(VSYPP/VSXPP)+pi;
        end
        if((VSXPP<0.0)&&(VSYPP<0.0))
            DIRL=pi+atan(VSYPP/VSXPP);
        end
        Frame(NIY*(i-1)+j,9)=x(i+1);
        Frame(NIY*(i-1)+j,10)=y(j+1);
        Frame(NIY*(i-1)+j,11)=VSPP;
        Frame(NIY*(i-1)+j,12)=DIRS;
        
        %--------------------- CPX --------------------------
        VSXPP=0.5*(VSXP.CPX(j+1,i)+VSXP.CPX(j+1,i+1));%x-axis solid velocity [m/sec]
        VSYPP=0.5*(VSYP.CPX(j,i+1)+VSYP.CPX(j+1,i+1));%y-axis solid velocity [m/sec]
        VSPP=sqrt(VSXPP^2+VSYPP^2);%absolute magnitute of solid velocity [m/sec]
        if((abs(VSYPP)<=1.0e-9)&&(abs(VSXPP)>1.0e-9))%VSYPP=0.0, VSXPP!=0.0
            if(VSXPP>0.0)
                DIRL=0.0;
            else
                DIRL=pi;
            end
        end
        if((abs(VSXPP)<=1.0e-9)&&(abs(VSYPP)>1.0e-9))%VSXPP=0.0, VSYPP!=0.0
            if(VSYPP>0.0)
                DIRL=0.5*pi;
            else
                DIRL=1.5*pi;
            end
        end
        if((abs(VSXPP)<=1.0e-9)&&(abs(VSYPP)<=1.0e-9))
            DIRL=0.0;
        end
        if((VSXPP>0.0)&&(VSYPP>0.0))
        DIRL=atan(VSYPP/VSXPP);%direction of liquid velocity assuming 0.0 is right horizon arrow [radians]
        end
        if((VSXPP>0.0)&&(VSYPP<0.0))
            DIRL=atan(VSYPP/VXPP)+2.0*pi;
        end
        if((VSXPP<0.0)&&(VSYPP>0.0))
            DIRL=atan(VSYPP/VSXPP)+pi;
        end
        if((VSXPP<0.0)&&(VSYPP<0.0))
            DIRL=pi+atan(VSYPP/VSXPP);
        end
        Frame(NIY*(i-1)+j,13)=x(i+1);
        Frame(NIY*(i-1)+j,14)=y(j+1);
        Frame(NIY*(i-1)+j,15)=VSPP;
        Frame(NIY*(i-1)+j,16)=DIRS;
        
        %--------------------- PL --------------------------
        VSXPP=0.5*(VSXP.PL(j+1,i)+VSXP.PL(j+1,i+1));%x-axis solid velocity [m/sec]
        VSYPP=0.5*(VSYP.PL(j,i+1)+VSYP.PL(j+1,i+1));%y-axis solid velocity [m/sec]
        VSPP=sqrt(VSXPP^2+VSYPP^2);%absolute magnitute of solid velocity [m/sec]
        if((abs(VSYPP)<=1.0e-9)&&(abs(VSXPP)>1.0e-9))%VSYPP=0.0, VSXPP!=0.0
            if(VSXPP>0.0)
                DIRL=0.0;
            else
                DIRL=pi;
            end
        end
        if((abs(VSXPP)<=1.0e-9)&&(abs(VSYPP)>1.0e-9))%VSXPP=0.0, VSYPP!=0.0
            if(VSYPP>0.0)
                DIRL=0.5*pi;
            else
                DIRL=1.5*pi;
            end
        end
        if((abs(VSXPP)<=1.0e-9)&&(abs(VSYPP)<=1.0e-9))
            DIRL=0.0;
        end
        if((VSXPP>0.0)&&(VSYPP>0.0))
        DIRL=atan(VSYPP/VSXPP);%direction of liquid velocity assuming 0.0 is right horizon arrow [radians]
        end
        if((VSXPP>0.0)&&(VSYPP<0.0))
            DIRL=atan(VSYPP/VXPP)+2.0*pi;
        end
        if((VSXPP<0.0)&&(VSYPP>0.0))
            DIRL=atan(VSYPP/VSXPP)+pi;
        end
        if((VSXPP<0.0)&&(VSYPP<0.0))
            DIRL=pi+atan(VSYPP/VSXPP);
        end
        Frame(NIY*(i-1)+j,17)=x(i+1);
        Frame(NIY*(i-1)+j,18)=y(j+1);
        Frame(NIY*(i-1)+j,19)=VSPP;
        Frame(NIY*(i-1)+j,20)=DIRS;
        
        %--------------------- ILM --------------------------
        VSXPP=0.5*(VSXP.ILM(j+1,i)+VSXP.ILM(j+1,i+1));%x-axis solid velocity [m/sec]
        VSYPP=0.5*(VSYP.ILM(j,i+1)+VSYP.ILM(j+1,i+1));%y-axis solid velocity [m/sec]
        VSPP=sqrt(VSXPP^2+VSYPP^2);%absolute magnitute of solid velocity [m/sec]
        if((abs(VSYPP)<=1.0e-9)&&(abs(VSXPP)>1.0e-9))%VSYPP=0.0, VSXPP!=0.0
            if(VSXPP>0.0)
                DIRL=0.0;
            else
                DIRL=pi;
            end
        end
        if((abs(VSXPP)<=1.0e-9)&&(abs(VSYPP)>1.0e-9))%VSXPP=0.0, VSYPP!=0.0
            if(VSYPP>0.0)
                DIRL=0.5*pi;
            else
                DIRL=1.5*pi;
            end
        end
        if((abs(VSXPP)<=1.0e-9)&&(abs(VSYPP)<=1.0e-9))
            DIRL=0.0;
        end
        if((VSXPP>0.0)&&(VSYPP>0.0))
        DIRL=atan(VSYPP/VSXPP);%direction of liquid velocity assuming 0.0 is right horizon arrow [radians]
        end
        if((VSXPP>0.0)&&(VSYPP<0.0))
            DIRL=atan(VSYPP/VXPP)+2.0*pi;
        end
        if((VSXPP<0.0)&&(VSYPP>0.0))
            DIRL=atan(VSYPP/VSXPP)+pi;
        end
        if((VSXPP<0.0)&&(VSYPP<0.0))
            DIRL=pi+atan(VSYPP/VSXPP);
        end
        Frame(NIY*(i-1)+j,21)=x(i+1);
        Frame(NIY*(i-1)+j,22)=y(j+1);
        Frame(NIY*(i-1)+j,23)=VSPP;
        Frame(NIY*(i-1)+j,24)=DIRS;
        
    end
end
xlswrite('VSV.xlsx',Frame(:,1:4),'LiqMov');
xlswrite('VSV.xlsx',Frame(:,5:8),'OLMov');
xlswrite('VSV.xlsx',Frame(:,9:12),'OPXMov');
xlswrite('VSV.xlsx',Frame(:,13:16),'CPXMov');
xlswrite('VSV.xlsx',Frame(:,17:20),'PLMov');
xlswrite('VSV.xlsx',Frame(:,21:24),'ILMMov');

TPI=TP(1:NIY+2,1:NIX+2);
TPI(1:NIY+1,NIX+2)=TPI(1:NIY+1,NIX+1);
TPI(NIY+2,1:NIX+2)=TPI(NIY+1,1:NIX+2);
xlswrite('G:\PhD\Img\TFS.xlsx',TPI-273.15,'Zhang2021T');

FSI=FSP(1:NIY+2,1:NIX+2);
FSI(1:NIY+1,NIX+2)=FSI(1:NIY+1,NIX+1);
FSI(NIY+2,1:NIX+2)=FSI(NIY+1,1:NIX+2);
xlswrite('G:\PhD\Img\TFS.xlsx',FSI,'Zhang2021FS');

end

