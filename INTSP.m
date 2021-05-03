function [FRCPS,FRCS,FK,FR,RSCPS,RLCPL] = INTSP(FRCPSTM,FRCSTM,FKTM,FRTM,CPSTM,CPLTM,CSTM,KSTM,RSTM,RLTM,dFSTM)%FRH,FRHTM,HLTM,
%To integrate solid properties, FS*RS*CpS, FS*RS*CS, FS*KS, FS*RS, RS*CpS, RL*CpL; since
%KS, CpS, HL, RS are functions of CS which is not homogeneous in solidified
%portion of finite volume.
%Created 2019-10-17

%---------- INPUT: old quantity -------
%FRCPSTM = FS*RS*CpS [NIY+2,NIX+2]
%FRCSTM = FS*RS*CS [NIY+2,NIX+2]
%FKTM = FS*KS [NIY+2,NIX+2]
%FRTM = FS*RS [NIY+2,NIX+2]
%CPSTM [NIY+2,NIX+2]
%CPLTM [NIY+2,NIX+2]
%CSTM [NIY+2,NIX+2]
%KSTM [NIY+2,NIX+2]
%RSTM [NIY+2,NIX+2]
%RLTM [NIY+2,NIX+2]
%dFSTM [NIY+2,NIX+2]

%---------- OUTPUT: new quantity -------
%FRCPS: FS*RS*CpS [NIY+2,NIX+2]
%FRCS: FS*RS*CS [NIY+2,NIX+2]
%FK: FS*KS [NIY+2,NIX+2]
%FR: FS*RS [NIY+2,NIX+2]
%FRH: FS*RS*HL(CS) [NIY+2,NIX+2]
%RSCPS: RS*CpS [NIY+2,NIX+2]
%RLCPL: RL*CpL [NIY+2,NIX+2]

global NIX
global NIY

FRCPS=zeros(NIY+2,NIX+2);%FS*RS*CpS
FRCS=zeros(NIY+2,NIX+2);%FS*RS*CS
FK=zeros(NIY+2,NIX+2);%FS*KS
FR=zeros(NIY+2,NIX+2);%FS*RS
%FRH=zeros(NIY,NIX);%FS*RS*HL(CS)
RSCPS=zeros(NIY+2,NIX+2);%RS*CPS for current step
RLCPL=zeros(NIY+2,NIX+2);%RL*CPL for current step

%New quantity=Old quantity+delta(quantity)
%In fact, only one sentence is needed (matrix addition); but we unfold these fomulae for
%better understanding as following
for i=1:NIX+2
    for j=1:NIY+2
        FRCPS(j,i)=FRCPSTM(j,i)+dFSTM(j,i)*RSTM(j,i)*CPSTM(j,i);
        FRCS(j,i)=FRCSTM(j,i)+dFSTM(j,i)*RSTM(j,i)*CSTM(j,i);
        FK(j,i)=FKTM(j,i)+dFSTM(j,i)*KSTM(j,i);
        FR(j,i)=FRTM(j,i)+dFSTM(j,i)*RSTM(j,i);
        %FRH(j,i)=FRHTM(j,i)+dFSTM(j,i)*RSTM(j+1,i+1)*HLTM(j+1,i+1);
        RSCPS(j,i)=RSTM(j,i)*CPSTM(j,i);
        RLCPL(j,i)=RLTM(j,i)*CPLTM(j,i);
    end
end

end

