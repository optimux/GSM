function FindNSE(FSO,FSN,step)
%[NFSS,NFSE] = FindNSE(dFSO,dFSN,step)
%To find step numbers which show solidification starts and ends,
%repectively, for step>=2
%Created 2019-10-20

%NFSS: record step at which solid begins to form
%NFSE: record step at which no more solid to form
%NOTE: these two variables above are used to cope with inhomogeneous Cu
%distribution in solid

%dFSO: old dFS
%dFSN: new dFS

global NIX
global NIY
global NFSS
global NFSE

errFS=0.0001;

%======================= VERY IMPORTANT NOTE ==============================
%Sometimes, remelting may occur, thus dFS<0.0; and also somtimes during
%solidification process, there is no crystallization for some steps. So,
%the judgement below based on dFS is no better than the one based on FS!
%======================= VERY IMPORTANT NOTE ==============================

% for i=1:NIX
%     for j=1:NIY
%         %Situation One: no solidification yet
%         if((dFSO(j,i)<=0.0)&&(dFSN(j,i)<=0.0))
%             NFSS(j,i)=step;%solidification start step number should change until real solid shows up
%             NFSE(j,i)=step;%solidification end step number should >=NFSS
%         end
%         %Situation Two: solidification just begins right now
%         if((dFSO(j,i)<=0.0)&&(dFSN(j,i)>0.0))
%             NFSS(j,i)=step;%solidification start step number should change until real solid shows up
%             NFSE(j,i)=step;%solidification end step number should >=NFSS
%         end
%         %Situation Three: solidification keeps working
%         if((dFSO(j,i)>0.0)&&(dFSN(j,i)>0.0))
%             %NFSS(j,i)=step;%solidification start step number has been recorded
%             NFSE(j,i)=step;%solidification end step number should >=now
%         end
%         %Situation Four: solidification ends
%         if((dFSO(j,i)>0.0)&&(dFSN(j,i)<=0.0))
%             %NFSS(j,i)=step;%solidification start step number has been recorded
%             NFSE(j,i)=step;%solidification end step number should >=now
%         end
%
%     end
% end

for i=1:NIX+2
    for j=1:NIY+2
        %Situation One: no solidification yet
        if((abs(FSO(j,i))<=errFS)&&(abs(FSN(j,i))<=errFS))
            NFSS(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE(j,i)=step;%solidification end step number should >=NFSS
        end
        %Situation Two: solidification just begins right now
        if((abs(FSO(j,i))<=errFS)&&(abs(FSN(j,i))>errFS)&&(abs(FSN(j,i)-1.0)>=errFS))
            NFSS(j,i)=step;%solidification start step number should change until real solid shows up
            NFSE(j,i)=step;%solidification end step number should >=NFSS
        end
        %Situation Three: solidification keeps working
        if((abs(FSO(j,i))>errFS)&&(abs(FSO(j,i)-1.0)>=errFS)&&(abs(FSN(j,i))>errFS)&&(abs(FSN(j,i)-1.0)>=errFS))
            %NFSS(j,i)=step;%solidification start step number has been recorded
            NFSE(j,i)=step;%solidification end step number should >=now
        end
        %Situation Four: solidification ends
        if((abs(FSO(j,i))>errFS)&&(abs(FSO(j,i)-1.0)>=errFS)&&(abs(FSN(j,i)-1.0)<errFS))
            %NFSS(j,i)=step;%solidification start step number has been recorded
            NFSE(j,i)=step;%solidification end step number should >=now
        end
        
    end
end


end

