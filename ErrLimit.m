function err=ErrLimit(A)
%This script is used to calculate error limit based on machine error
%created 2020-3-27

AE=log10(abs(A));
[n,m]=size(A);
machine=A;
err=A;

for i=1:m
    for j=1:n
        if(AE(j,i)>0.0)%|A|>1.0
            if(A(j,i)>0.0)
                machine(j,i)=15.0-floor(AE(j,i))-2.0;%15 is significant figure of double-precision floting point,2 is decimal point and 1~10 leading number
            else
                machine(j,i)=15-floor(AE(j,i))-3.0;%15 is significant figure of double-precision floting point,3 is decimal point + 1~10 leading number + minus sign
            end
            
        else%|A|<=1.0
            if(A(j,i)>0.0)
                machine(j,i)=15.0+floor(AE(j,i))-1.0;%1 is 1~10 leading number
            else
                machine(j,i)=15+floor(AE(j,i))-1.0+1.0;%1 is minus sign;1 is 1~10 leading number
            end
            
        end
            err(j,i)=1.0*10^-machine(j,i);
    end
end


end

