function DFLP = GaussSeidel(A,B)
%under-relaxation Gauss-Siedel iteration to sovle [A][dP]=[RES], mentioned
%in paper 'Numerical simulation of heat, mass and momentum transport
% behaviours in directionally solidifying alloy castings under
% electromagnetic elds using an extended
% Direct-SIMPLE scheme', Daming Xu, 2004

%Created 2019-11-18

global lambda


[r,c]=size(A);
L=length(B);

if(r~=c);fprintf('Coefficient matrix should be a square\n');end
if(L~=c);fprintf('Matrix does not match RHS\n');end

err0=1.0e-6;
errmax=1.0;
err=zeros(L,1);
SumRow=zeros(L,1);
X=0.1*ones(L,1);%initial guess
iters=0;
while(errmax>=err0)
    for i=1:L
        X0=X(i);
        SumRow(i)=A(i,:)*X-A(i,i)*X(i);
       X(i)=(1.0-lambda)*X(i)+lambda*(B(i)-SumRow(i))/A(i,i);
       err(i)=abs((X0-X(i))/X0);
    end
    errmax=max(err);
    iters=iters+1;
end
fprintf('Under-relaxation Gauss-Seidel iterations: %5d\n',iters);
DFLP=X;
end

