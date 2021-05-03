function x = conjgrad(A, b)
x0=ones(length(b),1);
r = b - A * x0;
if(r<1.0e-10)
    x=x0;
    return
end
x=x0;

p = r;
rsold = r' * r;
err0=1.0e-6;
err=1.0;
iters=1;
while(err>err0)
    Ap = A * p;
    alpha = rsold / (p' * Ap);
    err=max(abs(alpha*p));
    x = x + alpha * p;
    r = r - alpha * Ap;
    if r < 1.0e-10
        break;
    end
    rsnew = r' * r;
    
    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
    iters=iters+1;
end
end


