function X=PardisoSolver(A,B)
% This function utilizes PARDISO to solve large sparse, real non-symmetric matrix.
%Modifies from PARDISO example exampleunsym.m

verbose = false;

% Initialize the PARDISO internal data structures. We've told PARDISO to
% handle real non-symmetric matrices using the sparse direct solver.
info = pardisoinit(11,0);

% Analyze the matrix and compute a symbolic factorization.
info = pardisoreorder(A,info,verbose);
%fprintf('The factors have %d nonzero entries.\n',info.iparm(18));

% Compute the numeric factorization.
info = pardisofactor(A,info,verbose);

% Compute the solutions X using the symbolic factorization.
[X info] = pardisosolve(A,B,info,verbose);
iters=info.iparm(7);

% Compute the residuals.
R = max(abs(A*X - B));
err = max(R(:));

fprintf('Max LAS error:%E\n',err);%LAS=large algebraic system
fprintf('Max LAS iterations:%3d\n',iters);

% Free the PARDISO data structures.
pardisofree(info);
% clear info

end

