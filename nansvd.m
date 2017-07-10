function out = nansvd(X)
    
% Compute the column average
clmn_mn = nanmean(X);

clmn_mn_mtrx = ones(length(X), 1)*clmn_mn;

X_filled = X;
X_filled(isnan(X)) = clmn_mn_mtrx(isnan(X));

brk = 0;
v_prev = 0;


while brk == 0
    % Do the SVD
    [U, P, D] = svd(X_filled - clmn_mn_mtrx);
    
    tmp = U*P*D + clmn_mn_mtrx;
    X_filled(isnan(X)) = tmp(isnan(X));
    
    clmn_mn = mean(X_filled);
    clmn_mn_mtrx = ones(length(X), 1)*clmn_mn;
    
    v = sum(diag(P))
    it = 0;
    maxit = 10;
    
    if (v - v_prev)/v_prev < 1e-5 || it == maxit
        brk = 1;
        it = it + 1;
    end
    
    v_prev = v;
    
end
keyboard
    
    out = 1;