function [C, c, t, p] = nancov_matrix(X, Y, method, div, remmn)
% nancov_matrix computes the column-wise covariance matrix between X and Y.
% NaNs are treated as missing elements. The user can chose between two
% methods:
% - pairwise: for each pair of columns, the NaNs are removed separately,
%             which means that the covariances are based on different 
%             (e.g. time) periods
% - row:      if a row in X or Y contains NaN-values, it is removed
%             completely. Thus, the covariances are based on a consistent 
%             period, but some data might be unused
% Furthermore, the covariances are either normalized by N (div = 0) or 
% N - 1 (div = 1, default) to compute either the second moment matrix of 
% the observations about ther mean (div = 0) or the best unbiased estimate 
% of the covariance matrix (div = 1) if the observations are from a normal
% distribution.

if nargin < 5, remmn = 1; end
if nargin < 4, div = 1; end
if nargin < 3, method = 'pairwise'; end

% Check the dimensions of the input matrices
[r1, c1] = size(X);
[r2, c2] = size(Y);

if r1 ~= r2 | c1 ~= c2
    error('The input matrices must be of the same size')
end


if strcmp(method, 'pairwise')
    for i = 1:c1
        for j = 1:c2
            pair = [X(:, i), Y(:, j)];
            
            % Mask contains ones at the positions where either X or Y has
            % NaN-elements.
            mask              = zeros(size(pair));
            mask(isnan(pair)) = 1;
            
            % Compute the row-wise sum of mask. The rows where 
            % sum(mask, 2) > 0 are removed pair.
            mask              = sum(mask, 2);
            pair(mask > 0, :) = [];
            
            % Compute the number of observations after removing the 
            % NaN-elements
            N = size(pair, 1);
            
            if remmn == 1
                % Compute the mean
                mn_pair = mean(pair);
            
                % Put the columns of pair in x and y and remove the
                % mean from both vectors.
                x = pair(:, 1) - mn_pair(1);
                y = pair(:, 2) - mn_pair(2);
            else
                x = pair(:, 1);
                y = pair(:, 2);
            end
           
            % Compute the empirical covariance matrix between X and Y
            if div == 1
                C(i, j) = (1/(N-1))*x'*y;
                
            elseif div == 0
                C(i, j) = (1/(N))*x'*y;
            end
            c(i, j) = C(i, j)/(std(x)*std(y));
            t(i, j) = c(i, j)*sqrt((N - 2)/(1 - c(i, j)^2));
            T = tinv(0.95, N - 2);
            if abs(t(i, j)) > T || i == j
                p(i, j) = 1;
            else
                p(i, j) = 0;
            end
            
        end
    end
    
elseif strcmp(method, 'row')
    % Mask contains ones at the positions where either X or Y has
    % NaN-elements.
    mask           = zeros(size(X));
    mask(isnan(X)) = 1;
    mask(isnan(Y)) = 1;
    
    % Compute the row-wise sum of mask. The rows where sum(mask, 2) > 0 are
    % removed from both X and Y.
    mask = sum(mask, 2);
    X(mask > 0, :) = [];
    Y(mask > 0, :) = [];
    
    % Compute the number of observations after removing the NaN-elements
    N = size(X, 1);
    
    if remmn == 1
        % Compute the mean of X and Y
        mn_X = mean(X, 1);
        mn_Y = mean(Y, 1);
    
        % Remove the mean from the input matrices
        X = X - repmat(mn_X, N, 1);
        Y = Y - repmat(mn_Y, N, 1);
    end
    
    % Compute the empirical covariance matrix between X and Y
    if div == 1
        C = (1/(N-1))*(X'*Y);
    elseif div == 0
        C = (1/(N))*(X'*Y);
    end
end
    
