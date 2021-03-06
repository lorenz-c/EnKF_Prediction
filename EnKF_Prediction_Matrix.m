function [B, Q] = EnKF_Prediction_Matrix(settings)
% This function computes the prediction function for the EnKF-assimilation.


% Load some data
P = EnKF_Prepare_Data(settings.pred.data, ...
             'prec', ...
             'prec', ...
             settings.id_var, ...
             settings.region_ids, ...
             settings.pred.sdte, ...
             settings.pred.edte, ...
             settings.pred.remsc, ...
             settings.pred.filter_prec, ...
             false);
%               
E = EnKF_Prepare_Data(settings.pred.data, ...
             'evap', ...
             'evap', ...
             settings.id_var, ...
             settings.region_ids, ...
             settings.pred.sdte, ...
             settings.pred.edte, ...
             settings.pred.remsc, ...
             settings.pred.filter_evap, ...
             false);
         
R = EnKF_Prepare_Data(settings.pred.data, ...
             'runoff', ...
             'runoff', ...
             settings.id_var, ...
             settings.region_ids, ...
             settings.pred.sdte, ...
             settings.pred.edte, ...
             settings.pred.remsc, ...
             settings.pred.filter_runoff, ...
             false);
         
DM = EnKF_Prepare_Data(settings.pred.data, ...
             'twsc', ...
             'twsc', ...
             settings.id_var, ...
             settings.region_ids, ...
             settings.pred.sdte, ...
             settings.pred.edte, ...
             settings.pred.remsc, ...
             settings.pred.filter_twsc, ...
             false);

nr_P  = size(P, 3);
nr_E  = size(E, 3);
nr_R  = size(R, 3);
nr_DM = size(DM, 3);

sig   = zeros(4*settings.nr_regions, 4*settings.nr_regions);
sig_d = zeros(4*settings.nr_regions, 4*settings.nr_regions);

% Compute the big matrix which holds all the data
if strcmp(settings.pred.avg, 'mean')
    m = 1;
    for i = 1:size(P, 3)
        for j = 1:size(E, 3)
            for k = 1:size(R, 3)
                for l = 1:size(DM, 3)
        
                    Bigmat = [squeeze(P(:, :, i)) squeeze(E(:, :, j)) ...
                                 squeeze(R(:, :, k)) squeeze(DM(:, :, l))];

                    [sig_tmp, crr, t, p_tmp] = ...
                           nancov_matrix(Bigmat, Bigmat, 'pairwise', 1, 1);
                    [sig_d_tmp, crr, t, p_d_tmp] = ...
                           nancov_matrix(Bigmat(2:end, :), ...
                                     Bigmat(1:end-1, :), 'pairwise', 1, 1);
                    


                    if settings.pred.filter_sig == true
                        sig_tmp(p_tmp == 0) = 0;
                        sig_d_tmp(p_d_tmp == 0) = 0;
                    end
                    
                    
                    sig   = sig + sig_tmp;
                    sig_d = sig_d + sig_d_tmp;
                    
                    m = m + 1;
                end
            end
        end
    end
  
    sig = sig/m;
    sig_d = sig_d/m;

% 
%     sig   = mean(sig, 3);
%     sig_d = mean(sig_d, 3);    
else
    error('Covariance averaging method not known!')
end

            
if settings.pred.method == 1     
% -------------------------------------------------------------------------
%             Correlations between catchments and variables
% -------------------------------------------------------------------------
    B     = sig_d*inv(sig);
    Q     = sig - sig_d*inv(sig)*sig_d';   % --> always symmetric!
     
elseif settings.pred.method == 2
% -------------------------------------------------------------------------
%         Correlations between catchments but not between variables
% -------------------------------------------------------------------------   
    I = ones(settings.nr_regions, settings.nr_regions);
    O = zeros(settings.nr_regions, settings.nr_regions);


    dummy = [I O O O;
             O I O O;
             O O I O;
             O O O I];
        
    sig   = sig.*dummy;
    sig_d = sig_d.*dummy;
   
    B     = sig_d*inv(sig);
    Q     = sig - sig_d*inv(sig)*sig_d'; 

elseif settings.pred.method == 3
% -------------------------------------------------------------------------
%           No Correlations between catchments BUT between variables
% -------------------------------------------------------------------------    
    
    dummy = eye(length(settings.region_ids), length(settings.region_ids));
    dummy = repmat(dummy, 4, 4);

    sig   = sig.*dummy;
    sig_d = sig_d.*dummy;
    
    B     = sig_d*inv(sig);
    Q     = sig - sig_d*inv(sig)*sig_d'; 
    
elseif settings.pred.method == 4
    % A localization matrix is provided
  
    sig   = sig.*settings.pred.loc_matrix;
    sig_d = sig_d.*settings.pred.loc_matrix;
    
    B     = sig_d*inv(sig);
    Q     = sig - sig_d*inv(sig)*sig_d'; 
    
else
    error('Unknown prediction matrix method!')
end


if settings.pred.makesymmetric == true
    % Q might not be 100% symmetric (due to rounding errors...). So let's
    % fix that. (This does not change the correlation structure in Q, but 
    % corrects for numerical (e.g. rounding) errors).

    if iscell(Q)
        for i = 1:length(Q)
            Q{i} = 1/2*(Q + Q');
        end
    else
        Q = 1/2*(Q + Q');
    end
end

if settings.pred.fix_cov == true 
    % If Q is not positive-semi-definite, we can try to repair the matrix.
    % However, depending on the "level" of non-positive-definiteness, this 
    % can seriously damange the covariance structure and thus produce
    % unreasonable results... 

    EPS = 10^-6;
    ZERO = 10^-10;

    if iscell(Q)
        for i = 1:length(Q)
            [V, D]       = eig(Q{i});
            D(D <= ZERO) = EPS;
            Q{i}         = V*diag(diag(D))*V';
        end
    else
        [V, D]       = eig(Q);
        D(D <= ZERO) = EPS;
        Q            = V*diag(diag(D))*V';
    end
end

