function [TS_flt_c, TS_flt_h, TS_flt_s, TS_smth_c, TS_smth_h, TS_smth_s, final_stats, final_stats_anom] = EnKF_Runfile(settings)
%%
% -------------------------------------------------------------------------
%                  MAIN CODE OF THE ENKF-ASSIMILATION
% -------------------------------------------------------------------------

% Functions within the main file
% - EnKF_pred_mtrx.m       Computation of the least-squares prediction
%                          matrix and it's covariance
% - EnKF_cntrl_inpt.m      Computation of the cyclostationary control input
% - EnKF_obsvec.m          Computation of the observation vector
% - EnKF_obscov.m          Computation of the observation covariance
% - EnKF_WB_constraints.m  Computation of the waterbudget constraints
% - EnKF_reshape.m         Put the estimated paramter in a readable format
% - EnKF_val.m             Validation function

% Other functions (from the Matlab-GeoTS-toolbox)
% - netcdf2datastruct.m    Load a CF-formatted netcdf-file
% - dtevec.m               Compute a date vector from a given start- and
%                          end-date
% - assignvar.m            Load some .mat-file and assing the contents to a
%                          variable
% - trunc_TS.m             Truncate (or extend) a time-series to a given
%                          time-period
% - cov2corr.m             Transform a covariance- to a correlation matrix
% - structerrs.m           Compute errors between two data-structures
% - 


%% ------------------------------------------------------------------------
%                 1.   Compute the prediction function
%  ------------------------------------------------------------------------


%  1.1  Compute the least-squares prediction matrix and it's covariance
%       - b      Prediction matrix
%       - q_p    Prediction covariance
disp('Computing the prediction function.......')
[b, q_p] = EnKF_Prediction_Matrix(settings);
disp('Done!')

%  1.2  Compute the cyclostat. control input
%       - U      Control input
disp('Computing the control input.......')
U        = EnKF_Control_Input(settings);
disp('Done!')



%% ------------------------------------------------------------------------
%                 2.   Compute the observation function
%  ------------------------------------------------------------------------

%  2.1  Compute the observation vector
%       - Y      Observation vector
%       - H      Observation relation function
%       - Y_pred Observation vector where predicted elements are set to NaN
disp('Computing the observation vector.......')
[Y, H, Y_pred] = EnKF_Observation_Vector(settings);
disp('Done!')

%  2.2  Compute the observation errors
%       - Q_obs  Observation covariance
disp('Computing the observation covariance.......')
Q_obs          = EnKF_Observation_Covariance(settings);
disp('Done!')

%  2.3  Compute the waterbudget-constraints
%       - WB_var Waterbudget constraints
%       - Relation function for the constraints
disp('Computing the WB-constraints......')
[WB_var, G]    = EnKF_WB_Constraints(settings);
disp('Done!')


%% ------------------------------------------------------------------------
%                       3.   Initialize the EnKF
%  ------------------------------------------------------------------------
%  3.1  Choose the first element of the observation vector as start values
x_0             = Y(1:4*settings.nr_regions, 1);
x_0(isnan(x_0)) = 0;

%  3.2  Create empty cells for the covariance matrices
P_p             = cell(settings.nts, settings.nr_assim_regions);
P_c             = cell(settings.nts, settings.nr_assim_regions);
P_t1            = cell(settings.nts, settings.nr_assim_regions);
P_t2            = cell(settings.nts, settings.nr_assim_regions);

P_s             = cell(settings.nts, settings.nr_assim_regions);
P_s_t1          = cell(settings.nts, settings.nr_assim_regions);
P_s_t2          = cell(settings.nts, settings.nr_assim_regions);

%  3.3  Create empty cells for the state vectors
x_p             = cell(settings.nts, settings.nr_assim_regions);
x_c             = cell(settings.nts, settings.nr_assim_regions);
x_t1            = cell(settings.nts, settings.nr_assim_regions);
x_t2            = cell(settings.nts, settings.nr_assim_regions);

x_s             = cell(settings.nts, settings.nr_assim_regions);
x_s_t1          = cell(settings.nts, settings.nr_assim_regions);
x_s_t2          = cell(settings.nts, settings.nr_assim_regions);

x_c_mn          = cell(1, settings.nr_assim_regions);
x_t1_mn         = cell(1, settings.nr_assim_regions);
x_t2_mn         = cell(1, settings.nr_assim_regions);


%% 4.   Start filtering...
br    = waitbar(0, 'Filtering....');
for i = 2:settings.nts
    
    % Prepare the prediction matrix
    a = [-b eye(4*settings.nr_regions, 4*settings.nr_regions)] * U(:, i);
    a = repmat(a, 1, settings.assim.enssize);
    
    fprintf('Timestep %g; \n', i)   

    % Draw prediction errors
    errs_pred = mvnrnd(zeros(settings.assim.enssize, ...
                                             4*settings.nr_regions), q_p)';

    % If we're at the first time-step, we have to create the initial 
    % ensemble. This is the same for all predictions.
    if i == 2
        x_flt_c{1} = x_0*ones(1, settings.assim.enssize) + errs_pred;
        x_flt_h{1} = x_0*ones(1, settings.assim.enssize) + errs_pred;
        x_flt_s{1} = x_0*ones(1, settings.assim.enssize) + errs_pred;
        
        x_flt_p_anom{1} = x_0*ones(1, settings.assim.enssize) + errs_pred;
    end

	% Now, we have to prepare the observations and their errors
	h     = H;
	y     = Y_pred(:, i);
	q_o   = Q_obs(:, :, i);
        
 	% Check for missing values in the observation vector
	mvals = isnan(y);
        
  	if max(mvals) > 0
        % Remove the rows and columns from q_o which correspond to the
        % missing values
        q_o(mvals == 1, :) = [];
        q_o(:, mvals == 1) = [];
            
        % Remove the corresponding elements in the observation relation
        % matrix
        h(mvals == 1, :)   = [];

        % Finally, remove the missing observations from the observation
        % vector
        y(mvals == 1, :)   = [];

        % Sum up the number of missing values
        nr_mvals           = sum(mvals);
    else
        nr_mvals = 0;
    end  

    % 1. Draw observation errors
    errs_obs  = mvnrnd(zeros(settings.assim.enssize, length(y)), q_o)';    
    
    
    % ---------------------------------------------------------------------
    %                          CONSTRAINTS-SECTION
    % ---------------------------------------------------------------------
    % 2. Augment the observation relation matrix
    h_aug     = [h; G];
    
    % 3. Augment the observation covariance (hard constraints)
    q_h_aug   = [q_o zeros(size(q_o, 1), settings.nr_regions);
                 zeros(settings.nr_regions, size(q_o, 2) + ...
                                                     settings.nr_regions)];
              
    % 4. Augment the observation covariance (soft constraints)
    q_s_aug   = [q_o              zeros(size(q_o, 1), settings.nr_regions);
                 zeros(settings.nr_regions, size(q_o, 2)) ...
                                                       diag(WB_var(i, :))];
               
    y_aug     = [y; zeros(settings.nr_regions, 1)];

    % PERTUBE THE OBSERVATIONS
    y_o       = repmat(y, 1, settings.assim.enssize)     + errs_obs;
    y_o_aug   = repmat(y_aug, 1, settings.assim.enssize) + [errs_obs; ...
                       zeros(settings.nr_regions, settings.assim.enssize)];
        
    % PREDICT THE STATE 
    x_flt_p{i}   = b*x_flt_c{i-1} + a + errs_pred;
    x_flt_h_p{i} = b*x_flt_h{i-1} + a + errs_pred;
    x_flt_s_p{i} = b*x_flt_s{i-1} + a + errs_pred;
    
    x_flt_p_anom{i} = b*x_flt_p_anom{i-1} + errs_pred;
    
    % Compute the predicted covariance
    P_flt_p{i}   = cov(x_flt_p{i}');
    P_flt_h_p{i} = cov(x_flt_h_p{i}');
    P_flt_s_p{i} = cov(x_flt_s_p{i}');
        
    % Compute the Kalman Gain
    s       = h*P_flt_p{i}*h' + q_o;
    s_h_aug = h_aug*P_flt_h_p{i}*h_aug' + q_h_aug;
    s_s_aug = h_aug*P_flt_s_p{i}*h_aug' + q_s_aug;

    K       = P_flt_p{i}*h'*inv(s);
    K_h_aug = P_flt_h_p{i}*h_aug'*inv(s_h_aug);
    K_s_aug = P_flt_s_p{i}*h_aug'*inv(s_s_aug);
            
    % Correct the state
    x_flt_c{i} = x_flt_p{i} + K*(y_o - h*x_flt_p{i});
    x_flt_h{i} = x_flt_h_p{i} + K_h_aug*(y_o_aug - h_aug*x_flt_h_p{i});
    x_flt_s{i} = x_flt_s_p{i} + K_s_aug*(y_o_aug - h_aug*x_flt_s_p{i});
        
    % Compute the ensemble mean of the corrected and predicted state
    x_flt_c_mn(:, i) = mean(x_flt_c{i}, 2);
    x_flt_h_mn(:, i) = mean(x_flt_h{i}, 2);
    x_flt_s_mn(:, i) = mean(x_flt_s{i}, 2);
        
    % Compute the ensemble standard deviation of the corrected and 
    % predicted state (which is an empirical approximation of the errors in
    % the estimated parameters).
    std_flt_c(:, i)  = std(x_flt_c{i}, [], 2);
    std_flt_h(:, i)  = std(x_flt_h{i}, [], 2);
    std_flt_s(:, i)  = std(x_flt_s{i}, [], 2);

    % Compute the covariance matrices
    P_flt_c{i} = cov(x_flt_c{i}');
    P_flt_h{i} = cov(x_flt_h{i}');
    P_flt_s{i} = cov(x_flt_s{i}');

    % Compute the errors
    errs_flt_c(:, i) = sqrt(diag(P_flt_c{i}));
    errs_flt_h(:, i) = sqrt(diag(P_flt_h{i}));
    errs_flt_s(:, i) = sqrt(diag(P_flt_s{i}));
    
    waitbar(i/settings.nts)
end
close(br)


%% Apply the RTS-Algorithm -> Ensemble Kalman Smoother
% Initialize the smoothed solutions
x_smth_c(i) = x_flt_c(i);
x_smth_h(i) = x_flt_h(i);
x_smth_s(i) = x_flt_s(i);

P_smth_c(i) = P_flt_c(i);
P_smth_h(i) = P_flt_h(i);
P_smth_s(i) = P_flt_s(i);

x_smth_c_mn = x_flt_c_mn;
x_smth_h_mn = x_flt_h_mn;
x_smth_s_mn = x_flt_s_mn;

std_smth_c  = std_flt_c;
std_smth_h  = std_flt_h;
std_smth_s  = std_flt_s;

errs_smth_c = errs_flt_c;
errs_smth_h = errs_flt_h;
errs_smth_s = errs_flt_s;
        
br = waitbar(0, 'Smoothing....');
for i = settings.nts - 1:-1:2

    k_s            = P_flt_c{i}*b'*inv(P_flt_p{i+1});
    k_s_h_aug      = P_flt_h{i}*b'*inv(P_flt_h_p{i+1});
    k_s_s_aug      = P_flt_s{i}*b'*inv(P_flt_s_p{i+1});

    P_smth_c{i} = P_flt_c{i} + k_s*(P_smth_c{i+1} - P_flt_p{i+1})*k_s';
    P_smth_h{i} = P_flt_h{i} + k_s_h_aug*(P_smth_h{i+1} - P_flt_h_p{i+1})*k_s_h_aug';
    P_smth_s{i} = P_flt_s{i} + k_s_s_aug*(P_smth_s{i+1} - P_flt_s_p{i+1})*k_s_s_aug';
        
    x_smth_c{i} = x_flt_c{i} + k_s*(x_smth_c{i+1} - x_flt_p{i+1});
    x_smth_h{i} = x_flt_h{i} + k_s_h_aug*(x_smth_h{i+1} - x_flt_h_p{i+1});
    x_smth_s{i} = x_flt_s{i} + k_s_h_aug*(x_smth_s{i+1} - x_flt_s_p{i+1});
    
    x_smth_c_mn(:, i) = mean(x_smth_c{i}, 2);
    x_smth_h_mn(:, i) = mean(x_smth_h{i}, 2);
    x_smth_s_mn(:, i) = mean(x_smth_s{i}, 2);
        
    errs_smth_c(:, i) = sqrt(diag(P_smth_c{i}));
    errs_smth_h(:, i) = sqrt(diag(P_smth_h{i}));
    errs_smth_s(:, i) = sqrt(diag(P_smth_s{i}));
    
    std_smth_c(:, i)  = std(x_smth_c{i}, [], 2);
    std_smth_h(:, i)  = std(x_smth_h{i}, [], 2);
    std_smth_s(:, i)  = std(x_smth_s{i}, [], 2);
    
    waitbar((settings.nts-i)/settings.nts)
end
close(br);


% %% 6.   Re-arrange the observations and states in a data-structure
TS_flt_c  = EnKF_Reshape(x_flt_c_mn, std_flt_c, errs_flt_c, settings, ...
                                                      settings.sources{1});
                     
TS_flt_h  = EnKF_Reshape(x_flt_h_mn, std_flt_h, errs_flt_h, settings, ...
                                                      settings.sources{2});
TS_flt_s  = EnKF_Reshape(x_flt_s_mn, std_flt_s, errs_flt_s, settings, ...
                                                      settings.sources{3});

TS_smth_c = EnKF_Reshape(x_smth_c_mn, std_smth_c, errs_smth_c, ...
                                            settings, settings.sources{4});
TS_smth_h = EnKF_Reshape(x_smth_h_mn, std_smth_h, errs_smth_h, ...
                                            settings, settings.sources{5});
TS_smth_s = EnKF_Reshape(x_smth_s_mn, std_smth_s, errs_smth_s, ...
                                            settings, settings.sources{6});

                                    
[final_stats, final_stats_anom, P_val, E_val, R_val, TWSC_Val] = ...
       EnKF_Validation(settings, TS_flt_c, TS_flt_h, TS_flt_s, TS_smth_c, TS_smth_h, TS_smth_s);

EnKF_Show_Output(settings, final_stats, final_stats_anom)



