function [settings] = EnKF_Settings

%--------------------------------------------------------------------------
%                            General parameter
%--------------------------------------------------------------------------

% Output directory and filename
settings.outdir = '/Users/lorenz-c/R_estimation/Results/'; 
settings.outnme = [settings.outdir, 'Run1', '.txt'];

% Temporal resolution of the input data
settings.temp_res = 'monthly';

% Name of the variable, which holds the catchment ids
settings.id_var = 'regions';

% Load the region IDs
tmp = assignvar('Region_IDs/IDs_1980_2002_ge_240_2005_2010_ge_48.mat');

% Ids of the areas of interest
%settings.region_ids   = [193 22 24 26 40 50 51 74 112 134 183 238 247 375 29];
settings.region_ids   = [193 24 26 51 112 134 183 247 375 29];
%settings.region_ids = [4 6 10 13 14 16 17 21 22 23 24 26 29 42 45 51 52 57 74 78 81 91 100 102 103]; 
settings.region_ids   = [29]; 
%settings.region_ids = [4 6 10 13 14 16 17 21 22 23 24 26 29];
% settings.region_names = tmp.Names_pred;

%--------------------------------------------------------------------------
%                 Parameter for the prediction function
%--------------------------------------------------------------------------
% List of datasets for deriving long-term statistics for the prediction
% function. The .mat-file must be a structure with the follwing fields:
% - data.prec       
% - data.evap       
% - data.runoff     
% - data.twsc
% - data.datadir
% The four variable-structures must contain cell-arrays with the filenames
% of the different datasets. The datadir-string contains the directory,
% where the data is stored.
settings.pred.data = assignvar('Data_subsets/long_term_data.mat');

% Start- and end-date for the statistics
settings.pred.sdte = [1980 01];
settings.pred.edte = [2002 12];

% Averaging of the covariance matrices
settings.pred.avg = 'mean';

% Select the prediction matrix approach
settings.pred.method = 1;

% Remove seasonal cycle for the least-squares prediction?
settings.pred.remsc = true;

% Apply the [1/4 1/2 1/4]-filter to the data?
settings.pred.filter_prec   = false;
settings.pred.filter_evap   = false;
settings.pred.filter_runoff = false;
settings.pred.filter_twsc   = false;

% Switch for removing insignificant covariances
settings.pred.filter_sig = false;

% Switch for "repairing" not positive-definite matrices
settings.pred.fix_cov = true;

% Enforce symmetry of the covariance matrices
settings.pred.makesymmetric = true;

I = ones(length(settings.region_ids), length(settings.region_ids));
E = eye(length(settings.region_ids), length(settings.region_ids));

settings.pred.loc_matrix = [E E E E;
                            E E E E;
                            E E I E;
                            E E E E];
%--------------------------------------------------------------------------
%                 Parameter for the control input
%--------------------------------------------------------------------------

% List of datasets for the control input
settings.cntrl.data = assignvar('Data_subsets/long_term_data.mat');

% Start- and end-date for the control input
settings.cntrl.sdte = [1980 01];
settings.cntrl.edte = [2002 12];

% Apply the [1/4 1/2 1/4]-filter to the data?
settings.cntrl.filter  = false;

%--------------------------------------------------------------------------
%                 Parameter for the observation function
%--------------------------------------------------------------------------
settings.obs.data   = assignvar('Data_subsets/obs_data.mat');
settings.obs.remsc  = false;
settings.obs.filter = false;




%--------------------------------------------------------------------------
%              Parameter for the observation covariance
%--------------------------------------------------------------------------
settings.obscov.data_back  = assignvar('Data_subsets/long_term_data.mat');
settings.obscov.sdte_back  = [1980 01];
settings.obscov.edte_back  = [2002 12];
settings.obscov.flt_back   = false;
settings.obscov.remsc_back = true;

settings.obscov.data_scle  = assignvar('Data_subsets/obs_data.mat');
settings.obscov.sdte_scle  = [2003 01];
settings.obscov.edte_scle  = [2013 12];
settings.obscov.flt_scle   = false;
settings.obscov.remsc_scle = true;

% Select a method for deriving errors 
settings.obscov.errs_P.method  = 'ensemble';

settings.obscov.errs_E.method  = 'ensemble';

settings.obscov.errs_R.method  = 'percentage';
settings.obscov.errs_R.data    = 'http://imk-ifu-thred1.imk-ifu.kit.edu:8080/thredds/dodsC/Area_averages/GRDC_Runoff_2014.nc';
settings.obscov.errs_R.varnme  = 'runoff';
settings.obscov.errs_R.value   = 0.1;
settings.obscov.errs_R.remsc   = false;
settings.obscov.errs_R.filter  = false;


settings.obscov.errs_DM.method = 'constant';
settings.obscov.errs_DM.value  = 20; 

settings.obscov.draw_errs = 0;


%--------------------------------------------------------------------------
%                     Parameter for additional inputs
%--------------------------------------------------------------------------
% Additional input can be provided by adding a dataset via e.g. the
% settings.obs.p_alt parameter. In this case, the parameter contains the
% filename of an additional dataset, which is loaded and added to the
% observation vector.
% If additional data is provided, one must also provide some errors for
% this dataset. Similar to the observation errors, it can be chosen between
% different methods for computing the errors:

% - settings.obscov.errs_p_alt.method = 'data'
%   --> The data from the file is used as errors
%   Mandatory fields: .data   (filename)
%                     .varnme (name of the variable in the netcdf-file)  
%                     .remsc  (remove the annual cycle)
%                     .filter (apply 1/4-1/2-1/4-filter)
%
% - settings.obscov.errs_p_alt.method = 'percentage'
%   --> The errors are the squared products between the data and a 
%       defined value
%   Mandatory fields: .data   (filename)
%                     .varnme (name of the variable in the netcdf-file)  
%                     .remsc  (remove the annual cycle)
%                     .filter (apply 1/4-1/2-1/4-filter)
%                     .value  (scalar with which the errors are multiplied)
%                
% - settings.obscov.errs_p_alt.method = 'constant'
%   --> The errors are set to constant values, where 
%       settings.obscov.errs_p_alt.value can be a single value or a vector
%       (i.e. one value for each region)
%   Mantatory fields: .value  (scalar or vector with errors)



settings.obs.p_alt                = [];
settings.obscov.errs_p_alt        = [];

settings.obs.e_alt                = [];
settings.obscov.errs_e_alt        = [];

% settings.obs.r_alt                = 'http://imk-ifu-thred1.imk-ifu.kit.edu:8080/thredds/dodsC/Area_averages/Altimetry_R_GRDC_Basins.nc';
% settings.obscov.r_alt_errs.method = 'data';
% settings.obscov.r_alt_errs.data   = 'http://imk-ifu-thred1.imk-ifu.kit.edu:8080/thredds/dodsC/Area_averages/Altimetry_R_GRDC_Basins.nc';
% settings.obscov.r_alt_errs.varnme = 'runofferror';
% settings.obscov.r_alt_errs.remsc  = false;
% settings.obscov.r_alt_errs.filter = false;

settings.obs.r_alt               = [];
settings.obscov.errs_r_alt       = [];


settings.obs.dm_alt               = [];
settings.obscov.errs_dm_alt       = [];


%--------------------------------------------------------------------------
%                    Parameter for the constraints
%--------------------------------------------------------------------------

settings.constraints.method  = 'percentage';
settings.constraints.data    = 'http://imk-ifu-thred1.imk-ifu.kit.edu:8080/thredds/dodsC/Area_averages/GRDC_Runoff_2014.nc';
settings.constraints.varnme  = 'runoff';
settings.constraints.value   = 0.5;
settings.constraints.remsc   = false;
settings.constraints.filter  = false;
settings.constraints.ssnl    = true;

% settings.constraints.method  = 'constant';
% settings.constraints.value   = [10 10 1];

%--------------------------------------------------------------------------
%                 Parameter for the assimilation
%--------------------------------------------------------------------------
% Start- and end-date of the assimilation period
settings.assim.sdte = [2003 01];
settings.assim.edte = [2013 12];

% Select a catchment for prediction
% settings.assim.pred_regions = 1:50;
settings.assim.pred_regions = 29;
% settings.assim.pred_names   = tmp.Names_val;

% Select a variable for prediction
% 1: Precipitation; 2: Evapotranspiration; 3: Runoff; 4: Water storage
settings.assim.pred_vars  = 3;

% Select the number of months for the training period
settings.assim.spinup = 0;

% Size of the EnKF-ensemble
settings.assim.enssize = 5000;

% Number of assimilation runs
nr_runs = 1;


%--------------------------------------------------------------------------
%                 Parameter for the validation
%--------------------------------------------------------------------------


settings.val.prec.data     = [];
settings.val.prec.varnme   = [];

settings.val.evap.data     = [];
settings.val.evap.varnme   = [];

settings.val.runoff.data   = 'http://imk-ifu-thred1.imk-ifu.kit.edu:8080/thredds/dodsC/Area_averages/GRDC_Runoff_2014.nc';
settings.val.runoff.varnme = 'runoff';

settings.val.twsc.data       = [];
settings.val.twsc.varnme     = [];

settings.val.perfmeasures  = {'nse', 'corr', 'rmse'};

settings.val.sdte          = [2003 01];
settings.val.edte          = [2012 12];

%--------------------------------------------------------------------------
%                DON'T CHANGE ANYTHING BELOW THIS LINE
%--------------------------------------------------------------------------
% Perform checks for the set variables
settings.refvec = dtevec(settings.assim.sdte, ...
                                           settings.assim.edte, 'monthly');
settings.nts               = size(settings.refvec, 1);                                       
settings.nr_yrs            = size(settings.refvec, 1)/12;                                      
settings.nr_regions        = length(settings.region_ids);                                  
settings.nr_assim_regions  = length(settings.assim.pred_regions);


if ~isempty(settings.assim.pred_regions) || ~isempty(settings.assim.pred_vars)
    if length(settings.assim.pred_vars) > 1 && ...
                                             settings.nr_assim_regions == 1
        settings.assim.pred_regions = ...
                            ones(length(settings.assim.pred_vars), 1) * ...
                                               settings.assim.pred_regions;
     
    elseif settings.nr_assim_regions > 1 && ...
                                      length(settings.assim.pred_vars) == 1
        settings.assim.pred_vars = ...
                            ones(settings.nr_assim_regions, 1) * ...
                                                  settings.assim.pred_vars;
        
    end

    for i = 1:settings.nr_assim_regions
        settings.assim.c_indx(i) = find(settings.region_ids == ...
                                           settings.assim.pred_regions(i));   
    end
    
    if size(settings.assim.pred_vars, 1) > 1
        settings.assim.pred_vars = settings.assim.pred_vars';
    end
    
    if size(settings.assim.pred_regions, 1) > 1
        settings.assim.pred_regions = settings.assim.pred_regions';
    end
end




settings.sources = {'EnKF_no_constraints'; ...
                    'EnKF_hard_constraints'; ...
                    'EnKF_soft_constraints'; ...
                    'EnKS_no_constraints'; ...
                    'EnKS_hard_constraints'; ...
                    'EnKS_soft_constraints'};













