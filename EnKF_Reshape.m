function TS = EnKF_Reshape(x_in, std_in, errs_in, settings, source)

P      = x_in(1:settings.nr_regions, :)';
E      = x_in(settings.nr_regions+1:2*settings.nr_regions, :)';
R      = x_in(2*settings.nr_regions+1:3*settings.nr_regions, :)';
DM     = x_in(3*settings.nr_regions+1:4*settings.nr_regions, :)';

P_std  = std_in(1:settings.nr_regions, :)';
E_std  = std_in(settings.nr_regions+1:2*settings.nr_regions, :)';
R_std  = std_in(2*settings.nr_regions+1:3*settings.nr_regions, :)';
DM_std = std_in(3*settings.nr_regions+1:4*settings.nr_regions, :)';

P_err  = errs_in(1:settings.nr_regions, :)';
E_err  = errs_in(settings.nr_regions+1:2*settings.nr_regions, :)';
R_err  = errs_in(2*settings.nr_regions+1:3*settings.nr_regions, :)';
DM_err = errs_in(3*settings.nr_regions+1:4*settings.nr_regions, :)';



TS = prepare_output(settings, source);

TS.Data.prec         = P;
TS.Data.prec_std     = P_std;
TS.Data.prec_error   = P_err;

TS.Data.evap         = E;
TS.Data.evap_std     = E_std;
TS.Data.evap_error   = E_err;

TS.Data.runoff       = R;
TS.Data.runoff_std   = R_std;
TS.Data.runoff_error = R_err;

TS.Data.twsc         = DM;
TS.Data.twsc_std     = DM_std;   
TS.Data.twsc_error   = DM_err;

end


function TS_out = prepare_output(settings, source);


TS_out = create_datastruct({'prec'; 'prec_std'; 'prec_error'; 'evap'; 'evap_std'; 'evap_error'; 'runoff'; 'runoff_std'; 'runoff_error'; 'twsc'; 'twsc_std'; 'twsc_error'; }, '2d_regions', 'double');


TS_out.DataInfo.title       = 'EnKF-based WB estimates';
TS_out.DataInfo.source      = source;
TS_out.DataInfo.references  = 'Lorenz, C., M. J. Tourian, B. Devaraju, N. Sneeuw, and H. Kunstmann (2015), Basin-scale runoff prediction: An Ensemble Kalman Filter framework based on global hydrometeorological data sets, Water Resour. Res., 51, 8450?8475, doi:10.1002/2014WR016794.';
TS_out.DataInfo.comment     = ['Predicted regions: ', num2str(settings.assim.pred_regions), '; predicted variables: ', num2str(settings.assim.pred_vars), ': Spinup: ', num2str(settings.assim.spinup)];
TS_out.DataInfo.institution = 'Institute for Meteorology and Climate Research (IMK-IFU), Karlsruhe Institute of Technology (KIT)';

TS_out.Variables.prec.long_name       = 'Precipitation';
TS_out.Variables.prec.standard_name   = 'precipitation_flux';
TS_out.Variables.prec.units           = 'mm/month';

TS_out.Variables.prec_std.long_name   = 'Ensemble standard deviation of precipitation';
TS_out.Variables.prec_std.units       = 'mm/month';

TS_out.Variables.prec_error.long_name  = 'Errors of precipitation estimates';
TS_out.Variables.prec_error.units       = 'mm/month';

TS_out.Variables.evap.long_name       = 'Evapotranspiration';
TS_out.Variables.evap.standard_name   = 'evapotranspiration_flux';
TS_out.Variables.evap.units           = 'mm/month';

TS_out.Variables.evap_std.long_name   = 'Ensemble standard deviation of evapotranspiration';
TS_out.Variables.evap_std.units       = 'mm/month';

TS_out.Variables.evap_error.long_name  = 'Errors of evapotranspiration estimates';
TS_out.Variables.evap_error.units       = 'mm/month';


TS_out.Variables.runoff.long_name     = 'Runoff';
TS_out.Variables.runoff.standard_name = 'runoff_flux';
TS_out.Variables.runoff.units         = 'mm/month';

TS_out.Variables.runoff_std.long_name = 'Ensemble standard deviation of runoff';
TS_out.Variables.runoff_std.units     = 'mm/month';

TS_out.Variables.runoff_error.long_name  = 'Errors of runoff estimates';
TS_out.Variables.runoff_error.units       = 'mm/month';

TS_out.Variables.twsc.long_name       = 'Total water storage chanes';
TS_out.Variables.twsc.standard_name   = 'twsc';
TS_out.Variables.twsc.units           = 'mm/month';

TS_out.Variables.twsc_std.long_name   = 'Ensemble standard deviation of total water storage changes';
TS_out.Variables.twsc_std.units       = 'mm/month';

TS_out.Variables.twsc_error.long_name  = 'Errors of TWSC estimates';
TS_out.Variables.twsc_error.units       = 'mm/month';

TS_out.Data.time    = settings.refvec;
TS_out.TimeStamp    = datenum(settings.refvec);
TS_out.Data.regions = settings.region_ids';

TS_out.Dimensions.regions = settings.nr_regions;
TS_out.Dimensions.time    = Inf;

end








