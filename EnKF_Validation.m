
function [stats, stats_anom, varargout] = EnKF_Validation(settings, varargin)

if ~isempty(settings.val.sdte)
    sdte = settings.val.sdte;
else
    if isempty(settings.assim.spinup)
        sdte = settings.refvec(1, :);
    else
        sdte = settings.refvec(settings.assim.spinup+1, :);
    end
end

if ~isempty(settings.val.edte)
    edte = settings.val.edte;
else
    edte = settings.refvec(end, :);
end

refvec = dtevec(sdte, edte, settings.temp_res);

TS_in = varargin;

if ~isempty(settings.val.prec.data)
    [stats.prec, stats_anom.prec, varargout{1}] = val_dtasets(settings, 'prec', refvec, TS_in);
else
    varargout{1} = [];
end

if ~isempty(settings.val.evap.data)
    [stats.evap, stats_anom.evap, varargout{2}] = val_dtasets(settings, 'evap', refvec, TS_in);
else
    varargout{2} = [];
end

if ~isempty(settings.val.runoff.data)
    [stats.runoff, stats_anom.runoff, varargout{3}]= val_dtasets(settings, 'runoff', refvec, TS_in);
else
    varargout{3} = [];
end

if ~isempty(settings.val.twsc.data)
    [stats.twsc, stats_anom.twsc, varargout{4}] = val_dtasets(settings, 'twsc', refvec, TS_in);
else
    varargout{4} = [];
end



end



function [stats_out, stats_anom_out, TS_out] = val_dtasets(settings, varnme, refvec, TS_in)

    if isstr(settings.val.(varnme).data)
        Val_dta{1}  = netcdf2datastruct(settings.val.(varnme).data);
        Val_dta{1}  = findregions(Val_dta{1}, settings.region_ids);
    elseif iscell(settings.val.(varnme).data)
        for i = 1:length(settings.val.(varnme).data)
            Val_dta{i} = netcdf2datastruct(settings.val.(varnme).data{i});
            Val_dta{i} = findregions(Val_dta{i}, settings.region_ids);
        end 
    end
    
    if length(Val_dta) == 1
        for j = 1:length(TS_in)
            stats_out(j, :) = structerrs(Val_dta{1}, TS_in{j}, ...
                                {settings.val.(varnme).varnme; varnme}, ...
                                     settings.val.perfmeasures, 'time', ...
                                       [], [refvec(1, :); refvec(end, :)]);
                                   
            stats_anom_out(j, :) = structerrs(Val_dta{1}, TS_in{j}, ...
                                {settings.val.(varnme).varnme; varnme}, ...
                                     settings.val.perfmeasures, 'time', ...
                                 [], [refvec(1, :); refvec(end, :)], true);                        
        end
    else
        for i = 1:length(Val_dta)
            for j = 1:length(TS_in)
                stats_out{i}(j, :) = structerrs(Val_dta{i}, ...
                   TS_in{j}, {settings.val.(varnme).varnme; varnme}, ...
                                     settings.val.perfmeasures, 'time', ...
                                       [], [refvec(1, :); refvec(end, :)]);
                stats_anom_out{i}(j, :) = structerrs(Val_dta{i}, ...
                   TS_in{j}, {settings.val.(varnme).varnme; varnme}, ...
                                     settings.val.perfmeasures, 'time', ...
                                 [], [refvec(1, :); refvec(end, :)], true);                   
            end
        end
    end
   
    TS_out = Val_dta;
%
end



