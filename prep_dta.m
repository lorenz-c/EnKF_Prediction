
function data_out = prep_dta(fnmes, variable, data_variable, id_var, ...
                    region_ids, sdte, edte, remove_ssnl, filter, ssnl_otpt)
   

if isstruct(fnmes)
    if isfield(fnmes, 'datadir')
        datadir = fnmes.datadir;
    else
        datadir = '';
    end
        
    if isstr(fnmes.(variable))
        fnmes_in{1} = [datadir, fnmes.(variable)];
    elseif iscell(fnmes.(variable))
        for i = 1:length(fnmes.(variable))
            fnmes_in{i, 1} = [datadir, fnmes.(variable){i}];
        end
    end
elseif iscell(fnmes)
    for i = 1:length(fnmes)
        fnmes_in{i, 1} = fnmes{i};
    end
elseif isstr(fnmes)
    fnmes_in{1} = fnmes;
end


for i = 1:length(fnmes_in)
    % Load a dataset
    data_in{i} = netcdf2datastruct(fnmes_in{i});

    % Truncate (or extend) the dataset to the desired time period
    data_in{i} = trunc_TS(data_in{i}, sdte, edte);
    
    % Remove the seasonal cycle --> Anomalies
    if length(remove_ssnl) > 1
        if remove_ssnl(i) == true
            data_in{i} = remsc(data_in{i});
        end
    else
        if remove_ssnl == true
            data_in{i} = remsc(data_in{i});
        end
    end
    
    % Compute the annual cycle
    if ssnl_otpt == true
        if remove_ssnl == true
            warning('prep_data.m: Annual cycle from anomalies!')
        end
        data_in{i} = ts_average(data_in{i}, 'monthly_lt');
    end
    
    % Create an empty array for the data
    data_out(:, :, i) = NaN(length(data_in{i}.TimeStamp), ...
                                                       length(region_ids));                                          
    % Check for the different regions
    for j = 1:length(region_ids)
        indx = find(data_in{i}.Data.(id_var) == region_ids(j));
        
        if ~isempty(indx)
            data_out(:, j, i) = data_in{i}.Data.(data_variable)(:, indx);
        end 
    end
    
    if length(filter) > 1
        if filter(i) == true
            tmp      = [data_out(1, :, i); ...
                        data_out(1:end, :, i); ...
                        data_out(end, :, i)];
                    
            data_out(:, :, i) = 1/4*tmp(1:end-2, :) + ...
                       1/2*tmp(2:end-1, :) + ...
                       1/4*tmp(3:end, :);
        end
    else
        if filter == true
            tmp      = [data_out(1, :, i); ...
                        data_out(1:end, :, i); ...
                        data_out(end, :, i)];
            
            data_out(:, :, i) =  1/4*tmp(1:end-2, :) + ...
                        1/2*tmp(2:end-1, :) + ...
                        1/4*tmp(3:end, :);
        end
    end
end



    
    


