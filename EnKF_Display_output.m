
function [] = EnKF_display_output(settings, stats, stats_anom)
    
if ~isempty(settings.assim.c_indx)
    % Format for the performance metrics
    frmt = '%5s%1s%4.2f%1s%4.2f%1s';
    % --> Found some predicted regions
    tmp  = fieldnames(stats);
    % --> Cell array with the names of the performance measures
    
    for i = 1:length(settings.assim.c_indx)
        % Loop over the predicted regions
        disp(' ')
        disp(['Region ID: ', num2str(settings.region_ids(settings.assim.c_indx(i)))])
        
        for j = 1:length(tmp)
            % Loop over the predicted variables
            disp(['Statistics for ', tmp{j}])
            
            for k = 1:length(settings.sources)
                % Loop over the different configurations
                
                for l = 1:length(settings.val.perfmeasures)
                    if l == 1
                        strout{1, 1} = sprintf('%-23s', [settings.sources{k}, ': ']);
                    end
                    
                    out            = {settings.val.perfmeasures{l}, ': ', stats.(tmp{j})(k).(settings.val.perfmeasures{l})(i), ' (', stats_anom.(tmp{j})(k).(settings.val.perfmeasures{l})(i), ')'};
                    strout{1, l+1} = sprintf(frmt, out{:});
                end
                disp(strjoin(strout))
            end
        end
    end
end


% 
%             
%             
%             for k = 1:length(settings.sources)
%                 for l = 1:length(settings.val.perfmeasures)         
%                     perf_out{1, l} = strcat(settings.val.perfmeasures{l}, ': ', num2str(final_stats.(tmp{j})(k).(settings.val.perfmeasures{l}), '%4.2f'), ' (', 
%                 end
%                 fprintf([settings.sources{k}, ': ', strjoin(perf_out) '\n'])
%                 clear perf_out
%             end
%             
%             fprintf(['Statistics for ', tmp{j} ' (anomalies) \n'])
%             
%             for k = 1:length(settings.sources)
%                 for l = 1:length(settings.val.perfmeasures)         
%                     perf_out{1, l} = strcat(settings.val.perfmeasures{l}, ': ', num2str(final_stats_anom.(tmp{j})(k).(settings.val.perfmeasures{l}), '%10.5f'), '; ');
%                 end
%                 fprintf([settings.sources{k}, ': ', strjoin(perf_out) '\n'])
%                 clear perf_out
%             end       
%         end
%     end
% end