      
function [data] = max_trial_dist (data, max_ISI)
 
            textdata= data.time ;
            time_file = char(textdata(:,:));
            hour = str2num(time_file(:,1:2));
            minutes = str2num(time_file(:,4:5));
            sec = str2num(time_file(:,7:8));
           
            hour_diff_trial = diff (hour) ;
            min_diff_trial = diff (minutes) ;
            sec_diff_trial = diff(sec);
            
            for L = 1:length(sec_diff_trial)
                if hour_diff_trial(L) < 0 
                   hour_diff_trial(L) = dist(hour_diff_trial(end),-24) ;
                end
                
                if min_diff_trial(L) < 0 
                   min_diff_trial(L) = dist(min_diff_trial(end),-60) ;
                end
                
                if sec_diff_trial(L) < 0 
                   sec_diff_trial(L) = dist(sec_diff_trial(end),-60) ;
                end
            end 
            
            for L = 1:length(sec_diff_trial)
                if hour_diff_trial(L) > 1
                    min_diff_trial(L) = min_diff_trial(L)* (hour_diff_trial(L)*60);
                end
            end
            
            for L = 1:length(sec_diff_trial)
                if min_diff_trial(L) > 1
                    sec_diff_trial(L) = sec_diff_trial(L)* (min_diff_trial(L)*60);
                end
            end
        
        index = find (sec_diff_trial <= max_ISI);
        data = data(index,:);
               
end 