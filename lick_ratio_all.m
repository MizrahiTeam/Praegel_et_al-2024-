       
function [lick_all, lick_learned, stderr_lick_all, mean_lick_all, stderr_lick_learned, mean_lick_learned] = lick_ratio_all (catch_freqs, learned_freqs,mice, stim,data)



for  m=1:length(mice)
    
            mouse = data(data.mouse_num== mice(m),:);

            
             th = int16(length(mouse.level)/3);
            mouse = mouse(th*2:th*3-1,:);
            [lick_ratio_per_stim] =  stim_lick_ratios (learned_freqs, catch_freqs, mouse, stim);
            lick_all(mice(m),:) = lick_ratio_per_stim;
end

lick_all(isnan(lick_all))= 0;


[mean_lick_all, stderr_lick_all] = mean_stderr (lick_all);

        lick_learned = lick_all (:,1:max(learned_freqs)) ;
        
[mean_lick_learned, stderr_lick_learned] = mean_stderr (lick_learned);
end 