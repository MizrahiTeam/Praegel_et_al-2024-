
function [lick_all, stderr_lick_all, mean_lick_all, lick_learned,stderr_lick_learned, mean_lick_learned] = sort_freqs (level, l_all, stderr_l_all, mean_l_all)
if level == 4
    lick_all = [l_all(:,1) l_all(:,6) l_all(:,2) l_all(:,7) l_all(:,8)...
        l_all(:,9) l_all(:,3) l_all(:,10) l_all(:,4)];
    
    stderr_lick_all  = [stderr_l_all(:,1) stderr_l_all(:,6) stderr_l_all(:,2) stderr_l_all(:,7) stderr_l_all(:,8)...
        stderr_l_all(:,9) stderr_l_all(:,3) stderr_l_all(:,10) stderr_l_all(:,4)];
    
    mean_lick_all = [mean_l_all(:,1) mean_l_all(:,6) mean_l_all(:,2) mean_l_all(:,7) mean_l_all(:,8)...
        mean_l_all(:,9) mean_l_all(:,3) mean_l_all(:,10) mean_l_all(:,4)];
    
    lick_learned = [l_all(:,1)  l_all(:,2) l_all(:,3) l_all(:,4)];
    
    stderr_lick_learned  = [stderr_l_all(:,1) stderr_l_all(:,2) stderr_l_all(:,3) stderr_l_all(:,4)];
        
    mean_lick_learned  = [mean_l_all(:,1) mean_l_all(:,2)  mean_l_all(:,3) mean_l_all(:,4)];
    
elseif level == 3

lick_all = [l_all(:,1) l_all(:,4)  l_all(:,5) l_all(:,6)...
    l_all(:,7) l_all(:,8) l_all(:,2)];

stderr_lick_all  = [stderr_l_all(:,1) stderr_l_all(:,4)  stderr_l_all(:,5) stderr_l_all(:,6)...
    stderr_l_all(:,7) stderr_l_all(:,8) stderr_l_all(:,2)];

mean_lick_all  = [mean_l_all(:,1) mean_l_all(:,4)  mean_l_all(:,5) mean_l_all(:,6)...
    mean_l_all(:,7) mean_l_all(:,8) mean_l_all(:,2)];

lick_learned = [l_all(:,1)  l_all(:,2)];

stderr_lick_learned  = [stderr_l_all(:,1) stderr_l_all(:,2)];
    
mean_lick_learned  = [mean_l_all(:,1) mean_l_all(:,2)];

end
end