function [lick_ratio_per_stim] =  stim_lick_ratios (learned_freqs, catch_freqs, mouse, stim)

                        lick_ratio_per_stim = [];
                        lick_ratio_all = [];

                        for stim = [learned_freqs catch_freqs];                                          
                            stim_trials = mouse(mouse.stimID == stim ,: );                     
                            total_stim = size(stim_trials, 1);
                            
                            if stim <= max(learned_freqs)/2
                                lick_ratio = sum(stim_trials.score == 0) ./ total_stim;          
                            elseif stim <= max(learned_freqs)
                                lick_ratio = sum(stim_trials.score == 1) ./ total_stim;
                            elseif stim >= max(learned_freqs)
                                 lick_ratio = sum(stim_trials.score == 6 ) ./ total_stim;
                            end
                                lick_ratio_per_stim(stim) = lick_ratio ;         
                        end                 

end 