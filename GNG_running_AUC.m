function [r_AUC, r_AUC_abs, r_AUC_shuf, r_AUC_shuf_abs]...
    = GNG_running_AUC (GNG_rec_all_cell,a,b, A_k_ge, B_k_ge, trial_samples, it_size, run_window, window_length_ms, n_shuffle)
%Input:
% a = index number of trial type a e.g. hit = 1 
% b = index number of trial type b e.g FA = 2
% A_k_ge = Fr or binARRAY of trial type e.g hit
% B_k_ge =Fr or bin ARRAY of trial type e.g FA
% trial_samples = sample of trials if trials are unequal in size
% it_size = iterations of running window
% run_window = length of a single running window
% window_length_ms = = total length of running window
% n_shuffle = number of shuffled iteration

%Output:
% r_AUC = running AUC 0-1
% r_AUC_abs = running AUC 0.5-1
% r_AUC_shuf = shuffled AUC 0-1
% r_AUC_shuf_abs = shuffled AUC 0.5 -1 

% allocate the cell array
r_AUC= cell(numel(GNG_rec_all_cell),length(a)) ;
r_AUC_abs = cell(numel(GNG_rec_all_cell),length(a)) ;
r_AUC_shuf = cell(numel(GNG_rec_all_cell),length(a)) ;
r_AUC_shuf_abs = cell(numel(GNG_rec_all_cell),length(a)) ;

 tic
 for g   = 1:numel(GNG_rec_all_cell) % run per group
     for e = 1:length(a)

         % eallocate the matrix
        r_AUC{g,e} = nan(trial_samples,size(A_k_ge{e,g},2),length(1:it_size:( window_length_ms- run_window))) ;
        r_AUC_abs{g,e} = nan(trial_samples,size(A_k_ge{e,g},2),length(1:it_size:( window_length_ms- run_window))) ;
        r_AUC_shuf{g,e} = nan(trial_samples,size(A_k_ge{e,g},2),length(1:it_size:( window_length_ms- run_window))) ;
        r_AUC_shuf_abs{g,e} = nan(trial_samples,size(A_k_ge{e,g},2),length(1:it_size:( window_length_ms- run_window))) ;
         
         for c = 1:size(A_k_ge{e,g},2)
             for k = 1:trial_samples
                 
                 % retrieve the smoothed fr per trials
                 fr_trials_A =  squeeze(A_k_ge{e,g}{c}(k,:,:)) ;
                 fr_trials_B =  squeeze(B_k_ge{e,g}{c}(k,:,:)) ;
                 
                 for bin = 1:it_size:(length(fr_trials_A)-run_window) % run per bin per window
                     bin_pos = 1 ;
                     
                     if bin > 1
                         bin = bin - 1 ;
                         bin_pos = (bin / it_size) +1  ;
                     end
                     
                     A = [] ;
                     B = [] ;
                     
                     A = mean(fr_trials_A(:,bin:bin+run_window-1)') ;
                     B = mean(fr_trials_B(:,bin:bin+run_window-1)') ;
                     
                     % cat vectors
                     C = [A' ;B'] ;
                     
                     % assign ID
                     index = [ones(1,numel(A)),zeros(1,numel(B))];
                     %compare AUC
                     [~,~,T,AUC] = perfcurve(logical(index),C,'true');
                     
                     r_AUC{g,e}(k,c,bin_pos) = AUC ; % 0-1
                     
                     r_AUC_abs{g,e}(k,c,bin_pos) = 0.5 + (abs(AUC - 0.5)) ; %absolute 0.5-1
                     
                     AUC_shuffled_trail = zeros(1,n_shuffle);
                     index = [ones(1,numel(A)),zeros(1,numel(B))]; % assign ID
                     
                     for kk = 1:n_shuffle  % iteration of shuffle
                         shuffled_index = index(randperm(length(index))); % permute shufflede ID
                         [~,~,~,AUC_shuffled_trail(kk)] = perfcurve(logical(shuffled_index), C,'true');
                     end
                     
                     AUC_h_shuff = mean(AUC_shuffled_trail); % extract mean of all shuffeled iterations
                     
                    r_AUC_shuf{g,e}(k,c,bin_pos) = AUC_h_shuff ;  % 0-1
                    r_AUC_shuf_abs{g,e}(k,c,bin_pos) = 0.5 + (abs(AUC_h_shuff - 0.5)) ;  %absolute 0.5-1
                     
                 end
             end
         end
     end
 end
end 