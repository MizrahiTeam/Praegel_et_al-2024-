

function [A_k_ge, B_k_ge, idx_k_ge, T_k_ge, deasy_k_ge, dhard_k_ge, rec_idx, neuron_idx, mouse_idx]...
    = GNG_t_per_n_AUC (FR_array, GNG_rec_all_cell, behavior, areas, a, b,  trial_samples, min_n_trials)
% Input: 
% FR_array = output of previous analyis
% GNG_rec_all_cell = neuronal and behavioral ouput of  all recordings 
% behavior =   = ouput GNG_neuro_behavior
% areas= define according to rec location
% a= trial types e.g. hit compared to FA  (a = [1 2]) 
% b = trial types e.g. hit compared to FA  (a = [1 2]) 
% see GNG_rec_all_cell.eventTimes 
% trial_samples = how many times you want to take trial samples per
% neuron
% min_n_trials = minimal number of trials per neuron 

% Output:
% A_k_ge = fr or binarray of trial type A per trial
% B_k_g e =fr or binarray of trial type B per trial 
% idx_k_ge = idx of trials in rec
% T_k_ge = frarray per trial type o fhit,miss, fa and cr
% deasy_k_ge = behavir d' easy
% dhard_k_ge= behavir d' hard
% beasy_k_ge= behavir lick bias  easy
% bhard_k_ge= behavir lick bias hard

%% random perm of at least n: min_trials  of trial type A and trial type B per E(comparsion)
for g   = 1:numel(GNG_rec_all_cell) % run per group
    Recs = 1:numel(GNG_rec_all_cell{1,g}) ;
    for i = 1:length(Recs) % run per recording

        dprime_easy_rec = behavior.dprime_easy(g,i) ;
        dprime_hard_rec = behavior.dprime_hard(g,i); 
        cbias_easy_rec = behavior.cbias_easy(g,i) ;
        cbias_hard_rec = behavior.cbias_hard(g,i) ;


        for e = 1:length(a)
            % extract trial type a index per rec
            idx_a =  GNG_rec_all_cell{1, g}(i).eventTimes.tr_resp_idx(a(e),:);
            idx_a(isnan(idx_a)) = [] ;
            idx_a = idx_a(1:end-1) ;
            % extract trial type b index per rec
            idx_b =  GNG_rec_all_cell{1, g}(i).eventTimes.tr_resp_idx(b(e),:);
            idx_b(isnan(idx_b)) = [] ;
            idx_b = idx_b(1:end-1) ;
            
            % maximal shared number of trials
            size_max =  min([size(idx_a,2) size(idx_b,2)]) ;
            
            for area = 1:length(areas)
                
                for c = 1:size(FR_array{1,g},1)
                    % if c exists 
                    if ~isempty(FR_array{1,g}{c,i,area})
                        
                            idx_a_k = [] ;
                            idx_b_k = [];
                            
                            % take data sample of the maximal number of trials
                            % betwen a and b
                            idx_a_k = sort(datasample(idx_a,size_max)) ;
                            idx_b_k = sort(datasample(idx_b,size_max)) ;
                            
                            % if a trial type larger than min trial types
                            if length(idx_a_k) >= min_n_trials
                                % if btrial type larger than min trial types
                                if  length(idx_b_k) >= min_n_trials
                                    
                                    for k = 1:trial_samples
                                        %extract FRArray according to data sample
                                        A_frArray_k = FR_array{1,g}{c,i,area}(idx_a_k,:) ;
                                        B_frArray_k = FR_array{1,g}{c,i,area}(idx_b_k,:) ;
                                        
                                        %add to 3D mat of all trial type samples
                                        A_k_all{g,area}{i}{c,e}(k,:,:) = A_frArray_k(:,:) ;
                                        B_k_all{g,area}{i}{c,e}(k,:,:) = B_frArray_k(:,:) ;
                                    end
                                    
                                    % create heatmap of spike in differen trial types for each neuron
                                    %allocate
                                    idx_k_all{g,area}{i}{c,e} = nan(4,size(GNG_rec_all_cell{1, g}(i).eventTimes.tr_resp_idx,2)) ;
                                    for tv = [13 14 15 16]
                                        idx_tv = [] ;
                                        idx_tv =  GNG_rec_all_cell{1, g}(i).eventTimes.tr_resp_idx(tv,:); % hit all
                                        idx_tv(isnan(idx_tv)) = [] ;
                                        idx_tv = idx_tv(1:end-1) ;
                                        idx_k_all{g,area}{i}{c,e}(tv,1:length(idx_tv)) = idx_tv ;
                                    end
                                    
                                    T_k_all{g,area}{i}{c,e} = FR_array{1,g}{c,i,area};
                                    deasy_k_all{g,area}{i}{c,e} =  dprime_easy_rec  ;
                                    dhard_k_all{g,area}{i}{c,e} = dprime_hard_rec  ;
                                    rec_idx (g,area,i, c,e) = i ;
                                    mouse_idx(g,area, i, c, e) = GNG_rec_all_cell{1, g}(i).Mouse ; 
                                    neuron_idx (g,area, i, c, e) = c ;

                                    
                                elseif length(idx_a_k) < min_n_trials
                                    A_k_all{g,area}{i}{c,e} = [];
                                    B_k_all{g,area}{i}{c,e} = [];
                                    idx_k_all{g,area}{i}{c,e} = [];
                                    T_k_all{g,area}{i}{c,e} = [];
                                    deasy_k_all{g,area}{i}{c,e} = [] ;
                                    dhard_k_all{g,area}{i}{c,e} = [] ;
                                end
                                
                            elseif length(idx_b_k) < min_n_trials
                                A_k_all{g,area}{i}{c,e} = [];
                                B_k_all{g,area}{i}{c,e}= [];
                                idx_k_all{g,area}{i}{c,e} = [];
                                T_k_all{g,area}{i}{c,e} = [];
                                deasy_k_all{g,area}{i}{c,e} = [] ;
                                dhard_k_all{g,area}{i}{c,e} = [] ;
                            end
                    elseif isempty(FR_array{1,g}{c,i,area})
                        A_k_all{g,area}{i}{c,e} = [];
                        B_k_all{g,area}{i}{c,e} = [];
                        idx_k_all{g,area}{i}{c,e} = [];
                        T_k_all{g,area}{i}{c,e} = [];
                        deasy_k_all{g,area}{i}{c,e} = [] ;
                        dhard_k_all{g,area}{i}{c,e} = [] ;
                    end
                end
            end
        end
    end
end

%% remove recs and areas at do not have cells with efficient trials or were not modulated
for g   = 1:numel(GNG_rec_all_cell) % run per group

    Recs = 1:numel(GNG_rec_all_cell{1,g}) ;
    for area = 1:length(areas)
        A_k_all{g,area} = A_k_all{g,area}(~cellfun('isempty', A_k_all{g,area})) ;
        B_k_all{g,area} = B_k_all{g,area}(~cellfun('isempty', B_k_all{g,area})) ;
        idx_k_all{g,area}{i}{c,e} = idx_k_all{g,area}(~cellfun('isempty', idx_k_all{g,area})) ;
        T_k_all{g,area}{i}{c,e} = T_k_all{g,area}(~cellfun('isempty', T_k_all{g,area})) ;
        deasy_k_all{g,area}{i}{c,e} = deasy_k_all{g,area}(~cellfun('isempty', deasy_k_all{g,area})) ;
        dhard_k_all{g,area}{i}{c,e} = dhard_k_all{g,area}(~cellfun('isempty', dhard_k_all{g,area})) ;

        
        
        A_k_all{g,area} = reshape( A_k_all{g,area},[],1) ;
        B_k_all{g,area} = reshape( B_k_all{g,area},[],1) ;
        idx_k_all{g,area} = reshape( idx_k_all{g,area},[],1) ;
        T_k_all{g,area} = reshape( T_k_all{g,area},[],1) ;

        deasy_k_all{g,area} = reshape( deasy_k_all{g,area},[],1) ;
        dhard_k_all{g,area} = reshape( dhard_k_all{g,area},[],1) ;
    end
end

%% concatenate all recordings
for g   = 1:numel(GNG_rec_all_cell) % run per group
    for area = 1:length(areas)
        A_k_fr{g,area} = cat(1,A_k_all{g,area}{:}) ;
        B_k_fr{g,area} = cat(1,B_k_all{g,area}{:}) ;

        idx_k_fr{g,area} = cat(1,idx_k_all{g,area}{:}) ;
        T_k_fr{g,area} = cat(1,T_k_all{g,area}{:}) ;

        deasy_k_fr{g,area} = cat(1,deasy_k_all{g,area}{:}) ;
        dhard_k_fr{g,area} = cat(1,dhard_k_all{g,area}{:}) ;

    end
end

%% concatenate all brain areas
for g   = 1:numel(GNG_rec_all_cell) % run per group
    A_k_g{1,g} = cat(1,A_k_fr{g:numel(GNG_rec_all_cell):numel(GNG_rec_all_cell)*numel(GNG_rec_all_cell)})  ;
    B_k_g{1,g} = cat(1,B_k_fr{g:numel(GNG_rec_all_cell):numel(GNG_rec_all_cell)*numel(GNG_rec_all_cell)})  ;

    idx_k_g{1,g} = cat(1,idx_k_fr{g:numel(GNG_rec_all_cell):numel(GNG_rec_all_cell)*numel(GNG_rec_all_cell)})  ;
    T_k_g{1,g} = cat(1,T_k_fr{g:numel(GNG_rec_all_cell):numel(GNG_rec_all_cell)*numel(GNG_rec_all_cell)})  ;

    deasy_k_g{1,g} = cat(1,deasy_k_fr{g:numel(GNG_rec_all_cell):numel(GNG_rec_all_cell)*numel(GNG_rec_all_cell)})  ;
    dhard_k_g{1,g} = cat(1,dhard_k_fr{g:numel(GNG_rec_all_cell):numel(GNG_rec_all_cell)*numel(GNG_rec_all_cell)})  ;

end

%% reassign all comparisons (e) per group and remove empty cells
for g   = 1:numel(GNG_rec_all_cell) % run per group
    for e = 1:length(a)
        for c = 1:size(A_k_g{1,g},1)
dhard_k_ge{e,g}{c} = [] ;
  deasy_k_ge{e,g}{c} = [] ;
            A_k_ge{e,g}{c} =   A_k_g{1,g}{c,e};
            B_k_ge{e,g}{c} =    B_k_g{1,g}{c,e};

            idx_k_ge{e,g}{c} =   idx_k_g{1,g}{c,e};
            T_k_ge{e,g}{c} =    T_k_g{1,g}{c,e};

            deasy_k_ge{e,g}{c} =    deasy_k_g{1,g}{c,e};
            dhard_k_ge{e,g}{c} =    dhard_k_g{1,g}{c,e};

        end

        A_k_ge{e,g} =  A_k_ge{e,g}(~cellfun('isempty', A_k_ge{e,g})) ;
        B_k_ge{e,g} =  B_k_ge{e,g}(~cellfun('isempty', B_k_ge{e,g})) ;

        idx_k_ge{e,g} =  idx_k_ge{e,g}(~cellfun('isempty', idx_k_ge{e,g})) ;
        T_k_ge{e,g} =  T_k_ge{e,g}(~cellfun('isempty', T_k_ge{e,g})) ;

        deasy_k_ge{e,g} =  deasy_k_ge{e,g}(~cellfun('isempty', deasy_k_ge{e,g})) ;
        dhard_k_ge{e,g} =  dhard_k_ge{e,g}(~cellfun('isempty', dhard_k_ge{e,g})) ;
    end
end

end 