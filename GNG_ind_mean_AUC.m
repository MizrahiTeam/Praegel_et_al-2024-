
function [ind_AUC, ind_AUC_shuf, ind_AUC_all,ind_AUC_shuf_all,std_ind_AUC_shuf_all, mean_vec, stde_vec,mean_vec_shuf, std_vec_shuf ]...
    = GNG_ind_mean_AUC (GNG_rec_all_cell, r_AUC_abs, r_AUC_shuf_abs, a, trial_samples, run_window, stim_onset_bin, window, it_size, per_mean)




if per_mean  == 0

    for   g = 1:numel(GNG_rec_all_cell) % run per group
        for e = a
            for c = 1:size(r_AUC_abs{g,e},2)
                for k = 1:trial_samples
                    % retrieve the abosolute AUC
                    ind_AUC{e,g}(k,c,:)  = r_AUC_abs{g,e} (k,c,:) ;
                    ind_AUC_shuf{e,g}(k,c,:)   = r_AUC_shuf_abs{g,e}(k,c,:)  ;

                    % retrieve mean baseline and mean shuffle (since baseline is because of shufflled bins)
                    AUC_base = mean(ind_AUC{e,g}(k,c,1:stim_onset_bin)) ;
                    AUC_base_shuf = mean(ind_AUC_shuf{e,g}(k,c,:)) ;

                    % adjust the shift in baseline
                    ind_AUC{e,g}(k,c,:) = ((ind_AUC{e,g}(k,c,:) - AUC_base) + 0.5) ;
                    ind_AUC_shuf{e,g}(k,c,:) =  ((ind_AUC_shuf{e,g}(k,c,:) -   AUC_base_shuf) +  0.5) ;
                end

                % mean for all trial samples
                ind_AUC_all{e,g}(c,:) = squeeze(mean (ind_AUC{e,g}(:,c,:))) ;
                ind_AUC_shuf_all{e,g}(c,:) =  squeeze(mean (ind_AUC_shuf{e,g}(:,c,:)))' ;
                std_ind_AUC_shuf_all{e,g}(c,:) =  squeeze(std (ind_AUC_shuf{a(e),g}(:,c,:)))'*3 ;

            end

            for bin = it_size:it_size:((window(2)*1000) - run_window) % run per bin per window
                bin_pos = bin / it_size ;

                mean_vec{e,g}(1,bin_pos)   = nanmean(ind_AUC_all{e,g}(:,bin_pos)) ;
                stde_vec{e,g} (1,bin_pos)  = nanstd(ind_AUC_all{e,g}(:,bin_pos))./sqrt(size(ind_AUC_all{e,g},1)) ;

                mean_vec_shuf{e,g}(1,bin_pos)   = nanmean(ind_AUC_shuf_all{e,g}(:,bin_pos)) ;
                std_vec_shuf{e,g} (1,bin_pos)  = (nanstd(ind_AUC_shuf_all{e,g}(:,bin_pos))*3) ;

            end
        end
    end


elseif per_mean == 1


    size_g_e = nan (numel(GNG_rec_all_cell), max(a)) ;
    for  g = 1:numel(GNG_rec_all_cell)
        for e = a
            size_g_e(g,e) = size (r_AUC_abs{g,e},2) ;
        end
    end

    for   g = 1:numel(GNG_rec_all_cell) % run per group
        for e = 1:length(a)
            for c = 1:(min(size_g_e(g,:))) % change to 4 for naive
                for k = 1:trial_samples
                    % retrieve the abosolute AUC

                    ind_AUC_eh_mean = [] ;
                    ind_AUC_shuf_eh_mean = [] ;

                    for bin = 1: size(r_AUC_abs{g,a(e)}(k,c,:),3)

                        ind_AUC_eh = [] ;
                        ind_AUC_shuf_eh = [] ;

                        for ch = 1:length(a)

                            ind_AUC_eh(ch) =r_AUC_abs{g,a(ch)} (k,c,bin) ;
                            ind_AUC_shuf_eh(ch) =r_AUC_shuf_abs{g,a(ch)}(k,c,bin)  ;
                        end
                        ind_AUC_eh_mean(bin) = mean(ind_AUC_eh) ;
                        ind_AUC_shuf_eh_mean(bin) = mean(ind_AUC_shuf_eh) ;

                    end

                    ind_AUC{a(e),g}(k,c,:) =  ind_AUC_eh_mean ;
                    ind_AUC_shuf{a(e),g}(k,c,:)  = ind_AUC_shuf_eh_mean  ;

                    %retrieve mean baseline and mean shuffle (since baseline is because of shufflled bins)
                    AUC_base = mean(ind_AUC{a(e),g}(k,c,1:stim_onset_bin)) ;
                    AUC_base_shuf = mean(ind_AUC_shuf{a(e),g}(k,c,:)) ;

                    % adjust the shift in baseline
                    ind_AUC{a(e),g}(k,c,:) = ((ind_AUC{a(e),g}(k,c,:) - AUC_base) + 0.5) ;
                    ind_AUC_shuf{a(e),g}(k,c,:) =  ((ind_AUC_shuf{a(e),g}(k,c,:) -   AUC_base_shuf) +  0.5) ;

                end

                % mean for all trial samples
                ind_AUC_all{a(e),g}(c,:) = squeeze(mean (ind_AUC{a(e),g}(:,c,:))) ;
                ind_AUC_shuf_all{a(e),g}(c,:) =  squeeze(mean (ind_AUC_shuf{a(e),g}(:,c,:)))' ;
                std_ind_AUC_shuf_all{a(e),g}(c,:) =  squeeze(std (ind_AUC_shuf{a(e),g}(:,c,:)))'*3 ;

            end

            for bin = it_size:it_size:(800 - run_window) % run per bin per window
                bin_pos = bin / it_size ;

                mean_vec{a(e),g}(1,bin_pos)   = nanmean(ind_AUC_all{a(e),g}(:,bin_pos)) ;
                stde_vec{a(e),g} (1,bin_pos)  = nanstd(ind_AUC_all{a(e),g}(:,bin_pos))./sqrt(size(ind_AUC_all{a(e),g},1)) ;

                mean_vec_shuf{a(e),g}(1,bin_pos)   = nanmean(ind_AUC_shuf_all{a(e),g}(:,bin_pos)) ;
                std_vec_shuf{a(e),g} (1,bin_pos)  = (nanstd(ind_AUC_shuf_all{a(e),g}(:,bin_pos))*3) ;

            end
        end
    end

end
end