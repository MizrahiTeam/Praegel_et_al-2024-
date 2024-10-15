
function [latency_peak_AUC, onset_latency_AUC, width_AUC, max_AUC, rec_idx_th, mouse_idx_th]...
    = GNG_max_onset_width_AUC (GNG_rec_all_cell, ind_AUC_all, ind_AUC_shuf_all, std_ind_AUC_shuf_all, a, stim_onset_bin,it_size,rec_idx, mouse_idx)

for e = 1:length(a) % stimulus code


    for   g = 1:numel(GNG_rec_all_cell) % run per group
        rec_idx_ge = [] ;
        mouse_idx_ge = [] ;


        rec_idx_ge = reshape(rec_idx(g,:,:,:,a(e)),[],1)  ;
        rec_idx_ge(rec_idx_ge == 0) = nan ; 
        rec_idx_ge = rec_idx_ge(~isnan(rec_idx_ge)) ;

        mouse_idx_ge = reshape(mouse_idx(g,:,:,:,a(e)),[],1) ;
        mouse_idx_ge(mouse_idx_ge == 0) = nan ; 
        mouse_idx_ge = mouse_idx_ge(~isnan(mouse_idx_ge)) ;


        latency_peak_AUC{a(e),g} = nan(700,1) ;
        onset_latency_AUC{a(e),g} = nan(700,1) ;
        width_AUC {a(e),g} = nan(700,1) ;
        max_AUC{a(e),g} = nan(700,1) ;
        rec_idx_th{a(e),g} = nan(700,1) ;
        mouse_idx_th{a(e),g} = nan(700,1) ;

        for c =1:size(ind_AUC_all{a(e),g},1)
            idx_onset_latency = find(ind_AUC_all{a(e),g}(c,stim_onset_bin:end)...
                > (ind_AUC_shuf_all{a(e),g}(c,stim_onset_bin:end) + std_ind_AUC_shuf_all{a(e),g}(c,stim_onset_bin:end))) ;

            if ~isempty(idx_onset_latency)

                onset_latency_AUC{a(e),g}(c,:) = (idx_onset_latency(1) * it_size);

                width_AUC{a(e),g}(c,:) = length(idx_onset_latency) * it_size ;

                AUC_peak =  max(ind_AUC_all{a(e),g}(c,stim_onset_bin:end)) ;
                latency_peak_cell =  find(ind_AUC_all{a(e),g}(c,stim_onset_bin:end) == AUC_peak(1)) ;

                max_AUC{a(e),g}(c,:) = max(AUC_peak) ;
                latency_peak_AUC{a(e),g}(c,:) = latency_peak_cell(end)* it_size ;

                rec_idx_th{a(e),g}(c,:) = rec_idx_ge(c,:) ;
                mouse_idx_th{a(e),g}(c,:) = mouse_idx_ge(c,:) ;
            end
        end


    end
end
end