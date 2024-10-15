function [behavior] = GNG_neuro_behavior (GNG_rec_all_cell) 

behavior.dprime = nan(numel(GNG_rec_all_cell),15) ;
behavior.cbias = nan(numel(GNG_rec_all_cell),15) ;
behavior.dprime_easy = nan(numel(GNG_rec_all_cell),15) ;
behavior.cbias_easy = nan(numel(GNG_rec_all_cell),15) ;
behavior.dprime_hard = nan(numel(GNG_rec_all_cell),15) ;
behavior.cbias_hard = nan(numel(GNG_rec_all_cell),15) ;

for g = 1:numel(GNG_rec_all_cell)
    for i = 1:numel (GNG_rec_all_cell{1,g})
        if ~isempty(GNG_rec_all_cell{1,g}(i))

            hit_easy  = size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_hit_easy,2) ;
            if isempty(hit_easy)
                hit_easy = 0 ;
            end
            FA_easy =  size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_FA_easy,2) ;
            if isempty(FA_easy)
                FA_easy = 0 ;
            end
            miss_easy =  size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_miss_easy,2) ;
            if isempty(miss_easy)
                miss_easy = 0 ;
            end
            CR_easy =  size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_CR_easy,2) ;
            if isempty(CR_easy)
                CR_easy = 0 ;
            end
            hit_hard  = size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_hit_hard,2) ;
            if isempty(hit_hard)
                hit_hard = 0 ;
            end
            FA_hard =  size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_FA_hard,2) ;
            if isempty(FA_hard)
                FA_hard = 0 ;
            end
            miss_hard =  size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_miss_hard,2) ;
            if isempty(miss_hard)
                miss_hard = 0 ;
            end
            CR_hard=  size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_CR_hard,2) ;
            if isempty(CR_hard)
                CR_hard = 0 ;
            end

           hit  = hit_easy + hit_hard ;
           FA = FA_easy + FA_hard  ;
           miss = miss_easy + miss_hard ;
           CR = CR_easy + CR_hard ;

           go =  size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_go,2) ;
           no_go =  size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_no_go,2) ;
           go_easy = size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_go_easy,2) ;
           go_hard =  size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_go_hard,2) ;
           no_go_easy = size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_no_go_easy,2) ;
           no_go_hard = size(GNG_rec_all_cell{1,g}(i).eventTimes.times_vec_no_go_hard,2) ;

           lick_rate(g,i) = (hit + FA) / (hit + FA + miss + CR)  ;
           if isnan(lick_rate(g,i)) | isinf(lick_rate(g,i))
               lick_rate(g,i) = 0 ;
           end

           go_lick_rate(g,i)  = hit/ (hit + miss)  ;
           if isnan(go_lick_rate(g,i) ) | isinf(go_lick_rate(g,i) )
               go_lick_rate(g,i)  = 0 ;
           end

           no_go_lick_rate(g,i) = (FA) / (FA + CR)  ;
           if isnan(no_go_lick_rate(g,i) ) | isinf(no_go_lick_rate(g,i) )
               no_go_lick_rate(g,i)  = 0 ;
           end

           easy_lick_rate(g,i) = (hit_easy + FA_easy) / ...
               (hit_easy + FA_easy + miss_easy + CR_easy)   ;
           if isnan(easy_lick_rate(g,i) ) | isinf(easy_lick_rate(g,i) )
               easy_lick_rate(g,i)  = 0 ;
           end

           hard_lick_rate(g,i) = (hit_hard + FA_hard) / ...
               (hit_hard + FA_hard + miss_hard + CR_hard)   ;
           if isnan(hard_lick_rate(g,i) ) | isinf(hard_lick_rate(g,i) )
               hard_lick_rate(g,i)  = 0 ;
           end


           easy_go_lick_rate(g,i) = hit_easy / (hit_easy + miss_easy);
           if isnan(easy_go_lick_rate(g,i) ) | isinf(easy_go_lick_rate(g,i) )
               easy_go_lick_rate(g,i)  = 0 ;
           end

           hard_go_lick_rate(g,i) = hit_hard / (hit_hard + miss_hard);
           if isnan(hard_go_lick_rate(g,i) ) | isinf(hard_go_lick_rate(g,i) )
               hard_go_lick_rate(g,i)  = 0 ;
           end

           easy_nogo_lick_rate(g,i) = FA_easy / (FA_easy + CR_easy);
           if isnan(easy_nogo_lick_rate(g,i) ) | isinf(easy_nogo_lick_rate(g,i) )
               easy_nogo_lick_rate(g,i)  = 0 ;
           end

           hard_nogo_lick_rate(g,i) = FA_hard / (FA_hard + CR_hard);
           if isnan(hard_nogo_lick_rate(g,i) ) | isinf(hard_nogo_lick_rate(g,i) )
               hard_nogo_lick_rate(g,i)  = 0 ;
           end


           zHit = norminv(go_lick_rate(g,i)) ;

           if isnan(zHit ) | isinf(zHit)
               zHit = 0 ;
           end
           zFA = norminv(no_go_lick_rate(g,i)) ;

           if isnan(zFA) | isinf(zFA )
               zFA = 0 ;
           end
           behavior.dprime(g,i) = zHit - zFA ;
           behavior.cbias(g,i) = (zHit + zFA)./2;

           zHit_easy = norminv( easy_go_lick_rate(g,i)) ;
             if isnan(zHit_easy ) | isinf(zHit_easy)
               zHit_easy = 0 ;
             end

           zFA_easy = norminv(easy_nogo_lick_rate(g,i)) ;
           if isnan(zFA_easy ) | isinf(zFA_easy)
               zFA_easy = 0 ;
           end
           behavior.dprime_easy(g,i) = zHit_easy - zFA_easy ;
           behavior.cbias_easy(g,i) = (zHit_easy + zFA_easy)./2;

           zHit_hard = norminv(hard_go_lick_rate(g,i)) ;
             if isnan(zHit_hard ) | isinf(zHit_hard)
               zHit_hard = 0 ;
           end
           zFA_hard = norminv(hard_nogo_lick_rate(g,i)) ;
            if isnan(zFA_hard ) | isinf(zFA_hard)
               zFA_hard = 0 ;
           end
           behavior.dprime_hard(g,i) = zHit_hard - zFA_hard ;
           behavior.cbias_hard(g,i) = (zHit_hard + zFA_hard)./2;


           if isnan(behavior.dprime(g,i) ) | isinf(behavior.dprime(g,i) )
               behavior.dprime(g,i) = 0 ;
           end

           if isnan(behavior.cbias(g,i) ) | isinf(behavior.cbias(g,i) )
               behavior.cbias(g,i) = 0 ;
           end

           if isnan( behavior.dprime_hard(g,i)  ) | isinf( behavior.dprime_hard(g,i)  )
               behavior.dprime_hard(g,i) = 0 ;
           end

           if isnan(behavior.cbias_hard(g,i) ) | isinf(behavior.cbias_hard(g,i) )
               behavior.cbias_hard(g,i) = 0 ;
           end

           if isnan(behavior.dprime_easy(g,i) ) | isinf(behavior.dprime_easy(g,i) )
               behavior.dprime_easy(g,i) = 0 ;
           end

           if isnan(behavior.cbias_easy(g,i) ) | isinf(behavior.cbias_easy(g,i) )
               behavior.cbias_easy(g,i) = 0 ;
           end
        end
    end
end
end