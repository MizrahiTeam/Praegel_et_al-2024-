clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shuffeled datasamples will lead to different distributions than in paper
%figures. Yet these do not affect the significance and effect sizes of the
%results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load GNG_rec_all_cell & Fr_array
addpath(genpath('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_MATLABR2023b_scripts'))

% select recording sessions
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , 'Select GNG_rec_all_cell ');
addpath(path)
load (file)

% select recording sessions
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , 'Select Fr_array ');
addpath(path)
load (file)
cd (path)

%% areas of recording
area_str ={'AUDd','AUDp','AUDv','TEa'};
areas = 1:length(area_str) ;
group_str = {'adolescent','adult'};

% trial types for comparison
% stimulus + choice
% hit vs. fa easy ; hit vs fa hard; cr vs. fa easy; cr vs. fa hard ;
a = [1 7 2 8] ;
b = [2 8 4 10] ;

% AUC parameters
run_window = 50 ; % ms
it_size = 25 ; % ms
n_shuffle = 10 ; % n shuffles for AUC
trial_samples = 10 ;
min_n_trials = 15 ;

%time parameters
window = [-0.2; 0.6]; % time window according to stimulus onset
window_length_ms = dist(window(1),window(2))*1000 ;
n_bins = round( (( (length (1:window_length_ms)) - run_window )/it_size),0) ;

start_base = 1 ;
startStim = 0; % stimulus onset ms
stopStim = 0.1; % stimulus offset  ms
stim_start_ms = dist(window(1), startStim)*1000 ;
stim_stop_ms = dist(window(1), stopStim)*1000 ;
stim_onset_bin = 7 ;


% colors
color_eh= {[0.3010 0.7450 0.9330], [0 0.4470 0.7410], [0 0 1]};
color_eh_fade = {[0.3010 0.7450 0.9330 0.3], [0 0.4470 0.7410 0.3], [0 0 1 0.3]};
color_eh_patch = {[0 0.4470 0.7410],[0.3010 0.7450 0.9330], [0 0.4470 0.7410],[0.3010 0.7450 0.9330]};

L = {['--'], ['-'],['--'], ['-']};
M = {['o'],['o'],['o'], ['o']};

%% only consider experts
GNG_rec_all_cell_exp{1,1} = GNG_rec_all_cell{1,1}; 
GNG_rec_all_cell_exp{1,2} = GNG_rec_all_cell{1,2}; 

%% extract behavior per recording

clc
[behavior] = GNG_neuro_behavior (GNG_rec_all_cell_exp) ;

%% concatenate all cells per recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shuffeled datasamples will lead to different distributions than in paper
%figures. Yet these do not affect the significance and effect sizes of the
%results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
tic
[A_k_ge, B_k_ge, idx_k_ge, T_k_ge, deasy_k_ge, dhard_k_ge, rec_idx, neuron_idx, mouse_idx]...
    = GNG_t_per_n_AUC (FR_array, GNG_rec_all_cell_exp, behavior, areas, a, b,  trial_samples, min_n_trials) ;
toc

%% run the AUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shuffeled datasamples will lead to different distributions than in paper
%figures. Yet these do not affect the significance and effect sizes of the
%results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[r_AUC, r_AUC_abs, r_AUC_shuf, r_AUC_shuf_abs]...
    = GNG_running_AUC (GNG_rec_all_cell_exp,a,b, A_k_ge, B_k_ge, trial_samples, it_size, run_window, window_length_ms, n_shuffle)
toc
GNG_analysis.r_AUC = r_AUC ;
GNG_analysis.r_AUC_abs = r_AUC_abs ;
GNG_analysis.r_AUC_shuf = r_AUC_shuf ;
GNG_analysis.r_AUC_shuf_abs =  r_AUC_shuf_abs ;


%% mean of stimulus activity  of all trial samples
clc
per_mean = false ;
a = [1 2] ;
stim_onset_bin = abs((window(1)* 1000)/ run_window) ;
[ind_AUC, ind_AUC_shuf, ind_AUC_all,ind_AUC_shuf_all,std_ind_AUC_shuf_all, mean_vec, stde_vec,mean_vec_shuf, std_vec_shuf ]...
    = GNG_ind_mean_AUC (GNG_rec_all_cell_exp, GNG_analysis.r_AUC_abs, GNG_analysis.r_AUC_shuf_abs,...
    a, trial_samples, run_window, stim_onset_bin, window, it_size, per_mean) ;

GNG_analysis.ind_AUC = ind_AUC ;
GNG_analysis.ind_AUC_shuf  = ind_AUC_shuf ;
GNG_analysis.ind_AUC_all  =ind_AUC_all ;
GNG_analysis.ind_AUC_shuf_all  = ind_AUC_shuf_all ;
GNG_analysis.std_ind_AUC_shuf_all  =std_ind_AUC_shuf_all ;
GNG_analysis.mean_vec =mean_vec ;
GNG_analysis.stde_vec  = stde_vec ;
GNG_analysis.mean_vec_shuf  =mean_vec_shuf ;
GNG_analysis.std_vec_shuf  = std_vec_shuf;

%% mean choice activity of all trial samples
clc
% concatenate e = 3 and e = 4 (fa vs. CR in easy and hard condition
a = [3 4] ; % fa vs. cr mean of easy and hard
per_mean = true; % calculate the mean choice for easy and hard choice together 

[ind_AUC, ind_AUC_shuf, ind_AUC_all,ind_AUC_shuf_all,std_ind_AUC_shuf_all, mean_vec, stde_vec,mean_vec_shuf, std_vec_shuf ]...
    = GNG_ind_mean_AUC (GNG_rec_all_cell_exp, GNG_analysis.r_AUC_abs, GNG_analysis.r_AUC_shuf_abs,...
    a, trial_samples, run_window, stim_onset_bin, window, it_size, per_mean) ;

for e = a(1)
    for g = 1:numel(GNG_rec_all_cell_exp)
        GNG_analysis.ind_AUC{e,g} = ind_AUC{e,g} ;
        GNG_analysis.ind_AUC_shuf{e,g}  = ind_AUC_shuf{e,g} ;
        GNG_analysis.ind_AUC_all{e,g}  =ind_AUC_all{e,g} ;
        GNG_analysis.ind_AUC_shuf_all{e,g}  = ind_AUC_shuf_all{e,g} ;
        GNG_analysis.std_ind_AUC_shuf_all{e,g}  =std_ind_AUC_shuf_all{e,g} ;
        GNG_analysis.mean_vec{e,g}  =mean_vec{e,g} ;
        GNG_analysis.stde_vec{e,g}   = stde_vec{e,g} ;
        GNG_analysis.mean_vec_shuf{e,g}   =mean_vec_shuf{e,g} ;
        GNG_analysis.std_vec_shuf {e,g}  = std_vec_shuf{e,g};
    end
end

%% Max AUC & AUC Latency Caclulation
clc

a = [1 2 3] ; % stimulus easy ^ hard , choice mean 
[latency_peak_AUC, onset_latency_AUC, width_AUC, max_AUC, rec_idx_th, mouse_idx_th]...
    = GNG_max_onset_width_AUC (GNG_rec_all_cell_exp, GNG_analysis.ind_AUC_all, GNG_analysis.ind_AUC_shuf_all,...
    GNG_analysis.std_ind_AUC_shuf_all, a, stim_onset_bin,it_size,rec_idx, mouse_idx) ;

for e = a
    for g = 1:numel(GNG_rec_all_cell_exp)
        GNG_analysis.max_AUC{e,g} = max_AUC{e,g} ;
        GNG_analysis.width_AUC{e,g}  = width_AUC{e,g} ;
        GNG_analysis.onset_latency_AUC{e,g}  =onset_latency_AUC{e,g} ;
    end
end

 %% Figure 4 E & F  plot a single neuron example 
% plot hit vs fa easy and hard and choice fa vs cr
a = [ 1 2 3] ; 
ac = [1 2 3] ; % color code
close all
% plot mean (optional include shuffled mean)
    for   g = 2 % expert adult example 
        
        for c = 271 %example neuron 
             figure(c)
         
            
            for e = 1:length(a)
                subplot(1,3,e)
                plot( GNG_analysis.ind_AUC_all{e,g}(c,:),'Color',color_eh{ac(e)},...
                    'linestyle',L{g},'linewidth',4);
                hold on

                plot( GNG_analysis.ind_AUC_shuf_all{e,g}(c,:),'Color',[0.5 0.5 0.5],...
                    'linestyle',L{g},'linewidth',4)
                hold on
                xline(6,'--','Color','k');
                xticks([1 6 18 ])
                xticklabels([ -200 0 400])
                ylim([0.45 1])
                yticks([ 0.5 0.75 1])
                xlabel('time (ms)')
                ylabel('Discrimination (abs. AUC)')
                ax = gca;
                ax.XAxis.FontSize = 20;
                ax.YAxis.FontSize = 20;
                movegui('east');
                box off;

            end
     
        hit_idx = idx_k_ge{2,g}{c}(13,:)  ;
        fa_idx =  idx_k_ge{2,g}{c} (14,:) ;
        miss_idx =  idx_k_ge{2,g}{c} (15,:) ;
        cr_idx  =  idx_k_ge{2,g}{c} (16,:) ;

          hit_idx(hit_idx == 0) =nan ;
         hit_idx =  hit_idx(~isnan(hit_idx)) ;

          fa_idx(fa_idx == 0) =nan ;
         fa_idx =  fa_idx(~isnan(fa_idx)) ;

          miss_idx(miss_idx == 0) =nan ;
         miss_idx =  miss_idx(~isnan(miss_idx)) ;

          cr_idx(cr_idx == 0) =nan ;
         cr_idx =  cr_idx(~isnan(cr_idx)) ;

         color_mean = {[0.4660 0.6740 0.1880],[0.9290 0.6940 0.1250], [0.8500 0.3250 0.0980], [0.4940 0.1840 0.5560]};

         size_tr = [length(hit_idx),length(miss_idx), length(fa_idx),length(cr_idx)] ;
         max_size = max(size_tr) ;
        trial_types_raster = cat(1,T_k_ge{2,g}{c}(hit_idx,:),T_k_ge{2,g}{c}(miss_idx,:), T_k_ge{2,g}{c}(fa_idx,:),T_k_ge{2,g}{c}(cr_idx,:)) ;

        figure(c+1) 
     imagesc(mat2gray(trial_types_raster))

     hold on
     rectangle( 'Position' , [0 0 50 size_tr(1,1)+1],'Facecolor',color_mean{1},'EdgeColor','none')
     hold on
     yline(size_tr(1,1)+1,'-w','linewidth',3)
     hold on

     rectangle( 'Position' , [0 (size_tr(1,1)+1) 50 size_tr(1,2)+1],'Facecolor',color_mean{3},'EdgeColor','none')
     hold on
     yline((size_tr(1,1)+size_tr(1,2)),'-w','linewidth',3)
     hold on

     rectangle( 'Position' , [0 (size_tr(1,1)+size_tr(1,2)+1) 50 size_tr(1,3)],'Facecolor',color_mean{2},'EdgeColor','none')
     hold on
     yline((size_tr(1,1)+size_tr(1,2)+size_tr(1,3)+1),'-w','linewidth',3)
     hold on

     rectangle( 'Position' , [0 (size_tr(1,1)+size_tr(1,2)+ size_tr(1,3)+1) 50 size_tr(1,4)],'Facecolor',color_mean{4},'EdgeColor','none')
     hold on

     xticks([0:200:800])
     xticklabels([-200:200:600])
     ylim([0.5 sum(size_tr)])
     yticks([1 sum(size_tr)/2 sum(size_tr) ])
     xlabel('time(ms)')
     ylabel('trial number')

     c = colorbar 
c.FontSize = 20 
ax = gca;
ax.XAxis.FontSize =20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;
ylabel(c,'normalized FR (Hz)')

        end
    end
    
    
    %% Figure 4 G / plot stimulus related & choice related AUC mean 
    clc
%close all 
a = [1 2 3] ;
ac = a ;  % color code
% plot mean (optional include shuffled mean)
for e = 1:length(a) 
    figure(1)
    for   g = 1:numel(GNG_rec_all_cell_exp) % run per group
subplot(1,3,e)

         length_bins = size(  GNG_analysis.mean_vec{a(e),g},2) ;
        plot([1, 1:length_bins-1],  GNG_analysis.mean_vec{a(e),g},'Color',color_eh{ac(e)},...
            'linestyle',L{g},'linewidth',2);
        hold on
        
        patch([[1,1:length_bins-1] flip([1,1:length_bins-1])] , [  GNG_analysis.mean_vec{a(e),g} + ....
              GNG_analysis.stde_vec{a(e),g} flip(  GNG_analysis.mean_vec{a(e),g} -  GNG_analysis.stde_vec{a(e),g})], color_eh{ac(e)} ,...
            'facealpha' , 0.2, 'EdgeColor','none')
        hold on
        
        plot([1, 1:length_bins-1],   GNG_analysis.mean_vec_shuf{a(e),g}+ ....
             GNG_analysis.std_vec_shuf{a(e),g},'Color',[0.5 0.5 0.5],...
            'linestyle',L{g},'linewidth',2);
        hold on
        

        xline(6,'--','Color','k');
        xticks([1 6  18])
        xticklabels([-200 0  400])
        xlim([1 (length_bins)])
        ylim([0.5 0.6])
        xlabel('time (ms)')
        ylabel(' Discrimination (abs. AUC)')
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');
        box off;
    end
end
%% Figure 4H violinplot AUC onset latency
close all
clc
clear mouse rec var

var = GNG_analysis.onset_latency_AUC ;
mouse = mouse_idx_th ;
rec = rec_idx_th ;

for e = 1:3
    for g = 1:numel(GNG_rec_all_cell_exp)
        var{e,g} = var{e,g}(~isnan(var{e,g})) ;
        mouse{e,g} = mouse{e,g}(~isnan(mouse{e,g})) ;
        rec{e,g} = rec{e,g}(~isnan(rec{e,g})) ;
        mean_onset(e,g) = mean (var{e,g}) ; 
        ste_onset(e,g) = std (var{e,g}) / sqrt (size(var{e,g},1)) ;
    end
end

figure
for e = a
    subplot(1,3,e)
    violinplot([GNG_analysis.onset_latency_AUC{e,1},  GNG_analysis.onset_latency_AUC{e,2}]...
        ,{' adolescenct','adult'},"ViolinColor",...
        {[  color_eh{e};  color_eh{e}]} ) ;
    ylim ([0 1000])
    ylabel('time(ms)')

    box off;
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;

    LME_AUC_onset = table;
    LME_AUC_onset.onset = [  var{e,1}'  , var{e,2}']' ;
    LME_AUC_onset.mouse = [  mouse{e,1}'  , mouse{e,2}']' ;
    LME_AUC_onset.rec = [ rec{e,1}'  , rec{e,2}']' ;
    LME_AUC_onset.age = [ ones(size(var{e,1})) ;  ones(size(var{e,2}))+ 1] ;

    % Convert appropriate columns to categorical variables
    LME_AUC_onset.mouse = categorical(LME_AUC_onset.mouse);
    LME_AUC_onset.rec = categorical(LME_AUC_onset.rec);
    LME_AUC_onset.age = categorical(LME_AUC_onset.age);

    formula = 'onset ~ age  + (1|mouse) + (1|rec)';
    lme_onset{e} = fitlme(LME_AUC_onset, formula);
    % Extract the fixed effects table
    fixedEffectsTable = lme_onset{e}.Coefficients;

    resultstable = round([lme_onset{e}.Coefficients.Estimate,lme_onset{e}.Coefficients.SE,lme_onset{e}.Coefficients.tStat,lme_onset{e}.Coefficients.pValue],3) ; 

writematrix(resultstable,'LME_fig4_onset.xlsx','Sheet',e)

    % Extract p-values from the fixed effects table
    p = fixedEffectsTable.pValue(2)
   
        if  p(1,1) < 0.0005
            txt{1,1} = '***';
        elseif  p(1,1) < 0.005
            txt{1,1} = '**';
        elseif  p(1,1) < 0.05
            txt{1,1} = '*';
        elseif p(s,1) > 0.05
            txt{1,1} = 'n.s.';
        end
        
        line( linspace(1.1, 1.9), linspace(700 ,700),'color','k')
        hold on
        text(mean([1.5 1.5]), 750 ,txt{1,1},'FontSize',20)
        hold on
end
% %% LME for onset delay of discrimination 
% % Prepare the data
% LME_AUC_width = table; 
% LME_AUC_width.onset = [var{1,1}' , var{1,2}' , var{2,1}' , var{2,2}' , var{3,1}', var{3,2}']';
% LME_AUC_width.mouse = [mouse{1,1}' , mouse{1,2}' , mouse{2,1}' , mouse{2,2}' , mouse{3,1}', mouse{3,2}']';
% LME_AUC_width.rec = [rec{1,1}' , rec{1,2}' , rec{2,1}' , rec{2,2}' , rec{3,1}', rec{3,2}']';
% LME_AUC_width.age = [ones(size(var{1,1})) ; ones(size(var{1,2}))+ 1 ; ones(size(var{2,1})) ; ones(size(var{2,2}))+ 1 ; ones(size(var{3,1})) ; ones(size(var{3,2}))+ 1];
% LME_AUC_width.condition = [ones(size(var{1,1})) ; ones(size(var{1,2})) ; ones(size(var{2,1}))+1 ; ones(size(var{2,2}))+ 1 ; ones(size(var{3,1}))+ 2; ones(size(var{3,2}))+ 2];
% 
% % Convert appropriate columns to categorical variables
% LME_AUC_width.mouse = categorical(LME_AUC_width.mouse);
% LME_AUC_width.rec = categorical(LME_AUC_width.rec);
% LME_AUC_width.age = categorical(LME_AUC_width.age);
% LME_AUC_width.condition = categorical(LME_AUC_width.condition);
% 
% % Define the formula and fit the LME model
% formula = 'onset ~ age * condition + (1|mouse) + (1|rec)';
% lme = fitlme(LME_AUC_width, formula);
% 
% % Extract coefficients and covariance matrix from the model
% coefficients = lme.Coefficients;
% covarianceMatrix = lme.CoefficientCovariance;
% 
% % Find indices for coefficients
% interceptIdx = strcmp(coefficients.Name, '(Intercept)');
% condition2Idx = strcmp(coefficients.Name, 'condition_2');
% condition3Idx = strcmp(coefficients.Name, 'condition_3');
% age2Idx = strcmp(coefficients.Name, 'age_2');
% age2_condition2Idx = strcmp(coefficients.Name, 'age_2:condition_2');
% age2_condition3Idx = strcmp(coefficients.Name, 'age_2:condition_3');
% 
% 
% % Calculate estimate1: difference between condition3 and condition2 for age group 1
%     b_condition2 = coefficients.Estimate(condition2Idx);
%     b_condition3 = coefficients.Estimate(condition3Idx);
%     estimate1 = b_condition3 - b_condition2;
% 
%     % Variance for estimate1
%     var_estimate1 = covarianceMatrix(condition3Idx, condition3Idx) + ...
%                     covarianceMatrix(condition2Idx, condition2Idx) - ...
%                     2 * covarianceMatrix(condition2Idx, condition3Idx);
%     ste1 = sqrt(var_estimate1);
%     tStat1 = estimate1 / ste1;
%     pValue1 = 2 * (1 - tcdf(abs(tStat1), lme.DFE)); % Two-tailed test
% 
% 
%     b_age2 = coefficients.Estimate(age2Idx);
%     b_age2_condition2 = coefficients.Estimate(age2_condition2Idx);
%     b_age2_condition3 = coefficients.Estimate(age2_condition3Idx);
%     estimate2 = (b_age2_condition3 - b_age2) - (b_age2_condition2 - b_age2);
% 
%     % Variance for estimate2
%     var_age2 = covarianceMatrix(age2Idx, age2Idx);
%     var_age2_condition2 = covarianceMatrix(age2_condition2Idx, age2_condition2Idx);
%     var_age2_condition3 = covarianceMatrix(age2_condition3Idx, age2_condition3Idx);
%     cov_age2_age2_condition2 = covarianceMatrix(age2Idx, age2_condition2Idx);
%     cov_age2_age2_condition3 = covarianceMatrix(age2Idx, age2_condition3Idx);
%     cov_age2_condition2_age2_condition3 = covarianceMatrix(age2_condition2Idx, age2_condition3Idx);
% 
%     var_estimate2 = var_age2_condition3 + var_age2 - 2 * cov_age2_age2_condition3 + ...
%                     var_age2_condition2 + var_age2 - 2 * cov_age2_age2_condition2 + ...
%                     2 * cov_age2_condition2_age2_condition3;
%     ste2 = sqrt(var_estimate2);
%     tStat2 = estimate2 / ste2;
%     pValue2 = 2 * (1 - tcdf(abs(tStat2), lme.DFE)); % Two-tailed test
% 
% 
% resultstable = [lme.Coefficients.Estimate,lme.Coefficients.SE,lme.Coefficients.tStat,lme.Coefficients.pValue] ; 
% 
% 
% 
% % Create a new results table including the additional estimates and their statistics
% fixedeffects = table(round([resultstable(1,:);resultstable(2,:); resultstable(3,:); resultstable(4,:); [estimate1 ste1 tStat1 pValue1];...
% resultstable(5,:); resultstable(6,:);[estimate2 ste2 tStat2 pValue2]] ,3), ...
% 'RowNames',{'Intercept';'age';  'stim. easy vs. stim. hard';'stim. easy vs. choice';'stim hard vs. choice'; ...
%     'interaction stim. easy vs. stim. hard';' interaction stim. easy vs. choice';'interaction stim hard vs. choice'});
% 
% % Append to the existing results table
% disp(fixedeffects)
% writetable(fixedeffects,'LME_AUC_onset.xlsx','Sheet',1)
%% Figure 4I violinplot AUC width 
clc
clear mouse rec var
var = GNG_analysis.width_AUC ;
mouse = mouse_idx_th ;
rec = rec_idx_th ;

for e = 1:3
    for g = 1:numel(GNG_rec_all_cell_exp)
        var{e,g} = var{e,g}(~isnan(var{e,g})) ;
        mouse{e,g} = mouse{e,g}(~isnan(mouse{e,g})) ;
        rec{e,g} = rec{e,g}(~isnan(rec{e,g})) ;
         mean_width(e,g) = mean (var{e,g}) ; 
        ste_width(e,g) = std (var{e,g}) / sqrt (size(var{e,g},1)) ;
    end
end

figure
for e = a
    subplot(1,3,e)
    violinplot([GNG_analysis.width_AUC{e,1},  GNG_analysis.width_AUC{e,2}]...
        ,{' adolescenct','adult'},"ViolinColor",...
        {[  color_eh{e};  color_eh{e}]} ) ; 
    ylim ([0 1000])
    ylabel('time(ms)')
    
    box off;
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;


    LME_AUC_width = table; 
LME_AUC_width.width = [  var{e,1}'  , var{e,2}']' ;
LME_AUC_width.mouse = [  mouse{e,1}'  , mouse{e,2}']' ;
LME_AUC_width.rec = [ rec{e,1}'  , rec{e,2}']' ;
LME_AUC_width.age = [ ones(size(var{e,1})) ;  ones(size(var{e,2}))+ 1] ;

% Convert appropriate columns to categorical variables
LME_AUC_width.mouse = categorical(LME_AUC_width.mouse);
LME_AUC_width.rec = categorical(LME_AUC_width.rec);
LME_AUC_width.age = categorical(LME_AUC_width.age);

formula = 'width ~ age  + (1|mouse) + (1|rec)';
lme_width{e} = fitlme(LME_AUC_width, formula);
% Extract the fixed effects table
fixedEffectsTable = lme_width{e}.Coefficients;



resultstable = round([lme_width{e}.Coefficients.Estimate,lme_width{e}.Coefficients.SE,lme_width{e}.Coefficients.tStat,lme_width{e}.Coefficients.pValue],3) ; 

writematrix(resultstable,'LME_fig4_width.xlsx','Sheet',e)

% Extract p-values from the fixed effects table
p = fixedEffectsTable.pValue(2)   
        if  p(1,1) < 0.0005
            txt{1,1} = '***';
        elseif  p(1,1) < 0.005
            txt{1,1} = '**';
        elseif  p(1,1) < 0.05
            txt{1,1} = '*';
        elseif p(s,1) > 0.05
            txt{1,1} = 'n.s.';
        end
        
        line( linspace(1.1, 1.9), linspace(700 ,700),'color','k')
        hold on
        text(mean([1.5 1.5]), 750 ,txt{1,1},'FontSize',20)
        hold on
end
%% maximal discrimination (AUC)

clear mouse rec var
var = GNG_analysis.max_AUC ;
mouse = mouse_idx_th ;
rec = rec_idx_th ;

for e = 1:3
    for g = 1:numel(GNG_rec_all_cell_exp)
        var{e,g} = var{e,g}(~isnan(var{e,g})) ;
        mouse{e,g} = mouse{e,g}(~isnan(mouse{e,g})) ;
        rec{e,g} = rec{e,g}(~isnan(rec{e,g})) ;
              mean_max(e,g) = mean (var{e,g}) ; 
        ste_max(e,g) = std (var{e,g}) / sqrt (size(var{e,g},1)) ;
    end
end




figure
for e = a
    subplot(1,3,e)
    violinplot([GNG_analysis.max_AUC{e,1},  GNG_analysis.max_AUC{e,2}]...
        ,{' adolescenct','adult'},"ViolinColor",...
        {[  color_eh{e};  color_eh{e}]} ) ; 
    ylim ([0.5 1.1])
    ylabel('max discrimination (AUC)')
    
    box off;
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;


    LME_AUC_max = table; 
LME_AUC_max.maxi = [  var{e,1}'  , var{e,2}']' ;
LME_AUC_max.mouse = [  mouse{e,1}'  , mouse{e,2}']' ;
LME_AUC_max.rec = [ rec{e,1}'  , rec{e,2}']' ;
LME_AUC_max.age = [ ones(size(var{e,1})) ;  ones(size(var{e,2}))+ 1] ;

% Convert appropriate columns to categorical variables
LME_AUC_max.mouse = categorical(LME_AUC_max.mouse);
LME_AUC_max.rec = categorical(LME_AUC_max.rec);
LME_AUC_max.age = categorical(LME_AUC_max.age);

formula = 'maxi ~ age  + (1|mouse) + (1|rec)';
lme_max{e} = fitlme(LME_AUC_max, formula);
% Extract the fixed effects table
fixedEffectsTable = lme_max{e}.Coefficients;

% Extract p-values from the fixed effects table
p = fixedEffectsTable.pValue(2)


resultstable = round([lme_max{e}.Coefficients.Estimate,lme_max{e}.Coefficients.SE,lme_max{e}.Coefficients.tStat,lme_max{e}.Coefficients.pValue],3) ; 

writematrix(resultstable,'LME_fig4_max.xlsx','Sheet',e)
    
        
        if  p(1,1) < 0.0005
            txt{1,1} = '***';
        elseif  p(1,1) < 0.005
            txt{1,1} = '**';
        elseif  p(1,1) < 0.05
            txt{1,1} = '*';
        elseif p(s,1) > 0.05
            txt{1,1} = 'n.s.';
        end
        
        line( linspace(1.1, 1.9), linspace(1 ,1),'color','k')
        hold on
        text(mean([1.5 1.5]), 1.05 ,txt{1,1},'FontSize',20)
        hold on
end
clc