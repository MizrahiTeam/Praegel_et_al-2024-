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
a = [5 11] ;
b = [6 12] ;

% AUC parameters
run_window = 50 ; % ms
it_size = 25 ; % ms
n_shuffle = 10 ; % n shuffles for AUC
trial_samples = 10 ;
nan_units = 200 ; % to preallocate
min_n_trials = 25 ;

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

% modulation type

% colors
color_eh= {[0.3010 0.7450 0.9330], [0 0.4470 0.7410], [0 0 1]};
color_eh_fade = {[0.3010 0.7450 0.9330 0.3], [0 0.4470 0.7410 0.3], [0 0 1 0.3]};
color_eh_patch = {[0 0.4470 0.7410],[0.3010 0.7450 0.9330], [0 0.4470 0.7410],[0.3010 0.7450 0.9330]};

L = {['--'], ['-'],['--'], ['-']};
M = {['o'],['o'],['o'], ['o']};

%% extract behavior per recording

[behavior] = GNG_neuro_behavior (GNG_rec_all_cell) ;

%% extract trial samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shuffeled datasamples will lead to different distributions than in paper
%figures. Yet these do not affect the significance and effect sizes of the
%results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
tic
[A_k_ge, B_k_ge, idx_k_ge, T_k_ge, deasy_k_ge, dhard_k_ge, rec_idx, neuron_idx, mouse_idx]...
    = GNG_t_per_n_AUC (FR_array, GNG_rec_all_cell, behavior, areas, a, b,  trial_samples, min_n_trials) ;
toc

%% run the AUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%shuffeled datasamples will lead to different distributions than in paper
%figures. Yet these do not affect the significance and effect sizes of the
%results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
tic
[r_AUC, r_AUC_abs, r_AUC_shuf, r_AUC_shuf_abs]...
    = GNG_running_AUC (GNG_rec_all_cell,a,b, A_k_ge, B_k_ge, trial_samples, it_size, run_window, window_length_ms, n_shuffle)
toc

GNG_analysis.r_AUC = r_AUC ;
GNG_analysis.r_AUC_abs = r_AUC_abs ;
GNG_analysis.r_AUC_shuf = r_AUC_shuf ;
GNG_analysis.r_AUC_shuf_abs =  r_AUC_shuf_abs ;

%% normalize AUC to baseline AUC
clc
clc
per_mean = false ;
a = [1 2] ;%  easy  hard go and no go
stim_onset_bin = 7 ;
[ind_AUC, ind_AUC_shuf, ind_AUC_all,ind_AUC_shuf_all,std_ind_AUC_shuf_all, mean_vec, stde_vec,mean_vec_shuf, std_vec_shuf ]...
    = GNG_ind_mean_AUC (GNG_rec_all_cell, GNG_analysis.r_AUC_abs, GNG_analysis.r_AUC_shuf_abs,...
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

%% latency, duration and max discriminability calculation
clc

a = [1 2 ] ; %  easy  hard go and no go
[latency_peak_AUC, onset_latency_AUC, width_AUC, max_AUC, rec_idx_th, mouse_idx_th]...
    = GNG_max_onset_width_AUC (GNG_rec_all_cell, GNG_analysis.ind_AUC_all, GNG_analysis.ind_AUC_shuf_all,...
    GNG_analysis.std_ind_AUC_shuf_all, a, stim_onset_bin,it_size,rec_idx, mouse_idx) ;

for e = a
    for g = 1:numel(GNG_rec_all_cell)
        GNG_analysis.max_AUC{e,g} = max_AUC{e,g} ;
        GNG_analysis.width_AUC{e,g}  = width_AUC{e,g} ;
        GNG_analysis.onset_latency_AUC{e,g}  =onset_latency_AUC{e,g} ;
    end
end
%% Figure 4 E & F  plot a single neuron example
a = [ 1 2 ] ; % easy and hard task
ac = [1 2 ] ; % color code
close all

% one example per group #cell
example = [208 10 43 123] ;

for   g = 1:4
    figure
    for e = 1:length(a)

        plot( ind_AUC_all{e,g}(example(g),:),'Color',color_eh{ac(e)},...
            'linestyle',L{g},'linewidth',4);
        hold on

        plot( ind_AUC_shuf_all{e,g}(example(g),:),'Color',[0.5 0.5 0.5],...
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

    hit_idx = idx_k_ge{2,g}{example(g)}(13,:)  ;
    fa_idx =  idx_k_ge{2,g}{example(g)} (14,:) ;
    miss_idx =  idx_k_ge{2,g}{example(g)} (15,:) ;
    cr_idx  =  idx_k_ge{2,g}{example(g)} (16,:) ;

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
    trial_types_raster = cat(1,T_k_ge{2,g}{example(g)}(hit_idx,:),T_k_ge{2,g}{example(g)}(miss_idx,:), T_k_ge{2,g}{example(g)}(fa_idx,:),T_k_ge{2,g}{example(g)}(cr_idx,:)) ;

    figure
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
    yticks([1  sum(size_tr) ])
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



%% plot the AUC mean
close all
ac = a ;  % color code
% plot mean (optional include shuffled mean)
for e = 1:length(a)
    %figure% (1)
    for   g = 1:numel(GNG_rec_all_cell) % run per group
        figure(g)
        length_bins = size(  mean_vec{a(e),g},2) ;
        % subplot(2,2,g)
        plot([1, 1:length_bins-1], mean_vec{a(e),g},'Color',color_eh{ac(e)},...
            'linestyle',L{g},'linewidth',2);
        hold on

        patch([[1,1:length_bins-1] flip([1,1:length_bins-1])] , [  mean_vec{a(e),g} + ....
            stde_vec{a(e),g} flip(  mean_vec{a(e),g} -  stde_vec{a(e),g})], color_eh{ac(e)} ,...
            'facealpha' , 0.2, 'EdgeColor','none')
        hold on

        plot([1, 1:length_bins-1],   mean_vec_shuf{a(e),g}+ ....
            std_vec_shuf{a(e),g},'Color',[0.5 0.5 0.5],...
            'linestyle',L{g},'linewidth',2);
        hold on


        xline(6,'--','Color','k');
        xticks([1 6  18])
        xticklabels([-200 0  400])
        xlim([1 (length_bins)])
        yticks([0.5:0.05:0.6])

        ylim([0.48 0.6])
        xlabel('time (ms)')
        ylabel(' discrimination (AUC)')
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');
        box off;
    end
end
%%  statistics AUC onset latency

close all
clc
clear mouse rec var

var = GNG_analysis.onset_latency_AUC ;
mouse = mouse_idx_th ;
rec = rec_idx_th ;

for e = a
    for g = 1:numel(GNG_rec_all_cell)
        var{e,g} = var{e,g}(~isnan(var{e,g})) ;
        mouse{e,g} = mouse{e,g}(~isnan(mouse{e,g})) ;
        rec{e,g} = rec{e,g}(~isnan(rec{e,g})) ;
    end
end

LME_AUC_onset = table;
LME_AUC_onset.onset = [  var{1,1}'  , var{1,2}' ,  var{1,3}' ,  var{1,4}' , var{2,1}'  , var{2,2}' ,  var{2,3}' ,  var{2,4}' ]' ;
LME_AUC_onset.mouse = [  mouse{1,1}'  , mouse{1,2}' ,  mouse{1,3}' ,  mouse{1,4}' , mouse{2,1}'  , mouse{2,2}' ,  mouse{2,3}' ,  mouse{2,4}' ]' ;
LME_AUC_onset.rec =  [  rec{1,1}'  , rec{1,2}' ,  rec{1,3}' ,  rec{1,4}' , rec{2,1}'  , rec{2,2}' ,  rec{2,3}' ,  rec{2,4}' ]' ;
LME_AUC_onset.age = [ ones(size(var{1,1})) ;  ones(size(var{1,2}))+ 1 ;ones(size(var{1,3})) ;  ones(size(var{1,4}))+ 1 ;....
    ones(size(var{2,1})) ;  ones(size(var{2,2}))+ 1 ;ones(size(var{2,3})) ;  ones(size(var{2,4}))+ 1] ;
LME_AUC_onset.learning =  [ ones(size(var{1,1})) ;  ones(size(var{1,2})) ; ones(size(var{1,3})) + 1 ;  ones(size(var{1,4}))+ 1 ;....
    ones(size(var{2,1})) ;  ones(size(var{2,2})) ;ones(size(var{2,3})) + 1 ;  ones(size(var{2,4}))+ 1] ;
LME_AUC_onset.difficulty = [ ones(size(var{1,1})) ;  ones(size(var{1,2})) ; ones(size(var{1,3}))  ;  ones(size(var{1,4})) ;....
    ones(size(var{2,1})) + 1 ;  ones(size(var{2,2})) + 1 ;ones(size(var{2,3})) + 1 ;  ones(size(var{2,4}))+ 1] ;
% Convert appropriate columns to categorical variables
LME_AUC_onset.mouse = categorical(LME_AUC_onset.mouse);
LME_AUC_onset.rec = categorical(LME_AUC_onset.rec);
LME_AUC_onset.age = categorical(LME_AUC_onset.age);
LME_AUC_onset.learning = categorical(LME_AUC_onset.learning);
LME_AUC_onset.difficulty = categorical(LME_AUC_onset.difficulty);

formula = 'onset ~ age* learning * difficulty + (1|mouse)  + (1|rec)';

lme_onset = fitlme(LME_AUC_onset, formula);
disp(lme_onset)

resultstable = round([lme_onset.Coefficients.Estimate,lme_onset.Coefficients.SE,lme_onset.Coefficients.tStat,lme_onset.Coefficients.pValue],3) ;

writematrix(resultstable,'LME_fig6_onset.xlsx','Sheet',1)

%% Figure 6 C - plot the onset latency
clc
close all
var = GNG_analysis.onset_latency_AUC ;
pos_text = [1.35 3.35 1.85 2.85] ;
height_text = [750 750  850 950] ;

pose_left = [1.1 3.1 1.1 1.9] ;
pos_right = [1.9 3.9 2.9 3.9] ;
height_line = [700 700 800 900] ;

s1 = [6 1 2 5] ;

for e = a
    stat = [] ;
    c = [] ;
    [~,~,stats] = kruskalwallis([var{e,1} var{e,2} var{e,3} var{e,4}],[],'off');
    c = multcompare(stats,"Display","off");

    figure
    violinplot([ var{e,3},  var{e,4},var{e,1},  var{e,2}]...
        ,{'naive adolescent','naive adult','expert adolescent','expert adult'},"ViolinColor",...
        {[color_eh{e}; color_eh{e};  color_eh{e};  color_eh{e}]} ) ;
    ylim ([0 1000])
    ylabel('time(ms)')

    box off;
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;


    for s = 1:4
        p_e(s,e) = c(s1(s),6) ;

        if  p_e(s,e) < 0.0005
            txt = '***';
        elseif  p_e(s,e) < 0.005
            txt = '**';
        elseif  p_e(s,e) < 0.05
            txt= '*';
        elseif p_e(s,e) > 0.05
            txt = 'n.s.';
        end

        line( linspace(pose_left(s),pos_right(s)), linspace(height_line(s),height_line(s)),'color','k')
        hold on
        text(mean([pos_text(s) pos_text(s)]), height_text(s) ,txt,'FontSize',20)
        hold on

    end
end
disp(p_e)


%%  statistics AUC width
clc
clear mouse rec var

var = GNG_analysis.width_AUC ;
mouse = mouse_idx_th ;
rec = rec_idx_th ;

for e = a
    for g = 1:numel(GNG_rec_all_cell)
        var{e,g} = var{e,g}(~isnan(var{e,g})) ;
        mouse{e,g} = mouse{e,g}(~isnan(mouse{e,g})) ;
        rec{e,g} = rec{e,g}(~isnan(rec{e,g})) ;
    end
end

LME_AU_width = table;
LME_AU_width.width = [  var{1,1}'  , var{1,2}' ,  var{1,3}' ,  var{1,4}' , var{2,1}'  , var{2,2}' ,  var{2,3}' ,  var{2,4}' ]' ;
LME_AU_width.mouse = [  mouse{1,1}'  , mouse{1,2}' ,  mouse{1,3}' ,  mouse{1,4}' , mouse{2,1}'  , mouse{2,2}' ,  mouse{2,3}' ,  mouse{2,4}' ]' ;
LME_AU_width.rec =  [  rec{1,1}'  , rec{1,2}' ,  rec{1,3}' ,  rec{1,4}' , rec{2,1}'  , rec{2,2}' ,  rec{2,3}' ,  rec{2,4}' ]' ;
LME_AU_width.age = [ ones(size(var{1,1})) ;  ones(size(var{1,2}))+ 1 ;ones(size(var{1,3})) ;  ones(size(var{1,4}))+ 1 ;....
    ones(size(var{2,1})) ;  ones(size(var{2,2}))+ 1 ;ones(size(var{2,3})) ;  ones(size(var{2,4}))+ 1] ;
LME_AU_width.learning =  [ ones(size(var{1,1})) ;  ones(size(var{1,2})) ; ones(size(var{1,3})) + 1 ;  ones(size(var{1,4}))+ 1 ;....
    ones(size(var{2,1})) ;  ones(size(var{2,2})) ;ones(size(var{2,3})) + 1 ;  ones(size(var{2,4}))+ 1] ;
LME_AU_width.difficulty = [ ones(size(var{1,1})) ;  ones(size(var{1,2})) ; ones(size(var{1,3}))  ;  ones(size(var{1,4})) ;....
    ones(size(var{2,1})) + 1 ;  ones(size(var{2,2})) + 1 ;ones(size(var{2,3})) + 1 ;  ones(size(var{2,4}))+ 1] ;
% Convert appropriate columns to categorical variables
LME_AU_width.mouse = categorical(LME_AU_width.mouse);
LME_AU_width.rec = categorical(LME_AU_width.rec);
LME_AU_width.age = categorical(LME_AU_width.age);
LME_AU_width.learning = categorical(LME_AU_width.learning);
LME_AU_width.difficulty = categorical(LME_AU_width.difficulty);


formula = 'width ~ age* learning * difficulty + (1|mouse)  + (1|rec)';

lme_width = fitlme(LME_AU_width, formula);
disp(lme_width)


resultstable = round([lme_width.Coefficients.Estimate,lme_width.Coefficients.SE,lme_width.Coefficients.tStat,lme_width.Coefficients.pValue],3) ;

writematrix(resultstable,'LME_fig6_width.xlsx','Sheet',1)

%% Figure 6 D - plot the AUC width
var = GNG_analysis.width_AUC ;

pos_text = [1.35 3.35 1.85 2.85] ;
height_text = [750 750  850 950] ;
pose_left = [1.1 3.1 1.1 1.9] ;
pos_right = [1.9 3.9 2.9 3.9] ;
height_line = [700 700 800 900] ;
s1 = [6 1 2 5] ;


for e = a

    [~,~,stats] = kruskalwallis([var{e,1} var{e,2} var{e,3} var{e,4}],[],'off');
    c = multcompare(stats,"Display","off");

    figure
    violinplot([ var{e,3},  var{e,4},var{e,1},  var{e,2}]...
        ,{'naive adolescent','naive adult','expert adolescent','expert adult'},"ViolinColor",...
        {[color_eh{e}; color_eh{e};  color_eh{e};  color_eh{e}]} ) ;
    ylim ([0 1000])
    ylabel('time(ms)')

    box off;
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;

    for s = 1:4

        p_e(s,e) = c(s1(s),6) ;

        if   p_e(s,e) < 0.0005
            txt = '***';
        elseif   p_e(s,e) < 0.005
            txt = '**';
        elseif   p_e(s,e) < 0.05
            txt= '*';
        elseif  p_e(s,e) > 0.05
            txt = 'n.s.';
        end

        line( linspace(pose_left(s),pos_right(s)), linspace(height_line(s),height_line(s)),'color','k')
        hold on
        text(mean([pos_text(s) pos_text(s)]), height_text(s) ,txt,'FontSize',20)
        hold on

    end
end
disp( p_e)
%%  statistics max AUC

clc
clear mouse rec var

var = GNG_analysis.max_AUC ;
mouse = mouse_idx_th ;
rec = rec_idx_th ;

for e = a
    for g = 1:numel(GNG_rec_all_cell)
        var{e,g} = var{e,g}(~isnan(var{e,g})) ;
        mouse{e,g} = mouse{e,g}(~isnan(mouse{e,g})) ;
        rec{e,g} = rec{e,g}(~isnan(rec{e,g})) ;
    end
end
LME_AU_max = table;
LME_AU_max.max = [  var{1,1}'  , var{1,2}' ,  var{1,3}' ,  var{1,4}' , var{2,1}'  , var{2,2}' ,  var{2,3}' ,  var{2,4}' ]' ;
LME_AU_max.mouse = [  mouse{1,1}'  , mouse{1,2}' ,  mouse{1,3}' ,  mouse{1,4}' , mouse{2,1}'  , mouse{2,2}' ,  mouse{2,3}' ,  mouse{2,4}' ]' ;
LME_AU_max.rec =  [  rec{1,1}'  , rec{1,2}' ,  rec{1,3}' ,  rec{1,4}' , rec{2,1}'  , rec{2,2}' ,  rec{2,3}' ,  rec{2,4}' ]' ;
LME_AU_max.age = [ ones(size(var{1,1})) ;  ones(size(var{1,2}))+ 1 ;ones(size(var{1,3})) ;  ones(size(var{1,4}))+ 1 ;....
    ones(size(var{2,1})) ;  ones(size(var{2,2}))+ 1 ;ones(size(var{2,3})) ;  ones(size(var{2,4}))+ 1] ;
LME_AU_max.learning =  [ ones(size(var{1,1})) ;  ones(size(var{1,2})) ; ones(size(var{1,3})) + 1 ;  ones(size(var{1,4}))+ 1 ;....
    ones(size(var{2,1})) ;  ones(size(var{2,2})) ;ones(size(var{2,3})) + 1 ;  ones(size(var{2,4}))+ 1] ;
LME_AU_max.difficulty = [ ones(size(var{1,1})) ;  ones(size(var{1,2})) ; ones(size(var{1,3}))  ;  ones(size(var{1,4})) ;....
    ones(size(var{2,1})) + 1 ;  ones(size(var{2,2})) + 1 ;ones(size(var{2,3})) + 1 ;  ones(size(var{2,4}))+ 1] ;
% Convert appropriate columns to categorical variables
LME_AU_max.mouse = categorical(LME_AU_max.mouse);
LME_AU_max.rec = categorical(LME_AU_max.rec);
LME_AU_max.age = categorical(LME_AU_max.age);
LME_AU_max.learning = categorical(LME_AU_max.learning);
LME_AU_max.difficulty = categorical(LME_AU_max.difficulty);


formula = 'max ~ age* learning * difficulty + (1|mouse)  + (1|rec)';

lme_max = fitlme(LME_AU_max, formula);
disp(lme_max)


resultstable = round([lme_max.Coefficients.Estimate,lme_max.Coefficients.SE,lme_max.Coefficients.tStat,lme_max.Coefficients.pValue],3) ;

writematrix(resultstable,'LME_fig6_max.xlsx','Sheet',1)

%% Figure 6 E - plot the AUC max
var = GNG_analysis.max_AUC ;

pos_text = [1.35 3.35 1.85 2.85] ;
height_text = [1.05 1.05  1.15 1.25] ;

pose_left = [1.1 3.1 1.1 1.9] ;
pos_right = [1.9 3.9 2.9 3.9] ;
height_line = [1 1 1.1 1.2] ;

s1 = [6 1 2 5] ;


for e = a

    [~,~,stats] = kruskalwallis([var{e,1} var{e,2} var{e,3} var{e,4}],[],'off');
    c = multcompare(stats,"Display","off");

    figure
    violinplot([ var{e,3},  var{e,4},var{e,1},  var{e,2}]...
        ,{'naive adolescent','naive adult','expert adolescent','expert adult'},"ViolinColor",...
        {[color_eh{e}; color_eh{e};  color_eh{e};  color_eh{e}]} ) ;
    ylim ([0.5 1.3])
    yticks ([0.5:0.25:1])

    ylabel('max. discrimination')

    box off;
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;

    for s = 1:4
        p_e(s,e) = c(s1(s),6) ;

        if   p_e(s,e) < 0.0005
            txt = '***';
        elseif   p_e(s,e) < 0.005
            txt = '**';
        elseif   p_e(s,e) < 0.05
            txt= '*';
        elseif  p_e(s,e) > 0.05
            txt = 'n.s.';
        end


        line( linspace(pose_left(s),pos_right(s)), linspace(height_line(s),height_line(s)),'color','k')
        hold on
        text(mean([pos_text(s) pos_text(s)]), height_text(s) ,txt,'FontSize',20)
        hold on

    end
end
disp( p_e)

%%
close all

var = GNG_analysis.max_AUC ;
mouse = mouse_idx_th ;
rec = rec_idx_th ;
for e = a
    for g = 1:numel(GNG_rec_all_cell)
        var{e,g} = var{e,g}(~isnan(var{e,g})) ;
        rec{e,g} = rec{e,g}(~isnan(rec{e,g})) ;
    end
end
AUC_mean = nan (2,4, 15) ;
dprime = nan (2,4, 15) ;
for e = a
    for g = 1:numel(GNG_rec_all_cell)
        recs_i = [] ;
        recs_i = unique(rec{e,g}) ;
        for i = 1:length(unique(rec{e,g}))
            idx_c = [] ;
            idx_c =  find (rec{e,g} == recs_i(i)) ;
            AUC_mean(e,g,i) = mean(var{e,g}(idx_c)) ;

            if e == 1
                dprime(e,g,i) = behavior.dprime_easy(g,i) ;
            elseif e == 2
                dprime(e,g,i) = behavior.dprime_hard(g,i) ;
            end

        end
    end
end



for e = a
    for g = [1 2]

        dprime_e_g = [ squeeze(dprime(e,g,:)); squeeze(dprime(e,g+2,:))]  ;
        dprime_e_g = dprime_e_g(~isnan(dprime_e_g)) ;
        AUC_e_g =   [ squeeze(AUC_mean(e,g,:));  squeeze(AUC_mean(e,g+2,:))] ;
        AUC_e_g = AUC_e_g(~isnan(AUC_e_g)) ;

        figure
        scatter(dprime_e_g, AUC_e_g)
        hold on
        mdl = fitlm(dprime_e_g ,AUC_e_g)
        h =  plot(mdl)
        hold on
        scatter(squeeze(dprime(e,g,:)),  squeeze(AUC_mean(e,g,:)),'Marker','o','MarkerFaceColor',color_eh{e},'MarkerEdgeColor',color_eh{e})
        hold on
        scatter(squeeze(dprime(e,g+2,:)),  squeeze(AUC_mean(e,g+2,:)),'Marker','o','MarkerFaceColor','none','MarkerEdgeColor',color_eh{e})

        [R,P] = corr(dprime_e_g, AUC_e_g) ;
        txt =  ['R = ', num2str(R)];
        text(0, 0.67 ,txt,'FontSize',12)
        txt2 =  ['p  = ',num2str(P)];
        text(0, 0.68 ,txt2,'FontSize',12)

        hold on
        delete(h(1))
        h(2).Color = color_eh{e} ;
        h(3).Color = color_eh{e} ;
        h(2).LineWidth = 5 ;
        h(3).LineWidth =  5 ;
        ylim([0.55 0.7])
        xlim([-0.5 3.5])

        lgd = findobj('type', 'legend')
        delete(lgd)
        title('')
        ylabel('neuronal disc. (AUC)')
        xlabel("behavioral disc. (d')")
        box off
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');
        box off;

    end
end

