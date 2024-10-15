clear all
clc
close all
%% load the cell array of the behavioral session and the recording table
%add the Praegel_et_al_MATLABR2023b_scripts directory here
addpath(genpath('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_MATLABR2023b_scripts'))
% selection training session
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , ['Select training sessions: GNGtable_cell ']);
addpath(path)
load (file)
% select recording session
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , ['Select recording sessions: GNG_rec_all_cell ']);
addpath(path)
load (file)
cd(path)
% load excel file for stimulus identity
[file, path] = uigetfile('*.xlsx', ['Select the stimulus file: Disc_Rec_GNG ']);
tfreq = readtable(fullfile(path, file));

%% Parameters
tone_dur = 0.1 ; % sec
response_window = 2 ; % sec
dprime_threshold = 1 ;
min_binsize = 30;
startRange = 0; %beginning of baseline
stopRange = 0.1; %end of the offset
startReinf = 0.6 ; % start reinforcement
length_base = 200 ; %ms
reinforcement = 500 + length_base ; %ms
max_trial = 1000 ; % n trials
trial_cut_off = 100 ; % n tirals
startRange_resp = - 0.1 ;% sec
stopRange_resp = 2 ; %sec
binSize = 0.001 ; %sec
smoothSize = 5 ; %ms

%stimulus parameters
length_stim = (abs(startRange_resp) + stopRange_resp)*1000 ;
stim = [1 2 3 4 5 6 7];
freqs =sort(tfreq.Frequency(tfreq.StimulusID(stim(1:4))))';
catch_freqs =sort(tfreq.Frequency(tfreq.StimulusID(stim)))';

binBorders = startRange_resp:binSize:stopRange_resp ; %define the time points
numBins = length(binBorders)-1 ;


% Figure Parameters
color_eh = {[0.3010 0.7450 0.9330], [0 0.4470 0.7410]} ;
color_eh_mean = {[0.3010 0.7450 0.9330], [0 0.4470 0.7410]} ;
color_eh_plot = {[0.3010 0.7450 0.9330 0.4], [0 0.4470 0.7410 0.4]} ;
color_eh_all = {[0 0 1]} ;
color_eh_all_fade = {[0 0 1 0.4]} ;
color_mean_marker = {[0.5 0.5 0.5]} ;
L = {['-'], ['-']};
L2 = {['--'], ['-'],['--'], ['-']};


%% calculate the d' per training session per mouse
for g   = 1:numel(GNGtable_cell)

    for  n = 1:height (GNGtable_cell{1,g});
        for s = 1:length(GNGtable_cell{1,g});
            if ~isempty(GNGtable_cell{1,g}{n,s})

                stim_types = GNGtable_cell{1,g}{n,s}.stim_types;
                stim_ids   = GNGtable_cell{1,g}{n,s}.stim_ids;
                trial_responses = GNGtable_cell{1,g}{n,s}.trial_responses;
                lick_times = GNGtable_cell{1,g}{n,s}.lick_times;
                stim_times = GNGtable_cell{1,g}{n,s}.stim_times;

                % size of the training session
                size_session(g,n,s) = length (stim_ids) ;

                if max(stim_ids) == 2
                    % easy task level 3
                    stim = [1 2];

                    index_stim = find (stim_ids == stim(1) | stim_ids == stim(2));

                    [~,~,dprimes_easy{g,1}(n,s),~,~,~] = GNG_lick_rt_dprime (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;


                elseif max(stim_ids) == 4
                    stim = [1 4 2 3];
                    index_stim = find (stim_ids == stim(1) | stim_ids == stim(2));
                    % easy task level 4
                    [~,~,dprimes_easy{g,1}(n,s),~,~,~] = GNG_lick_rt_dprime(index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;
                    index_stim = find (stim_ids == stim(3) | stim_ids == stim(4));
                    %hard task level 4
                    [~,~,dprimes_hard{g,1}(n,s),~,~,~] = GNG_lick_rt_dprime (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;
                end

            end
        end
        % d' threshold
        idx_th =  find (dprimes_easy{g,1}(n,:) >= dprime_threshold) ;
        dprime_th(g,n) = idx_th(1) ;
    end
end
%% stats of sessions per group
% size of the session
size_session (size_session == 0) = nan ;
for i = 1:2
    h(i) = kstest(reshape(size_session(1,:,:),[],1)) ;
end
% test according to the distribution of h
if sum(h) > 0
    [p,~,~] = ranksum (reshape(size_session(1,:,:),[],1),reshape(size_session(2,:,:),[],1),'alpha',0.05,'tail','both') ;

elseif sum(h) == 0
    [~,p,~] =  ttest2 (reshape(size_session(1,:,:),[],1),reshape(size_session(2,:,:),[],1),'alpha',0.05,'tail','both') ;;
end
disp(p)
disp(nanmean(reshape(size_session(1,:,:),[],1)))
disp(nanstd(reshape(size_session(1,:,:),[],1))/sqrt(size(reshape(size_session(1,:,:),[],1),1)))
disp(nanmean(reshape(size_session(2,:,:),[],1)))
disp(nanstd(reshape(size_session(2,:,:),[],1))/sqrt(size(reshape(size_session(2,:,:),[],1),1)))

size_session (size_session == 0) = nan ;
for i = 1:2
    h(i) = kstest(reshape(dprime_th,[],1)) ;
end


% procedural learning
dprime_th(dprime_th == 0 ) = nan ;
% test according to the distribution of h
h = kstest(reshape(dprime_th,[],1)) ;
if h > 0
    [p,~,~] = ranksum (dprime_th(1,:), dprime_th (2,:),'alpha',0.05,'tail','both') ;

elseif h == 0
    [~,p,~] =  ttest2 (dprime_th(1,:), dprime_th (2,:),'alpha',0.05,'tail','both') ;;
end
disp(p)
%% FIGURE 2D: dprime  individual per session
clc
close all
for g   = 1:numel(GNGtable_cell)
    dprimes_easy{g,1} (dprimes_easy{g,1}== 0) = nan ;
    dprimes_hard{g,1} (dprimes_hard{g,1}== 0) = nan ;

    figure(101+g)
    for  n = 1:height (GNGtable_cell{1,g});
        plot(dprimes_easy{g,1}(n,:)','Color',color_eh_mean{1},'linewidth',2,'linestyle',L{g})
        hold on
        plot (dprimes_hard{g,1}(n,:)','Color',color_eh_mean{2},'linewidth',2,'linestyle',L{g})
    end

    ylim([-1 4]);
    yticks([-1:1:4])
    xlim([0 12])
    yline(1,'--','Color','k');
    xlabel('session')
    ylabel("discrimination (d')")
    box off;
    makepretty;
    hold on
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    yticks([-1:2:4])
    xticks([1:4:12])
    movegui('east');
end
clc
%% Figure 2E left: Comparison of first to last session dprime easy
close all

for g   = 1:numel(GNGtable_cell)
    dprimes_easy{g,1}(:,size(dprimes_easy{g,1},2)+1) = 0 ;
    dprimes_easy{g,1}(  dprimes_easy{g,1}== 0) = nan;

    for n = 1:height(dprimes_easy{g,1})

        all = find(~isnan(dprimes_easy{g,1}(n,:)));
        dprime_startend{g,n} = [dprimes_easy{g,1}(n,min(all)) dprimes_easy{g,1}(n,max(all))] ;
        dprime_start(g,n) = dprimes_easy{g,1}(n,min(all)) ;
        dprime_end(g,n) =  dprimes_easy{g,1}(n,max(all)) ;
    end
    mean_dprime_startend(g,:) = [mean(dprime_start(g,:)'), mean(dprime_end(g,:)')]
    ste_dprime_startend(g,:) = [std(dprime_start(g,:)')/sqrt(size(dprime_start(g,:)',1)), std(dprime_end(g,:)')/sqrt(size(dprime_end(g,:)',1))]
end



figure
for n = 1:height(dprimes_easy{1,1})
    plot([1 2],dprime_startend{1,n},'Color',[0.5 0.5 0.5 0.4],'linewidth',2,'linestyle',L{g},...
        'Marker','.','MarkerFaceColor',color_eh_mean{1},'MarkerEdgeColor',color_eh_mean{1})
    hold on
end
for n = 1:height(dprimes_easy{2,1})
    plot([3 4],dprime_startend{2,n},'Color',[0.5 0.5 0.5 0.4],'linewidth',2,'linestyle',L{g},...
        'Marker','.','MarkerFaceColor',color_eh_mean{1},'MarkerEdgeColor',color_eh_mean{1})
    hold on
end

errorbar( mean_dprime_startend(1,:), ste_dprime_startend(1,:) ,'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)
errorbar( [nan nan mean_dprime_startend(2,:)], [nan nan ste_dprime_startend(2,:)] ,'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)

ylim([-1 5]);
yticks([-1:1:5])
xlim([0 5])
xticks([1 2 3 4])
yline(1,'--','Color','k');
ylabel("d'")
xticklabels({'adolescent naive','adolescent exp.','adult naive','adult exp.'})
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');



h = kstest (dprime_start)
if h > 0
    p = [] ;
    [p(:,3),~,~] = ranksum(dprime_start(1,:),dprime_start(2,:),'alpha',0.05,'tail','both')
    [p(:,4),~,~] = ranksum(dprime_end(1,:),dprime_end(2,:),'alpha',0.05,'tail','both')
    [~,p(:,3:4),~]= bonferroni_holm(p(:,3:4),0.05) ;
    [p(:,1),~,~] = signrank(dprime_start(1,:),dprime_end(1,:),'alpha',0.05,'tail','both')
    [p(:,2),~,~] = signrank(dprime_start(2,:),dprime_end(2,:),'alpha',0.05,'tail','both')
    [~,p(:,1:2),~]= bonferroni_holm(p(:,1:2),0.05) ;

elseif h == 0
    p = [] ;
    [~,p(:,3),~] = ttest2(dprime_start(1,:),dprime_start(2,:),'alpha',0.05,'tail','both')
    [~,p(:,4),~] = ttest2(dprime_end(1,:),dprime_end(2,:),'alpha',0.05,'tail','both')
    [~,p(:,3:4),~]= bonferroni_holm(p(:,3:4),0.05) ;
    [~,p(:,1),~] = ttest(dprime_start(1,:),dprime_end(1,:),'alpha',0.05,'tail','both')
    [~,p(:,2),~] = ttest(dprime_start(2,:),dprime_end(2,:),'alpha',0.05,'tail','both')
    [~,p(:,1:2),~]= bonferroni_holm(p(:,1:2),0.05) ;
end
disp(p)
if  p(:,1) < 0.0005
    txt = '***';
elseif  p(:,1) < 0.005
    txt = '**';
elseif  p(:,1) < 0.05
    txt = '*';
elseif p(:,1)> 0.05
    txt = 'n.s.';
end
line( linspace(1.05,1.95), linspace(4,4),'color','k')
hold on
text(mean([1.5 1.5]), 4.2 ,txt,'FontSize',20)

if  p(:,2) < 0.0005
    txt = '***';
elseif  p(:,2) < 0.005
    txt = '**';
elseif  p(:,2) < 0.05
    txt = '*';
elseif p(:,2) > 0.05
    txt = 'n.s.';
end
line( linspace(3.05,3.95), linspace(4,4),'color','k')
hold on
text(mean([3.5 3.5]), 4.2 ,txt,'FontSize',20)

if  p(:,3) < 0.0005
    txt = '***';
elseif  p(:,3) < 0.005
    txt = '**';
elseif  p(:,3) < 0.05
    txt = '*';
elseif p(:,3) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,3.05), linspace(4.4,4.4),'color','k')
hold on
text(mean([2.05 2.05]), 4.6 ,txt,'FontSize',20)

if  p(:,4) < 0.0005
    txt = '***';
elseif  p(:,4) < 0.005
    txt = '**';
elseif  p(:,4) < 0.05
    txt = '*';
elseif p(:,4) > 0.05
    txt = 'n.s.';
end
line( linspace(1.95,3.95), linspace(4.8,4.8),'color','k')
hold on
text(mean([3 3]), 5 ,txt,'FontSize',20)


%% Figure 2E right: Comparison of first to last session dprime hard
for g   = 1:numel(GNGtable_cell)
    dprimes_hard{g,1}(:,size(dprimes_hard{g,1},2)+1) = 0 ;
    dprimes_hard{g,1}(isnan(dprimes_hard{g,1})) = 0 ;

    for n = 1:height(dprimes_hard{g,1})
        all = find(dprimes_hard{g,1}(n,:)~= 0);
        dprime_startend{g,n} = [dprimes_hard{g,1}(n,min(all)) dprimes_hard{g,1}(n,max(all))] ;
        dprime_start(g,n) = dprimes_hard{g,1}(n,min(all)) ;
        dprime_end(g,n) =  dprimes_hard{g,1}(n,max(all)) ;
    end

    mean_dprime_startend(g,:) = [mean(dprime_start(g,:)'), mean(dprime_end(g,:)')]
    ste_dprime_startend(g,:) = [std(dprime_start(g,:)')/sqrt(size(dprime_start(g,:)',1)), std(dprime_end(g,:)')/sqrt(size(dprime_end(g,:)',1))]

end

figure
for n = 1:height(dprimes_easy{1,1})
    plot([1 2],dprime_startend{1,n},'Color',[0.5 0.5 0.5 0.4],'linewidth',2,'linestyle',L{g},...
        'Marker','.','MarkerFaceColor',color_eh_mean{2},'MarkerEdgeColor',color_eh_mean{2})
    hold on
end
for n = 1:height(dprimes_easy{2,1})
    plot([3 4],dprime_startend{2,n},'Color',[0.5 0.5 0.5 0.4],'linewidth',2,'linestyle',L{g},...
        'Marker','.','MarkerFaceColor',color_eh_mean{2},'MarkerEdgeColor',color_eh_mean{2})
    hold on
end

errorbar( mean_dprime_startend(1,:), ste_dprime_startend(1,:) ,'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)
errorbar( [nan nan mean_dprime_startend(2,:)], [nan nan ste_dprime_startend(2,:)] ,'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)

ylim([-1 5]);
yticks([-1:1:5])
xlim([0 5])
xticks([1 2 3 4])
yline(1,'--','Color','k');
ylabel("d'")
xticklabels({'adolescent naive','adolescent exp.','adult naive','adult exp.'})
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');


%statistics
h = kstest (dprime_start)
if h > 0
    p = [] ;
    [p(:,3),p,~] = ranksum(dprime_start(1,:),dprime_start(2,:),'alpha',0.05,'tail','both')
    [p(:,4),~,~] = ranksum(dprime_end(1,:),dprime_end(2,:),'alpha',0.05,'tail','both')
    [~,p(:,3:4),~]= bonferroni_holm(p(:,3:4),0.05) ;
    [p(:,1),~,~] = signrank(dprime_start(1,:),dprime_end(1,:),'alpha',0.05,'tail','both')
    [p(:,2),~,~] = signrank(dprime_start(2,:),dprime_end(2,:),'alpha',0.05,'tail','both')
    [~,p(:,1:2),~]= bonferroni_holm(p(:,1:2),0.05) ;

elseif h == 0
    p = [] ;
    [~,p(:,3),~] = ttest2(dprime_start(1,:),dprime_start(2,:),'alpha',0.05,'tail','both')
    [~,p(:,4),~] = ttest2(dprime_end(1,:),dprime_end(2,:),'alpha',0.05,'tail','both')
    [~,p(:,3:4),~]= bonferroni_holm(p(:,3:4),0.05) ;
    [~,p(:,1),~] = ttest(dprime_start(1,:),dprime_end(1,:),'alpha',0.05,'tail','both')
    [~,p(:,2),~] = ttest(dprime_start(2,:),dprime_end(2,:),'alpha',0.05,'tail','both')
    [~,p(:,1:2),~]= bonferroni_holm(p(:,1:2),0.05) ;
end
disp(p)
if  p(:,1) < 0.0005
    txt = '***';
elseif  p(:,1) < 0.005
    txt = '**';
elseif  p(:,1) < 0.05
    txt = '*';
elseif p(:,1)> 0.05
    txt = 'n.s.';
end
line( linspace(1.05,1.95), linspace(4,4),'color','k')
hold on
text(mean([1.5 1.5]), 4.2 ,txt,'FontSize',20)

if  p(:,2) < 0.0005
    txt = '***';
elseif  p(:,2) < 0.005
    txt = '**';
elseif  p(:,2) < 0.05
    txt = '*';
elseif p(:,2) > 0.05
    txt = 'n.s.';
end
line( linspace(3.05,3.95), linspace(4,4),'color','k')
hold on
text(mean([3.5 3.5]), 4.2 ,txt,'FontSize',20)

if  p(:,3) < 0.0005
    txt = '***';
elseif  p(:,3) < 0.005
    txt = '**';
elseif  p(:,3) < 0.05
    txt = '*';
elseif p(:,3) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,3.05), linspace(4.4,4.4),'color','k')
hold on
text(mean([2.05 2.05]), 4.6 ,txt,'FontSize',20)

if  p(:,4) < 0.0005
    txt = '***';
elseif  p(:,4) < 0.005
    txt = '**';
elseif  p(:,4) < 0.05
    txt = '*';
elseif p(:,4) > 0.05
    txt = 'n.s.';
end
line( linspace(1.95,3.95), linspace(4.8,4.8),'color','k')
hold on
text(mean([3 3]), 5 ,txt,'FontSize',20)

%% Behavior during the Recording
close all
clc
%% calculate the dprime across the recording

for g = 1:numel(GNG_rec_all_cell)

    % d' per recording per stimulus difficulty
    for i   = 1:size(GNG_rec_all_cell{1,g},1)

        % get trial stamps in each recording
        stim_types = GNG_rec_all_cell{1,g}(i).Behavior.stim_types;
        stim_ids =  GNG_rec_all_cell{1,g}(i).Behavior.stim_ids;
        trial_responses =  GNG_rec_all_cell{1,g}(i).Behavior.trial_responses;
        lick_times =  GNG_rec_all_cell{1,g}(i).Behavior.lick_times;
        stim_times =  GNG_rec_all_cell{1,g}(i).Behavior.stim_times;

        % retrieve trial identities
        hit_ID{g,i} = find (stim_types == 1 & trial_responses == 1);
        FA_ID{g,i} = find (stim_types == -1 & trial_responses == 1);
        miss_ID{g,i} = find (stim_types == 1 & trial_responses == 0);
        CR_ID{g,i} = find (stim_types == -1 & trial_responses == 0);

        size_trials (g,i)  = size(stim_ids,2) ;

        stim = [1 4 2 3];

        index_stim = find (stim_ids == stim(1) | stim_ids == stim(2));

        [~,~,dprimes_easy_rec(g,i),~,~,~] = GNG_lick_rt_dprime(index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;
        index_stim = find (stim_ids == stim(3) | stim_ids == stim(4));

        [~,~,dprimes_hard_rec(g,i),~,~,~] = GNG_lick_rt_dprime (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

        % % extract mouse and rec ID
        mice(g,i) = GNG_rec_all_cell{1,g}(i).Mouse ;
        recs (g,i) = GNG_rec_all_cell{1,g}(i).Rec ;


        % dprime curve per trial
        for  bin = 10:1:size(trial_responses ,2)

            index_stim = find (stim_ids(1:bin) == stim(1) | stim_ids(1:bin) == stim(2) ...
                | stim_ids(1:bin) == stim(3) | stim_ids(1:bin) == stim(4)); ;

            [~,~,dprimes_cum(g,i,bin),~,~,~] =...
                GNG_lick_rt_dprime ...
                (index_stim, lick_times, stim_times(1:bin),stim_types(1:bin), ...
                trial_responses(1:bin), stim_ids(1:bin), tone_dur , response_window,stim) ;
        end

        % lick raster of recording session
        eventTimes_all = GNG_rec_all_cell{1, g}(i).eventTimes.all;
        eventTimes_licks = GNG_rec_all_cell{1, g}(i).eventTimes.times_vec_licks ;

        % creater raster matrix of all lick times across all trials
        for it = 1:length(eventTimes_all) ;
            licks_after_stim = lick_times (lick_times > stim_times(it)) ;
            licks_per_trial = licks_after_stim(licks_after_stim < stim_times(it)+tone_dur+response_window) ;
            licks_from_stim_onset = licks_per_trial - stim_times(it) ;
            licks_raster{g,i}(it,1:length(licks_from_stim_onset)) = licks_from_stim_onset ;
        end

        eventTimes_all = GNG_rec_all_cell{1, g}(i).eventTimes.all;
        eventTimes_licks = GNG_rec_all_cell{1, g}(i).eventTimes.times_vec_licks ;

        lick_time{g,i} = cell (length(eventTimes_all),1) ;
        for E = 1:length(eventTimes_all)

            idx_lick{E,1} =  find(eventTimes_licks > (eventTimes_all(1,E) + startRange)...
                & eventTimes_licks < (eventTimes_all(1,E) + stopRange)) ;
            if ~isempty(idx_lick{E,1})
                lick_time{g,i}{E,1} =  eventTimes_licks(idx_lick{E,1}) - eventTimes_all(1,E) ;
            end
        end
    end


    % % average dprime per mouse per recording
    for m = unique(mice(g,mice(g,:)>0))
        for r = unique(recs(g,recs(g,:)>0))
            idx =  find (mice(g,:)== m & recs(g,:) == r) ;
            if ~ isempty (idx)
                dprime_mouse(g,m,r) = nanmean([dprimes_easy_rec(g,idx) dprimes_hard_rec(g,idx)]) ;
            end
        end
    end

end
dprime_mouse(dprime_mouse == 0) = nan ;
toc
clc
close all
%% Supplementary Figure 2.1A plot the dprime criterion per recording and stimulus difficulty for experts
clc
close all

dprimes_easy_rec(dprimes_easy_rec == 0)= nan ;
dprimes_hard_rec(dprimes_hard_rec == 0)= nan ;

dprime_adolescent = [dprimes_easy_rec(1,:)', dprimes_hard_rec(1,:)'] ;
dprime_adult = [dprimes_easy_rec(2,:)', dprimes_hard_rec(2,:)'] ;

figure
% adolescent recs
plot([1 2],dprime_adolescent','Color',[0.5 0.5 0.5 0.2],'Marker','o','MarkerFaceColor',color_eh{1},'MarkerEdgeColor',color_eh{1})
hold on
scatter(ones(size(dprimes_hard_rec(1,~isnan(dprimes_hard_rec(1,:))))).*(1 + ones(size(dprimes_hard_rec(1,~isnan(dprimes_hard_rec(1,:)))))),dprimes_hard_rec(1,~isnan(dprimes_hard_rec(1,:)))...
    ,'Marker','o','MarkerFaceColor',color_eh{2},'MarkerEdgeColor',color_eh{2})
hold on
errorbar(nanmean(dprime_adolescent), nanstd(dprime_adolescent)/ sqrt(height(dprime_adolescent))  ,'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)

% adult recs
plot([3 4],dprime_adult','Color',[0.5 0.5 0.5 0.2],'Marker','o','MarkerFaceColor',color_eh{1},'MarkerEdgeColor',color_eh{1})
hold on
scatter(ones(size(dprimes_hard_rec(2,~isnan(dprimes_hard_rec(2,:))))).*(3 + ones(size(dprimes_hard_rec(2,~isnan(dprimes_hard_rec(2,:)))))),dprimes_hard_rec(2,~isnan(dprimes_hard_rec(2,:)))...
    ,'Marker','o','MarkerFaceColor',color_eh{2},'MarkerEdgeColor',color_eh{2})
hold on
errorbar([nan nan nanmean( dprime_adolescent)], [nan nan nanstd(dprime_adolescent)/ sqrt(height(dprime_adolescent))]  ,'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)

xlim([0 5])
xticks([1.5 3.5])
ylim([-1 5])
yticks([-1:2:5])
ylabel("discrimination (d') ")
yline(1,'linestyle','--','linewidth',2)
xticklabels({'adolescent','adult'})
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;


h = kstest (dprime_adolescent)
if h > 0
    p = [] ;
    [p(:,1),~,~] = signrank(dprime_adolescent(:,1),dprime_adolescent(:,2),'alpha',0.05,'tail','both')
    [p(:,2),~,~] = signrank(dprime_adult(:,1)',dprime_adult(:,2)','alpha',0.05,'tail','both')
    [~,p(:,1:2),~]= bonferroni_holm(p(:,1:2),0.05) ;
    [p(:,3),~,~] = ranksum(dprime_adolescent(:,1),dprime_adult(:,1),'alpha',0.05,'tail','both')
    [p(:,4),~,~] = ranksum(dprime_adolescent(:,2),dprime_adult(:,2),'alpha',0.05,'tail','both')
    [~,p(:,3:4),~]= bonferroni_holm(p(:,3:4),0.05) ;

elseif h == 0
    p = [] ;
    [~,p(:,1),~] = ttest(dprime_adolescent(:,1),dprime_adolescent(:,2),'alpha',0.05,'tail','both')
    [~,p(:,2),~] = ttest(dprime_adult(:,1),dprime_adult(:,2),'alpha',0.05,'tail','both')
    [~,p(:,1:2),~]= bonferroni_holm(p(:,1:2),0.05) ;
    [~,p(:,3),~] = ttest2(dprime_adolescent(:,1),dprime_adult(:,1),'alpha',0.05,'tail','both')
    [~,p(:,4),~] = ttest2(dprime_adolescent(:,2),dprime_adult(:,2),'alpha',0.05,'tail','both')
    [~,p(:,3:4),~]= bonferroni_holm(p(:,3:4),0.05) ;
end
disp(p)
if  p(:,1) < 0.0005
    txt = '***';
elseif  p(:,1) < 0.005
    txt = '**';
elseif  p(:,1) < 0.05
    txt = '*';
elseif p(:,1)> 0.05
    txt = 'n.s.';
end
line( linspace(1.05,1.95), linspace(4,4),'color','k')
hold on
text(mean([1.5 1.5]), 4.2 ,txt,'FontSize',20)

if  p(:,2) < 0.0005
    txt = '***';
elseif  p(:,2) < 0.005
    txt = '**';
elseif  p(:,2) < 0.05
    txt = '*';
elseif p(:,2) > 0.05
    txt = 'n.s.';
end
line( linspace(3.05,3.95), linspace(4,4),'color','k')
hold on
text(mean([3.5 3.5]), 4.2 ,txt,'FontSize',20)

if  p(:,3) < 0.0005
    txt = '***';
elseif  p(:,3) < 0.005
    txt = '**';
elseif  p(:,3) < 0.05
    txt = '*';
elseif p(:,3) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,3.05), linspace(4.4,4.4),'color','k')
hold on
text(mean([2.05 2.05]), 4.6 ,txt,'FontSize',20)

if  p(:,4) < 0.0005
    txt = '***';
elseif  p(:,4) < 0.005
    txt = '**';
elseif  p(:,4) < 0.05
    txt = '*';
elseif p(:,4) > 0.05
    txt = 'n.s.';
end
line( linspace(1.95,3.95), linspace(4.8,4.8),'color','k')
hold on
text(mean([3 3]), 5 ,txt,'FontSize',20)

%% Supplementary Figure 6.1 A dprime per recording and stimulus difficulty in naive mice
clc
close all

dprimes_easy_rec(dprimes_easy_rec == 0)= nan ;
dprimes_hard_rec(dprimes_hard_rec == 0)= nan ;

dprime_adolescent = [dprimes_easy_rec(3,:)', dprimes_hard_rec(3,:)'] ;
dprime_adult = [dprimes_easy_rec(4,:)', dprimes_hard_rec(4,:)'] ;

figure
% adolescent recs
plot([1 2],dprime_adolescent','Color',[0.5 0.5 0.5 0.2],'Marker','o','MarkerFaceColor',color_eh{1},'MarkerEdgeColor',color_eh{1})
hold on
scatter(ones(size(dprimes_hard_rec(3,~isnan(dprimes_hard_rec(1,:))))).*(1 + ones(size(dprimes_hard_rec(3,~isnan(dprimes_hard_rec(1,:)))))),dprimes_hard_rec(3,~isnan(dprimes_hard_rec(1,:)))...
    ,'Marker','o','MarkerFaceColor',color_eh{2},'MarkerEdgeColor',color_eh{2})
hold on
errorbar(nanmean(dprime_adolescent), nanstd(dprime_adolescent)/ sqrt(height(dprime_adolescent))  ,'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)

% adult recs
plot([3 4],dprime_adult','Color',[0.5 0.5 0.5 0.2],'Marker','o','MarkerFaceColor',color_eh{1},'MarkerEdgeColor',color_eh{1})
hold on
scatter(ones(size(dprimes_hard_rec(4,~isnan(dprimes_hard_rec(2,:))))).*(3 + ones(size(dprimes_hard_rec(4,~isnan(dprimes_hard_rec(2,:)))))),dprimes_hard_rec(4,~isnan(dprimes_hard_rec(2,:)))...
    ,'Marker','o','MarkerFaceColor',color_eh{2},'MarkerEdgeColor',color_eh{2})
hold on
errorbar([nan nan nanmean( dprime_adolescent)], [nan nan nanstd(dprime_adolescent)/ sqrt(height(dprime_adolescent))]  ,'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)

xlim([0 5])
xticks([1.5 3.5])
ylim([-1 5])
yticks([-1:2:5])
ylabel("discrimination (d') ")
yline(1,'linestyle','--','linewidth',2)
xticklabels({'adolescent','adult'})
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;


h = kstest (dprime_adolescent)
if h > 0
    p = [] ;
    [p(:,1),~,~] = signrank(dprime_adolescent(:,1),dprime_adolescent(:,2),'alpha',0.05,'tail','both')
    [p(:,2),~,~] = signrank(dprime_adult(:,1)',dprime_adult(:,2)','alpha',0.05,'tail','both')
    [~,p(:,1:2),~]= bonferroni_holm(p(:,1:2),0.05) ;
    [p(:,3),~,~] = ranksum(dprime_adolescent(:,1),dprime_adult(:,1),'alpha',0.05,'tail','both')
    [p(:,4),~,~] = ranksum(dprime_adolescent(:,2),dprime_adult(:,2),'alpha',0.05,'tail','both')
    [~,p(:,3:4),~]= bonferroni_holm(p(:,3:4),0.05) ;

elseif h == 0
    p = [] ;
    [~,p(:,1),~] = ttest(dprime_adolescent(:,1),dprime_adolescent(:,2),'alpha',0.05,'tail','both')
    [~,p(:,2),~] = ttest(dprime_adult(:,1),dprime_adult(:,2),'alpha',0.05,'tail','both')
    [~,p(:,1:2),~]= bonferroni_holm(p(:,1:2),0.05) ;
    [~,p(:,3),~] = ttest2(dprime_adolescent(:,1),dprime_adult(:,1),'alpha',0.05,'tail','both')
    [~,p(:,4),~] = ttest2(dprime_adolescent(:,2),dprime_adult(:,2),'alpha',0.05,'tail','both')
    [~,p(:,3:4),~]= bonferroni_holm(p(:,3:4),0.05) ;
end
disp(p)
if  p(:,1) < 0.0005
    txt = '***';
elseif  p(:,1) < 0.005
    txt = '**';
elseif  p(:,1) < 0.05
    txt = '*';
elseif p(:,1)> 0.05
    txt = 'n.s.';
end
line( linspace(1.05,1.95), linspace(4,4),'color','k')
hold on
text(mean([1.5 1.5]), 4.2 ,txt,'FontSize',20)

if  p(:,2) < 0.0005
    txt = '***';
elseif  p(:,2) < 0.005
    txt = '**';
elseif  p(:,2) < 0.05
    txt = '*';
elseif p(:,2) > 0.05
    txt = 'n.s.';
end
line( linspace(3.05,3.95), linspace(4,4),'color','k')
hold on
text(mean([3.5 3.5]), 4.2 ,txt,'FontSize',20)

if  p(:,3) < 0.0005
    txt = '***';
elseif  p(:,3) < 0.005
    txt = '**';
elseif  p(:,3) < 0.05
    txt = '*';
elseif p(:,3) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,3.05), linspace(4.4,4.4),'color','k')
hold on
text(mean([2.05 2.05]), 4.6 ,txt,'FontSize',20)

if  p(:,4) < 0.0005
    txt = '***';
elseif  p(:,4) < 0.005
    txt = '**';
elseif  p(:,4) < 0.05
    txt = '*';
elseif p(:,4) > 0.05
    txt = 'n.s.';
end
line( linspace(1.95,3.95), linspace(4.8,4.8),'color','k')
hold on
text(mean([3 3]), 5 ,txt,'FontSize',20)
%% Supplementary Figure 2.1 b:  expert recordings average d' perfromance per mouse across recordings
close all
clc
for g = 1:numel(GNG_rec_all_cell)-2
    dprime_recs = [] ;
    for m = unique(mice(g,mice(g,:)>0))
        figure(g)
        plot(1:length(reshape(dprime_mouse(g,m,:),[],1))',reshape(dprime_mouse(g,m,:),[],1)',...
            'Color',color_mean_marker{1},'Marker','o','MarkerFaceColor',[0.3 0.3 0.3],'MarkerEdgeColor',[0.3 0.3 0.3])
        hold on
        xlim([0 4])
        xticks([])
        ylim([-1 5])
        yticks([-1:2:5])
        ylabel("discrimination (d')")
        xticks([0:1:4])
        xticklabels({'','1st rec.','2nd rec.','3rd rec.',''})
        yline(1,'linestyle','--','linewidth',2)

        dprime_recs(m,:)= dprime_mouse(g,m,:) ;
        ax = gca;
        ax.XAxis.FontSize =20;
        ax.YAxis.FontSize = 20;
        movegui('east');
        box off;

    end

    [p(:,1),~,~] = signrank(dprime_recs(:,1),dprime_recs(:,2),'alpha',0.05,'tail','both') ;
    [p(:,2),~,~] = signrank(dprime_recs(:,1)',dprime_recs(:,3)','alpha',0.05,'tail','both') ;
    [p(:,3),~,~] = signrank(dprime_recs(:,2),dprime_recs(:,3),'alpha',0.05,'tail','both') ;
    [~,p,~]= bonferroni_holm(p,0.05) ;
    disp(p)
    if  p(:,1) < 0.0005
        txt = '***';
    elseif  p(:,1) < 0.005
        txt = '**';
    elseif  p(:,1) < 0.05
        txt = '*';
    elseif p(:,1)> 0.05
        txt = 'n.s.';
    end
    line( linspace(1.05,1.95), linspace(4.4,4.4),'color','k')
    hold on
    text(mean([1.5 1.5]), 4.6 ,txt,'FontSize',20)

    if  p(:,2) < 0.0005
        txt = '***';
    elseif  p(:,2) < 0.005
        txt = '**';
    elseif  p(:,2) < 0.05
        txt = '*';
    elseif p(:,2) > 0.05
        txt = 'n.s.';
    end
    line( linspace(1.05,2.95), linspace(4,4),'color','k')
    hold on
    text(mean([2 2]), 4.2 ,txt,'FontSize',20)

    if  p(:,3) < 0.0005
        txt = '***';
    elseif  p(:,3) < 0.005
        txt = '**';
    elseif  p(:,3) < 0.05
        txt = '*';
    elseif p(:,3) > 0.05
        txt = 'n.s.';
    end
    line( linspace(2.05,2.95), linspace(4.4,4.4),'color','k')
    hold on
    text(mean([2.5 2.5]), 4.6 ,txt,'FontSize',20)
end
%% supplementary Figure 6.1 B:   naive recordings average d' perfromance per mouse across recordings
close all
clc
figure
dprime_recs = [] ;
% adolescent
for m = unique(mice(3,mice(3,:)>0))
    d_naive = reshape(dprime_mouse(3,m,:),[],1)' ;
    plot([1 2],d_naive(~isnan(d_naive)),...
        'Color',color_mean_marker{1},'linestyle','-','Marker','o','MarkerFaceColor',color_eh_all{1},'MarkerEdgeColor',color_eh_all{1})
    hold on

    dprime_recs(m,:)= dprime_mouse(3,m,:) ;
    p = [] ;
    [p(:,1),~,~] = signrank(dprime_recs(:,1),dprime_recs(:,2),'alpha',0.05,'tail','both') ;

    disp(p)
    if  p(:,1) < 0.0005
        txt = '***';
    elseif  p(:,1) < 0.005
        txt = '**';
    elseif  p(:,1) < 0.05
        txt = '*';
    elseif p(:,1)> 0.05
        txt = 'n.s.';
    end
    line( linspace(1.05,1.95), linspace(1.5,1.5),'color','k')
    hold on
    text(mean([1.5 1.5]), 1.7 ,txt,'FontSize',20)
end

% adult
for m = unique(mice(4,mice(4,:)>0))
    d_naive = reshape(dprime_mouse(4,m,:),[],1)' ;

    plot([3 4 ],d_naive(~isnan(d_naive)),...
        'Color',color_mean_marker{1},'linestyle','-','Marker','o','MarkerFaceColor',color_eh_all{1},'MarkerEdgeColor',color_eh_all{1})
    hold on
    dprime_recs(m,:)= dprime_mouse(4,m,:) ;


    % adult
    [p(:,1),~,~] = signrank(dprime_recs(:,1),dprime_recs(:,2),'alpha',0.05,'tail','both') ;
    disp(p)
    if  p(:,1) < 0.0005
        txt = '***';
    elseif  p(:,1) < 0.005
        txt = '**';
    elseif  p(:,1) < 0.05
        txt = '*';
    elseif p(:,1)> 0.05
        txt = 'n.s.';
    end
    line( linspace(3.05,3.95), linspace(1.5,1.5),'color','k')
    hold on
    text(mean([3.5 3.5]), 1.7 ,txt,'FontSize',20)
end
hold on
xlim([0 5])
xticks([])
ylim([-1 2])
ylabel("discrimination (d')")
xticks([0:1:4])
xticklabels({'', '1st rec adolescent','2nd rec. adolescent','1st rec adult','2nd rec. adult','',})
yline(1,'linestyle','--','linewidth',2)

dprime_recs(m,:)= dprime_mouse(3,m,:) ;
ax = gca;
ax.XAxis.FontSize =20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;


%% supplementary Figure 2.2:  Comparison between Educage 'experts' and head-fixed 'experts'

% select behavior file from Educage
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\Educage'...
    , ['Select behavior file: dprime_edu ']);
addpath(path)
load (file)

for g = 1:numel(GNG_rec_all_cell)-2 % expert recordings only
    % concatenate the dprime per age-group
    all_dprimes = nan (15,4) ;
    all_dprimes (1:size(dprime_edu.dprime_easy(dprime_edu.dprime_easy(:,g) >0),1),1) = dprime_edu.dprime_easy(dprime_edu.dprime_easy(:,g) >0 ,g) ;
    all_dprimes(1:size( dprimes_easy_rec(g, dprimes_easy_rec(g,:)>0),2),2) =  dprimes_easy_rec(g, dprimes_easy_rec(g,:)>0)' ;
    all_dprimes (1:size(dprime_edu.dprime_hard(dprime_edu.dprime_hard(:,g) >0),1),3) = dprime_edu.dprime_hard(dprime_edu.dprime_hard(:,g) >0 ,g) ;
    all_dprimes(1:size( dprimes_hard_rec(g, dprimes_hard_rec(g,:)>0),2),4) =  dprimes_hard_rec(g, dprimes_hard_rec(g,:)>0)' ;

    figure
    violinplot( [all_dprimes], {'Educage 1 oct.','set 1 oct.','Educage 0.25 oct.','set 0.25 oct.'},"ViolinColor",...
        {[0.3010 0.7450 0.9330; .3010 0.7450 0.9330; 0 0.4470 0.7410; 0 0.4470 0.7410]})
    xlim([0 5])
    ylim([-1 4.5])
    yticks([-1:1:4])
    ylabel("d'")
    ax = gca;
    ax.XAxis.FontSize =20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;


    hold on;
    h = kstest (all_dprimes) ;
    if h == 1
        [p(1),~,~] = ranksum(all_dprimes(:,1),all_dprimes(:,2),'alpha',0.05,'tail','both') ;
        [p(2),~,~] = ranksum(all_dprimes(:,3),all_dprimes(:,4),'alpha',0.05,'tail','both') ;
        [~,p,~]= bonferroni_holm(p,0.05) ;

    elseif h == 0
        [~,p(1),~] = ttest2(all_dprimes(:,1),all_dprimes(:,2),'alpha',0.05,'tail','both') ;
        [~,p(2),~] = ttest2(all_dprimes(:,3),all_dprimes(:,4),'alpha',0.05,'tail','both') ;
        [~,p,~]= bonferroni_holm(p,0.05) ;

    end
    disp(p)
    if  p(1) < 0.0005
        txt = '***';
    elseif  p(1) < 0.005
        txt = '**';
    elseif  p(1) < 0.05
        txt = '*';
    elseif p(1) > 0.05
        txt = 'n.s.';
    end

    line( linspace(1.1,1.9), linspace(4.1,4.1),'color','k')
    hold on
    text(mean([1.5 1.5]), 4.3 ,txt,'FontSize',15)

    if  p(2) < 0.0005
        txt = '***';
    elseif  p(2) < 0.005
        txt = '**';
    elseif  p(2) < 0.05
        txt = '*';
    elseif p(2) > 0.05
        txt = 'n.s.';
    end

    line( linspace(3.1,3.9), linspace(4.1,4.1),'color','k')
    hold on
    text(mean([3.5 3.5]), 4.3 ,txt,'FontSize',15)
end
%% stats of trial size during set recordings
for g = 1:numel(GNG_rec_all_cell)
    for i   = 1:size(GNG_rec_all_cell{1,g},1)
        stim_ids =  GNG_rec_all_cell{1,g}(i).Behavior.stim_ids;
        size_trials (g,i)  = size(stim_ids,2) ;
    end
end

% size of the session
size_trials (size_trials == 0) = nan ;
for i = 1:2
    h(i) = kstest(size_trials(1,:)) ;
end
% test according to the distribution of h
if sum(h) > 0
    [p,~,~] = ranksum (size_trials(1,:),size_trials(2,:),'alpha',0.05,'tail','both') ;

elseif sum(h) == 0
    [~,p,~] =  ttest2 (size_trials(1,:),size_trials(2,:),'alpha',0.05,'tail','both') ;
end
disp(p)
disp(nanmean(size_trials(1,:)))
disp(nanstd(size_trials(1,:))/sqrt(size(size_trials(1,:),2)))
disp(nanmean(size_trials(2,:)))
disp(nanstd(size_trials(2,:))/sqrt(size(size_trials(2,:),2)))

%% Figure 2C: lick raster of one example recording
close all
tic
% choose example mouse
g = 2 % adult
i = 8 % recording

size_E =  size(GNG_rec_all_cell{1,g}(i).eventTimes.all,2) ;
licks_raster{g,i} (licks_raster{g,i} == 0) = nan ;
figure
for it = 1:size_E
    if ismember (it, hit_ID{g,i})
        rectangle( 'Position' , [0 it tone_dur it+1],'Facecolor',[0.4660 0.6740 0.1880],'EdgeColor','none')
        hold on
        scatter(licks_raster{g,i}(it,:) , it + 0.5, 'MarkerFaceColor' ,[0.5 0.5 0.5],'Marker','.','MarkerEdgeColor',[0.5 0.5 0.5])
        hold on
        %   scatter (reward_raster{g,i}(it,:) + 0.02, it + 0.5, 'MarkerFaceColor',[0 0 0],'Marker','o','MarkerEdgeColor','none')
        if  ~isempty(lick_time{g,i}{it,1})
            scatter(lick_time{g,i}{it,1}(:,:), it + 0.5, 'MarkerFaceColor' ,[0.5 0.5 0.5],'Marker','.','MarkerEdgeColor',[0.5 0.5 0.5])
        end
    elseif ismember (it, FA_ID{g,i})
        rectangle( 'Position' , [0 it tone_dur it+1],'Facecolor',[0.9290 0.6940 0.1250],'EdgeColor','none')
        scatter(licks_raster{g,i}(it,:) , it + 0.5 , 'MarkerFaceColor',[0.5 0.5 0.5],'Marker','.','MarkerEdgeColor',[0.5 0.5 0.5])
        % scatter(pun_raster{g,i}(it,:)+ 0.02 , it + 0.5 ,'MarkerFaceColor',[0 0 0],'Marker','o','MarkerEdgeColor','none')

        if  ~isempty(lick_time{g,i}{it,1})
            scatter(lick_time{g,i}{it,1}(:,:), it + 0.5, 'MarkerFaceColor' ,[0.5 0.5 0.5],'Marker','.','MarkerEdgeColor',[0.5 0.5 0.5])
        end

    elseif ismember (it, miss_ID{g,i})
        rectangle( 'Position' , [0 it tone_dur it+1],'Facecolor',[0.8500 0.3250 0.0980],'EdgeColor','none')
    elseif ismember (it, CR_ID{g,i})
        rectangle( 'Position' , [0 it tone_dur it+1],'Facecolor',[0 0.4470 0.7410],'EdgeColor','none')
    end
end
rectangle( 'Position' , [0 (size_E+1) tone_dur (size_E+1)+1],'Facecolor',[1 1 1],'EdgeColor','none')

clear xlim
ylim ([0  1000])
xlim ([0 tone_dur+response_window])
xticks ([0 0.1 0.6 2])
yticks ([0:250:1000])
xlabel ('time(sec)')
ylabel ('trials')
xline(0.6,'-.k')

ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
movegui('east');

%% Figure 2F: dprime trajectory for all recordings

dprimes_cum(dprimes_cum == 0) = nan ;
close all

for g = 1:numel(GNG_rec_all_cell)-2 % only for experts
    figure
    for i   = 1:numel(GNG_rec_all_cell{1,g})
        plot_cum_dprime =  squeeze(dprimes_cum (g,i,:))' ;
        plot ( smoothdata(plot_cum_dprime),'Color',[0.3 0.3 0.3],...
            'linestyle',L2{g},'linewidth',2) ;
        hold on
    end

    xlim([ 0 850])
    xticks([1 150 850])
    ylim([0 4])
    yticks([0:1:4])
    yline(1,'--k')
    ylabel("d'")
    xlabel('trials ')
    ax = gca;
    ax.XAxis.FontSize =20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;
    ax = gca;
    ax.XAxis.FontSize =20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;

end
%% Analyze the minimal number of shared trials between all mice
size_trials = nan (numel(GNG_rec_all_cell),numel(GNG_rec_all_cell{1,2})) ;
% determine the minimal number of shared trials between all recs in both groups
for g = 1:numel(GNG_rec_all_cell)
    Recs = 1:numel(GNG_rec_all_cell{1,g});
    for i = 1:length(Recs)
        size_trials (g,i)  = size(GNG_rec_all_cell{1, g}(i).Behavior.stim_ids,2) ;
    end
    min_g(g,1) = min(size_trials(g,:)) ;
end
min_trials = min(min_g) ;

% dprime, lick psth and lick parameters
tic
for g = 1:numel(GNG_rec_all_cell)
    Recs = 1:numel(GNG_rec_all_cell{1,g});
    for i = 1:length(Recs)
        if ~isempty(GNG_rec_all_cell{1,g})

            % extract trial stamps
            stim_types = GNG_rec_all_cell{1, g}(i).Behavior.stim_types(:,1:min_trials) ;
            stim_ids  = GNG_rec_all_cell{1, g}(i).Behavior.stim_ids (:,1:min_trials) ;
            trial_responses = GNG_rec_all_cell{1, g}(i).Behavior.trial_responses(:,1:min_trials) ;
            lick_times = GNG_rec_all_cell{1, g}(i).Behavior.lick_times ;
            stim_times = GNG_rec_all_cell{1, g}(i).Behavior.stim_times(:,1:min_trials) ;
            reward_times =GNG_rec_all_cell{1, g}(i).Behavior.reward_times  ;
            punishment_times = GNG_rec_all_cell{1, g}(i).Behavior.punishment_times ;

            % dprime curve per trial
            for  bin = 10:1:size(trial_responses ,2)
                index_stim = find (stim_ids(1:bin) == stim(1) | stim_ids(1:bin) == stim(2) ...
                    | stim_ids(1:bin) == stim(3) | stim_ids(1:bin) == stim(4)); ;

                [~,~,dprimes_cum(g,i,bin),~,~,~,cbias_cum(g,i,bin)] = GNG_lick_rt_dprime(index_stim, lick_times, stim_times(1:bin),stim_types(1:bin), ...
                    trial_responses(1:bin), stim_ids(1:bin), tone_dur , response_window,stim) ;
            end
            binArray = zeros(length(stim_times), numBins) ;

            % create a binary array of the licks during the reinforcment delay
            for r = 1:length(stim_times)
                [n,binCenters] = histdiff(lick_times, stim_times(r), binBorders);
                binArray(r,:) = n;
            end
            binArray_lick_all{g,i} = binArray ;

            % Lick Latency, lick iti and lick accumulation / expectancy
            lick_latency_trial = nan(size(binArray_lick_all{g,i},1),1) ;
            lick_expec = nan(size(binArray_lick_all{g,i},1),1) ;

            for t = 1:size(binArray_lick_all{g,i},1)
                lick_times_trial = find (binArray_lick_all{g,i}(t,:) ~= 0) ;

                lick_latency_trial(t,:) = NaN ;
                lick_expec(t,:) = NaN ;


                if  isempty(lick_times_trial) ;
                    lick_latency_trial(t,:) = NaN ;
                    lick_expec(t,:) = NaN ;

                elseif ~isempty(lick_times_trial) ;
                    lick_times_trial = lick_times_trial (lick_times_trial > length_base) ;

                    if ~isempty(lick_times_trial) ;
                        lick_latency_1 = lick_times_trial(1) ;
                        lick_latency_trial(t,:) = (lick_latency_1) - length_base  ;


                        if  lick_latency_1 < reinforcement
                            lick_expec(t,:) = length(find (binArray_lick_all{g,i}(t,lick_latency_1:(lick_latency_1 + dist(lick_latency_1,reinforcement))) ~= 0)) ;
                        elseif lick_latency_1 >= reinforcement
                            lick_expec(t,:) = NaN ;
                        end

                    elseif isempty(lick_times_trial) ;
                        lick_latency_trial(t,:) = NaN ;
                        lick_expec(t,:) = NaN ;
                    end
                end
            end

            lick_latency_diff{g,i} = lick_latency_trial ;
            lick_expec_diff{g,i} = lick_expec ;
            rec_length(g,i) = length(stim_times) ;
            stim_ids_all{g,i} = stim_ids ;
        end
    end

    lick_latency_all = nan(2, length(Recs), min_trials);
    lick_expec_all =  nan(2, length(Recs), min_trials);

    for i   = 1:numel(GNG_rec_all_cell{1,g})
        lick_latency_diff{g,i}(lick_latency_diff{g,i} == 0) = nan;
        lick_expec_diff{g,i}(lick_expec_diff{g,i} == 0) = nan;

        lick_latency_all(g,i,1:length(lick_latency_diff{g,i})) = lick_latency_diff{g,i};
        lick_expec_all(g,i,1:length(lick_expec_diff{g,i})) = lick_expec_diff{g,i};

        lick_latency_m(g,i) = nanmean(squeeze(lick_latency_all(g,i,:))) ;
        lick_expec_m(g,i) = nanmean(squeeze(lick_expec_all(g,i,:))) ;
    end

    for t = 1:min_trials
        lick_latency_group (g,t) = nanmean(lick_latency_all(g,:,t)) ;
        lick_expec_group (g,t) = nanmean(lick_expec_all(g,:,t)) ;
    end

    lick_latency_group(g,lick_latency_group(g,:) == 0) = nan;
    lick_expec_group(g,lick_expec_group(g,:) == 0) = nan;

    idx_nan = [] ;
    idx_nan = find(~isnan(lick_latency_group(g,:))) ;

    idx_nan = find(~isnan(lick_latency_group(g,:))) ;
    lick_latency_cell{g,1}  =  lick_latency_group(g,idx_nan) ;
    lick_latency_cut{g,1}(:,1:min_trials) = lick_latency_cell{g,1}(:,1:min_trials) ;

    lick_expec_cell{g,1}  =  lick_expec_group(g,idx_nan) ;
    lick_expec_cut{g,1}(:,1:min_trials) = lick_expec_cell{g,1}(:,1:min_trials) ;

end
toc
%% Figure 2 G: population dprime trajectory for minimal number of shared trials
dprime_cum_all = nan(2,min_trials)

figure
for g = 1:2
    for  bin = 10:1:min_trials  ;

        mean_dprime_cum_all (g,bin) =   nanmean(squeeze(dprimes_cum(g,:,bin))) ;
        stde_dprime_cum_all (g,bin) =   nanstd(squeeze(dprimes_cum(g,:,bin)))./ sqrt(size(dprimes_cum,2)) ;
    end
end

for g = 1:2

    plot( mean_dprime_cum_all (g,1:min_trials),'Color',[0.3 0.3 0.3],...
        'linestyle',L2{g},'linewidth',2)
    hold on

    patch([[1,1:min_trials-1] flip([1,1:min_trials-1])] , [ mean_dprime_cum_all(g,1:min_trials) + ....
        stde_dprime_cum_all(g,1:min_trials) flip( mean_dprime_cum_all(g,1:min_trials)...
        -  stde_dprime_cum_all(g,1:min_trials))], [0.3 0.3 0.3],...
        'facealpha' , 0.2, 'EdgeColor','none')
    xlim([10 150])
    ylim([0 2.5])
    yticks([0:0.5:2.5])
    ylabel("d'")
    xlabel('n trials')
    ax = gca;
    ax.XAxis.FontSize =20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;
end


h = kstest( mean_dprime_cum_all )
if h == 1
    h = []
    for bin = 10:size(mean_dprime_cum_all,2)
        [p(bin,:),~,~] = ranksum(squeeze(dprimes_cum(1,:,bin)),...
            squeeze(dprimes_cum(2,:,bin)),'alpha',0.05,'tail','both')
    end

elseif h == 0
    h = []
    for bin = 10:size(mean_dprime_cum_all,2)
        [~,p(bin,:),~] = ttest2(squeeze(dprimes_cum(1,:,bin)),...
            squeeze(dprimes_cum(2,:,bin)),'alpha',0.05,'tail','both')
    end
end
idx_sig  = find (p < 0.05) ;
disp(idx_sig(end)+10)
p = p(~isnan(p))
p = - p ; % change to minus for visualization of p-value
figure
imagesc(p')
xticks([])
yticks([])
h = colorbar
h.FontSize = 20
ax = gca;
ax.XAxis.FontSize =20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;


%%  Figure 2I: cumulative lick curve adults and adolescents
clc
close all
for g = 1:numel(GNG_rec_all_cell)
    Recs = 1:numel(GNG_rec_all_cell{1,g});
    for i = 1:length(Recs)
        for t = 1:size(binArray_lick_all{g,i},1)
            binArray_cum(t,:) = cumsum(binArray_lick_all{g,i}(t,:)) ;
        end
        idx_0 =  find(binArray_cum(:,end) >0) ;
        binArray_cum = binArray_cum(idx_0,:) ;

        cum_t(g,i,:) = mean(binArray_cum) ;
        binArray = binArray_lick_all{g,i}(idx_0,:) ;

    end

    stderr_cum(g,:) = (nanstd(cum_t(g,:,:))./ (sqrt (size(cum_t(g,:,:),2)))) ;
    mean_cum(g,:) = nanmean(cum_t(g,:,:)) ;

end

for g = 1:numel(GNG_rec_all_cell)-2
    plot([1,1:length_stim-1], mean_cum(g,:),'Color',[0.3 0.3 0.3],...
        'linestyle',L2{g},'linewidth',2);
    hold on

    patch([[1,1:length_stim-1] flip([1,1:length_stim-1])] , [mean_cum(g,:) + ....
        stderr_cum(g,:) flip(mean_cum(g,:) - stderr_cum(g,:))],[0.3 0.3 0.3] ,...
        'facealpha' , 0.2, 'EdgeColor','none')
    hold on

    xticks([100 350 600 2100])
    xticklabels([0 250 500 2000])
    xlim([0 700])
    yline(1,'--','Color','k','linewidth',2);

    xline(100,'--','Color','k','linewidth',2);
    xline(200,'--','Color','k','linewidth',2);
    ylim([0 4])
    yticks([1 2 3 4])
    xlabel('time(sec)')
    ylabel('cumulative licks')
    box off;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    box off;
    hold on;
end
%% LME stats lick latency, lick iti and lick accumulation

lick_latency_m(lick_latency_m == 0) = nan  ;
lick_expec_m(lick_expec_m == 0) = nan ;

age_r = nan(size(vertcat(lick_latency_diff{1:2,:}),1),1)
mice_r = nan(size(vertcat(lick_latency_diff{1:2,:}),1),1)
recs_r = nan(size(vertcat(lick_latency_diff{1:2,:}),1),1)

age_r(1:size(vertcat(lick_latency_diff{1,:}),1),1) = 1 ;
age_r(size(vertcat(lick_latency_diff{1,:}),1)+1:size(vertcat(lick_latency_diff{1,:}),1) + size(vertcat(lick_latency_diff{2,:}),1),1) = 2 ;

for g = 1:numel(GNG_rec_all_cell)-2 % only experts
    Recs = 1:numel(GNG_rec_all_cell{1,g});

    for i = 1:length(Recs)
        if i == 1
            mice_r (1:size(vertcat(lick_latency_diff{1,:}),1),1) = mice(g,i) ;
            recs_r (1:size(vertcat(lick_latency_diff{1,:}),1),1) = recs(g,i) ;
            %d_r (1:min_trials,1) = squeeze(dprimes_cum (g,i,1:min_trials)) ;

        elseif i > 1
            mice_r (size(vertcat(lick_latency_diff{g,1:i-1}),1)+1:size(vertcat(lick_latency_diff{g,1:i}),1),1) = mice(g,i) ;
            recs_r (size(vertcat(lick_latency_diff{g,1:i-1}),1)+1:size(vertcat(lick_latency_diff{g,1:i}),1),1) = recs(g,i) ;
            % d_r ((min_trials*(i-1))+1:min_trials*i,1) = squeeze(dprimes_cum (g,i,1:min_trials)) ;
        end
    end
    % if g == 1
    %   dr_all(1:size(d_r,1),1) = d_r ;
    % elseif g== 2
    %   dr_all(size(dr_all,1)+ 1:size(dr_all,1) + size(d_r,1),1) = d_r ;
    % end
end

LME_licks = table;
LME_licks.Group  = age_r  ;%age
LME_licks.MouseID =  mice_r ;% mouse
LME_licks.RecordingID  = recs_r ;% rec
LME_licks.TrialID = (1:size(vertcat(lick_latency_diff{1:2,:}),1))' ;% trial ID
LME_licks.LickLatency = vertcat(lick_latency_diff{1:2,:}) ;
LME_licks.LickCount = vertcat(lick_expec_diff{1:2,:}) ;
% Convert categorical variables
LME_licks.Group = categorical(LME_licks.Group);
LME_licks.MouseID = categorical(LME_licks.MouseID);
LME_licks.RecordingID = categorical(LME_licks.RecordingID);


% Fit mixed-effects model for Lick Count
model_LickCount = fitlme(LME_licks, 'LickCount ~ Group * LickLatency + (1|MouseID) + (1|MouseID:RecordingID)');
disp(model_LickCount);

% Perform ANOVA on the model to get p-values for fixed effects
anovaTable = anova(model_LickCount);

% Display the ANOVA table
disp(anovaTable);

% Extract p-values for each fixed effect
pValues = anovaTable.pValue;

% Display p-values for each fixed effect
disp('P-values:');
disp(['Group (Adult): ', num2str(pValues(2))]);
disp(['Lick Latency: ', num2str(pValues(3))]);
disp(['Group (Adult): Lick Latency: ', num2str(pValues(4))]);


%% table 2
% Extracting necessary information from the LinearMixedModel object
fixedEffects = model_LickCount.Coefficients;

% Create a table for fixed effects
fixedEffectsTable = table(round(fixedEffects.Estimate,3), round(fixedEffects.SE,3), round(fixedEffects.tStat,3), fixedEffects.pValue, ...
    'VariableNames', {'Estimate', 'StandardError', 'tStatistic', 'pValue'}, ...
    'RowNames', fixedEffects.Name);

writetable(fixedEffectsTable,'LME_lick_behavior.xlsx','Sheet',1)

%%  Figure 2J display lick lantecy

lick_latency_adol = nan (size(vertcat(lick_latency_diff{1:2,:}),1),1)
lick_latency_adult = nan (size(vertcat(lick_latency_diff{1:2,:}),1),1)

lick_latency_adol(1:size(vertcat(lick_latency_diff{1,:}),1),1) = vertcat(lick_latency_diff{1,:}) ;
lick_latency_adult(1:size(vertcat(lick_latency_diff{2,:}),1),1) = vertcat(lick_latency_diff{2,:}) ;
figure
violinplot([lick_latency_adol lick_latency_adult],{'adolescent','adult'},"ViolinColor",{[0.5 0.5 0.5],[0.5 0.5 0.5]})

xlim([0 3])
ylim([1 2000])
yticks([0:500:2000])
ylabel('time(ms)')

disp(model_LickCount);

ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;
%% Figure 2K display lick expectation
lick_expec_adol = nan (size(vertcat(lick_expec_diff{1:2,:}),1),1)
lick_expec_adult = nan (size(vertcat(lick_expec_diff{1:2,:}),1),1)

disp(nanmean(vertcat(lick_expec_diff{1,:}),1)) % adolescents
disp(nanmean(vertcat(lick_expec_diff{2,:}),1)) % adolescents
disp(nanstd(vertcat(lick_expec_diff{1,:}),1)./ sqrt(size(vertcat(lick_expec_diff{1,:}),1))) % adolescents
disp(nanstd(vertcat(lick_expec_diff{2,:}),1)./sqrt(size(vertcat(lick_expec_diff{2,:}),1)))  % adolescents

lick_expec_adol(1:size(vertcat(lick_expec_diff{1,:}),1),1) = vertcat(lick_expec_diff{1,:}) ;
lick_expec_adult(1:size(vertcat(lick_expec_diff{2,:}),1),1) = vertcat(lick_expec_diff{2,:}) ;
figure
violinplot([lick_expec_adol lick_expec_adult],{'adolescent','adult'},"ViolinColor",{[0.5 0.5 0.5],[0.5 0.5 0.5]})
hold on
xlim([0 3])
ylim([0 10])
yticks([0:2:10])
ylabel('licks')

disp(model_LickCount);

ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;
