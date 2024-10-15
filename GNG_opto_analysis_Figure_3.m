clear all
close all
clc
%% load the the opto file
%add the Praegel_et_al_MATLABR2023b_scripts directory here 

addpath(genpath('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_MATLABR2023b_scripts'))

[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , ['Select opto sessions: GNG_opto_cell ']);
addpath(path)
 load (file)

%% Parameters
tone_dur = 0.1 ; % sec
response_window = 2 ; % sec
stim = [1:16] ; % ID
freqs = [7.07 9.17 10.95 14.14] % frequencies of stim IDs

% 
color_eh= {[0.3010 0.7450 0.9330], [0 0.4470 0.7410]} ;
color_opto = {[0.6350 0.0780 0.1840]} ;
color_opto_fade = {[0.6350 0.0780 0.1840 0.5]} ;
color_mean_marker = {[0.5 0.5 0.5 0.5]} ;
L = {['-'], ['-']} ;

%% analyse dprime lick rates and reaction times 

% g = 1 GtACR2
% g= 2 CAMKII
% preallocate the variables
go_licks_all_easy = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;
ngo_licks_all_easy = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;
dprimes_easy = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;

go_licks_all_easy_opto =nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;
ngo_licks_all_easy_opto =nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;
dprimes_easy_opto = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;

go_licks_all_hard = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;
ngo_licks_all_hard = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;
dprimes_hard = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;

go_licks_all_hard_opto = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;
ngo_licks_all_hard_opto = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;
dprimes_hard_opto = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;


mouse_ID = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;
session_ID = nan(size(GNG_opto_all_cell,1),size(GNG_opto_all_cell,2),5) ;

for g = 1:size(GNG_opto_all_cell,1)
    for i   = 1:size(GNG_opto_all_cell,2) % mouse
        for s = 1:numel(GNG_opto_all_cell{g,i}) % session

        stim_types = GNG_opto_all_cell{g,i}(s).Behavior.stim_types;
        stim_ids =GNG_opto_all_cell{g,i}(s).Behavior.stim_ids;
        trial_responses = GNG_opto_all_cell{g,i}(s).Behavior.trial_responses;
        lick_times = GNG_opto_all_cell{g,i}(s).Behavior.lick_times;
        stim_times = GNG_opto_all_cell{g,i}(s).Behavior.stim_times;

        %% easy
        index_stim = find (stim_ids == stim(1) | stim_ids == stim(2));
        [go_licks_all_easy(g,i,s),ngo_licks_all_easy(g,i,s),dprimes_easy(g,i,s),~,~,~] =....
            GNG_lick_rt_dprime ...
            (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;
        
        %% hard
        index_stim = find (stim_ids == stim(3) | stim_ids == stim(4));
        [go_licks_all_hard(g,i,s),ngo_licks_all_hard(g,i,s),dprimes_hard(g,i,s),~,~,~] =...
            GNG_lick_rt_dprime ...
            (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;
        
        %% easy opto
             index_stim = find (stim_ids == stim(8) | stim_ids == stim(9));
        [go_licks_all_easy_opto(g,i,s),ngo_licks_all_easy_opto(g,i,s),dprimes_easy_opto(g,i,s),~,~,~] =....
            GNG_lick_rt_dprime ...
            (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;
        
        %% hard opto
        index_stim = find (stim_ids == stim(10) | stim_ids == stim(11));
        [go_licks_all_hard_opto(g,i,s),ngo_licks_all_hard_opto(g,i,s),dprimes_hard_opto(g,i,s),~,~,~] =...
            GNG_lick_rt_dprime ...
            (index_stim, lick_times, stim_times,stim_types, trial_responses, stim_ids, tone_dur , response_window,stim) ;

        mouse_ID(g,i,s) = i ;
        session_ID(g,i,s) = s ;

        end
    end
end

 %% Supplementary Table 1 LME model of optogenetic effect on performance
clc
%close all

dprimes_easy_g = reshape(squeeze(dprimes_easy(1,:,:)),[],1) ;
dprimes_easy_g = dprimes_easy_g(~isnan(dprimes_easy_g)); 

dprimes_easy_opto_g = reshape(squeeze(dprimes_easy_opto(1,:,:)),[],1) ;
dprimes_easy_opto_g = dprimes_easy_opto_g(~isnan(dprimes_easy_opto_g)); 

dprimes_hard_g = reshape(squeeze(dprimes_hard(1,:,:)),[],1) ;
dprimes_hard_g = dprimes_hard_g(~isnan(dprimes_hard_g)); 

dprimes_hard_opto_g = reshape(squeeze(dprimes_hard_opto(1,:,:)),[],1) ;
dprimes_hard_opto_g = dprimes_hard_opto_g(~isnan(dprimes_hard_opto_g)); 

mouse_ID_g = reshape(mouse_ID(1,:,:),[],1) ;
mouse_ID_g = mouse_ID_g(~isnan(mouse_ID_g)) ;

session_ID_g = reshape(session_ID(1,:,:),[],1) ;
session_ID_g = session_ID_g(~isnan(session_ID_g)) ;

easy = nan(length(dprimes_easy_g),1);
easy(1:end,1) = 1 ;
hard = nan(length(dprimes_easy_g),1);
hard(1:end,1) = 2 ;
light_on = nan(length(dprimes_easy_g),1);
light_on(1:end,1) = 1 ;
light_off =  nan(length(dprimes_easy_g),1);
light_off(1:end,1) = 2 ;


LME_opto = table;
LME_opto.MouseID = reshape([mouse_ID_g mouse_ID_g mouse_ID_g mouse_ID_g],[],1) ;
LME_opto.sessionID = reshape([session_ID_g session_ID_g session_ID_g session_ID_g],[],1) ;
LME_opto.dprime = reshape([dprimes_easy_g dprimes_easy_opto_g dprimes_hard_g dprimes_hard_opto_g],[],1) ;
LME_opto.difficulty = reshape([easy easy hard hard],[],1) ;
LME_opto.light = reshape([light_on light_off light_on light_off],[],1) ;


% Convert categorical variables
LME_opto.MouseID = categorical(LME_opto.MouseID);
LME_opto.sessionID = categorical(LME_opto.sessionID);
LME_opto.difficulty = categorical(LME_opto.difficulty);
LME_opto.difficulty = categorical(LME_opto.difficulty);

% Fit mixed-effects model for Lick Latency
model_opto = fitlme(LME_opto, 'dprime ~ light * difficulty + (1|MouseID) + (1|sessionID)')
disp(model_opto);
%% supplementary Table 1
% Extracting necessary information from the LinearMixedModel object
fixedEffects = model_opto.Coefficients;

% Create a table for fixed effects
fixedEffectsTable = table(round(fixedEffects.Estimate,3), round(fixedEffects.SE,3), round(fixedEffects.tStat,3), fixedEffects.pValue, ...
    'VariableNames', {'Estimate', 'StandardError', 'tStatistic', 'pValue'}, ...
    'RowNames', fixedEffects.Name);

writetable(fixedEffectsTable,'LME_opto.xlsx','Sheet',1)

 %% Figure 3D: plot the dprime per experiment and per session 
 clc
for g = 1:size(GNG_opto_all_cell,1)

dprimes_easy_g = reshape(squeeze(dprimes_easy(g,:,:)),[],1) ;
dprimes_easy_g = dprimes_easy_g(~isnan(dprimes_easy_g)); 

dprimes_easy_opto_g = reshape(squeeze(dprimes_easy_opto(g,:,:)),[],1) ;
dprimes_easy_opto_g = dprimes_easy_opto_g(~isnan(dprimes_easy_opto_g)); 

dprimes_hard_g = reshape(squeeze(dprimes_hard(g,:,:)),[],1) ;
dprimes_hard_g = dprimes_hard_g(~isnan(dprimes_hard_g)); 

dprimes_hard_opto_g = reshape(squeeze(dprimes_hard_opto(g,:,:)),[],1) ;
dprimes_hard_opto_g = dprimes_hard_opto_g(~isnan(dprimes_hard_opto_g)); 


 figure
 plot([1 2],[dprimes_easy_g dprimes_easy_opto_g],'Color',[0.5 0.5 0.5 0.2],'Marker','o','MarkerFaceColor',color_eh{1},'MarkerEdgeColor',color_eh{1})
 hold on

scatter(ones(size(dprimes_easy_opto_g)).*(1 + ones(size(dprimes_easy_opto_g))),dprimes_easy_opto_g...
     ,'Marker','o','MarkerFaceColor',color_opto{1},'MarkerEdgeColor',color_opto{1})
 hold on

 errorbar( nanmean([dprimes_easy_g dprimes_easy_opto_g] ),std([dprimes_easy_g dprimes_easy_opto_g])/...
    sqrt(size([dprimes_easy_g dprimes_easy_opto_g],1)),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)
hold on

 plot([3 4],[dprimes_hard_g dprimes_hard_opto_g],'Color',[0.5 0.5 0.5 0.2],'Marker','o','MarkerFaceColor',color_eh{2},'MarkerEdgeColor',color_eh{2})
 hold on

scatter(ones(size(dprimes_hard_opto_g)).*(3 + ones(size(dprimes_hard_opto_g))),dprimes_hard_opto_g...
     ,'Marker','o','MarkerFaceColor',color_opto{1},'MarkerEdgeColor',color_opto{1})
 hold on

 errorbar([nan nan nanmean([dprimes_hard_g dprimes_hard_opto_g])], [nan nan std([dprimes_hard_g dprimes_hard_opto_g])]/...
    sqrt(size([dprimes_hard_g dprimes_hard_opto_g],1)),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)
hold on

 ylim([-2 5]);
yticks([-1:1:4])
yline(1,'--','Color','k');
xlim([0 5])
xticks([1.5 3.5])
 yline(1,'linestyle','--','linewidth',2)
 xticklabels({'easy','hard'})
ylabel("discrimination (d')");
hold on    
box off
ax = gca;
ax.XAxis.FontSize = 20 ;
ax.YAxis.FontSize = 20 ;
ax.Title.FontSize = 20 ;
ax.Title.FontSize = 20 ;
movegui('east');

clear h
for i = 1:2
 h(i) = kstest(dprimes_easy_g) ;
end 
if sum(h) > 0
    p = [] ;
[p(1,1),~,~] = signrank(dprimes_easy_g,dprimes_easy_opto_g,'alpha',0.05,'tail','both') ;
[p(2,1),~,~] = signrank(dprimes_hard_g,dprimes_hard_opto_g,'alpha',0.05,'tail','both') ;
[p(3,1),~,~] = signrank(dprimes_easy_g,dprimes_hard_g,'alpha',0.05,'tail','both') ;
[p(4,1),~,~] = signrank(dprimes_easy_opto_g,dprimes_hard_opto_g,'alpha',0.05,'tail','both') ;
[~,p(1:4,1),~]= bonferroni_holm(p(1:4,1),0.05) ;

elseif sum(h) == 0 
    p = [] ;
[~,p(1,1),~] = ttest(dprimes_easy_g,dprimes_easy_opto_g,'alpha',0.05,'tail','both') ;
[~,p(2,1),~] = ttest(dprimes_hard_g,dprimes_hard_opto_g,'alpha',0.05,'tail','both') ;
[~,p(3,1),~] = ttest(dprimes_easy_g,dprimes_hard_g,'alpha',0.05,'tail','both')  ;
[~,p(4,1),~] = ttest(dprimes_easy_opto_g,dprimes_hard_opto_g,'alpha',0.05,'tail','both') ;
[~,p(1:4,1),~]= bonferroni_holm(p(1:4,1),0.05) ;
end 
disp(p)

if  p(1,:) < 0.0005
    txt = '***';
elseif  p(1,:) < 0.005
    txt = '**';
elseif  p(1,:) < 0.05
    txt = '*';
elseif p(1,:)> 0.05
    txt = 'n.s.';
end
text(mean([1.5 1.5]), 4.15 ,txt,'FontSize',15)
hold on
line( linspace(1.05,1.95), linspace(4,4),'color','k')

if  p(2,:) < 0.0005
    txt = '***';
elseif  p(2,:) < 0.005
    txt = '**';
elseif  p(2,:) < 0.05
    txt = '*';
elseif p(2,:)> 0.05
    txt = 'n.s.';
end
text(mean([3.5 3.5]), 4.15 ,txt,'FontSize',15)
hold on
line( linspace(3.05,3.95), linspace(4,4),'color','k')

if  p(3,:) < 0.0005
    txt = '***';
elseif  p(3,:) < 0.005
    txt = '**';
elseif  p(3,:) < 0.05
    txt = '*';
elseif p(3,:)> 0.05
    txt = 'n.s.';
end
text(mean([1.05 3.05]), 4.45 ,txt,'FontSize',15)
hold on
line( linspace(1.05, 3.05), linspace(4.3,4.3),'color','k')

if  p(4,:) < 0.0005
    txt = '***';
elseif  p(4,:) < 0.005
    txt = '**';
elseif  p(4,:) < 0.05
    txt = '*';
elseif p(4,:)> 0.05
    txt = 'n.s.';
end
text(mean([1.95 3.95]), 4.75 ,txt,'FontSize',15)
hold on
line( linspace(1.95, 3.95), linspace(4.6,4.6),'color','k')

end


 %% Figure 3E:  plot the lick rate per experiment for every session 
clc
% close all
for g = 1:size(GNG_opto_all_cell,1)
    figure
    for i   = 1:size(GNG_opto_all_cell,2) % mouse
        for s = 1:numel(GNG_opto_all_cell{g,i}) % session

            plot(freqs, [go_licks_all_easy_opto(g,i,s) go_licks_all_hard_opto(g,i,s) ngo_licks_all_hard_opto(g,i,s) ngo_licks_all_easy(g,i,s)]...
                ,'Color',color_opto_fade{1},'linewidth',2,'linestyle','-')
            hold on
            plot(freqs, [go_licks_all_easy(g,i,s) go_licks_all_hard(g,i,s) ngo_licks_all_hard(g,i,s) ngo_licks_all_easy(g,i,s)]...
                ,'Color',color_mean_marker{1},'linewidth',2,'linestyle','-')
        end



        mean_licks_mice(i,:) =  [nanmean(go_licks_all_easy(g,i,:)) nanmean(go_licks_all_hard(g,i,:)) nanmean(ngo_licks_all_hard(g,i,:)) nanmean(ngo_licks_all_easy(g,i,:))] ;
        std_licks_mice(i,:) =  [nanstd(go_licks_all_easy(g,i,:)) nanstd(go_licks_all_hard(g,i,:)) nanstd(ngo_licks_all_hard(g,i,:)) nanstd(ngo_licks_all_easy(g,i,:))] ;

        mean_licks_mice_opto(i,:) =  [nanmean(go_licks_all_easy_opto(g,i,:)) nanmean(go_licks_all_hard_opto(g,i,:)) nanmean(ngo_licks_all_hard_opto(g,i,:)) nanmean(ngo_licks_all_easy_opto(g,i,:))] ;
        std_licks_mice_opto(i,:) =  [nanstd(go_licks_all_easy_opto(g,i,:)) nanstd(go_licks_all_hard_opto(g,i,:)) nanstd(ngo_licks_all_hard_opto(g,i,:)) nanstd(ngo_licks_all_easy_opto(g,i,:))] ;

    end

    mean_licks = mean(mean_licks_mice) ;
    std_licks = mean(std_licks_mice) ;

    mean_licks_opto = mean(mean_licks_mice_opto) ;
    std_licks_opto = mean(std_licks_mice_opto) ;
    session_ID_g = reshape(session_ID(g,:,:),[],1) ;

    errorbar(freqs, mean_licks_opto ,std_licks_opto/...
        sqrt(length(session_ID_g)),'Color',color_opto_fade{1},'linestyle',L{1},'linewidth',5)
    hold on

    errorbar(freqs, mean_licks ,std_licks/...
        sqrt(length(session_ID_g)),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)
    hold on

 

ylim([0 1]);
yticks([0:0.2:1])
yline(0.5,'--','Color','k','linewidth',2);
xline(10,'--','Color','k','linewidth',2);
xlim([7 14.1])
xticks([7 9.2 10.9 14.1])
xlabel('freq (kHz)');
ylabel('lick rate');
set(gca, 'XDir','reverse','XScale','log');
box off;
hold on    
ax = gca;
ax.XAxis.FontSize = 20 ;
ax.YAxis.FontSize = 20 ;
ax.Title.FontSize = 20 ;
ax.Title.FontSize = 20 ;


h = kstest(reshape(go_licks_all_easy_opto(g,:,:),[],1)) ;

if h == 1
    [p_opto(1,1),~,~] = signrank( reshape( [reshape(go_licks_all_easy_opto(g,:,:),[],1) reshape(go_licks_all_hard_opto(g,:,:),[],1)] ,[],1),...
      reshape( [reshape(go_licks_all_easy(g,:,:),[],1) reshape(go_licks_all_hard(g,:,:),[],1)] ,[],1),'alpha',0.05,'tail','both') ;

    [p_opto(2,1),~,~] = signrank( reshape( [reshape(ngo_licks_all_easy_opto(g,:,:),[],1) reshape(ngo_licks_all_hard_opto(g,:,:),[],1)] ,[],1),...
      reshape( [reshape(ngo_licks_all_easy(g,:,:),[],1) reshape(ngo_licks_all_hard(g,:,:),[],1)] ,[],1),'alpha',0.05,'tail','both') ;
elseif h == 0 
    [~, p_opto(1,1),~] = ttest( reshape( [reshape(go_licks_all_easy_opto(g,:,:),[],1) reshape(go_licks_all_hard_opto(g,:,:),[],1)] ,[],1),...
      reshape( [reshape(go_licks_all_easy(g,:,:),[],1) reshape(go_licks_all_hard(g,:,:),[],1)] ,[],1),'alpha',0.05,'tail','both') ;

    [~,p_opto(2,1),~] = ttest( reshape( [reshape(ngo_licks_all_easy_opto(g,:,:),[],1) reshape(ngo_licks_all_hard_opto(g,:,:),[],1)] ,[],1),...
      reshape( [reshape(ngo_licks_all_easy(g,:,:),[],1) reshape(ngo_licks_all_hard(g,:,:),[],1)] ,[],1),'alpha',0.05,'tail','both') ;

end

[~,p_opto,~]= bonferroni_holm(p_opto,0.05) ;
disp(p_opto)

end

