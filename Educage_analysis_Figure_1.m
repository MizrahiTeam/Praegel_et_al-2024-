clear all
close all
clc
%% Open the cell array of the adolescent and adult table
%add the Praegel_et_al_MATLABR2023b_scripts directory here
addpath(genpath('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_MATLABR2023b_scripts'))
%add the Praegel_et_al_data directory here
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\Educage\EducageTable_cell.mat'...
    , ['Select table (EducageTable_cell) ']);
cd(path)
addpath(path)
load(file)
%% load stimulus excel file
% select the Disc_hard file here
[file, path] = uigetfile('*.xlsx',['Select stimulus file (Disc_hard_Educage) ']);
t = readtable(fullfile(path, file)) ;

%% Parameters
% Discrimination parameters
binsize = 100; % binsize for dprime trajectory
dprime_threshold = 1; % trials to threshold
numtrials = 1000; % naive and expert trial number
max_ISI = 60; % maximal inter trial interval
smoothing_window = 20 ; % smoothing of the trajectory

% colors
color_eh = {[0.3010 0.7450 0.9330], [0 0.4470 0.7410]} ;
color_eh_alpha = {[0.3010 0.7450 0.9330 0.4 ], [0 0.4470 0.7410 0.4 ]} ;
color_lick = [0.4660 0.6740 0.1880] ;
color_no_lick = [0.6350 0.0780 0.1840] ;
color_lick_back = [0.4660 0.6740 0.1880 0.1] ;
color_no_lick_back = [0.6350 0.0780 0.1840 0.1] ;
color_mean_marker = {[0.5 0.5 0.5 0.8]} ;
color_plot_marker = {[0.5 0.5 0.5 0.2]} ;

% line styles
L = {['-'], ['-']} ;
L2 = {['--'], ['-']};
M = {['o'],['o']} ;

% determine stimulus IDs for psychometric cruve
level = 4;
learned_freqs = [1:4];
catch_freqs = [5:11];
stim = [learned_freqs catch_freqs];
disc_freqs =t.Frequency(t.StimulusID(stim(1:4)))';
freqs =sort(t.Frequency(t.StimulusID(stim)))';
freqs = freqs(2:end-1);
go_freqs = freqs(1:5) ;
ngo_freqs = freqs(5:9) ;


%% calculate the trial outcomes and the dprime + cbias
clc
tic
for i   = 1:numel(EducageTable_cell)

    EducageTable = cell2table(EducageTable_cell(i));
    EducageTable = EducageTable.Var1{1,1};
    mice = unique(EducageTable.mouse_num);

    for n=1:length(mice)

        mouse_table = EducageTable(EducageTable.mouse_num == mice(n),:);

        %exclude trials above max ISI
        [mouse_table] = max_trial_dist (mouse_table, max_ISI);

        % trial outcomes throughout the task discrimination (1 oct. and 0.25 oct. together)
        index_all = find (mouse_table.level == 3 | mouse_table.level == 4);
        [hit_all{n,i}, fa_all{n,i},~, cr_all{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_all, mouse_table);
        [~, ~, dPrime_all{n,i}, cbias_all{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);

        % trial outcomes in the easy level:  1 oct. discrimination
        index_easy_level =find (mouse_table.level == 3);
        [hit_easy_level{n,i}, fa_easy_level{n,i}, ~, cr_easy_level{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_easy_level, mouse_table);
        [~, ~, dPrime_easy_level{n,i}, cbias_easy_level{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);

        % trial outcomes in the hard level + 1 oct. + 0.25 oct. discrimination
        index_easy = find (mouse_table.level == 4 & (mouse_table.stimID ==1 | mouse_table.stimID == 4));
        index_hard = find (mouse_table.level == 4 & (mouse_table.stimID ==2  | mouse_table.stimID == 3));

        % 1 oc.t discrimination
        [hit_easy{n,i}, fa_easy{n,i}, ~, cr_easy{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_easy, mouse_table);
        [~, ~, dPrime_easy{n,i}, cbias_easy{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);
        % 0.25 oct. discrimination
        [hit_hard{n,i}, fa_hard{n,i}, ~, cr_hard{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_hard, mouse_table);
        [~, ~, dPrime_hard{n,i}, cbias_hard{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);

        %Coefficient of variation
        cv_dPrime (n,i) = std(dPrime_all{n,i},0,'all') / mean(dPrime_all{n,i}) ;
        cv_cbias (n,i) = std(cbias_all{n,i},0,'all') / mean(cbias_all{n,i}) ;

        % trial to threshold - procedural learning
        trial_criterion_easy_level = find((dPrime_easy_level{n,i}(2:end-1)>= dprime_threshold))+1;
        if ~isempty(trial_criterion_easy_level);
            trial_criterion_easy_all_level(n,i) = trial_criterion_easy_level (1);
        elseif isempty (trial_criterion_easy_level)
            trial_criterion_easy_all_level(n,i) = NaN;
        end

        % naive performance
        Naive_dPrime_easy(n,i) = mean(dPrime_easy_level{n,i}(1:1:(numtrials/binsize))) ;
        Naive_dPrime_hard(n,i) = mean(dPrime_hard{n,i}(1:1:(numtrials/binsize))) ;

        %expert performance
        Expert_dPrime_easy_level(n,i) = mean(dPrime_easy_level{n,i}(size(dPrime_easy_level{n,i},2):-1:size(dPrime_easy_level{n,i},2)...
            -(numtrials/binsize)));
        Expert_dPrime_easy(n,i) = mean(dPrime_easy{n,i}(size(dPrime_easy{n,i},2):-1:size(dPrime_easy{n,i},2) -(numtrials/binsize)));
        Expert_dPrime_hard(n,i) = mean(dPrime_hard{n,i}(size(dPrime_hard{n,i},2):-1:size(dPrime_hard{n,i},2) -(numtrials/binsize)));

        cage_ID (n,i) = unique(mouse_table.table_num) ;


        % task plotting
        data_plot_easy{1,i}(n,:) = [Naive_dPrime_easy(n,i) Expert_dPrime_easy(n,i)];
        data_plot_hard{1,i}(n,:) =  [Naive_dPrime_hard(n,i) Expert_dPrime_hard(n,i)];
    end

end
toc
%% figure 1 d: plot an adolescent and adult example mouse
clc
close all
figure
for i = [1 2]
    if i == 1
        n = 2 ; % adolescent mouse
    elseif i == 2
        n = 5 ; % adult mouse
    end

    nan_easy_level = [] ;
    d_easy = [dPrime_easy_level{n,i} dPrime_easy{n,i}]  ;
    nan_easy_level = nan(1,size(dPrime_easy_level{n,i},2)) ;

    d_hard = [nan_easy_level dPrime_hard{n,i}] ;
    %smooth d' trakectory
    d_easy= smoothdata(d_easy,"gaussian",20) ;
    d_hard (size(nan_easy_level,2):end)= smoothdata(d_hard(size(nan_easy_level,2):end),"gaussian",5) ;

    subplot(1,2,i);
    plot(d_easy,'Color',color_eh{1},'linestyle',L{i},'linewidth',3)
    hold on
    plot(d_hard,'Color',color_eh{2},'linestyle',L{i},'linewidth',3)
    hold on
    ylim([-1 5])
    xlim([0 150])
    xline(size(nan_easy_level,2),'--','linewidth',3)
    yline(1,'--','linewidth',3)
    xticks([0:50:150])
    xticklabels ([0:5000:15000])
    yticks([-1:2:5])
    yticklabels ([-1:2:5])
    xlabel('trials')
    ylabel("d'")
    box off;
    hold on
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
end
%% Figure 1 e - trials to threshold - learning criterion

clc
figure
trial_criterion_easy_all_level(trial_criterion_easy_all_level == 0) = nan ;
violinplot([ trial_criterion_easy_all_level(:,1)*100,trial_criterion_easy_all_level(:,2)*100],{'adol.','adult'},"ViolinColor",{[0.5 0.5 0.5;0.5 0.5 0.5]} ) ;
xticks ([1 2])
xlim([0 3]);
ylim([0 6000]);
yticks([0:2000:6000])
yticklabels([0:2000:6000])
xticklabels({'adolescent','adult'})
ylabel(' trials to threshold')
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');


%statistics
for i = 1:2
    h(i) = kstest(data_plot_easy{1,i}) ;
end
if sum(h) > 0
    h = [] ;
    [p,~,~] = ranksum(trial_criterion_easy_all_level(:,1),trial_criterion_easy_all_level(:,2),'alpha',0.05,'tail','both')
elseif sum(h) == 0
    h = [] ;
    [~,p,~] = ttest2(trial_criterion_easy_all_level(:,1),trial_criterion_easy_all_level(:,2),'alpha',0.05,'tail','both') ;
end
if p < 0.0005
    txt = '***';
elseif p < 0.005
    txt = '**';
elseif  p < 0.05
    txt = '*';
elseif p > 0.05
    txt = 'n.s.';
end
hold on
text(mean([1.4 1.4]), 5700 ,txt,'FontSize',20)
%% Figure 1f 1 oct. d' naive and expert   (also supplementary Fiure 1.2 A left)
clc
close all

% plot d' per mean and per mouse
figure
errorbar( nanmean(data_plot_easy{1,1}),std(data_plot_easy{1,1})/...
    sqrt(size(data_plot_easy{1,1},1)),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)
hold on
plot([1 2],data_plot_easy{1,1},'Color',color_plot_marker{1},'linestyle',L{1},...
    'Marker','.','MarkerFaceColor',color_eh{1},'MarkerEdgeColor',color_eh{1}) %adolescent

errorbar( [nan nan nanmean(data_plot_easy{1,2})], [nan nan std(data_plot_easy{1,2})]/...
    sqrt(size(data_plot_easy{1,2},1)),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)
plot([3 4], data_plot_easy{1,2},'Color',color_plot_marker{1},'linestyle',L{1},...
    'Marker','.','MarkerFaceColor',color_eh{1},'MarkerEdgeColor',color_eh{1}) %adult

hold on
xticks ([1 2 3 4 ])
xlim([0 5]);
yticks ([-1:2:5])
ylim([-1 5.2]);
yline(1,'linestyle','--','linewidth',2)
xticklabels({'adolescent naive','adolescent exp.','adult naive','adult exp.'})
ylabel("d'")

box off;
makepretty;
hold on
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
for i = 1:2
    h(i) = kstest(data_plot_easy{1,i}) ;
end
if sum(h) > 0
    h = [] ;
    [p_easy(3,1),p,~] = ranksum(data_plot_easy{1,1}(:,1),data_plot_easy{1,2}(:,1),'alpha',0.05,'tail','both') ;
    [p_easy(4,1),p,~] = ranksum(data_plot_easy{1,1}(:,2),data_plot_easy{1,2}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_easy(3:4,1),~]= bonferroni_holm(p_easy(3:4,1),0.05) ;
    [p_easy(1,1),p,~] = signrank(data_plot_easy{1,1}(:,1),data_plot_easy{1,1}(:,2),'alpha',0.05,'tail','both') ;
    [p_easy(2,1),p,~] = signrank(data_plot_easy{1,2}(:,1),data_plot_easy{1,2}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_easy(1:2,1),~]= bonferroni_holm(p_easy(1:2,1),0.05) ;

elseif sum(h) == 0
    h = [] ;
    [~,p_easy(3,1),~] = ttest2(data_plot_easy{1,1}(:,1),data_plot_easy{1,2}(:,1),'alpha',0.05,'tail','both') ;
    [~,p_easy(4,1),~] = ttest2(data_plot_easy{1,1}(:,2),data_plot_easy{1,2}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_easy(3:4,1),~]= bonferroni_holm(h(3:4,1),0.05) ;
    [~,p_easy(1,1),~] = ttest(data_plot_easy{1,1}(:,1),data_plot_easy{1,1}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_easy(2,1),~] = ttest(data_plot_easy{1,2}(:,1),data_plot_easy{1,2}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_easy(2:1,1),~]= bonferroni_holm(p_easy(2:1,1),0.05) ;

end
disp(p_easy)

if  p_easy(1,1) < 0.0005
    txt = '***';
elseif  p_easy(1,1) < 0.005
    txt = '**';
elseif  p_easy(1,1) < 0.05
    txt = '*';
elseif p_easy(1,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,1.95), linspace(4,4),'color','k')
hold on
text(mean([1.5 1.5]), 4.2 ,txt,'FontSize',20)

if  p_easy(2,1) < 0.0005
    txt = '***';
elseif  p_easy(2,1) < 0.005
    txt = '**';
elseif  p_easy(2,1) < 0.05
    txt = '*';
elseif p_easy(2,1) > 0.05
    txt = 'n.s.';
end
line( linspace(3.05,3.95), linspace(4,4),'color','k')
hold on
text(mean([3.5 3.5]), 4.2 ,txt,'FontSize',20)

if  p_easy(3,1) < 0.0005
    txt = '***';
elseif  p_easy(3,1) < 0.005
    txt = '**';
elseif  p_easy(3,1) < 0.05
    txt = '*';
elseif p_easy(3,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,3.05), linspace(4.4,4.4),'color','k')
hold on
text(mean([2.05 2.05]), 4.6 ,txt,'FontSize',20)

if  p_easy(4,1) < 0.0005
    txt = '***';
elseif  p_easy(4,1) < 0.005
    txt = '**';
elseif  p_easy(4,1) < 0.05
    txt = '*';
elseif p_easy(4,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.95,3.95), linspace(4.8,4.8),'color','k')
hold on
text(mean([3 3]), 5 ,txt,'FontSize',20)

%% Figure 1g 0.25 oct naive and expert d'  (also supplementary Fiure 1.2 A right)

% plot d' per mean and per mouse
figure
errorbar( nanmean(data_plot_hard{1,1}),std(data_plot_hard{1,1})/...
    sqrt(size(data_plot_hard{1,1},1)),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)
hold on
plot([1 2],data_plot_hard{1,1},'Color',color_plot_marker{1},'linestyle',L{1},...
    'Marker','.','MarkerFaceColor',color_eh{2},'MarkerEdgeColor',color_eh{2}) %adolescent

errorbar( [nan nan nanmean(data_plot_hard{1,2})], [nan nan std(data_plot_hard{1,2})]/...
    sqrt(size(data_plot_hard{1,2},1)),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5)
plot([3 4], data_plot_hard{1,2},'Color',color_plot_marker{1},'linestyle',L{1},...
    'Marker','.','MarkerFaceColor',color_eh{2},'MarkerEdgeColor',color_eh{2}) %adult

hold on
xticks ([1 2 3 4])
xlim([0 5]);
yticks ([-1:2:5])
ylim([-1 5.2]);
yline(1,'linestyle','--','linewidth',2)
xticklabels({'adolescent naive','adolescent exp.','adult naive','adult exp.'})
ylabel("d'")
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');


%statistics
for i = 1:2
    h(i) = kstest(data_plot_hard{1,i}) ;
end
if sum(h) > 0
    h = [] ;

    [p_hard(3,1),p,~] = ranksum(data_plot_hard{1,1}(:,1),data_plot_hard{1,2}(:,1),'alpha',0.05,'tail','both') ;
    [p_hard(4,1),p,~] = ranksum(data_plot_hard{1,1}(:,2),data_plot_hard{1,2}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_hard(3:4,1),~]= bonferroni_holm(p_hard(3:4,1),0.05) ;
    [p_hard(1,1),p,~] = signrank(data_plot_hard{1,1}(:,1),data_plot_hard{1,1}(:,2),'alpha',0.05,'tail','both') ;
    [p_hard(2,1),p,~] = signrank(data_plot_hard{1,2}(:,1),data_plot_hard{1,2}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_hard(1:2,1),~]= bonferroni_holm(p_hard(1:2,1),0.05) ;

elseif sum(h) == 0
    h = [] ;

    [~,p_hard(3,1),~] = ttest2(data_plot_hard{1,1}(:,1),data_plot_hard{1,2}(:,1),'alpha',0.05,'tail','both') ;
    [~,p_hard(4,1),~] = ttest2(data_plot_hard{1,1}(:,2),data_plot_hard{1,2}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_hard(3:4,1),~]= bonferroni_holm(p_hard(3:4,1),0.05) ;
    [~,p_hard(1,1),~] = ttest(data_plot_hard{1,1}(:,1),data_plot_hard{1,1}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_hard(2,1),~] = ttest(data_plot_hard{1,2}(:,1),data_plot_hard{1,2}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_hard(1:2,1),~]= bonferroni_holm(p_hard(1:2,1),0.05) ;

end

disp(p_hard)
if  p_hard(1,1) < 0.0005
    txt = '***';
elseif  p_hard(1,1) < 0.005
    txt = '**';
elseif  p_hard(1,1) < 0.05
    txt = '*';
elseif p_hard(1,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,1.95), linspace(4,4),'color','k')
hold on
text(mean([1.5 1.5]), 4.2 ,txt,'FontSize',20)

if  p_hard(2,1) < 0.0005
    txt = '***';
elseif  p_hard(2,1) < 0.005
    txt = '**';
elseif  p_hard(2,1) < 0.05
    txt = '*';
elseif p_hard(2,1) > 0.05
    txt = 'n.s.';
end
line( linspace(3.05,3.95), linspace(4,4),'color','k')
hold on
text(mean([3.5 3.5]), 4.2 ,txt,'FontSize',20)

if  p_hard(3,1) < 0.0005
    txt = '***';
elseif  p_hard(3,1) < 0.005
    txt = '**';
elseif  p_hard(3,1) < 0.05
    txt = '*';
elseif p_hard(3,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,3.05), linspace(4.4,4.4),'color','k')
hold on
text(mean([2.05 2.05]), 4.6 ,txt,'FontSize',20)

if  p_hard(4,1) < 0.0005
    txt = '***';
elseif  p_hard(4,1) < 0.005
    txt = '**';
elseif  p_hard(4,1) < 0.05
    txt = '*';
elseif p_hard(4,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.95,3.95), linspace(4.8,4.8),'color','k')
hold on
text(mean([3 3]), 5 ,txt,'FontSize',20)
%% LME effect of age,  group housing and sex on dprime
clc
for i   = 1:numel(EducageTable_cell) % perform per group

    dprime_average(:,i) = mean ([Expert_dPrime_easy(:,i)'; Expert_dPrime_hard(:,i)']) ;

    idx_f = find(EducageTable_cell{1,i}.sex == 'female') ;
    idx_m = find(EducageTable_cell{1,i}.sex == 'male') ;

    female_mice = unique(EducageTable_cell{1,i}.mouse_num(idx_f,:)) ;
    male_mice =unique(EducageTable_cell{1,i}.mouse_num(idx_m,:)) ;

    male_ID = nan(height(female_mice),1) ;
    female_ID = nan(height(male_mice),1) ;
    male_ID(1:end,:) = 2 ;
    female_ID(1:end,:) = 1 ;

    age_sex_ID(:,i) = [male_ID' female_ID']
end

%LME test if the effect is relevent per group
adult_ID = nan(height(Expert_dPrime_easy),1) ;
adolescent_ID = nan(height(Expert_dPrime_easy),1) ;
adult_ID(1:end,:) = 2 ;
adolescent_ID(1:end,:) = 1 ;

LME_edu = table ;
LME_edu.group = [adolescent_ID' adult_ID']' ;
LME_edu.dprime = [dprime_average(:,1)' dprime_average(:,2)']' ;
LME_edu.sex = [age_sex_ID(:,1)' age_sex_ID(:,2)']' ;
LME_edu.cage = [cage_ID(:,1)' (cage_ID(:,2)' + max(cage_ID(:,1)))]' ;

% Convert categorical variables
LME_edu.group = categorical(LME_edu.group);
LME_edu.sex = categorical(LME_edu.sex);
LME_edu.cage = categorical(LME_edu.cage);


model_edu = fitlme(LME_edu, 'dprime ~ group + sex + (1|cage)') ;
disp(model_edu);

% Extracting necessary information from the LinearMixedModel object
fixedEffects = model_edu.Coefficients;
randomEffects = model_edu.randomEffects;

% Create a table for fixed effects
fixedEffectsTable = table(fixedEffects.Estimate, fixedEffects.SE, fixedEffects.tStat, fixedEffects.pValue, ...
    'VariableNames', {'Estimate', 'StandardError', 'tStatistic', 'pValue'}, ...
    'RowNames', fixedEffects.Name);

writetable(fixedEffectsTable,'LME_edu.xlsx','Sheet',1)

%% psychometric curve calculation
clc
close all
%% calculate lick rates
if ~any(strcmp(EducageTable_cell,'1'))
    for i   = 1:numel(EducageTable_cell)
        EducageTable = cell2table(EducageTable_cell(i));
        EducageTable = EducageTable.Var1{1,1};
        data = EducageTable( EducageTable.level == level, :);

        [data] = max_trial_dist (data, max_ISI);

        mice = unique(EducageTable.mouse_num);

        % extract the lick rates per mouse
        [l_all, l_learned, stderr_l_all, mean_l_all, stderr_l_learend, mean_l_learned] = lick_ratio_all (catch_freqs, learned_freqs,mice, stim,data);
        % extract the mean lick rate per group
        [lick_all{i,1}, stderr_lick_all{i,1}, mean_lick_all{i,1}, lick_learned{i,1}, stderr_lick_learned{i,1}, mean_lick_learned{i,1}] = sort_freqs (level, l_all, stderr_l_all, mean_l_all);
        % fit a sigmoid function and normalize the curve
        [fitted_curve_all{i,1},mean_fitted_curve{i,1}, stderr_fitted_curve{i,1}] = fit_norm_sigmoid (freqs,lick_all{i,1});

    end
end
%% Fig.1 h : plot raw psychometric curves   (also supplementary Fiure 1.2 B)
clc
close all
for i   = 1:numel(EducageTable_cell)

    figure(3);
    plot(go_freqs,nanmean(lick_all{i,1}(:,1:5)),'linewidth',2,'linestyle',L2{i},'Color',color_lick)
    hold on
    plot(ngo_freqs,nanmean(lick_all{i,1}(:,5:9)),'linewidth',2,'linestyle',L2{i},'Color',color_no_lick)
    hold on

    patch([go_freqs flip(go_freqs)] , [(nanmean(lick_all{i,1}(:,1:5))) + (stderr_lick_all{i,1}(:,1:5))...
        flip((nanmean(lick_all{i,1}(:,1:5))) - (stderr_lick_all{i,1}(:,1:5)))],color_lick,...
        'facealpha' , 0.1, 'EdgeColor','none')

    hold on
    patch([ngo_freqs flip(ngo_freqs)] , [(nanmean(lick_all{i,1}(:,5:9))) + (stderr_lick_all{i,1}(:,5:9))...
        flip((nanmean(lick_all{i,1}(:,5:9))) - (stderr_lick_all{i,1}(:,5:9)))],color_no_lick,...
        'facealpha' , 0.1, 'EdgeColor','none')

    hold on
    ylim([0 1]);
    yticks([0:0.2:1]);
    xline(10,'--','Color','k','linewidth',1);
    yline(0.5,'--','Color','k','linewidth',1);

    xlim([min(t.Frequency) max(t.Frequency)])
    xlabel('freq (kHz)');
    ylabel('lick rate');
    set(gca, 'XDir','reverse','XScale','log');
    box off;
    makepretty;
    hold on
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');

end

%statistics
h = kstest(lick_all{1,1}) ;
if h == 1
    [p_freq(1,1),~,~] = ranksum(reshape(lick_all{1,1}(:,[1 3]),[],1),reshape(lick_all{2,1}(:,[1 3]),[],1),'alpha',0.05,'tail','both')
    [p_freq(2,1),~,~] = ranksum(reshape(lick_all{1,1}(:,[7 9]),[],1),reshape(lick_all{2,1}(:,[7 9]),[],1),'alpha',0.05,'tail','both')
elseif h == 0
    [~,p_freq(1,1),~] = ttest2(reshape(lick_all{1,1}(:,[1 3]),[],1),reshape(lick_all{1,1}(:,[1 3]),[],1),'alpha',0.05,'tail','both') ;
    [~,p_freq(2,1),~] = ttest2(reshape(lick_all{1,1}(:,[7 9]),[],1),reshape(lick_all{1,1}(:,[7 9]),[],1),'alpha',0.05,'tail','both') ;
end

[~,p_freq,~]= bonferroni_holm(p_freq,0.05) ;
disp(p_freq)
%% Fig. 1 h: normalized psychometric curve fitted to a sigmoid curve
clc
close all

for i   = 1:numel(EducageTable_cell)
    for n = 1:size(fitted_curve_all{i,1},1);
        perc = find (fitted_curve_all{i,1}(n,:) > 0.5);
        percile(n,i) = perc(end);
        percentile(n,i) = fitted_curve_all{1,1}(n,percile(n,i));
        freq_percentile(n,i) = freqs(:,percile(n,i));
    end

    figure(1)
    plot(go_freqs,mean_fitted_curve{i,1}(:,1:5),'linewidth',2,'linestyle',L2{i},'Color',color_lick)
    hold on
    plot(ngo_freqs,mean_fitted_curve{i,1}(:,5:9),'linewidth',2,'linestyle',L2{i},'Color',color_no_lick)
    hold on

    patch([go_freqs flip(go_freqs)] , [mean_fitted_curve{i,1}(:,1:5) + (stderr_fitted_curve{i,1}(:,1:5))...
        flip(mean_fitted_curve{i,1}(:,1:5) - (stderr_fitted_curve{i,1}(:,1:5)))],color_lick,...
        'facealpha' , 0.1, 'EdgeColor','none')

    hold on
    patch([ngo_freqs flip(ngo_freqs)] , [mean_fitted_curve{i,1}(:,5:9) + (stderr_fitted_curve{i,1}(:,5:9))...
        flip(mean_fitted_curve{i,1}(:,5:9) - (stderr_fitted_curve{i,1}(:,5:9)))],color_no_lick,...
        'facealpha' , 0.1, 'EdgeColor','none')

    hold on
    ylim([0 1]);
    yticks([0:0.2:1]);
    xline(10,'--','Color','k');
    yline(0.5,'--','Color','k');
    xlim([min(t.Frequency) max(t.Frequency)])
    set(gca, 'XDir','reverse','XScale','log');
    xlabel('freq  (kHz)');
    ylabel('normalized lick rate');
    box off;
    makepretty;
    hold on
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
end

%statistics
[p_freq_perc(1,1),~,~] = ranksum(freq_percentile(:,1),freq_percentile(:,2),'alpha',0.05,'tail','both')
[p_freq_perc(2,1),~,~] = ranksum(percentile(:,1),percentile(:,2),'alpha',0.05,'tail','both')
[~,p_freq_perc,~]= bonferroni_holm(p_freq_perc,0.05) ;
disp(p_freq_perc)
%% Fig 1 j: Calculate maximal c-bias
clc
tic

for i   = 1:numel(EducageTable_cell)
    EducageTable = cell2table(EducageTable_cell(i));
    EducageTable = EducageTable.Var1{1,1};

    mice = unique(EducageTable.mouse_num);

    for n=1:length(mice)

        mouse_table = EducageTable(EducageTable.mouse_num == mice(n),:);

        numtrials_b = int16(length(mouse_table.level(mouse_table.level == 4))/3);
        [mouse_table] = max_trial_dist (mouse_table, max_ISI);

        index_all = find ( mouse_table.level == 4);
        [~, ~, ~, ~, go_licks, ngo_licks] = trial_outcomes (binsize, index_all, mouse_table);
        [~, ~, dPrime_all{n,i}, cbias_all{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);

        Expert_bias_all{1,i}(n,:) =  min(cbias_all{n,i}(size(cbias_all{n,i},2):-1:size(cbias_all{n,i},2) -(numtrials_b/binsize)));
    end
end
toc

figure
violinplot([ Expert_bias_all{1,1},Expert_bias_all{1,2}],{'adolescent','adult'},"ViolinColor",{[0.5 0.5 0.5;0.5 0.5 0.5]} )
xticks ([1 2])
xlim([0 3]);
ylim([-3 0]);
yticks([-3:1:0])
yticklabels([-3:1:0])
xticklabels({'adolescent','adult'})
ylabel(' max. lick bias (c-bias)')
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');


%statistics
for i = 1:2
    h(i) = kstest(Expert_bias_all{1,i}) ;
end
if sum(h) > 0
    h = [] ;
    [p,~,~] = ranksum(Expert_bias_all{1,1},Expert_bias_all{1,2},'alpha',0.05,'tail','both') ;
elseif sum(h) == 0
    h = [] ;
    [~,p,~] = ttest2(Expert_bias_all{1,1},Expert_bias_all{1,2},'alpha',0.05,'tail','both') ;
end
disp(p)

if p < 0.0005
    txt = '***';
elseif p < 0.005
    txt = '**';
elseif  p < 0.05
    txt = '*';
elseif p > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,1.95), linspace(-0.2,-0.2),'color','k')
hold on
text(mean([1.5 1.5]), 0 ,txt,'FontSize',20)

%% supplmentary figure 1.3 A: CR rate throughout the experiment per mouse
clc
close all
easy_all_level = dPrime_easy_level ;

for i   = 1:numel(EducageTable_cell)
    n = max(cellfun(@length,easy_all_level));

    nan_easy_level = [] ;

    figure
    for n=1:max(mice)

        nan_easy_level = nan(1,size(easy_all_level{n,i},2)) ;

        nan_easy_level_mean(i,n) =    size(nan_easy_level,2) ;

        cr_all{n,i} = smoothdata(cr_all{n,i},"gaussian",20) ;

        plot(cr_all{n,i}(1:end-1),'Color',[0.4940 0.1840 0.5560 0.4],'linestyle',L{i},'LineWidth',3)
        hold on
        hold on
        xline(size(nan_easy_level,2),'--')

    end

    ylim([0 55])
    xlim([0 301])
    yticks([0:25:50])
    yticklabels ([0:25:50])
    xticks([0:270:270])
    xticklabels ([0:27000:27000])
    xlabel('trials')
    ylabel('CR Rate (%)')
    box off;
    hold on
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');

end

%% supplmentary figure 1.3 C: d' rate throughout the experiment per mouse


for i   = 1:numel(EducageTable_cell)
    n = max(cellfun(@length,dPrime_easy_level));

    figure
    for n=1:max(mice)
        nan_easy_level = [] ;
        d_easy = [dPrime_easy_level{n,i} dPrime_easy{n,i}]  ;

        nan_easy_level = nan(1,size(dPrime_easy_level{n,i},2)) ;

        nan_easy_level_mean(i,n) = size(nan_easy_level,2) ;

        d_hard = [nan_easy_level dPrime_hard{n,i}] ;

        d_easy= smoothdata(d_easy,"gaussian",20) ;
        d_hard= smoothdata(d_hard,"gaussian",20) ;

        plot(d_easy,'Color',color_eh_alpha{1},'linestyle',L{i},'LineWidth',3)
        hold on
        plot(d_hard,'Color',color_eh_alpha{2},'linestyle',L{i},'LineWidth',3)
        hold on
        xline(size(nan_easy_level,2),'--','Color',[0.5 0.5 0.5 0.1])
    end

    hold on
    ylim([-1 5])
    xlim([0 201])
    yline(1,'--')
    xticks([0:100:201])
    xticklabels ([0:10000:20100])
    yticks([-1:2:5])
    yticklabels ([-1:2:5])
    xlabel(' trials')
    ylabel("d'")
    box off;
    hold on
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
end
%% supplementartyFigure 1.3 B: coefficient of variation of the CR rate

for i   = 1:numel(EducageTable_cell)
    EducageTable = cell2table(EducageTable_cell(i));
    EducageTable = EducageTable.Var1{1,1};

    mice = unique(EducageTable.mouse_num);

    for n=1:length(mice)
        mouse_table = EducageTable(EducageTable.mouse_num == mice(n),:);
        numtrials_b = int16(length(mouse_table.level(mouse_table.level == 4))/3);
        [mouse_table] = max_trial_dist (mouse_table, max_ISI);
        index_all = find ( mouse_table.level == 3 | mouse_table.level == 4);
        %calculate across the whole task 
        [hit_all{n,i},  fa_all{n,i}, ~,cr_all{n,i}, go_licks, ngo_licks] = trial_outcomes (binsize, index_all, mouse_table);
        [~, ~, dPrime_all{n,i}, cbias_all{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);
        %calculate the CV
        cv_cbias (n,i) = abs( std( cr_all{n,i}) / mean( cr_all{n,i}) );
    end
end

figure
violinplot([ cv_cbias(:,1),cv_cbias(:,2)],{'adol.','adult'},"ViolinColor",{[0.4940 0.1840 0.5560; 0.4940 0.1840 0.5560]} )

xticks ([1 2])
xlim([0 3]);
ylim([0 2]);
yticks([0:1:2])
yticklabels([0:1:2])
xticklabels({'adol.','adult'})
ylabel('CV of CR Rate')
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
for i = 1:2
    h(i) = kstest(cv_cbias(:,1)) ;
end
if sum(h) > 0
    h = [] ;
    [p,~,~] = ranksum(cv_cbias(:,1),cv_cbias(:,2),'alpha',0.05,'tail','both')
elseif sum(h) == 0
    h = [] ;
    [~,p,~] = ttest2(cv_cbias(:,1),cv_cbias(:,2),'alpha',0.05,'tail','both') ;
end
disp(p)
if p < 0.0005
    txt = '***';
elseif p < 0.005
    txt = '**';
elseif  p < 0.05
    txt = '*';
elseif p > 0.05
    txt = 'n.s.';
end
hold on
text(mean([1.4 1.4]), 1.5 ,txt,'FontSize',20)

%% supplementary figre 1.3 D: coefficient of variation of dprime
clc
tic
for i   = 1:numel(EducageTable_cell)
    EducageTable = cell2table(EducageTable_cell(i));
    EducageTable = EducageTable.Var1{1,1};
    mice = unique(EducageTable.mouse_num);
    for n=1:length(mice)
        mouse_table = EducageTable(EducageTable.mouse_num == mice(n),:);
        numtrials_b = int16(length(mouse_table.level(mouse_table.level == 4))/3);
        [mouse_table] = max_trial_dist (mouse_table, max_ISI);
        % caclulate for level 4
        index_all = find ( mouse_table.level == 4);
        [~, ~, ~, ~, go_licks, ngo_licks] = trial_outcomes (binsize, index_all, mouse_table);
        [~, ~, dPrime_all{n,i}, cbias_all{n,i}] = dprime(binsize, go_licks, ngo_licks, mice);
        % calculate the CV
        cv_dPrime (n,i) = std(dPrime_all{n,i},0,'all') / mean(dPrime_all{n,i}) ;
    end
end

figure
violinplot([ cv_dPrime(:,1),cv_dPrime(:,2)],{'adol.','adult'},"ViolinColor",{[color_eh{1}; color_eh{1}]} )
xticks ([1 2])
xlim([0 3]);
ylim([0 2]);
yticks([0:0.5:2])
yticklabels([0:0.5:2])
xticklabels({'adol.','adult'})
ylabel('CV of dprime')
box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');

%statistics
for i = 1:2
    h(i) = kstest([cv_dPrime(:,1) cv_dPrime(:,2)]) ;
end
if sum(h) > 0
    h = [] ;
    [~,p,~] = ttest2(cv_dPrime(:,1),cv_dPrime(:,2),'alpha',0.05,'tail','both') ;
elseif sum(h) == 0
    h = [] ;
    [p,~,~] = ranksum(cv_dPrime(:,1),cv_dPrime(:,2),'alpha',0.05,'tail','both')
end
disp(p)
if p < 0.0005
    txt = '***';
elseif p < 0.005
    txt = '**';
elseif  p < 0.05
    txt = '*';
elseif p > 0.05
    txt = 'n.s.';
end
hold on
text(mean([1.4 1.4]), 1.9 ,txt,'FontSize',20)



%% save Educage d' output for comparison with the set
for i = 1:numel(EducageTable_cell)
    % sort out mice with d' easy <1 to ensure the same criterion as in the
    % recording
    dprime_edu.dprime_easy(1:size(Expert_dPrime_easy(Expert_dPrime_easy(:,i)>1,i),1),i)=  Expert_dPrime_easy(Expert_dPrime_easy(:,i)>1,i);
    dprime_edu.dprime_hard(1:size(Expert_dPrime_easy(Expert_dPrime_easy(:,i)>1,i),1),i) =  Expert_dPrime_hard(Expert_dPrime_easy(:,i)>1,i);
end
save('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\Educage\dprime_edu.mat', 'dprime_edu')

%% supplementary Figure 1.1 :  control experiment - hard task only
clc
close all
% choose EducageTable_control_cell.mat
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\Educage\EducageTable_cell.mat'...
    , ['Select table: Educage_control_cell ']);
addpath(path)
load (file)
%% calculate dprime in the control task
tic
for i   = 1:numel(EducageTable_control_cell)
    EducageTable = cell2table(EducageTable_control_cell(i));
    EducageTable = EducageTable.Var1{1,1};

    mice = unique(EducageTable.mouse_num);

    for n=1:length(mice)
        mouse_table = EducageTable(EducageTable.mouse_num == mice(n),:);

        [mouse_table] = max_trial_dist (mouse_table, max_ISI);

        index_hard = find (mouse_table.level == 3 & (mouse_table.stimID ==1  | mouse_table.stimID == 2));

        [~, ~, ~, ~, go_licks, ngo_licks] = trial_outcomes (binsize, index_hard, mouse_table);
        [~, ~, dPrime_hard, ~] = dprime(binsize, go_licks, ngo_licks, mice);

        dPrime_control_all{n,i} = dPrime_hard;

        Naive_dPrime_control(n,i) = mean(dPrime_hard(1:1:(numtrials/binsize))) ;

        Expert_dPrime_control(n,i) = mean(dPrime_hard(size(dPrime_hard,2):-1:size(dPrime_hard,2) -(numtrials/binsize)));

        data_plot_control{1,i}(n,:) =  [Naive_dPrime_control(n,i) Expert_dPrime_control(n,i)];

    end
end
toc
%% supplementary Figure 1.1 a: d' control trajectory

for i   = 1:numel(EducageTable_control_cell)
    n = max(cellfun(@length,dPrime_control_all));

    figure
    for n=1:max(mice)

        d_hard = [dPrime_control_all{n,i}] ;

        d_hard = smoothdata(d_hard,'gaussian',smoothing_window)

        plot(d_hard,'Color',color_eh{2},'linestyle',L{i},'linewidth',4)
        hold on
        ylim([-1 5])
        xlim([0 101])
        yline(1,'--')
        xticks([0:50:101])
        xticklabels ([0:5000:10000])
        yticks([-1:2:5])
        yticklabels ([-1:2:5])
        xlabel('trials')
        ylabel("d'")
        box off;
        makepretty;
        hold on
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        movegui('east');


    end
end
%% supplementary Figure 1.1 B: d' control naive a expert d'
clc
close all

figure
plot([1 2],nanmean(data_plot_control{1,1}),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5) %adolescent
hold on

plot([1 2],data_plot_control{1,1},'Color',color_plot_marker{1},'linestyle',L{1},...
    'Marker','.','MarkerFaceColor',color_eh{2},'MarkerEdgeColor',color_eh{2}) %adolescent
hold on

plot([3 4],nanmean(data_plot_control{1,2}),'Color',color_mean_marker{1},'linestyle',L{1},'linewidth',5) %adolescent
hold on

plot([3 4],data_plot_control{1,2},'Color',color_plot_marker{1},'linestyle',L{1},...
    'Marker','.','MarkerFaceColor',color_eh{2},'MarkerEdgeColor',color_eh{2}) %adult
hold on
xticks ([1 2 3 4])
xlim([0 5]);
yticks ([-1:5])
ylim([-1 5]);
ylabel("d'")
yline(1,'linestyle','--','linewidth',2)
xticklabels({'adolescent naive','adolescent exp.','adult naive','adult exp.',})

box off;
makepretty;
hold on
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');


%statistics
h = kstest(data_plot_control{1,1})
if h > 0
    [p_hard(3,1),~,~] = ranksum(data_plot_control{1,1}(:,1),data_plot_control{1,2}(:,1),'alpha',0.05,'tail','both') ;
    [p_hard(4,1),~,~] = ranksum(data_plot_control{1,1}(:,2),data_plot_control{1,2}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_hard(3:4,1),~]= bonferroni_holm(p_hard(3:4,1),0.05) ;
    [p_hard(1,1),~,~] = signrank(data_plot_control{1,1}(:,1),data_plot_control{1,1}(:,2),'alpha',0.05,'tail','both') ;
    [p_hard(2,1),~,~] = signrank(data_plot_control{1,2}(:,1),data_plot_control{1,2}(:,2),'alpha',0.05,'tail','both') ;

elseif h == 0
    [~,p_hard(3,1),~] = ttest2(data_plot_control{1,1}(:,1),data_plot_control{1,2}(:,1),'alpha',0.05,'tail','both') ;
    [~,p_hard(4,1),~] = ttest2(data_plot_control{1,1}(:,2),data_plot_control{1,2}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_hard(3:4,1),~]= bonferroni_holm(p_hard(3:4,1),0.05) ;
    [~,p_hard(1,1),~] = ttest(data_plot_control{1,1}(:,1),data_plot_control{1,1}(:,2),'alpha',0.05,'tail','both') ;
    [~,p_hard(2,1),~] = ttest(data_plot_control{1,2}(:,1),data_plot_control{1,2}(:,2),'alpha',0.05,'tail','both') ;
end
disp(p_hard)

if  p_hard(1,1) < 0.0005
    txt = '***';
elseif  p_hard(1,1) < 0.005
    txt = '**';
elseif  p_hard(1,1) < 0.05
    txt = '*';
elseif p_hard(1,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,1.95), linspace(4,4),'color','k')
hold on
text(mean([1.5 1.5]), 4.2 ,txt,'FontSize',20)

if  p_hard(2,1) < 0.0005
    txt = '***';
elseif  p_hard(2,1) < 0.005
    txt = '**';
elseif  p_hard(2,1) < 0.05
    txt = '*';
elseif p_hard(2,1) > 0.05
    txt = 'n.s.';
end
line( linspace(3.05,3.95), linspace(4,4),'color','k')
hold on
text(mean([3.5 3.5]), 4.2 ,txt,'FontSize',20)

if  p_hard(3,1) < 0.0005
    txt = '***';
elseif  p_hard(3,1) < 0.005
    txt = '**';
elseif  p_hard(3,1) < 0.05
    txt = '*';
elseif p_hard(3,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.05,3.05), linspace(4.4,4.4),'color','k')
hold on
text(mean([2.05 2.05]), 4.6 ,txt,'FontSize',20)

if  p_hard(4,1) < 0.0005
    txt = '***';
elseif  p_hard(4,1) < 0.005
    txt = '**';
elseif  p_hard(4,1) < 0.05
    txt = '*';
elseif p_hard(4,1) > 0.05
    txt = 'n.s.';
end
line( linspace(1.95,3.95), linspace(4.8,4.8),'color','k')
hold on
text(mean([3 3]), 5 ,txt,'FontSize',20)
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');