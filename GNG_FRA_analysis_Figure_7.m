%% clearvars -except GNG_rec_all_cell GNG_analysis
clear all
close all
clc
%% load GNG_rec_all_cell & Fr_array
%put in your directory
addpath(genpath('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_MATLABR2023b_scripts'))

% select recording sessions
[file, path] = uigetfile('Z:\Shared\Benne\Praegel_et_al_2024\Praegel_et_al_data\'...
    , 'Select GNG_rec_all_cell ');
cd(path)
load (file)

%% Parameters
% Area and Layer Parameters
area_str ={'AUDd','AUDp','AUDv','TEa'}; % all subregions included in the analysis
areas = 1:length(area_str) ;
cortex = ({'AUDd6a', 'AUDp6a', 'AUDv6a', 'TEa6a',...
    'AUDd5', 'AUDp5', 'AUDv5', 'TEa5'});
table_cat_areas = categorical({'AUDd'; 'AUDp'; 'AUDv'; 'TEa';})  ;
group_str = {'adol.','adult'};

% time parameters
baseline = [-0.2 -0.05]; % baseline according to stimulus onset
start_baseline = 1; % baseline onset ms
stop_baseline = dist(baseline(1), baseline(2))* 1000; % baseline  ms
startStim = dist(baseline(1), 0)* 1000 ; % stimulus onset
stopStim = dist(baseline(1), 0.1)* 1000 ; % stimulus offset + 20 ms

window = [-0.2; 0.6]; % time window according to stimulus onset
startRange = window(1); %beginning of baseline
stopRange = window(2); %last timepoint before the reward

% minimal requirements of Rec
min_spikes = 100;

% iteration parameters
binSize = 0.001; % bin = 1 ms
rasterScale = 1; % binsize of raster plot
smoothSize = 5; % binsize of PSTh smoothing

% FRA parameters
oct_freqs = log2 (40 / 4)./ 20 ;
rep_num = 80 ;
attenuation = 20 ; 
chosen_attenutation = 2 ;
easy_go = 7.07 ;
learned_stim = 7:11 ;


% Color parameters
Colors_area = {[.1 .3 .8], [.5 .4 .9], [0 .5 .6],[0.4940 0.1840 0.5560]};
colors_ado_adu = { [.5 .7 .2],[ .2 .4 .2],[.5 .7 .2],[ .2 .4 .2]} ; 
alpha_f = [1:-0.05:0.05];

L = {['--'], ['-']};
M = {['o'],['o']};


%%  extracting auditory excitatory units in recordings that were followed by an FRA protocol
tic
% extract all neurons modulated by pure tones in highest attenuation 
for g   = 1:numel(GNG_rec_all_cell) % run per group
        Recs =1:numel(GNG_rec_all_cell{1,g});
    GNG_analysis.clusterIDs_all_FRA {1,g} = nan (300,length(Recs), length(areas)) ;
    GNG_analysis.clusterIDs_mod_FRA{1,g} = nan(300,length(Recs), length(areas)) 
    for i = 1:length(Recs) % run per recording
        if ~isempty(GNG_rec_all_cell{1, g}(i).Protocols)

            event_freq{g,i}(1,:) = GNG_rec_all_cell{1, g}(i).Protocols(1).max_time + ...
                GNG_rec_all_cell{1, g}(i).Protocols(2).trig_starts  ; % time (ms)
            event_freq{g,i}(2,:) = GNG_rec_all_cell{1, g}(i).FRA.perStimAtten(1,:) ;  % frequency
            event_freq{g,i}(3,:) = GNG_rec_all_cell{1, g}(i).FRA.perStimAtten(2,:) ; % attentuation
            unique_freq = unique(event_freq{g,i}(2,:)) ;
            unique_att = unique(event_freq{g,i}(3,:)) ;
            category_boundary = interp1(unique_freq,[1:20],10) ;
            spikes = GNG_rec_all_cell{1,g}(i).Units; %extract spikes
            spikeTimes = spikes.st; %extract spike times
            clusters = spikes.clu; %extract cluster ids per spike
            clusterIDs = unique(GNG_rec_all_cell{1,g}(i).Units.Unit_table.clu); % extract the unique clusterIDs

            idx_att_highest = find(event_freq{g,i}(3,:) == attenuation) ;
            idx_freq_chosen_attenuation = event_freq{g,i}(2,idx_att_highest) ; 
            eventTimes_all =  event_freq{g,i}(1,idx_att_highest) ;
            eventTimes_all_freq = event_freq{g,i}(1,:) ;
            
            % run this for loop to determine significantly modulated units
            % across all events.
            for area = 1:length(areas)
                
            % choose per area and layer 
                    all_clusters = find(strcmp(GNG_rec_all_cell{1,g}(i).Units.Unit_table.area_acronym,char(cortex(area)))...
                        |strcmp(GNG_rec_all_cell{1,g}(i).Units.Unit_table.area_acronym,char(cortex(area+4)))) ;
                
               
                clusterIDs = GNG_rec_all_cell{1,g}(i).Units.Unit_table.clu(all_clusters) ;
                clusterIDs =  clusterIDs (clusterIDs >0) ;

                GNG_analysis.clusterIDs_all_FRA{1,g}(1:length(clusterIDs),i,area) = clusterIDs ;
                
                for c = 1:length(clusterIDs)  %run per cluster = per neuron

                    
                    st = spikeTimes(clusters == clusterIDs(c));
                    %sort out spike times outside of the event time range
                    spikeTime = st(st>min(eventTimes_all+startRange) & st<max(eventTimes_all+stopRange));
                    
                    if length(spikeTime) > min_spikes
                        
                        %extract the basics
                        [~, psth,binArray,binCenters] = GNG_binning_PSTH(startRange,stopRange, binSize, eventTimes_all, spikeTime);
                        
                        % extract the smoothed basics
                        [~,psth_Smoothed, fr_Smoothed,spikeCounts] ...
                            = GNG_smoothed_PSTH (startRange, stopRange, binCenters,binArray, binSize, smoothSize);
                        
                        baseline_psth = psth(start_baseline:stop_baseline);  % baseline
                        onset_psth = psth(startStim +1 :stopStim + 50) ; % onset + 50 ms after tone offset to account for delayed firing
                        
                        [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','right');
                        if  p < 0.05 % sort out non-excited units
                            GNG_analysis.clusterIDs_mod_FRA{1,g}(c,i,area) = clusterIDs(c);
                            GNG_analysis.psth_Smoothed_mod_FRA{1,g}{c,i,area} = psth_Smoothed ;
                             GNG_analysis.spiketimes_mod_FRA{1,g}{c,i,area} = st ;
                            GNG_analysis.et_mod_FRA{1,g}{c,i,area} = event_freq{g,i} ;
                            GNG_analysis.frArray_mod_FRA_all{1,g}{c,i,area} = fr_Smoothed  ;
                        GNG_analysis.binArray_mod_FRA_all{1,g}{c,i,area} = binArray ;
                        GNG_analysis.idx_freq{1,g}{c,i,area} = idx_freq_chosen_attenuation ;

                        idx_cell{1,g}{c,i,area} = i ; 

                        end
                    end
                end
            end
        end
    end
end
toc
%% create table of excited units in expert and naive passive listening FRA protocol 

for g   = 1:numel(GNG_rec_all_cell) % run per group
    
    Recs = 1:numel(GNG_rec_all_cell{1,g}) ; % number of recordings
   
    for i = 1:length(Recs) % run per recording
        if ~isempty(GNG_rec_all_cell{1, g}(i).Protocols)
            
            for area = 1:length(areas)
                idx_all_area = ~isnan(GNG_analysis.clusterIDs_all_FRA{1,g}(:,i,area))  ;
                sum_all_area(1,g,i,area) = size(idx_all_area(idx_all_area == 1),1) ;
            end
            
            idx_all = ~isnan(GNG_analysis.clusterIDs_all_FRA{1,g}(:,i,:)) ;
            sum_all(g,i) = size(idx_all(idx_all == 1),1) ;
            
            for area = 1:length(areas)
                idx = find(GNG_analysis.clusterIDs_mod_FRA{1,g}(:,i,area)== 0) ;
                GNG_analysis.clusterIDs_mod_FRA{1,g}(idx,i,area) = nan ;
                idx_mod_area = ~isnan(GNG_analysis.clusterIDs_mod_FRA{1,g}(:,i,area))  ;
                sum_mod_area(1,g,i,area) = size(idx_mod_area(idx_mod_area == 1),1) ;
                idx = [] ;
            end
            idx_mod = ~isnan(GNG_analysis.clusterIDs_mod_FRA{1,g}(:,i,:))  ;
            sum_mod(1,g,i) = size(idx_mod(idx_mod == 1),1) ;
        end
    end
end

 for g   = 1:numel(GNG_rec_all_cell) % run per group
     format short
     sum_mod_recs(g,:) = sum(squeeze(sum_mod(1,g,:))) ; % all modulated neurons in ACx
     sum_all_g(g,:) = sum(sum_all(g,:)); % all neurons in ACx
     frac_mod_all(g,:) =   round(sum_mod_recs(g,:)/ sum_all_g(g,:),2); % fraction of modulated units
     
     sum_area(g,:,:) = squeeze(sum_all_area(1,g,:,:)) ; % per area all neurons
     sum_mod_area_all(g,:,:) =  squeeze(sum_mod_area(1,g,:,:)) ; % per area only modulated neurons
     
     for area = 1:length(areas)
         sum_area_recs(area,g) = sum(squeeze(sum_area(g,:,area))); % all units per area
         sum_mod_area_recs(area,g) = sum(squeeze(sum_mod_area_all(g,:,area))); % all modulated units per area
         frac_mod_area(area,g) =   round(sum_mod_area_recs(area,g)/ sum_area_recs(area,g),2) ; % fraction of modulated units per area
     end
 end

 for g = [1 3]
GNG_SU_area_table = table(table_cat_areas,sum_area_recs(:,g), sum_mod_area_recs(:,g), frac_mod_area(:,g),...
    sum_area_recs(:,g+1), sum_mod_area_recs(:,g+1), frac_mod_area(:,g+1),'VariableNames',...
    {'areas','neurons adol.','exc. neurons adol.','% adol.',...
        'neurons adults','exc. neurons adults','% adults'}) ;
writetable(GNG_SU_area_table,'GNG_ecx_area_table_FRA.xlsx','Sheet',g) ;
 end 
 %% extract the FR per frequency and attenuation 
close all
clc
tic
for g = 1:numel(GNG_rec_all_cell)

  spiketimes_all_cell{1,g}  = GNG_analysis.spiketimes_mod_FRA{1,g}(~cellfun('isempty', GNG_analysis.spiketimes_mod_FRA{1,g})) ;
    et_all_cell{1,g} = GNG_analysis.et_mod_FRA {1,g}(~cellfun('isempty', GNG_analysis.et_mod_FRA{1,g})) ;
     idx_all_cell{1,g} = idx_cell{1,g}(~cellfun('isempty', idx_cell{1,g})); 


    for c = 1:length(spiketimes_all_cell{1,g})

        spiketime =  spiketimes_all_cell{1,g}{c,1} ;


        for f = 1:size(unique_freq,2) % per frequency

            eventTimes_freq = [] ;
            idx_freq = find (et_all_cell{1,g}{c,1}(2,:) == unique_freq(f)) ;
            eventTimes_freq = et_all_cell{1,g}{c,1}(1,idx_freq) ; % sort events according to frequency
            spiketime_freq = [] ;

            spiketime_freq = spiketime(spiketime>min(eventTimes_freq+startRange) ...
                & spiketime<max(eventTimes_freq+stopRange)) ; % sort spieks per frequency

            if size(spiketime_freq,1) > 1 % if you have spikes

                [~,psth, binArray, binCenters] = GNG_binning_PSTH...
                    (startRange,stopRange, binSize, eventTimes_freq, spiketime_freq);

                [~,psth_Smoothed,frArray_Smoothed,spikeCounts] = GNG_smoothed_PSTH ...
                    (startRange, stopRange, binCenters, binArray, binSize, smoothSize);

                baseline_psth = psth(start_baseline:stop_baseline);  % baseline
                onset_psth = psth(startStim +1 :stopStim + 50) ; % onset + 50 ms after tone offset to account for delayed firing
                [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','right');
                if p <= 0.05
                    sig_freq {1,g}(c,f) = 1 ; 
                    sig_freq_FR{1,g}(c,f) = mean(onset_psth-baseline_psth); 

                elseif p > 0.05
                    sig_freq {1,g}(c,f) = 0 ; 
                    sig_frqe_FR{1,g}(c,f) = nan; 
                end


                mean_fr_base{1,g}(c,f) = mean (psth_Smoothed(1,1:150)) ;
                mean_fr{1,g}(c,f) = mean (psth_Smoothed(1,201:350)) ;
                psth_freq{1,g}{c,f} = psth_Smoothed ;
                frArray_freq{1,g}{c,f} = frArray_Smoothed ;
                binArray_freq{1,g}{c,f} = binArray ;

            elseif size(spiketime_freq,1) < 1
                mean_fr_base{1,g}(c,f) = 0;
                mean_fr{1,g}(c,f) = 0;
                psth_freq{1,g}{c,f} = nan(1,800);
                frArray_freq{1,g}{c,f}  = nan(80,800);
                binArray_freq{1,g}{c,f}  = nan(80,800);
            end

            for a = 1:size(unique_att,2)
                eventTimes_att = [] ;
                idx_att = find (et_all_cell{1,g}{c,1}(3,:) == unique_att(a) ) ;

                idx_comb = ismember (idx_freq,idx_att);
                idx_freq_att = idx_freq(idx_comb == 1) ;

                eventTimes_freq_att = et_all_cell{1,g}{c,1}(1,idx_freq_att) ;

                spiketime_freq_att = [] ;

                spiketime_freq_att = spiketime(spiketime>min(eventTimes_freq_att+startRange) ...
                    & spiketime<max(eventTimes_freq_att+stopRange)) ;

                if size(spiketime_freq_att,1) > 1

                    [~, psth, binArray, binCenters] = GNG_binning_PSTH...
                        (startRange,stopRange, binSize, eventTimes_freq_att, spiketime_freq_att);

                    [~,psth_Smoothed,frArray_Smoothed,spikeCounts] = GNG_smoothed_PSTH ...
                        (startRange, stopRange, binCenters, binArray, binSize, smoothSize);

                    baseline_psth = psth(start_baseline:stop_baseline);  % baseline
                    onset_psth = psth(startStim +1 :stopStim + 50) ; % onset + 50 ms after tone offset to account for delayed firing

                    [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','right');
                     if p <= 0.05
                    sig_freq_att{1,g}(c,a,f) = 1 ; 
                    sig_fre_att_FR{1,g}(c,a,f) = mean(onset_psth-baseline_psth); 
                elseif p > 0.05
                    sig_freq_att{1,g}(c,a,f) = 0 ; 
                    sig_fre_att_FR{1,g}(c,a,f) = nan; 

                end

                    mean_fr_att_base{1,g}(c,a,f) = mean (psth_Smoothed(1,1:150)) ;
                    mean_fr_att{1,g}(c,a,f) = mean (psth_Smoothed(1,201:350)) ;
                    psth_freq_att{1,g}{c,a,f} = psth_Smoothed ;
                    frArray_freq_att{1,g}{c,a,f} = frArray_Smoothed ;
                    binArray_freq_att{1,g}{c,a,f} = binArray ;

                elseif size(spiketime_freq_att,1) < 1
                    mean_fr_att_base{1,g}(c,a,f) = 0;
                    mean_fr_att{1,g}(c,a,f) = 0;
                    psth_freq_att{1,g}{c,a,f} = nan(1,800);
                    frArray_freq_att{1,g}{c,a,f} = nan(16,800);
                    binArray_freq_att{1,g}{c,a,f} = nan(16,800);


                end
            end
        end

        FRA_max_att = [] ;
       % FRA_max_att =  mean(squeeze(mean_fr_att{1,g}(c,:,:))) ;
                FRA_max_att =  (squeeze(mean_fr_att{1,g}(c,2,:))) ;

        max_FR = max(max(FRA_max_att)) ;

        BF = find(max_FR == FRA_max_att) ;
        max_FR_all{1,g}(c,1) = max_FR (1) ;
        BF_all{1,g}(c,1) = BF(1) ;
    end
end
toc

%% Figure 7 B C D plot  FRAs and PSTHS per frequency per unit
close all
clc

for g =1: numel(GNG_rec_all_cell)-2 % expert examples
    if g == 1
        c = 4 ;% adolescent example
    elseif g == 2
        c = 26; % adult example neuron
    end


% Figure 7B PSTH
figure
for f = 1:size(unique_freq,2)

    plot( psth_freq{1,g}{c,f} ,'Color',[colors_ado_adu{g} alpha_f(f)])
    hold on
end
xlabel('time(ms)')
xline(200,'--','linewidth',2)
xticks([0:200:800])
xticklabels([-200:200:600])
ylabel('Fr(Hz)')
box off
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
            
    % Figure 7C RASTERPLOT
    rasterScale = 5
    for f = 1:length(unique_freq)
        [rasterX{f} rasterY{f}] = GNG_raster (binArray_freq{1,g}{c,f}, binSize, binCenters, rasterScale);
        if f>1
        rasterY{f} = rasterY{f} + (80*f) ;
        end 
    end
    rasterX_mat =  cell2mat(rasterX) ;
    rasterY_mat =  cell2mat(rasterY) ;

    figure
    plot(rasterX_mat,rasterY_mat,'Color',[0.3 0.3 0.3])
    yticks([80 category_boundary*80 1600])
    yticklabels([4 10 40])
    ylabel('# trials per frequency ')
    xticks([-0.2:0.2:0.6])
    xticklabels(-200:200:600)
    yline(category_boundary* 80,'--k','linewidth',2)
    xlabel('time(ms)')
    box off
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');
    set(gca,'LineWidth',1.5,'TickLength',[0.025 0.025]);


figure
    imagesc(squeeze(mean_fr_att{1,g}(c,:,:)) - squeeze( mean_fr_att_base{1,g}(c,:,:)))
hold on
xlabel('frequency')
xticks([1 category_boundary 20])
yticks([1:4:5])
xline(category_boundary,'--w','linewidth',3)
xticklabels([4 10 40])
ylabel('attenuation')
box off;
hold on;
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
end




 %%  Calculation of tuning bandwidth and  Population responsiveness
clc
tic
max_size = max([size(BF_all{1,1} ,1) size(BF_all{1,2},1) size(BF_all{1,3} ,1) size(BF_all{1,4} ,1)]) ;

for g = 1:numel(GNG_rec_all_cell)
     sum_mod_freq_FRA{1,g} = nan(max_size,1) ;
     oct_mod_freq_FRA{1,g} = nan(max_size,1) ; 
    for c = 1:size(frArray_freq{1,g},1 )      
        for f = 1:size(unique_freq,2) % per frequency
            baseline_psth = mean(frArray_freq{1,g}{c,f}(:,start_baseline:stop_baseline));  % baseline
            onset_psth = mean( frArray_freq{1,g}{c,f}(:,startStim +1 :stopStim + 50)) ; % onset + 50 ms after tone offset to account for delayed firing
            
            [p,~,~] = ranksum(onset_psth , baseline_psth,0.05,'tail','right');
            if  p < 0.05 % sort out non-excited units
                mod_freq_FRA{1,g}(c,f) = 1;
                
            end
        end
      sum_mod_freq_FRA{1,g}(c,1) = sum(mod_freq_FRA{1,g}(c,:));
        oct_mod_freq_FRA{1,g}(c,1) =  (sum_mod_freq_FRA{1,g}(c,1) - (length(unique_freq)/5)) * oct_freqs ;

    end
          sum_mod_freq_FRA{1,g} = (sum_mod_freq_FRA{1,g} - nanmin(sum_mod_freq_FRA{1,g})) / ( nanmax(sum_mod_freq_FRA{1,g}) - min(sum_mod_freq_FRA{1,g}) ) ;
end

%% Figure 7E tuning bandwidth in oct. 
clc
toc
close all

figure
violinplot([oct_mod_freq_FRA{1,3} oct_mod_freq_FRA{1,1} oct_mod_freq_FRA{1,4} oct_mod_freq_FRA{1,2}]...
    ,{'adolescent naive','adult naive','adolescent expert','adult expert'},"ViolinColor",{[.5 .7 .2; .5 .7 .2; .2 .4 .2; .2 .4 .2]} )
xlim([0 5])
ylim([0 4])
yticks([0:0.5:3])
ylabel('tuning bandwidth (oct.)')
    box off
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');

 [~,~,stats] = kruskalwallis([ oct_mod_freq_FRA{1,3} oct_mod_freq_FRA{1,4} oct_mod_freq_FRA{1,1} oct_mod_freq_FRA{1,2}],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6)

line_left = [1.05 3.05 1.05 1.95]
line_right =  [1.95 3.95 3.05 3.95]
line_height = [3.1 3.1 3.3 3.5]
text_mean = [1.5 3.5 2 3 ]
text_height = [3.2 3.2 3.4 3.6] 

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'n.s.';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end
   %% Figure 7F population responsiveness
toc

figure
violinplot([sum_mod_freq_FRA{1,3} sum_mod_freq_FRA{1,1} sum_mod_freq_FRA{1,4} sum_mod_freq_FRA{1,2}],...
    {'adolescent naive','adolescent expert','adult naive','adult expert'},"ViolinColor",{[.5 .7 .2; .5 .7 .2; .2 .4 .2; .2 .4 .2]} ) ;
xlim([0 5])
ylim([0 1.6])
yticks([0:0.25:1])
ylabel('population responsiveness')
    box off
    hold on;
    ax = gca;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    movegui('east');

 [~,~,stats] = kruskalwallis([ sum_mod_freq_FRA{1,3} sum_mod_freq_FRA{1,4} sum_mod_freq_FRA{1,1} sum_mod_freq_FRA{1,2}],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6) ;

line_left = [1.05 3.05 1.05 1.95] ;
line_right =  [1.95 3.95 3.05 3.95] ;
line_height = [1.1 1.1 1.3 1.5] ;
text_mean = [1.5 3.5 2 3 ] ;
text_height = [1.2 1.2 1.4 1.6] ;

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'n.s.';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end

  %% distance of BF from Go and No Go

close all

max_size = max([size(BF_all{1,1} ,1) size(BF_all{1,2},1) size(BF_all{1,3} ,1) size(BF_all{1,4} ,1)]) ;

for g   = 1:numel(GNG_rec_all_cell) % run per group
     
    dist_easy_go{g,1} = nan (max_size,1) ;

    for f = 1:length(unique_freq)
        idx_f = find(BF_all{1,g}(:,1) == f) ;
        freq_BF_all{g,1}(idx_f,1) = unique_freq(f) ;
    end
    %dist per octave
    dist_easy_go{g,1}(1:size(freq_BF_all{g,1},1)) = abs(log2(abs(freq_BF_all{g,1} / easy_go))) ;

end 
   %% Figure 7G distance of BF to Go stimuli
   clc
figure
violinplot([dist_easy_go{3,1} dist_easy_go{1,1} dist_easy_go{4,1} dist_easy_go{2,1}],...
    {'adolescent naive','adolescent expert','adult naive','adult expert'},"ViolinColor",{[.5 .7 .2; .5 .7 .2; .2 .4 .2; .2 .4 .2]} ) ;
xlim([0 5])
%ylim([0 3.5])
%yticks([0:1:3])
ylabel('BF to Go (oct.)')
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;

 [~,~,stats] = kruskalwallis([ dist_easy_go{3,1} dist_easy_go{4,1} dist_easy_go{1,1} dist_easy_go{2,1}],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6);

line_left = [1.05 3.05 1.05 1.95] ;
line_right =  [1.95 3.95 3.05 3.95] ;
line_height = [2.8 2.8 3 3.2] ;
text_mean = [1.5 3.5 2 3 ] ;
text_height = [2.9 2.9 3.1 3.5]  ;

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'n.s.';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end

%% Figure 7H neuronal dprime of the per single neuron 

clc
tic
for g   = 1:numel(GNG_rec_all_cell) % run per group
    for c = 1:length (BF_all{1,g}) %run per cluster = per neuron
        for f = 1:size(unique_freq,2) % per frequency
            avg_resp_f = [] ; 
            resp_f = [] ; 
            
            for f2 = 1:size(unique_freq,2) % per compared frequency
                resp_f2 = [] ; 
                avg_resp_f2 = [] ; 
                distance_ff2 = [] ; 

                avg_resp_f = mean_fr_att{1,g}(c,chosen_attenutation,f)-  mean_fr_att_base{1,g}(c,chosen_attenutation,f);

                avg_resp_f2 = mean_fr_att{1,g}(c,chosen_attenutation,f2) - mean_fr_att_base{1,g}(c,chosen_attenutation,f2) ;
                
            distance_ff2 = sqrt(sum((avg_resp_f-avg_resp_f2).^2));
            
            rep_num = size( frArray_freq_att{1,g}{c,chosen_attenutation,f2},1);
                
            avg_trial_f = [] ;
            avg_trial_f2 = [] ;
            for t = 1:rep_num 
                 avg_trial_f(t) = mean(frArray_freq_att{1,g}{c,chosen_attenutation,f}(t,startStim +1:stopStim  + 50)') ...
                    - mean(frArray_freq_att{1,g}{c,chosen_attenutation,f}(t,1:startStim - 50)');
                avg_trial_f2(t) = mean(frArray_freq_att{1,g}{c,chosen_attenutation,f2}(t,startStim +1:stopStim + 50 )') ...
                     - mean(frArray_freq_att{1,g}{c,chosen_attenutation,f2}(t,1:startStim - 50)');  

                                 resp_f(t)=sqrt(sum((avg_trial_f(t)-avg_resp_f).^2));
                                 resp_f2(t)=sqrt(sum((avg_trial_f2(t)-avg_resp_f2).^2));
            end 
            
            innerdist_f= mean(resp_f);
            innerdist_f2= mean(resp_f2); 
            
            d_prime{1,g}(f,f2,c)=distance_ff2/ mean([innerdist_f innerdist_f2]);
            
            if isnan(d_prime{1,g}(f,f2,c))
                d_prime{1,g}(f,f2,c) = 0 ;
            end
            if isinf(d_prime{1,g}(f,f2,c))
                d_prime{1,g}(f,f2,c) = 0 ;
            end
            
        end
        
        end
    end
        GNG_analysis.d_prime{1,g} = d_prime{1,g} ;
end 
toc 

 %% Figure 7H neuronal dprime for learned stim range

 for g   = 1:numel(GNG_rec_all_cell) % run per group
    mean_dprime_neuron{1,g} = nan(1,150) ;

    for c = 1:length (BF_all{1,g}) %run per cluster = per neuron
        d_prime{1,g}(d_prime{1,g} == 0) = nan ; 
            mean_dprime_neuron{1,g}(c)  = nanmean(nanmean(squeeze(d_prime{1,g}(learned_stim,learned_stim,c)))) ;

    end
end 

          figure
violinplot([mean_dprime_neuron{1,3}'  mean_dprime_neuron{1,1}' mean_dprime_neuron{1,4}'  mean_dprime_neuron{1,2}' ],...
    {'adolescent naive','adolescent expert','adult naive','adult expert'},"ViolinColor",{[.5 .7 .2; .5 .7 .2; .2 .4 .2; .2 .4 .2]} ) ;
xlim([0 5])
ylim([0 4.5])
yticks([0:1:4])
ylabel("freq. discrimination(d')")


ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
movegui('east');
box off;
 [~,~,stats] = kruskalwallis([ mean_dprime_neuron{1,3}' mean_dprime_neuron{1,4}' mean_dprime_neuron{1,1}' mean_dprime_neuron{1,2}'],[],'off');
   stat_multi = multcompare(stats,"Display","off");
   p =stat_multi(:,6)

line_left = [1.05 3.05 1.05 1.95] ;
line_right =  [1.95 3.95 3.05 3.95] ;
line_height = [3.8 3.8 4 4.2] ;
text_mean = [1.5 3.5 2 3 ] ;
text_height = [3.9 3.9 4.1 4.5]  ;

s1 = [2 5 1 6] ;
   for s = 1:length(s1)
       if  p(s1(s)) < 0.0005
           txt = '***';
       elseif  p(s1(s)) < 0.005
           txt = '**';
       elseif   p(s1(s)) < 0.05
           txt = '*';
       elseif p(s1(s)) > 0.05
           txt = 'n.s.';
       end

line( linspace(line_left(s),line_right(s)), linspace(line_height(s),line_height(s)),'color','k')
hold on
text(mean([text_mean(s) text_mean(s)]), text_height(s) ,txt,'FontSize',20)
   end
   
   
   %% Figure 7D BF overview
close all 
for g   = 1:numel(GNG_rec_all_cell) % run per group
figure(g);

[sorted_BF{1,g}, original_BF{1,g}] = sort(BF_all{1,g})


% Use tiled layout for better control over the layout
tiledlayout(numel(mean_fr_att{1,g})/100, 1, 'TileSpacing', 'none', 'Padding', 'none');  % 123 rows, 1 column   

for c = 1:length (BF_all{1,g})
    nexttile;
    imagesc(reshape(mean_fr_att{1,g}(original_BF{1,g}(c,1),2,:),[],20));  % Plot the heatmap for the row (1x20)
    set(gca, 'XTick', []);  % Remove X-axis ticks
    set(gca, 'YTick', []);  % Remove Y-axis ticks
        axis off;  % Turn off the axis to remove black outline
        xline(category_boundary,'w--','linewidth',3)
end



box off

end 
clc