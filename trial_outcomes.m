function [hit, miss, fa, cr, go_licks, ngo_licks] = trial_outcomes (binsize, chosen_trials, mouse_table);

bining=zeros(1,length(chosen_trials));
binind=1:binsize:length(bining);
for b=1:length(binind)-1;
    bining(binind(b):end)=bining(binind(b):end)+1;
end
bining=bining';
    

U = unique(bining,'stable');

for h= 1:length (U)
    hit(h)= sum (bining(:,1) == U(h) & mouse_table.score(chosen_trials)==0);
    miss(h)=sum(bining(:,1)== U(h) & mouse_table.score(chosen_trials)==2);
    fa(h)=sum(bining(:,1)== U(h) & mouse_table.score(chosen_trials)==1);
    cr(h)=sum(bining(:,1)== U(h) & mouse_table.score(chosen_trials)==3);
    go_licks(h) =  hit(h) / (hit(h)+ miss(h));
    ngo_licks(h) = fa(h) / (fa(h)+ cr(h));
    
end
end