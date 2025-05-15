
filename = 'MD2_2022-08-09-11-45_trialsummary.csv';
% 
opts = detectImportOptions(filename,'NumHeaderLines',0);
% 
T1 =readtable(filename,opts);
%%
lastrow = find(T1.trial==127,1,'last');
alltrs = T1.trial;
alltrs(lastrow+1:end) = alltrs(lastrow+1:end) + 128;

T1.trial = alltrs;
%%
writetable(T1, filename);
%%
% T2 = readtable('MD2_2022-08-09-12-17_framebyframe.csv',opts);
% 
% %
% T3= [T1;T2];
% %%
% lastrow = size(T1,1);
% firstrow = lastrow+1;
% lasttrial = T1.trial(lastrow);
% 
% wrongtrials= T2.trial;
% newtrials = wrongtrials+lasttrial+1;
% 
% alltrials = [T1.trial; newtrials];
% T3.trial=alltrials;
% 
% %%
% writetable(T3, filename);
% %% 
% %same for summary.
% 
% 
% filename = 'MD2_2022-08-09-11-45_trialsummary.csv';
% 
% opts = detectImportOptions(filename,'NumHeaderLines',0);
% 
% T1 =readtable('MD2_2022-08-09-11-45_trialsummary.csv',opts);
% T2 = readtable('MD2_2022-08-09-12-17_trialsummary.csv',opts);
% 
% %
% T3= [T1;T2];
% %%
% lastrow = size(T1,1);
% firstrow = lastrow+1;
% lasttrial = T1.trial(lastrow);
% 
% wrongtrials= T2.trial;
% newtrials = wrongtrials+lasttrial+1;
% 
% alltrials = [T1.trial; newtrials];
% T3.trial=alltrials;
% 
% %%
% writetable(T3, filename);
%%
% filename= 'LL1_2022-07-25-02-11_framebyframe.csv';
% filename='LL1_2022-07-25-02-11_trialsummary.csv';
% filename = 'LL2_2022-08-25-11-03_framebyframe.csv';
filename = 'LL2_2022-08-25-11-03_trialsummary.csv';

opts = detectImportOptions(filename,'NumHeaderLines',0);
T= readtable(filename, opts);

T.participant = repmat('LL2', [size(T,1),1]);
%%

writetable(T,'LL2_2022-08-25-11-03_trialsummary.csv');