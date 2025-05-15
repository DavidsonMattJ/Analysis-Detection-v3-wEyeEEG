% calc_cps_rawstride_correlations


% review figure, testing whether the differences in individual cps effects
% correlate with stride durations (in unit time).

% load the ppant level observed effects
% load the walk params.

% correlate stride duration with AMP at 2 cps,
% correlate stride duration with AMP at 4 cps,

% correlate stride duration with peak Freq
% correlate stride duration with peak Amp


%%
% cd('/Users/matthewdavidson/Documents/GitHub/Analysis-Detection-v3-wEyeEEG/data_Processed/GFX');
cd('C:\Users\mdav0285\Documents\GitHub\Analysis-Detection-v3-wEyeEEG\data_Processed\GFX');
load('GFX_Data_inGaits.mat','avWalkParams','GFX_TargPosData');
%%
load('GFX_Data_inGaits_FourierFits.mat','PFX_FourierNull', 'Hzspace');



clf
figure(1);
for ippant = 1:36
plot(Hzspace, PFX_FourierNull(ippant).TargetOns_doubgc_binned_Acc_fitsRsq_Obs, 'b');
hold on

end

%%
clf
DVs={'Acc', 'rts', 'counts'};
Onsets={'Target', 'Target', 'Response'};
range2= dsearchn(Hzspace', [1.8 2.2]');
range4= dsearchn(Hzspace', [3.8 4.2]');
titlesAre= {'Accuracy', 'RTs', 'Responses '}; 
ColsAre= {'b', 'r', 'm'};
AMP2=[];
AMP4=[];
peakF=[]; % this is the peak in Hz.

peakAmp=[];
for icond=1:3 % acc, rt, resp
figure(1)
    %all Amp at 2cps
    for ippant= 1:36
        AMP2(ippant) = max(PFX_FourierNull(ippant).([Onsets{icond} 'Ons_doubgc_binned_' DVs{icond} '_fitsRsq_Obs'])(range2));
        AMP4(ippant) = max(PFX_FourierNull(ippant).([Onsets{icond} 'Ons_doubgc_binned_' DVs{icond} '_fitsRsq_Obs'])(range4));
        
        [peakAmp(ippant), peakF(ippant)]= max(PFX_FourierNull(ippant).([Onsets{icond} 'Ons_doubgc_binned_' DVs{icond} '_fitsRsq_Obs']));
    
    end
    
    peakF = Hzspace(peakF);

    % as per rvwr 2, also plot peak in Hz vs stride duration:
    % cycles in duration to Hz:
    peakF_Hz = peakF'./avWalkParams(:,2);
    %plot all our correlations:
    subplot(4,3,icond);
    % plot 
    Xdata= avWalkParams(:,2); % this is the duration of each ppants stride.
    Ydata= {AMP2, AMP4, peakF, peakAmp};
    yLabels = {'R^2 @ 2 cps', 'R^2 @ 4 cps', 'Max frequency', 'Max R^2'};
    
    for iY = 1:3%
        subplot(3,3,icond + (3*(iY-1)));
        
        sc=scatter(Xdata, Ydata{iY},  'filled')
        sc.MarkerFaceColor= ColsAre{icond};
        
        
        
        % Fit a linear regression line
        coefficients = polyfit(Xdata, Ydata{iY}, 1);
        yfit = polyval(coefficients, Xdata);
        
        % Plot the regression line
        
        hold on;
        [r,p]= corrcoef(Xdata, Ydata{iY});
        
%         if p(1,2)<=.05
%         titletex= ['r= ' sprintf('%.3f', r(1,2)) ', p = ' sprintf('%.3f', p(1,2))];
%         plot(Xdata, yfit, 'r', 'LineWidth', 2);
%         else
            titletex= ['ns'];
        plot(Xdata, yfit, 'color',[.7 .7 .7], 'LineWidth', 2);
%         end
        
        title({[titlesAre{icond}];[titletex]})

        axis square
        box on
        ylabel(yLabels{iY})
        set(gca,'fontsize',12)
    end
    xlabel('Stride cycle duration (sec)')
    




figure(2); 
% first plot total scatteR:
subplot(3,3,1+3*(icond-1));
% now plot them all:
sc=scatter(Xdata,peakF_Hz);
sc.MarkerFaceColor= ColsAre{icond};
hold on
xlabel('Stride-cycle duration (second)');
ylabel('Best Fourier fit (Hz)')
% Fit a linear regression line
coefficients = polyfit(Xdata, peakF_Hz, 1);
yfit = polyval(coefficients, Xdata);
plot(Xdata, yfit, 'color',[.7 .7 .7], 'LineWidth', 2);
% Plot the regression line
hold on;
[r,p]= corrcoef(Xdata, peakF_Hz);

if p(1,2)<=.05
    titletex= ['r= ' sprintf('%.3f', r(1,2)) ', p = ' sprintf('%.3f', p(1,2))];
   pr= plot(Xdata, yfit, 'r', 'LineWidth', 2);
else
    titletex= ['ns'];
pr=plot(Xdata, yfit, 'color',[.7 .7 .7], 'LineWidth', 2);
end
box on
title({[titlesAre{icond} ' (N=36)']})
    legend(pr, titletex)

ylim([0 10])
xlim([1 1.4])
        set(gca,'fontsize',13)


%for those ppants with peaks ~ 2 cps, and 4 cps 
% first , find the ppants: 
ranges= [1.2,2.5; 2.5 4.5]

ts= {'~2 Hz cluster', '~4 Hz cluster'};
for irange= 1:2
    
a= find(peakF>ranges(irange,1));
b= find(peakF<ranges(irange,2));
c= intersect(a,b);
subplot(3,3, 2+ 3*(icond-1)+ (irange-1));
sc=scatter(Xdata(c),peakF_Hz(c));
sc.MarkerFaceColor= ColsAre{icond};
hold on
xlabel('Stride-cycle duration (second)');
ylabel('Best Fourier fit (Hz)')
% Fit a linear regression line
coefficients = polyfit(Xdata(c), peakF_Hz(c), 1);
yfit = polyval(coefficients, Xdata(c));
% Plot the regression line
hold on;
[r,p]= corrcoef(Xdata(c), peakF_Hz(c));

if p(1,2)<=.05
    titletex= ['r= ' sprintf('%.3f', r(1,2)) ', p = ' sprintf('%.3f', p(1,2))];
   pr= plot(Xdata(c), yfit, 'k', 'LineWidth', 2);
else
    titletex= ['ns'];
    pr=plot(Xdata(c), yfit, 'color',[.7 .7 .7], 'LineWidth', 2);
end
        box on
        title({[titlesAre{icond} ', ' ts{irange}]});
        legend(pr, titletex)
        set(gca,'fontsize',13)
xlim([1 1.4]);
ylim(ranges(irange,:))
end % irange
%
end 




%% 