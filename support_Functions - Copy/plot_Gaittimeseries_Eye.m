function plot_Gaittimeseries_Eye(GFX_headY, GFX_EyeD)

%% helper function to plot the (cleaned) time series of head, eye dir, and eye pos data.

% called from j6__
%prepare data:





% set up figure.
figure(1); clf; set(gcf,'color', 'w', 'units', 'normalized', 'position',[.1 .1 .8 .8]);






%prepare head Y.
headY=[];
for ippant = 1:length(GFX_headY)
headY(ippant,:) = GFX_headY(ippant).gc;
end


% % plot then hold :

mD = squeeze(nanmean(headY,1));
stE= CousineauSEM(headY);
tsAre={'Gaze direction', 'Gaze origin'};
for isub=1:2
    subplot(1,2,isub)
% sh=shadedErrorBar(1:length(mD), mD, stE, {'color','k', 'linewidth', 2},1);
% hHead = sh.mainLine;
title(tsAre{isub});
xlabel('% step completion'); 
ylabel('Position change [m]')
hold on
ylim([-.05 .05])
set(gca,'fontsize', 15)
end
%
% now plot the remainin fields:
plotfields = {'gc_dirY', 'gc_dirZ', 'gc_posY', 'gc_posZ'};

hEye=[]; % legend handles.
hCols = {'b',' r', 'b', 'r'};
hSts = {'-', '-', ':', ':'};
subP= [1,1,2,2];

for iplot= 1:length(plotfields )

    plotData=[];
    for ippant = 1:length(GFX_headY)
        plotData(ippant,:) = GFX_EyeD(ippant).([plotfields{iplot}]);
    end

    mD = squeeze(nanmean(plotData,1));
    stE= CousineauSEM(plotData);

    if iplot==3
        mD= mD+max(mD)
    elseif iplot==4
        mD= mD-mean(mD)
    end
        subplot(1,2,subP(iplot))
%         yyaxis right
    sh=shadedErrorBar(1:length(mD), mD, stE, {'color',hCols{iplot}, 'linestyle', hSts{iplot}, 'linewidth', 2},1);
    hEye(iplot) = sh.mainLine;
   

end

legend([hEye(1), hEye(2)], {'Y', 'X'})

% % 
% legend([hHead, hEye(1),hEye(2),hEye(3),hEye(4)], {'head', 'dirY', 'dirX', 'posY', 'posX'})
shg






end