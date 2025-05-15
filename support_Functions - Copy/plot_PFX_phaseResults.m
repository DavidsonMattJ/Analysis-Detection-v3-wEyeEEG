
function plot_PFX_phaseResults(cfg, testData, PFX_FourierNull);

nGaits_toPlot=2;
pc=1; % plot counter
barCols= {'b', 'r', 'm','k'};

% both this and the next use the same figure function:
iLR=3;
gaitfield = {'gc', 'doubgc'};
gaitprint = {'gait', 'stride'};
binfield = {'','_binned'};
usebin=1
legp=[]; % for legend
fntsize=12;

%%
%set up figure:

hDataAll=[];
figure(1); clf;
set(gcf, 'color', 'w', 'units','normalized', 'position', [.1 .1 .6 .8])

for id=1:length(cfg.DV)

    dataIN= testData{id};
    %which field of datastructure to plot?
    if strcmp(cfg.DV{id}, 'RT')
        usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];
    elseif strcmp(cfg.DV{id}, 'Accuracy')
        usefield = [gaitfield{nGaits_toPlot}  binfield{usebin+1} '_Acc'];

    elseif strcmp(cfg.DV{id}, 'Counts')
        usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_counts'];
            if strcmp(cfg.type{id},'Saccade')
        usefield = [gaitfield{nGaits_toPlot} '_sacc_all' binfield{usebin+1} '_counts'];
            end
    end
    %%
    ppantData_ObsFITS=[]; % ppant Data is the observed per freq. vs shuffled
    ppantBounds=[]; % shuffle bounds
    ppantNHST=[]; % result comparing observed to shuffled.
    ppantDV_Data=[];
    ppantSINE=[]; % perform sinefits per sig frequency (below).

    for ippant = 1:length(PFX_FourierNull)

        ppantData_ObsFITS(ippant,:) = PFX_FourierNull(ippant).([cfg.type{id} 'Ons_' usefield '_fitsRsq_Obs']);

        ppantBounds(ippant,:) =PFX_FourierNull(ippant).([cfg.type{id} 'Ons_' usefield '_fitsRsq_ShuffCV'])(3,:); % 3rd pos is 95%

        ppantDV_Data(ippant,:)= dataIN(ippant,iLR).(usefield);


        % binary NHST:
        ppantNHST(ippant,:) = ppantData_ObsFITS(ippant,:) > ppantBounds(ippant,:);



    end % ippant
    %%
    nsubs = size(ppantNHST,1);
    NHST = sum(ppantNHST,1);

    testF = [2,4]; % cpg.
    hzind = dsearchn(cfg.Hzspace', testF');
    %%
    xvec=linspace(1,100,41);
    xvec=xvec(1:40);
    %
%     clf
    % perform fit per ppant at this freq,
    for ifreq=1:2;
        %find the index in Hz space:

    subplot(length(cfg.DV),2,(ifreq+ 2*(id-1)));
        % now force a fit at this frequency, to the participant data to extract the
        % plot:

        % Declaring the type of fit.
        FitType = 'fourier1';
        % Creating and showing a table array to specify bounds.
        CoeffNames_F = coeffnames(fittype(FitType));
        %set bounds for w
        CoeffBounds_F = array2table([-Inf(1,length(CoeffNames_F));...
            Inf(1,length(CoeffNames_F))],'RowNames',...
            ["lower bound", "upper bound"],'VariableNames',CoeffNames_F);

        %specify range for fit, from 0 to 10 Hz:
        %
        %     testw = 2*pi*10/100;
        %
        %     CoeffBounds.w(1) = 0;
        %     CoeffBounds.w(2) = testw;

        testw = 2*pi*cfg.Hzspace(hzind(ifreq))/xvec(end);

        CoeffBounds_F.w(1) = testw;
        CoeffBounds_F.w(2) = testw;
        %update fit opts settings

        %Update Fit Options setting.
        FitOpts_F = fitoptions('Method','NonlinearLeastSquares','Lower',table2array(CoeffBounds_F(1,:)),...
            'Upper',table2array(CoeffBounds_F(2,:)));

% perform same for sine model:
 FitType = 'sin1';
        % Creating and showing a table array to specify bounds.
        CoeffNames_S = coeffnames(fittype(FitType));
        %set bounds for w
        CoeffBounds_S = array2table([-Inf(1,length(CoeffNames_S));...
            Inf(1,length(CoeffNames_S))],'RowNames',...
            ["lower bound", "upper bound"],'VariableNames',CoeffNames_S);

      

        testw = 2*pi*cfg.Hzspace(hzind(ifreq))/xvec(end);

        CoeffBounds_S.b1(1) = testw;
        CoeffBounds_S.b1(2) = testw;
        %update fit opts settings

        %Update Fit Options setting.
        FitOpts_S = fitoptions('Method','NonlinearLeastSquares','Lower',table2array(CoeffBounds_S(1,:)),...
            'Upper',table2array(CoeffBounds_S(2,:)));

allFs= nan(nsubs,1040);
allTheta= nan(nsubs,1);
allRho = nan(nsubs,1);
allTheta_ns= nan(nsubs,1);
        %% restrict fit to shorter range(?)
%         clf
        for ippant= 1:size(ppantDV_Data,1)


            % may need to replacenan.
            tmpD= ppantDV_Data(ippant,:);
            if any(isnan(tmpD))
                disp(['warning: replacing nan-> 0 for ippant, ' num2str(ippant)]);
            
                tmpD(isnan(tmpD))=0;
            end
            [f,gof]= fit(xvec',tmpD',  'fourier1',FitOpts_F);
          
            % plot if this ppant was sig.
            if ppantNHST(ippant,hzind(ifreq))
%                 [f,gof]= fit(xvec',ppantDV_Data(ippant,:)',  'sin1', FitOpts_S);

          
                hold on;

          
                h=plot(f, xvec, tmpD');%,
                h(1).Visible ='off'; % hide data points.
                h(2).Color = barCols{id};
%                 h(2).


                allFs(ippant,:) = h(2).YData;

                % the phase can be expressed as the arctangend of b1/a1,
                % the ratio of coefficients of the sine and cosine terms. 
                allTheta(ippant) = f.b1/f.a1;

                % RHO is the square rt of the sum of squares of the
                % coefficients of the sine and cosine terms.
                %rho = sqrt(a1^2 + b1^2)
                 allRho(ippant) = sqrt(f.a1^2 + f.b1^2);
%                 legend off
            else % store non sig?
                allTheta_ns(ippant)= sqrt(f.a1^2 + f.b1^2);
            end


        end % each subject
        hold on;
        %%
            legend([h(2)], [ num2str(sum(ppantNHST(:,hzind(ifreq)))) '/36'], 'AutoUpdate','off');
        
propSig = sum(ppantNHST(:,hzind(ifreq)));
        %% add average of fits:

        fitData = imresize(allFs, [nsubs, 100]);
        mD= squeeze(nanmean(fitData,1));
        errD = CousineauSEM(fitData);
        sh= shadedErrorBar(1:100, mD, errD, {'k', 'linew',2}, 1);
% tidy axes:
if id==3
xlabel('% gait completion');
else
    xlabel('')
end
ylabel(cfg.DV{id});

set(gca,'fontsize',14)

% ax= get(gca,'axes');

       %%
%        clf
%        subplot(2,2,ifreq);
%make space for phase plot:
if ifreq==1
ylts= get(gca, 'ylim');
end
ylim([ylts(1), (ylts(1) + diff(ylts)*1.3)]);

prvpos = get(gca,'position');
%
innerax = axes('Position', [prvpos(1)-.01 prvpos(2)+prvpos(4)*.7, prvpos(3)*.4 prvpos(4)*.4], 'units', 'normalized');
%%
       ph=polarhistogram(allTheta,20)
       ph.FaceColor= barCols{id};
       ph.FaceAlpha= .6;
       hold on;

thetaticks([0 90 180 270, 360])
thetaticklabels({'0','\pi/2','\pi','3\pi/2' })
ax=gca;
ax.RAxisLocation= 180;
% test for non-uniformity:
[p, r] = circ_rtest(allTheta(~isnan(allTheta)));
rl= rlim;
if p<.05
text(1, rl(2)*1.2, ['\itp\rm = ' sprintf('%.5f', p)]);
else
text(1, rl(2)*1.2, '\itns');
end
shg
set(gca,'fontsize',10);
       %%
    end % ifreq



    %%

end % id

cd(cfg.figdir);
cd('Participant phase results');
print('-dpng', ['DVs combined phase distribution']);




end % function