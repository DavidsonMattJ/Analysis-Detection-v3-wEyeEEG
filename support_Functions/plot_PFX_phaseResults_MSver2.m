
function plot_PFX_phaseResults_MSver2(cfg, testData, PFX_FourierNull);

nGaits_toPlot=2;
pc=1; % plot counter
barCols= {'b', 'r', 'm','k'};

% both this and the next use the same figure function:
iLR=3;
gaitfield = {'gc', 'doubgc'};
gaitprint = {'gait', 'stride'};
binfield = {'','_binned'};
usebin=1;
legp=[]; % for legend
fntsize=cfg.fntsize;
%%
%set up figure:

hDataAll=[];
figure(1); clf;
set(gcf, 'color', 'w', 'units','normalized', 'position', [.1 .1 .55 .8])

printy= {'Accuracy', 'Reaction times', 'Response counts','Saccade counts'};
for id=1:length(cfg.DV)

    dataIN= testData{id};
    %which field of datastructure to plot?
    if strcmp(cfg.DV{id}, 'RT')
        usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];
    elseif strcmp(cfg.DV{id}, 'Accuracy')
        usefield = [gaitfield{nGaits_toPlot}  binfield{usebin+1} '_Acc'];

    elseif strcmp(cfg.DV{id}, 'Counts')
               usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_counts'];

        if strcmp(cfg.type{id}, 'Saccade')
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

%     if id==2
        % subtract the mean from each ppants RT.
        meanper = nanmean(ppantDV_Data,2);
        divBy = repmat(meanper, [1, size(ppantDV_Data,2)]);
        %relative change:
        ppantDV_Data= (ppantDV_Data - divBy)./divBy;
%         ppantDV_Data=ppantDV_Data-1;
%     end

    nsubs = size(ppantNHST,1);
    NHST = sum(ppantNHST,1);

    testF = [1.5,2.5];
%         testF = [1.9,2.1];
        testF = [3.5,4.5]; 
    hzind = dsearchn(cfg.Hzspace', testF');
    %%
    xvec=linspace(1,100,41);
    xvec=xvec(1:40);
    %
%     clf
    % perform fit per ppant at this freq,
    
        %find the index in Hz space:

    subplot(length(cfg.DV),3,([1,2] +(3*(id-1))));
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

        testw1 = 2*pi*cfg.Hzspace(hzind(1))/xvec(end);
        testw2 = 2*pi*cfg.Hzspace(hzind(end))/xvec(end);

        CoeffBounds_F.w(1) = testw1;
        CoeffBounds_F.w(2) = testw2;
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

      

%         testw = 2*pi*cfg.Hzspace(hzind(ifreq))/xvec(end);

        CoeffBounds_S.b1(1) = testw1;
        CoeffBounds_S.b1(2) = testw2;
        %update fit opts settings

        %Update Fit Options setting.
        FitOpts_S = fitoptions('Method','NonlinearLeastSquares','Lower',table2array(CoeffBounds_S(1,:)),...
            'Upper',table2array(CoeffBounds_S(2,:)));

allFs= nan(nsubs,100);
allTheta= nan(nsubs,1);
allRho = nan(nsubs,1);
allTheta_index = nan(nsubs,1);
allTheta_ns= nan(nsubs,1);
        %% restrict fit to shorter range(?)
%         clf
        for ippant= 1:size(ppantDV_Data,1)


            tmpD=ppantDV_Data(ippant,:);
            if any(isnan(tmpD));
                disp(['warning: replacing nan -> 0 for ppant ' num2str(ippant)]);
                tmpD(isnan(tmpD))=0;
            end
            [f,gof]= fit(xvec',tmpD',  'fourier1',FitOpts_F);
          
            % plot if this ppant was sig.
            if any(ppantNHST(ippant,hzind(1):hzind(2)))
%                 [f,gof]= fit(xvec',ppantDV_Data(ippant,:)',  'sin1', FitOpts_S);

          
                hold on;

          
                h=plot(f, xvec,tmpD');%,
                h(1).Visible ='off'; % hide data points.
                h(2).Color = barCols{id};
                h(2).LineWidth=2;
                h(2).Color= [h(2).Color, .5];
               
%                 h(2).


                allFs(ippant,:) = imresize(h(2).YData,[1,100]);

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
           lg= legend([h(2)], [ num2str(length(find(sum(ppantNHST(:,hzind(1):hzind(2)),2)))) '/36'], 'AutoUpdate','off');
           lg.ItemTokenSize(1) = 8;
% legend off;

%% calculate phase (altnerate version).
%alternatively, calc phase at particular time:
nanp = find(isnan(mean(allFs,2)));
tmpF = allFs;
tmpF(nanp,:)=[];

signals=tmpF; %avoid bias (subzero amp.)
  
Pvals=[]; 
Zvals=[];
Varvals=[];
allTheta_index=[];
for index = 1:size(signals,2);
% index=
  % Calculate the complex representation and phase for each trial
complex_signals = hilbert(signals);
phase = unwrap(angle(complex_signals));

allTheta_index(:,index) = rad2deg(phase(:,index));
[Pvals(index), Zvals(index)] = circ_rtest(allTheta_index(:,index));
Varvals(index) = circ_var(allTheta_index(:,index));
end
% hold on; plot([index index], ylim, 'r-')



        %% add average of fits:

%         fitData = imresize(allFs, [nsubs, 100]);
%         mD= squeeze(nanmean(fitData,1));
%         errD = CousineauSEM(fitData);
%         sh= shadedErrorBar(1:100, mD, errD, {'k', 'linew',2}, 1);
% tidy axes:
if id<3
xlabel('Target onset (% stride-cycle)');
elseif id==3
    xlabel('Response onset (% stride-cycle)');
elseif id==4
    xlabel('Saccade onset (% stride-cycle)');

end
% ylabel([cfg.DV{id} ', (n='  num2str(sum(ppantNHST(:,hzind(ifreq)))) '/36)']);
ylabel({[printy{id}];['(relative change)']});

set(gca,'fontsize',cfg.fntsize)
if id<3
    ylim([-.15 .15]);
   
else
    ylim([-.4 .4]);
end
ys= get(gca,'ylim');
axis tight
[p,plotindex ]= max(Zvals); 
hold on;
% plot([plotindex-1, plotindex-1], [ys(1), ys(2)], 'k-');
% plot([plotindex+1, plotindex+1], [ys(1), ys(2)], 'k-');
% plot([plotindex-1, plotindex+1], [ys(1), ys(1)], 'k-')
% plot([plotindex-1, plotindex+1], [ys(2), ys(2)], 'k-')

%
xpch = [plotindex-1, plotindex-1, plotindex+1 plotindex+1];
ypch = [ys(1), ys(2), ys(2) ys(1)];
pch = patch(xpch, ypch, [.8 .8 .8],'FaceAlpha',.5);
subplot(length(cfg.DV),3,3+3*(id-1))
%%
cla
% [p,plotindex ]= max(Zvals); 
       ph=polarhistogram(allTheta_index(:,plotindex),20);
       ph.FaceColor= barCols{id};
       ph.FaceAlpha= .6;
       hold on;
%%
thetaticks([0 90 180 270, 360])
thetaticklabels({'0','\pi/2','\pi','3\pi/2' })
ax=gca;
ax.RAxisLocation= 180;
% test for non-uniformity:
[p, Z] = circ_rtest(allTheta_index(:,plotindex));
V = circ_var(allTheta_index(:,plotindex));
% rlim([0 5]);
rl= rlim;
disp(['DV ' num2str(id)]);
disp(['Z = ' sprintf('%.5f', Z)]);
 disp(['p = ' sprintf('%.5f', p)]);
  disp(['circ V = ' sprintf('%.5f', V)]);
if p<.05
%     text(0.65, rl(2)*1.3, ['Z\rm = ' sprintf('%.5f', Z)]);
%     text(0.5, rl(2)*1.3, ['\itp\rm = ' sprintf('%.5f', p)]);
% disp(['Z = ' sprintf('%.5f', Z)]);
%  disp(['p = ' sprintf('%.5f', p)]);
    
else
text(1, rl(2)*1.2, '\itns');
end
shg
% Example data


set(gca,'fontsize',cfg.fntsize);
       %%
%        legend off
  



    %%

end % id

cd(cfg.figdir);
cd('Participant phase results');
print('-dpng', ['DVs combined phase distribution']);




end % function