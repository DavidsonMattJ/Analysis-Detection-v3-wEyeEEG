% jXX_imageEyeandTargLocation%


% This script will attempt to capture error by target location, to see if
% there are hemifield effects


%%%%%% Quest detect v3 %%%%%%

setmydirs_detectv3;
cd(procdatadir)
% load:
pfols= dir([pwd  filesep '*summary_data.mat']);
nsubs= length(pfols);
%show ppant list:
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%% implementing a few different jobs.

% consider splitting the target location into a smaller grid.

job.calcTrial_targetEyelocation= 1;
% 
job.plotPFX=1;
% 
% job.concatGFX =0;
% 
% job.plotGFX =1;
%%
% calculate density per grid cell.

   [gridlimsX,gridlimsY]=deal( -.3:.01:.3);
%%
if job.calcTrial_targetEyelocation

    %keep a running tally per trial, to average over trial types.
    % densGrid = nan(length(gridlimsX), length(gridlimsY));
    for isub = 1%:nsubs
        cd(procdatadir )
        %%load data from import job.
        load(pfols(isub).name, 'EyeDir', 'EyePos' ,'TargPos','HeadPos','subjID','trial_summaryTable', 'gait_ts_gData');

        % we have error per sample (HandPos.
        savename = [subjID '_PFX_data.mat'];
        disp(['preparing whole trial target location error for ' subjID]);
        %%
        [PFX_densGrid_Target,PFX_densGrid_EyeLoc]=deal([]); % target location 
          %%
        [PFX_densGrid_Targtmp, PFX_densGrid_Eyetmp]=deal([]); % changes size on loops, so can't preallocate.
        trialTypesGrid=[];
        
        keept=1;
        [targlocY_hist,targlocZ_hist]=deal([]); % build distribution of all locs (to debug range)
        for  itrial=1:size(HeadPos,2)

            relvrows = find(trial_summaryTable.trial==itrial);
      
             isPrac = trial_summaryTable.isPrac( relvrows(1));
             isStationary= trial_summaryTable.isStationary( relvrows(1));
            %% skip bad trials:
            if isPrac || isStationary
                continue
            end
            % subj specific trial rejection
            skip=0;
            rejTrials_detectv3; %toggles 'skip' based on bad trial ID
            if skip==1
                continue
            end

            %per sample we want Eye and the target location  (Y-Z is relevant):

            %% debug
            %normalize target space, accounting for scaling:
            tZ= TargPos(itrial).Z-.03; % centre at zero (weirdly)
            tY= TargPos(itrial).Y- TargPos(itrial).Y(1);
            eZ= EyeDir(itrial).Z;
            eY= EyeDir(itrial).Y;

            %detect and remove outliers (discontinuities in eye data).
            oZ = isoutlier(eZ);
            oY= isoutlier(eY);
            eZ(oZ)=nan;
            eY(oZ)=nan;
%             %
%             clf; subplot(211);
%             plot(tZ); hold on;
%             plot(eZ);shg
%             subplot(212)
%             plot(tY); hold on;
%             plot(eY);
targlocZ_hist    = [targlocZ_hist, tZ(:)'];       
targlocY_hist= [targlocY_hist, tY(:)'];
            %%
            %for each sample, find closest X and Y location in our grid
            % %used for calculating density.
            [assignXTarg] = dsearchn(gridlimsX', tZ);
            [assignYTarg] = dsearchn(gridlimsY', tY);


             [assignXEye] = dsearchn(gridlimsX', eZ);
             [assignYEye] = dsearchn(gridlimsY', eY);

            %trial specific location density map:
            [densGridTarg,  densGridEye] = deal(zeros(length(gridlimsX), length(gridlimsY)));
            
            relvrows_g = find(gait_ts_gData.trialallocation==itrial);
             %restruct to friendly within trial window (first to last gait with target events).
            sampstart = gait_ts_gData.gaitStart(relvrows_g(1));
            sampend =  gait_ts_gData.gaitStart(relvrows_g(end)) + gait_ts_gData.gaitSamps(relvrows_g(end));

            
            for isamp = sampstart:sampend
                % for each, store the location in our density maps:

                %first store the location in our density map:
                rowCol = [assignXTarg(isamp), assignYTarg(isamp)];
                %increment by 1
                densGridTarg(rowCol(1), rowCol(2)) = densGridTarg(rowCol(1), rowCol(2)) + 1;
               
                rowCol = [assignXEye(isamp), assignYEye(isamp)];
                 %increment by 1
                densGridEye(rowCol(1), rowCol(2)) = densGridEye(rowCol(1), rowCol(2)) + 1;
               

            end


            %%
            PFX_densGrid_Targtmp(keept,:,:)= densGridTarg;
            
            PFX_densGrid_Eyetmp(keept,:,:)= densGridEye;

           
            %%

            trialTypesGrid(keept)= ~isStationary;
           
            keept=keept+1;
        end % all trials.
        
        % store PFX 
        for itype=1:2
          
            subset = find(trialTypesGrid==itype-1);           

% total density per type:
            PFX_densGrid_Target( itype, :,:) =  squeeze(sum(PFX_densGrid_Targtmp(subset,:,:),1,'omitnan'));
            
            PFX_densGrid_EyeLoc( itype, :,:) = squeeze(sum(PFX_densGrid_Eyetmp(subset,:,:),1,'omitnan'));
%            
            
        end

%         save(savename, 'PFX_densGrid_Target', 'PFX_densGrid_EyeLoc', '-append');
    end
end

%%
if job.plotPFX
    % errMap = cbrewer('seq', 'Reds', 5);
    colormap('Inferno');
    % adjust to show background colour.
    %% adjust background colour.
    ctmp = colormap;
    ctmp(1,:) = [.9 .9 .9];
    colormap(ctmp);

    cd([procdatadir ])
    pfiles= dir([pwd filesep '*_PFX_data.mat']);
    figure(1); clf; set(gcf, 'units','normalized', 'position', [.1 .1 .8 .8]);
    figure(2); clf; set(gcf, 'units','normalized', 'position', [.1 .1 .8 .8]);
    for isub = 1%:nsubs
        cd([procdatadir ])
%         load(pfiles(isub).name, 'PFX_densGrid', 'PFX_avErrorGrid', 'trialTypesGrid', 'subjID');

        %% plot Grand error:
        %PFX, average, and remove centre point (for plots)
        PFX_dens = squeeze(sum(PFX_densGrid_Target,1, 'omitnan'));
        [rws,cls]= size(PFX_dens);
       
        figure(1); clf;
        subplot(121)
        imagesc(gridlimsX, gridlimsY,PFX_dens);
        axis square;
        c=colorbar; title('All trials')
        ylabel(c, 'Target location ')
        set(gca,'fontsize',15)
%         xlim([-.1 .1]); ylim([-.1 .1])


        PFX_densEye= squeeze(mean(PFX_densGrid_EyeLoc,1, 'omitnan'));
        
        subplot(122);
        imagesc(gridlimsX, gridlimsY, PFX_densEye)
        axis square
        c=colorbar;title('All trials')
        ylabel(c, 'Eye direction');
        set(gca,'fontsize',15)
        
        %% print >>>
%         cd([figdir filesep 'Target_locationError']);
%         figure(1); gcf;
%         print('-dpng', [subjID(1:2) '_TargetLocation_alltrials'])

        
        %% repeat for within ppant x type:
        figure(2); clf;
        ic=1;
        plotorder= [3,1,4,2];
        for itype=1:4

            plottype = plotorder(itype);
            %plot
            subplot(2,4,itype)
            imagesc(gridlimsX, gridlimsY,squeeze(PFX_densGrid(plottype,:,:)));
            axis square;
            title(useLeg{plottype})
            caxis([0 75])

            xlabel('X [m]');
            c=colorbar;

            if itype==1
                ylabel('Y [m]');
            elseif itype==4
                ylabel(c, 'Frame count', 'fontsize',12)
            end
            set(gca,'fontsize',12)
            ic=ic+1;

            %plot

            subplot(2,4,itype + 4);
            imagesc(gridlimsX, gridlimsY, squeeze(PFX_avErrorGrid(plottype,:,:)));
            axis square

            xlabel('X [m]');
            caxis([0 .15])
            c=colorbar;

            if itype==1
                ylabel('Y [m]');
            elseif itype==4
                ylabel(c, 'Error', 'fontsize',12)
            end
            set(gca,'fontsize',12)
            ic=ic+1;
        end
        colormap(ctmp)
        %% print 
        figure(2); gcf;
        print('-dpng', [subjID(1:2) '_TargetLocation_bytrialTypes'])
    end

end

%%
if job.concatGFX
    cd([procdatadir filesep 'PFX'])
    pfiles= dir([pwd filesep '*_PFX_data.mat']);

    GFX_densGrid= zeros(nsubs,4,length(gridlimsY), length(gridlimsY));
    GFX_avErrGrid= zeros(nsubs,4,length(gridlimsY), length(gridlimsY));

    for isub = 1:nsubs
        cd([procdatadir filesep 'PFX'])
        load(pfiles(isub).name, 'PFX_densGrid', 'PFX_avErrorGrid');

        for itype=1:4

            subset = find(trialTypesGrid==itype);

            GFX_densGrid(isub, itype, :,:) =  squeeze(PFX_densGrid(itype,:,:));
            GFX_avErrGrid(isub, itype, :,:)=  squeeze(PFX_avErrorGrid(itype,:,:));

        end
    end

 % before saving, restrict the grid to locations with data for all
    % subjects (avoids large effects at target boundary).
    
    for itype=1:4
        %sum of data.
        for ir=1:length(gridlimsX)
            for ic=1:length(gridlimsY)

                allD = squeeze(GFX_avErrGrid(:,itype,ir,ic));

                if sum(~isnan(allD))~=nsubs
                    GFX_densGrid(:,itype,ir,ic)= nan;
                    GFX_avErrGrid(:,itype,ir,ic)= nan;

                end

            end
        end


    end





    cd([procdatadir filesep 'GFX']);
    save('GFX_targError_byLocation', 'GFX_avErrGrid', 'GFX_densGrid', 'gridlimsY', 'gridlimsX', '-append');
end
%%
if job.plotGFX
  
    cd([procdatadir filesep 'GFX'])
    load('GFX_targError_byLocation');

    figure(1); clf; set(gcf, 'units','normalized', 'position', [.1 .1 .8 .8], 'color', 'w');
    figure(2); clf; set(gcf, 'units','normalized', 'position', [.1 .1 .8 .8], 'color', 'w');
  % errMap = cbrewer('seq', 'Reds', 5);
    colormap('Inferno');
    % adjust to show background colour.
    %% adjust background colour.
 
   ctmp = colormap;
    ctmp(1,:) = [.9 .9 .9];
    colormap(ctmp);

    %% plot Grand error:
    %PFX, average, and remove centre point (for plots)
    GFX_denst = squeeze(mean(GFX_densGrid,1, 'omitnan'));
    %mean over trialtypes:
    GFX_dens= squeeze(mean(GFX_denst,1, 'omitnan'));
    [rws,cls]= size(GFX_dens);
%     GFX_dens(ceil(rws/2),ceil(cls/2))=0;

    figure(1); clf;
    subplot(121)
%     Z=inpaint_nans(GFX_dens,2);
    Z= inpaint_row_col(GFX_dens);
    imagesc(gridlimsX, gridlimsY,Z);
    axis square;
    c=colorbar; title('Target position')% (\itN\rm=26)')
    ylabel(c, 'Avg. Frame count', 'fontsize', 20)
    set(gca,'fontsize',20)
    caxis([0 60])
    xlabel('X [m]');
    ylabel('Y [m]')
    xlim([-.2 .2]);
    ylim([-.2 .2]);


    GFX_errt = squeeze(mean(GFX_avErrGrid,1, 'omitnan'));
    %mean over trialtypes:
    GFX_err= squeeze(mean(GFX_errt,1, 'omitnan'));
%     GFX_err(ceil(rws/2),ceil(cls/2))=0;

    subplot(122);
    Z=inpaint_row_col(GFX_err);

    imagesc(gridlimsX, gridlimsY, Z)
    axis square
    c=colorbar;
    title('Target error')
    ylabel(c, 'Avg. Distance [m]');
    set(gca,'fontsize',20)
%     caxis([0 .12])
     xlabel('X [m]');
%     ylabel('Y [m]')
    colormap(ctmp)
        xlim([-.2 .2]);
    ylim([-.2 .2]);
    %% >>> print

    cd([figdir filesep 'Target_locationError_GFX']);
    figure(1); gcf;
    print('-dpng', 'GFX_TargetLocation_alltrials')
    
    
    %% and surf map for ms:
   figure(10);
    Z=inpaint_row_col(GFX_err);

surf(gridlimsX,gridlimsY,Z,'FaceColor', 'interp', 'EdgeColor', 'interp');
colorbar;
colormap('inferno');
caxis([.095 .134])
set(gcf, 'color', 'w')
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Avg.Distance [m]');
view([40 20]);
set(gca,'fontsize', 16)
% axis square
xlim([-.3 .3]);
ylim([-.3 .3]);
% colormap(ctmp)
%%

figure(10); gcf;
print('-dpng', 'GFX_TargetLocation_surfError');

    
    
    
    
    
    
    %% repeat for each type:
    figure(2); clf;
%     ic=1;

typeorder= [3,1,4,2];

for itype=1:4

    ttype = typeorder(itype);

    %prep
    GFX_err = squeeze(mean(GFX_avErrGrid(:,ttype,:,:),1, 'omitnan'));

   Z= inpaint_row_col(GFX_err);
    %plot
    subplot(1,4,itype);

   %3D:
    surf(gridlimsX, gridlimsY, Z, 'FaceColor', 'interp', 'EdgeColor', 'interp');
 
    %%
    
    view([35 18])
    caxis([.09 .13])
    title(useLeg{ttype}, 'fontsize',15)

   
    ylabel('Y [m]')   
    xlabel('X [m]'); 
    zlim([0.09 .13]);
end
    colormap(ctmp)
    %%
  
    figure(2); gcf;
    print('-dpng', 'GFX_TargetLocation_bytrialTypes')

%% show collapsed across dims (X dim, Ydim);
figure(2); clf;
typeorder= [3,1,4,2];
for itype=1:4

    ttype = typeorder(itype);
       %prep
    GFX_err = squeeze(GFX_avErrGrid(:,ttype,:,:));
    
    GFX_inperr=[];
    for ippant=1:size(GFX_err,1)
    Z= inpaint_row_col(squeeze(GFX_err(ippant,:,:)));
    GFX_inperr(ippant,:,:)= Z;
    end
    
    pDRow = squeeze(mean(GFX_inperr,2, 'omitnan'));
    pDCol= squeeze(mean(GFX_inperr,3, 'omitnan'));
    

    %plot
    subplot(1,2,1); hold on; %xdim
    %2D:
    %shaderrorbar:
    mP= squeeze(mean(pDRow,1,'omitnan'));
    stE= CousineauSEM(pDRow);
    sh=shadedErrorBar(gridlimsX,mP, stE,[],1);
    sh.mainLine.Color = usecolsWalk{ttype};
    xlabel('Xposition');
    ylabel('mean error');
% plot(gridlimsX,squeeze(mean(GFX_err,1,'omitnan')));

 %plot
    subplot(1,2,2); hold on; %xdim
    %2D:
     mP= squeeze(mean(pDCol,1,'omitnan'));
    stE= CousineauSEM(pDCol);
    shadedErrorBar(gridlimsX,mP, stE);

    xlabel('Y position');
    ylabel('mean error');
    

end
%     colormap(ctmp)






end
%% %called above
% the following function hacks a 2D interpolation, between the bounds of
% real data in each row (first) and then column (second).

function gmat=inpaint_row_col(gmat)
    
    for ir=1:size(gmat,1)
        rdata = gmat(ir,:);
        % if more than a single data point, interp between first and last.
        if sum(~isnan(rdata))>=2
            %interp between index:
            tmpD = rdata;
            xps= find(~isnan(tmpD));
            [s,l]=bounds(xps);
            zD=tmpD;
            zt = inpaint_nans(tmpD(s:l),1);
            zD(s:l)=zt;
            gmat(ir,:) = zD;
        end

    end
    for ic=1:size(gmat,2)
        rdata = gmat(:,ic);
        % if more than a single data point, interp between first and last.
        if sum(~isnan(rdata))>=2
            %interp between index:
            tmpD = rdata;
            xps= find(~isnan(tmpD));
            [s,l]=bounds(xps);
            zD=tmpD;
            zt = inpaint_nans(tmpD(s:l),1);
            zD(s:l)=zt;
            gmat(:,ic) = zD;
        end

    end
 end
