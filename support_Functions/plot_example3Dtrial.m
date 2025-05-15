%plot_example3Dtrial 

% example walk?

clf
setmydirs_detectv3;
cd(procdatadir);
cd('GFX')
load('GFX_Data_inGaits.mat','GFX_headY');
headY=[];
for ippant = 1:35
    headY(ippant,:) = GFX_headY(ippant).doubgc;
end
cd(procdatadir);
subjID= 'BL';

% itrial = 23, 60,84
itrial=84;
cd(procdatadir)
lfile  = dir([pwd filesep subjID '*']);
load(lfile.name)
% fntsize= cfg.fntsize;
fntsize=16;
HeadPos(itrial).Z= HeadPos(itrial).Z- mean(HeadPos(itrial).Z); 
timevec= trialInfo(itrial).times;
% timevec= HeadPos(itrial).X +5;
%
figure(2); clf;
set(gcf,'units','normalized','position', [.02 .33 .9 .43]);
figure(1); clf
set(gcf,'units','normalized','position', [.02 .33 .9 .43]);


clf
figure(1)
% subplot(4,4,[1:3,5:7,9:11])
p3=plot3(timevec, HeadPos(itrial).Z,HeadPos(itrial).Y, ['b:'], 'LineWidth',2);
%trick to add shading.
x= timevec;
y= HeadPos(itrial).Z;
z= HeadPos(itrial).Y;
r = sqrt(x.^2 + y.^2 + z.^2);
g = patch('Vertices', [x(:), y(:),z(:); nan nan nan], 'Faces', (1:length(x)+1).', 'FaceVertexCData', [r(:); nan], 'EdgeColor', 'interp', 'Marker','o','linew',1)
%
% p.Alpha=.2
% axis tight;
backwall = .15;
floorwall = 1.6;
xlabel('Trial time (sec)');
ylabel('Head sway (m)', 'rotation', -35);
set(gca,'ytick',[])
ylim([-backwall backwall])
zlabel('Head height (m)');
zlim([floorwall floorwall+.15])
%
set(gca,'fontsize', fntsize);
hold on;
plot3(timevec, repmat(backwall, [1 length(trialInfo(itrial).times)]), HeadPos(itrial).Y, 'color', [.7 .7 .7], 'linew', 2)
% ylim([1.69 1.78])
hold on;
pS=plot3(timevec, HeadPos(itrial).Z, repmat(floorwall, [1 length(trialInfo(itrial).times)]), 'color', [.7 .7 .7], 'linew', 2)
% ylim([1.69 1.78])


grid on
view([-15 38]); % 
% view([-5.2 54])
% view([18 35]); % 
set(gca,'fontsize',fntsize, 'xtick', [0:1:10], 'xticklabels', ...
    {'0',...'',...
    '1',...'',...
    '2',...'',...
    '3',...''...
    '4',...''...
    '5',...''...
    '6',...''...
    '7',...''...
    '8',...''...
    '9',...''...
    '10'})
%

% add markers.
 
hold on;
% plot3([timevec(83) timevec(83)], [.05, HeadPos(itrial).Z(83)], [HeadPos(itrial).Y(83),HeadPos(itrial).Y(83)] )

% add target results?

% tOnsets = HeadPos(124).tr   
relT= find(trial_summaryTable.trial == itrial);
tOnsets = trial_summaryTable.targOnset(relT);
tResult = trial_summaryTable.targCor(relT);
tResult(4)=0;
hold on;
trialTimes = trialInfo(itrial).times;
 HeadData= HeadPos(itrial).Y;
txtHeight = 1.7;
 for itarg = 1:length(tOnsets)
   
     pAt = dsearchn(trialTimes, tOnsets(itarg)');

if tResult(itarg)==1% Hit.
    fillcol=  [.7 .7 .7];
     fillcol3=  'b';
     fillcolx='k';
else
    fillcol='w';
    fillcol3='w';
    fillcolx='w';
end

%add to back wall:
    tl=plot3(timevec(pAt), backwall, HeadData(pAt), 'o', 'color', [.7 .7 .7], 'MarkerFaceColor', fillcol,'linew',2, 'markersize', 10);    
% %and floor: 
plot3(timevec(pAt), HeadPos(itrial).Z(pAt), floorwall, 'o','color', [.7 .7 .7],'MarkerFaceColor',fillcol, 'linew',2, 'markersize', 10);

%
%and 3D trace?
tl=    plot3(timevec(pAt), HeadPos(itrial).Z(pAt), HeadPos(itrial).Y(pAt), 'o', 'color', 'b', 'MarkerFaceColor',fillcol3,'linew',2,'markersize',10);

    if tResult(itarg)==1 % store legend for hit and miss
        hl = tl;
    else
        ml=tl;
    end

%add to marker on time axis:
plot3([timevec(pAt)], [ -.15], [floorwall+.005], 'ko', 'linew',2,'MarkerFaceColor', fillcol3,'markersize',10)
% verical orientation:
plot3([timevec(pAt),timevec(pAt)], [-.15 -.15], [floorwall,floorwall+.025], 'k-', 'linew',2)
%perp orientation:
% plot3([timevec(pAt),timevec(pAt)], [-.15 -.125], [floorwall,floorwall], 'k-', 'linew',2)


% plot3([timevec(pAt)], [ -.15], [floorwall+.025], 'k-', 'linew',2)

% 
    if itarg==1
    
    %link to backwall:
            plot3([timevec(pAt) timevec(pAt)], [backwall, HeadPos(itrial).Z(pAt)], [HeadPos(itrial).Y(pAt),HeadPos(itrial).Y(pAt)],['-.'], 'color', ['k'])    
    %link to floor:
            plot3([timevec(pAt) timevec(pAt)], [ HeadPos(itrial).Z(pAt), HeadPos(itrial).Z(pAt)], [floorwall,HeadPos(itrial).Y(pAt)],['-.'], 'color', 'k')

    %link floor to time ax:
%     plot3([timevec(pAt) timevec(pAt)], [ - backwall HeadPos(itrial).Z(pAt)], [floorwall,floorwall],['-.'], 'color', 'k')
    
    %
    end
 end
 
 legend([p3, pS,hl,ml], {'3D Head position', '2D Head position', 'Target Hit', 'Target Miss'}, 'autoupdate','off', 'Location','NorthEast', 'fontsize', 16)
%  legend([p3,hl,ml], {'3D Head pos.','Target Hit', 'Target Miss'}, 'autoupdate','off', 'Location','NorthEast', 'fontsize', 16, 'NumColumns',1)


%  tt=text(0, -.15, 1.65, 'Target onsets:', 'fontsize', 20)
  tt=text(0, -.15, 1.625, 'Target onsets:', 'fontsize', 20)

% tt=text(3, backwall, 1.76, 'Walking direction', 'fontsize', 20,'BackgroundColor','w')
tt=text(1, backwall, 1.76, 'Walking direction', 'fontsize', 20,'BackgroundColor','w')

plot3([1 3], [backwall backwall], [1.72 1.72 ], 'k-', 'linew',4)
plot3([3 3], [backwall backwall], [1.72 1.72], 'k>', 'linew',4)
%  text(0.1, 1.71, {['Target'];['onset:']}, 'HorizontalAlignment', 'left', 'fontsize', fntsize-2);
 xlim([0 9])
 box on
 % enclosures?
%  plot3([timevec(1) timevec(1)], [backwall -backwall], [floorwall+.15 floorwall+.15], '-', 'color', [.5 .5 .5])
%   plot3([timevec(1) timevec(1)], [-backwall -backwall], [floorwall floorwall+.15], '-', 'color', [.5 .5 .5])
% 
%   plot3([timevec(1) timevec(end)], [-backwall -backwall], [floorwall+.15 floorwall+.15], '-', 'color', [.5 .5 .5])
%
  % shadow at floor?
   plot3([timevec(1) timevec(end)], [backwall backwall], [floorwall floorwall], '-', 'color',[.2 .2 .2 .1], 'linew',6)
   plot3([timevec(end) timevec(end)], [backwall backwall], [floorwall floorwall+.15], '-', 'color',[.2 .2 .2 .1], 'linew',6)
   plot3([timevec(end) timevec(end)], [backwall -backwall], [floorwall floorwall], '-', 'color',[.2 .2 .2 .1], 'linew',6)

%    set(gca,'PlotBoxAspectRatio', [1,.3 .4])

   
   %%
   % now plot the head height change
figure(1);
subplot(5,4,[4,12])
plot([1:200]./2, mean(headY,1), ['ko-'],'MarkerSize',10,'MarkerIndices',1:2:200);
ylim([-.02 .0275])
hold on
ht= plot(1, mean(headY(:,1),1), ['ko'],'MarkerSize',10,'MarkerIndices',1:2:200)
xlabel('');% stride-cycle');
legend(ht, 'Target onsets','fontsize', 18)
ylabel(['Head height (m)'],'fontsize',18);
set(gca,'ytick', [],'XTickLabel',[]);
box off

subplot(5,4,16);
plot([1:200]./2, repmat(1,[1,200]), 'w')
xlabel('stride-cycle (%)');
set(gca,'ytick',[], 'YColor', 'w')

set(gca,'fontsize',18);
% tightfig
 %%
% % add troughs?
% trs = HeadPos(itrial).Y_gait_troughs;
% 
% linktroughs= [1, 2,3];
% linktroughs= [];
% for itrough = 1:length(trs);
%     
%     plotsamp= trs(itrough);
%     %first add to 3D trace:
%     plot3(timevec(plotsamp), HeadPos(itrial).Z(plotsamp), HeadPos(itrial).Y(plotsamp), 'bo','linew',3);
% %    % plot , note the time will be fixed
% %    t3 = [timevec(trs(itrough)), timevec(trs(itrough))];
% %    % head pos
% %    Y3
% 
% % add to back wall projection:
% plot3(timevec(plotsamp), backwall, HeadPos(itrial).Y(plotsamp), 'ko','linew',3);
% 
% %and floor:
% 
% plot3(timevec(plotsamp), HeadPos(itrial).Z(plotsamp), floorwall, 'ko','linew',3);
% 
%     if ismember(itrough, linktroughs)
%         
%         %link to backwall:
%         plot3([timevec(plotsamp) timevec(plotsamp)], [backwall, HeadPos(itrial).Z(plotsamp)], [HeadPos(itrial).Y(plotsamp),HeadPos(itrial).Y(plotsamp)],['-.'], 'color', ['k'])
%         %link to floor:
%         plot3([timevec(plotsamp) timevec(plotsamp)], [ HeadPos(itrial).Z(plotsamp), HeadPos(itrial).Z(plotsamp)], [floorwall,HeadPos(itrial).Y(plotsamp)],['-.'], 'color', 'k')
%         
% 
%     end
%     
% end

%% make the plots for unity
