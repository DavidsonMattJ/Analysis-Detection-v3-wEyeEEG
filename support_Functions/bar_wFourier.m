function [bh, Fh]= bar_wFourier(xvec, ppantData)

% called to plot the nice bar charts with fourier overlayed. returns handle
% for bar (bh) and Fourier fit (Fh)

gM = squeeze(nanmean(ppantData));
stE = CousineauSEM(ppantData);
hold on;
% finely sampled bar, each gait "%" point.
bh=bar(xvec, gM);
hold on;
errorbar(xvec, gM, stE, ...
    'color', 'k',...
    'linestyle', 'none',...
    'linew', 1);
bh.FaceColor = [.9 .9 .9];
% bh.EdgeColor = barCols{id};

sdrange = max(gM) - min(gM);
ylim([min(gM)-.5*sdrange max(gM)+1*sdrange])
ylsat = get(gca, 'ylim');


midp=xvec(ceil(length(xvec)/2));

%% apply best fit (overlay)
[f,gof]= fit(xvec',gM',  'fourier1');

% plot:
hold on;
%             yyaxis right
Fh=plot(f, xvec, gM);%,
Fh(2).LineWidth = 2;
%treat max xvec as our full 'period'
% fitperiod = f.w;
%convert to period per samples.
% include period and Rsquared
%treat max xvec as our full 'period'
% Hzapp = xvec(end)/ (2*pi/(f.w));
% legdetails = [sprintf('%.2f', Hzapp) ' cycles per, R^2 = ' sprintf('%.2f', gof.rsquare) ];
% legdetails = [' R^2 = ' sprintf('%.2f', gof.rsquare) ];

% legend(Fh(2), legdetails, 'fontsize', 8, 'autoupdate', 'off', 'Location', 'NorthEast')


end