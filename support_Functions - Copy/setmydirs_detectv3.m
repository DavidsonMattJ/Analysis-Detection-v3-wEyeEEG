%setmydirs_detectv3
% should be relative to current data location. 
% PC
% datadir='C:\Users\User\Documents\matt\GitHub\Analysis-Tracking-v1\data_Raw';
% new laptop
datadir= [pwd filesep 'data_Raw'];

cd(datadir)
try cd ../../data_Processed
catch 
    cd ../../
    mkdir(data_Processed)
    cd(data_Processed)
end
procdatadir = pwd;

try cd ../Figures;
catch
    cd ../
    mkdir('Figures');
    cd('Figures');
end
figdir= pwd;
cd(procdatadir)

