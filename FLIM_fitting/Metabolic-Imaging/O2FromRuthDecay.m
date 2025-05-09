function [O2T, O2] = O2FromRuthDecay(decaypath,DropRng,RestRng,O2StartC)
% O2 experiments performed with alternating acquisitions with sample
% position and a ruthenium dye droplet. Dye sensitive to O2 concentration
% in the chamber. Expteriment flushes O2 out of chamber, waits
% for O2 to exchange in solution, then O2 is flushed back in.
% Inputs:
% -decaypath
% -DropRng: O2 drop start finish times
% -RestRng: O2 restore start finish times
% -O2StartC: Starting O2 concentration
% NOTE: all times are relative to the save time of the first frame

% clear all
% decaypath = 'C:\Users\Tim\Desktop\test\Ext_O2droprestore\fits_Pos1_GenMask.mat';
% DropRng = [0 1.8];
% RestRng = [2.75 5.5];
% O2StartC = 20;

timeunit = 3600; % 60 gives minutes, 3600 for hours, 86400 for days
switch timeunit
    case 1;
        tlab = 'Time (s)';
    case 60;
        tlab = 'Time (min)';
    case 3600;
        tlab = 'Time (h)';
    case 86400;
        tlab = 'Time (days)';
end
% Versions:



%% Group data in one cell. Also construct a cell that can be written in spreadsheet form
AllDecays = [];
labels = {};
tm = [];
chs = [];
load(decaypath);
% Get t and ch from deay structures instead of relying on nameinds at
% this point
decays_fits_struct(cellfun('isempty',decays_fits_struct)) = [];
for j = 1:size(decays_fits_struct,1)
    slashes = strfind(decaypath,'\');
    [T,ch,Z] = MultiD_Parse_FName(decays_fits_struct{j}.filename);
    t = decays_fits_struct{j}.timestp;
    
    tm = [tm; t];
    chs = [chs; {ch}];
    label = [decaypath(slashes(end)+1:end-4) '_' ch '_t' num2str(T) '_z' num2str(Z)];
    labels = [labels; label];
end
currL = length(decays_fits_struct);
AllDecays = [AllDecays; decays_fits_struct];

tm = (tm - min(tm))*86400/timeunit;

L = length(AllDecays);
for i = 1:L
    fitres = AllDecays{i}.fit_result;
    Chisqred(i) = AllDecays{i}.Chi_sq;
    AllParams{i} = fitres;
    irrs(i) = AllDecays{i}.irr;
    numphot(i) = sum(AllDecays{i}.decay);
end
ParamsArr = reshape(cell2mat(AllParams),5,3,length(AllDecays));

Nind = strcmp(chs,'NADH');
Find = strcmp(chs,'FAD');

%% Get FAD channel ruthenium info
ts = tm(Find);
RuthIrrs = irrs(Find);

tdiffs = abs(ts-DropRng(1)); ind1 = find(tdiffs==min(tdiffs));
tdiffs = abs(ts-DropRng(2)); ind2 = find(tdiffs==min(tdiffs));
DropTs = ts(ind1:ind2);
DropIrrs = RuthIrrs(ind1:ind2);
DropIrrs = -DropIrrs; DropIrrs = DropIrrs - min(DropIrrs); 

tdiffs = abs(ts-RestRng(1)); ind3 = find(tdiffs==min(tdiffs));
DropTs2 = ts(ind1:ind3-1); % Longer time frame going up to start of Restore
tdiffs = abs(ts-RestRng(2)); ind4 = find(tdiffs==min(tdiffs));
RestTs = ts(ind3:ind4);
RestTs = RestTs;
RestIrrs = RuthIrrs(ind3:ind4);
RestIrrs = RestIrrs - min(RestIrrs);  

% Fits - shift curves to start at t=0
[fitresult, gof] = ExpFitO2(DropTs-DropTs(1), DropIrrs'); 
coeffs = coeffvalues(fitresult);
DropFitT = (ts(ind1):0.005:ts(ind3)) - ts(ind1);
DropFit = O2StartC*exp(-DropFitT.*coeffs(2));
% Reshift time back to experimental frame
DropFitT = ts(ind1):0.005:ts(ind3);
DropDecTau = coeffs(2);

[fitresult, gof] = ExpFitO2(RestTs-RestTs(1), RestIrrs');
coeffs = coeffvalues(fitresult);
RestFitT = (ts(ind3):0.005:ts(ind4)) - ts(ind3);
RestFit = O2StartC-(O2StartC-DropFit(end))*exp(-RestFitT.*coeffs(2));
RestFitT = ts(ind3):0.005:ts(ind4);
RestDecTau = coeffs(2);

%% Stitch full O2 vs exp time (relative to 1st frame save time), and save
O2T = [DropFitT'; RestFitT'];
O2 = [DropFit'; RestFit'];

save([UpOneDir(decaypath) '\O2Info'],'O2T','O2','DropRng','RestRng','O2StartC','DropDecTau','RestDecTau')
