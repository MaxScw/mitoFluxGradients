% Script for reading the spreadsheet in the same path, containing what
% analysis to do on each sample. 

pth = pwd;
[num,txt,res]= xlsread([pth '\AcqAnalysisList.xlsx']);

% for i = 2:size(res,2)
i=2;

acqpath = [res{i,1} '\' res{i,2}];

Bools = [res{i,3:end}];

if Bools(1)
%     Run stuff;
end
if Bools(2)
%     Run stuff2;
end