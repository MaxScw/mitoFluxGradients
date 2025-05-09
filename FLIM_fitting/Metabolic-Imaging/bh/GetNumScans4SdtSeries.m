
close all;
clear all;
path = 'Z:\Lab\Marta\2016-11-17 ROS assay tests\2017-05-19 Batch 2\s1_a1\sorted_sdts\Pos0';

D = dir([path '\*.sdt']);

for i = 1:length(D)
    file = [path '\' D(i).name];
    
    sdt = bh_readsetup(file);
    meas = bh_getmeasdesc(sdt,1);
    numscans(i) = double(meas.hist_fida_points);

% ch = bh_getdatablock_v095(sdt,1);
% img = uint8(squeeze(sum(ch,1)));
% flim = ch;
end

plot(numscans,'o')
xlabel('frame'); ylabel('num scans')
