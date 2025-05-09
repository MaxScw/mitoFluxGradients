clear all
clc

temperature = 22;
rings=10;
lowoxy=0.2;
highoxy=0.3;

splitLow = strsplit(string(lowoxy), '.');
splitHigh = strsplit(string(highoxy), '.');
savepath = join(['temp_',string(temperature),'C/run_oxy_l',...
            splitLow(1), 'p', splitLow(2), '_h',...
            splitHigh(1), 'p', splitHigh(2)], '');
%%
for i=3:5
    i
get_intensity_decay_dist_in_rings(i,rings,1,1,0.7, savepath);
fit_data_in_rings(i,1,5,lowoxy,highoxy, savepath);
end

%%
for i=3:5
    i
average_int_and_dist(i,lowoxy,highoxy, savepath);
end

%%
for i=3:5
cellular_kn(i, savepath);
end

%%
for i=3:5
jox_cellular_kn(i, temperature, savepath);
end
%%
for i=3:5
jox_ring_kn(i, temperature, savepath);
end
