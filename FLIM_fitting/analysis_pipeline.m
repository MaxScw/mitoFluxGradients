clear all;

lowoxy=0.5;
highoxy=1.0;

for i=3:20
%get_intensity_decay_dist_in_rings(i,10,1,1,0.7);
fit_data_in_rings(i,1,5,lowoxy,highoxy);
end

%%
for i=3:20
average_int_and_dist(i,lowoxy,highoxy);
end

%%
for i=3:20
cellular_kn(i);
end
%%
for i=3:20
jox_cellular_kn(i);
end
%%
for i=3:20
jox_ring_kn(i);
end
%%
for i=3:20
average_lifetime_and_dist(i,lowoxy,highoxy)
end

%%
for i=3:20
mito_area_and_dist(i,lowoxy,highoxy)
end
