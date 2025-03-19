clear all;

lowoxy=0;
highoxy=0.5;

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
jox_ring_kn_with_eq_lifetime(i);
end

