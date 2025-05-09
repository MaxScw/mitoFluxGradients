function CopyName2MultiDIndices(acqpath)

cd(acqpath)

D = subdir('name_indexes.mat')

for i = 1:length(D)
    load(D(i).name);
    if isstr(nameinds{1,2})
        nameinds = nameinds(:,[1 3 2 4 5 6 7 8 9]);
    end
    save([D(i).folder '\multiD_indices.mat'],'nameinds')
    delete(D(i).name);
end     