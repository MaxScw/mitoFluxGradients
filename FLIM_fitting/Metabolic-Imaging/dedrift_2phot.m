function rdf=dedrift(lub, dr)
%re-written by Tim Sanchez 4-11-10
%to remove drift.  Much faster now.

for i = 1:length(lub(:,1))
     
    if (lub(i,3)>1)&(lub(i,3)<=length(dr)+1)
        lub(i,1) = lub(i,1)-dr(lub(i,3)-1,1);
        lub(i,2) = lub(i,2)-dr(lub(i,3)-1,2);
    end
    
end
rdf=lub;