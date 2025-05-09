
function peaks = findpeaksT(data,iter,rad)
% Found that 'FastPeakFind' and 'extrema2' from Matlab Central didn't work
% well. This is a makeshift iterative way of finding peaks.
% -Find global maximum
% -Delete around that peak with a radius, 'rad'
% -Find the next global maximum, and so forth for 'iter' iterations

for i = 1:iter
    % first find the biggest peak
    [y(i) x(i)] = find(data==max(data(:)));
    scores(i) = max(data(:));
    % Then delete a radius around that about equal to the egg dimension
    [X,Y] = meshgrid((1:size(data,2))-x(i),(1:size(data,1))-y(i));
    Z = sqrt(X.^2+Y.^2);
    pMask = find(Z<rad);
    
%     % Uncomment for diagnostic plots
%     subplot(1,2,1)
%     surf(data,'edgecolor','none'); hold on;
%     view([0,0,1]); xlim([0 512]); ylim([0 512])
%     plot(x(i),y(i),'rx','markersize',20,'linewidth',3); hold off;
%     subplot(1,2,2)
%     surf(Z,'edgecolor','none'); hold on;
%     plot(x(i),y(i),'rx','markersize',20,'linewidth',3); 
%     hold off; xlim([0 512]); ylim([0 512])
%     view([0,0,1]);
    
    data(pMask) = min(data(:));
    
end
peaks = x';
peaks(:,2) = y;
peaks(:,3) = scores;