function [SVMModel,ph] = SVM_3dPlot_plane_TS(meas,species,cost,BoxConstraint)
% Function for plotting the decision boundary plane resulting from an SVM
% fitting using 'fitcsvm'. Only works for 2 groups, but could adapt to
% more.
% INPUTS:
%  -meas: N-by-K matrix with N data points and K measurement parameters
%  -species: N-row array of strings that indicate which of the two groups
%            each data point belongs to.
% OUTPUTS:
%  -SVMModel: outputs the SVM model trained by the input meas
%  -Also creates a plot of the plane and data right over gcf.

% clear all; close all;

if ~exist('cost') cost = [0,1;1,0]; end
if ~exist('BoxConstraint') BoxConstraint = 100; end
labs = unique(species);
colors = 'rgb';
markers = 'osd';
% plot(meas(strcmp(species,'virginica'),1),meas(strcmp(species,'virginica'),2),'b.'); hold on;plot(meas(~strcmp(species,'virginica'),1),meas(~strcmp(species,'virginica'),2),'r.');hold off

% fitcsvm

Xrng = min(meas(:,1)):(max(meas(:,1))-min(meas(:,1)))/10:max(meas(:,1));
Yrng = min(meas(:,2)):(max(meas(:,2))-min(meas(:,2)))/10:max(meas(:,2));
Zrng = min(meas(:,3)):(max(meas(:,3))-min(meas(:,3)))/10:max(meas(:,3));
[x1Grid,x2Grid,x3Grid] = meshgrid(Xrng,Yrng,Zrng);
xGrid = [x1Grid(:),x2Grid(:),x3Grid(:)];

% Fit model from data
SVMModel = fitcsvm(meas,species,'Standardize',true,'BoxConstraint',BoxConstraint,'Cost',cost);

% Note, this hinges on Matlab's 'predict' function, outputing 'scores'
%    "Predicted class posterior probabilities, returned as a numeric matrix 
%    of size N-by-K. N is the number of observations (rows) in meas, and K is 
%    the number of classes (in Mdl.ClassNames). score(i,j) is the posterior probability that observation i in meas is of class j in Mdl.ClassNames."
% Other inputs are: [label,score,cost] = predict(...)
% In case of plane, the model will predict that all points 'above' the
% plain will have a positive probability (of being in group 1), and all the
% points 'below' the plane have negative probabilities. The scores thus
% define the regions, and so define the dividing plane.
[~,scores] = predict(SVMModel,xGrid);



% Plot the data and the decision boundary
% figure;
% for i = 1:length(labs)
%     data = meas(strcmp(species,labs{i}),:);
%     plot3(data(:,1), data(:,2), data(:,3), [colors(i) markers(i)]);
%     hold on;
% end
% 
% % Support vectors
% h(3) = plot3(meas(SVMModel.IsSupportVector,1),meas(SVMModel.IsSupportVector,2),meas(SVMModel.IsSupportVector,3),'k.');

% Optimal plane plot
% Use 'isosurface' to find the bounding surface of one of the two 'positive
% probability' regions - the region above the plane, for example. Just pick 
% either group, the outcome will be the same.
v0=reshape(scores(:,2),size(x1Grid));
fv = isosurface(x1Grid, x2Grid, x3Grid, v0, 0);
ph = patch(fv,'edgecolor', 'none', 'FaceAlpha', 0.4);


% legend([unique(species); 'Support Vectors']','location','best');
hold off





%% Working code for the fishery example


% % clear all; close all;
% % load data
% % load fisheriris
% % inds = ~strcmp(species,'versicolor');
% % inds = ~strcmp(species,'setosa');
% X = meas(inds,2:4);
% y = species(inds);
% labs = unique(y);
% colors = 'rgb';
% markers = 'osd';
% % plot(X(strcmp(y,'virginica'),1),X(strcmp(y,'virginica'),2),'b.'); hold on;plot(X(~strcmp(y,'virginica'),1),X(~strcmp(y,'virginica'),2),'r.');hold off
% 
% % fitcsvm
% 
% d = 0.02;
% [x1Grid,x2Grid,x3Grid] = meshgrid(min(X(:,1)):d:max(X(:,1)),...
%     min(X(:,2)):d:max(X(:,2)),min(X(:,3)):d:max(X(:,3)));
% xGrid = [x1Grid(:),x2Grid(:),x3Grid(:)];
% 
% % Fit model from data
% SVMModel = fitcsvm(X,y);
% 
% % Note, this hinges on Matlab's 'predict' function, outputing 'scores'
% %    "Predicted class posterior probabilities, returned as a numeric matrix 
% %    of size N-by-K. N is the number of observations (rows) in X, and K is 
% %    the number of classes (in Mdl.ClassNames). score(i,j) is the posterior probability that observation i in X is of class j in Mdl.ClassNames."
% % Other inputs are: [label,score,cost] = predict(...)
% % In case of plane, the model will predict that all points 'above' the
% % plain will have a positive probability (of being in group 1), and all the
% % points 'below' the plane have negative probabilities. The scores thus
% % define the regions, and so define the dividing plane.
% [~,scores] = predict(SVMModel,xGrid);
% 
% 
% 
% % Plot the data and the decision boundary
% % figure;
% % for i = 1:length(labs)
% %     data = X(strcmp(y,labs{i}),:);
% %     plot3(data(:,1), data(:,2), data(:,3), [colors(i) markers(i)]);
% %     hold on;
% % end
% 
% % Support vectors
% h(3) = plot3(X(SVMModel.IsSupportVector,1),X(SVMModel.IsSupportVector,2),X(SVMModel.IsSupportVector,3),'ko');
% 
% % Optimal plane plot
% % Use 'isosurface' to find the bounding surface of one of the two 'positive
% % probability' regions - the region above the plane, for example. Just pick 
% % either group, the outcome will be the same.
% v0=reshape(scores(:,2),size(x1Grid));
% fv = isosurface(x1Grid, x2Grid, x3Grid, v0, 0);
% p = patch(fv,'edgecolor', 'none', 'FaceAlpha', 0.5);
% 
% 
% legend([unique(y); 'Support Vectors']','location','best');
% hold off





%% Working code for 2D plot. Expand to 3D plane above.
% clear all; close all;
% % load data
% load fisheriris
% inds = ~strcmp(species,'versicolor');
% % inds = ~strcmp(species,'setosa');
% X = meas(inds,3:4);
% y = species(inds);
% % plot(X(strcmp(y,'virginica'),1),X(strcmp(y,'virginica'),2),'b.'); hold on;plot(X(~strcmp(y,'virginica'),1),X(~strcmp(y,'virginica'),2),'r.');hold off
% 
% % fitcsvm
% 
% d = 0.02;
% [x1Grid,x2Grid] = meshgrid(min(X(:,1)):d:max(X(:,1)),...
%     min(X(:,2)):d:max(X(:,2)));
% xGrid = [x1Grid(:),x2Grid(:)];
% 
% SVMModel = fitcsvm(X,y);
% [~,scores] = predict(SVMModel,xGrid);
% 
% 
% % Plot the data and the decision boundary
% figure;
% h(1:2) = gscatter(X(:,1),X(:,2),y,'rb','.');
% hold on
% h(3) = plot(X(SVMModel.IsSupportVector,1),X(SVMModel.IsSupportVector,2),'ko');
% contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');
% legend(h,[unique(y); 'Support Vectors']);
% hold off




