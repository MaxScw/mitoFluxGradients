
function errb = ErrorBars95Conf(Arr)

% inputs an array of data values, calculates the mean, standard deviation, and the
% confidence interval (error bars) for a 95% confidence.

% datamean = mean(Arr(Arr>0)); % I don't remember why I was only taking positive vals
% datastddev= std(Arr(Arr>0)); 
datamean = nanmean(Arr);
datastddev = nanstd(Arr);

%calculate the lower confidence limit and the upper confidence limit.  To
%do this, use the formula LCL=mu-Z*sigma/sqrt(N). where mu is the mean, sigma is the standard deviation, N
%is the number of data points, and Z is a value based on integrating a
%normal (gaussian) function.  Z is a standard value that changes depending
%on what confidence percentage you want.  In our case, 95% gives Z=1.95996
errb = 1.95996*datastddev/sqrt(length(Arr));
