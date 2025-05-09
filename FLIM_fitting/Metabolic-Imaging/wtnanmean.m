function y = wtnanmean(x,w,dim)
% 2018-08-31 Tim Sanchez function, adapted from Matlab Central functions 
% 'nanmean' and 'wmean' to allow for weighted means, ignoring NaNs.
%
% FORMAT: Y = NANMEAN(X,DIM)
% 
%    Weighted average or mean value ignoring NaNs
%
%    This function enhances the functionality of NANMEAN as distributed in
%    the MATLAB Statistics Toolbox and is meant as a replacement (hence the
%    identical name).  
%
%    NANMEAN(X,DIM) calculates the mean along any dimension of the N-D
%    array X ignoring NaNs.  If DIM is omitted NANMEAN averages along the
%    first non-singleton dimension of X.
%
%    Similar replacements exist for NANSTD, NANMEDIAN, NANMIN, NANMAX, and
%    NANSUM which are all part of the NaN-suite.
%
%    See also MEAN

% -------------------------------------------------------------------------
%    author:      Jan Gläscher
%    affiliation: Neuroimage Nord, University of Hamburg, Germany
%    email:       glaescher@uke.uni-hamburg.de
%    
%    $Revision: 1.1 $ $Date: 2004/07/15 22:42:13 $


if nargin<2
    error('Not enough input arguments.');
end

if length(x)==1 % If only one element present, exit function
    y = x;
    return;
end

if ~exist('dim') dim = find(size(x)==max(size(x))); end

% Check that dimensions of X match those of W.
if(~isequal(size(x), size(w)))
    error('Inputs x and w must be the same size.');
end

% Check that all of W are non-negative.
if (any(w(:)<0))
    error('All weights, W, must be non-negative.');
end

% Check that there is at least one non-zero weight.
if (all(w(:)==0))
    error('At least one weight must be non-zero.');
end

if isempty(x)
	y = NaN;
	return
end

if nargin < 2
	dim = min(find(size(x)~=1));
	if isempty(dim)
		dim = 1;
	end
end


y = nansum(w.*x,dim)./nansum(w,dim);

% $Id: nanmean.m,v 1.1 2004/07/15 22:42:13 glaescher Exp glaescher $
