function [invIm]=InvertIm(im, varargin)
% InvertIm generalized image inversion (grayscale reversal)
%    invIm=InvertIm(im) inverts an image array im in the sense of reversing
%    its grayscale (NOT matrix inversion).  InvertIm effects the inversion
%    by computing new image data rather than by reversing the colormap or
%    look-up-table (LUT), which is often, but not always, a better
%    approach.  im can be an array of any dimension composed of numeric or
%    binary elements (integers [uint8, int8, uint16, int16, etc.], floating
%    point [double or single], or logical).  The result invIm is the same
%    data class as im.
%    
%    For numeric arrays, the inverse is computed using 
%
%       invIm = (maxVal + minVal) - im
%
%    with appropriate treatment of the different data types.  For integer
%    arrays, minVal and maxVal are the minimum and maximum allowed
%    integers.  For floating point arrays, if the image data lie within the
%    range [0,1], the inversion is performed with default minVal and maxVal
%    of 0 and 1.  If the image data lie outside [0,1], then the minimum and
%    maxmimum values of the data are used for minVal and maxVal.
%
%    Invert(im, bitsStored) returns the inverse of the image im using a
%    maxVal of 2^bitsStored.  This option is convenient when inverting
%    images with a limited dynamic range, e.g., when you have a 16-bit
%    image that only uses the first 12 bits of data.  This option is only
%    applicable to integer data.
%
%    invIm=Invert(im, [], lowerLimit) returns the inverse of the image im
%    using a minVal of lowerLimit.  This option is useful, for example,
%    when inverting signed integer image data that is only positive.
%
%    invIm=Invert(im, [], [], upperLimit) returns the inverse of the image
%    im using a maxVal of upperLimit.
% 
%    All three optional arguments may be used together (e.g.,
%    invIm=Invert(im, bitsStored, lowerLimit, upperLimit)), but upperLimit
%    will supercede bitsStored.
%
%
%    Example:
%      im=imread('cameraman.tif');  % Brigg's field at MIT
%      fig=figure; colormap(gray)
%      subplot(2,2,1)
%      subimage(im)
%      axis off
%      subplot(2,2,2)
%      imhist(im)
%      ylabel('Number of Pixels')
%
%      invIm=InvertIm(im);
%      subplot(2,2,3)
%      subimage(invIm)
%      axis off
%      subplot(2,2,4)
%      imhist(invIm)
%      ylabel('Number of Pixels')
%
%
%    See also colormap

%    Written by Stead Kiger, Ph.D.
%               Department of Radiation Oncology
%               Beth Israel Deaconess Medical Center
%               Harvard Medical School
%               Boston, MA  02215
%
%    Last revised on April 11, 2007

% Initialize parameters
bitsStored=[];
lowerLimit=[];
upperLimit=[];

% Assign and check optional arguments for parameters
error(nargchk(1, 4, nargin, 'struct'))
% bitsStored
if nargin >= 2,
    bitsStored=varargin{1};
    if ~(isnumeric(bitsStored) || isempty(bitsStored)),
        error('Optional argument bitsStored must be either numeric or empty.')
    end
end

% lowerLimit
if nargin >= 3,
    lowerLimit=varargin{2};
    if ~(isnumeric(lowerLimit) || isempty(lowerLimit)),
        error('Optional argument lowerLimit must be either numeric or empty.')
    end
end

% upperLimit
if nargin == 4,
    upperLimit=varargin{3};
    if ~(isnumeric(upperLimit) || isempty(upperLimit)),
        error('Optional argument upperLimit must be either numeric or empty.')
    end
end

% Check that the image array is either numeric or logical
if ~(islogical(im) || isnumeric(im)),
    error('Input array must be logical or numeric, but is %s',class(im))
end


if islogical(im),
    % logical - just apply not
    invIm=~im;

else

    % determine the image class, i.e., data type; e.g., uint8, int16 etc.
    imClass=class(im);

    % create the function handle for casting the double precision image back
    % into its original class
    clsHdl=str2func(imClass);

    if isinteger(im),
        % get the minimum integer
        minVal=double(intmin(imClass));
        if ~isempty(lowerLimit),
            minVal=double(max(minVal, lowerLimit));
        end

        % get the maximum integer
        maxVal=double(intmax(imClass));
        if ~isempty(bitsStored),
            maxVal=double(min(maxVal, 2^bitsStored - 1));
        end
        if ~isempty(upperLimit),
            maxVal=min(maxVal, upperLimit);
        end

    elseif isfloat(im)
        imMin=min(im(:));
        imMax=max(im(:));
        if imMin >= 0 && imMax <= 1,
            % Intensity image or RGB
            minVal=0;
            maxVal=1;
        else
            minVal=imMin;
            maxVal=imMax;
        end

        if ~isempty(lowerLimit),
            minVal=lowerLimit;
        end
        if ~isempty(upperLimit),
            maxVal=upperLimit;
        end

    else
        error('Input array must be logical or numeric, but is %s',class(im))
    end

    % Compute the inverse and cast the array back into its original class
    invIm = feval(clsHdl, (maxVal + minVal) - double(im));

end