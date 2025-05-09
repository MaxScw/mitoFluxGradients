function StkWrite(img, filename, nameappend)
% wr_img(img, filename): write a struct array
% to a 3D image file in uncompressed tif format
% nameappend should be a string to append to the filename of the original
% file.  If you want to create a shorter stack, for example, you could
% enter the string, '_short', and a new file would be created with the name
% 'filename_short'

filename = filename(1:(strfind(filename, '.tif')-1));

if nargin ~= 3
    nameappend='';
else
    nameappend = ['_' nameappend];
end

wr_mode = 'overwrite';

% Check img type.
% Option 1: It's a struct with images in 'imstruct(SliceNumber).data'
if isstruct(img)
    for i=1:size(img,2)
        if ~isempty(img(i).data)
            imwrite(img(i).data,[filename nameappend '.tif'], 'tif', 'Compression',...
                'none','WriteMode', wr_mode);
            wr_mode = 'append';
        end
    end
elseif length(size(img))==3
    for i=1:size(img,3)
        imwrite(img(:,:,i),[filename nameappend '.tif'], 'tif', 'Compression',...
            'none','WriteMode', wr_mode);
        wr_mode = 'append';
    end
end
% Option 2: It's a 3D array with [x,y,SliceNumber]


