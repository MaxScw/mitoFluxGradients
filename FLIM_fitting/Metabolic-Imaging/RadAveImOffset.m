function Pf = RadAveImOffset(im,xoff,yoff)
% im = im2;
% xoff = 0; yoff = 0;

% Adapted from the Matlab Central function raPsd2d(img,res) ((C) E.
% Ruzanski, RCG, 2009)
% That program which took the fft2 of the image and averaged radially.  
% Because I already have the ffts, or I can fft2 outside of the program, 
% I just feed it whatever I want to radially average.  The program does 
% a fftshift of the zero point to the center, so this is designed for 
% data that has a peak at zero.
%

%% Process image size information
[N M] = size(im);

%% Adjust PSD size
dimDiff = abs(N-M);
dimMax = max(N,M);
% Make square
if N > M                                                                    % More rows than columns
    if ~mod(dimDiff,2)                                                      % Even difference
        im = [NaN(N,dimDiff/2) im NaN(N,dimDiff/2)];                  % Pad columns to match dimensions
    else                                                                    % Odd difference
        im = [NaN(N,floor(dimDiff/2)) im NaN(N,floor(dimDiff/2)+1)];
    end
elseif N < M                                                                % More columns than rows
    if ~mod(dimDiff,2)                                                      % Even difference
        im = [NaN(dimDiff/2,M); im; NaN(dimDiff/2,M)];                % Pad rows to match dimensions
    else
        im = [NaN(floor(dimDiff/2),M); im; NaN(floor(dimDiff/2)+1,M)];% Pad rows to match dimensions
    end
end

%% Radially average power spectrum
[X Y] = meshgrid(-dimMax/2:dimMax/2-1, -dimMax/2:dimMax/2-1);               % Make Cartesian grid
[theta rho] = cart2pol(X, Y);                                               % Convert to polar coordinate axes
rho = round(rho);
% Offsets - circshift center of rho, then set over-hanging elements on other
% side to -1 so they are excluded from any averaging below.
% NOTE: as always, y coord is reversed for images. So a positive y-off
% means the center of the radial average is BELOW the image center
rho = circshift(rho,yoff,1);
rho = circshift(rho,xoff,2);
if yoff>0
    rho(1:yoff,:) = -1;
else
    rho(end-yoff:end,:) = -1;
end
if xoff>0
    rho(:,1:xoff) = -1;
else
    rho(:,end-xoff:end) = -1;
end


i = cell(floor(dimMax/2) + 1, 1);
for r = 0:floor(dimMax/2)
    i{r + 1} = find(rho == r);
    Li(r+1) = length(i{r + 1});
end
Pf = zeros(1, floor(dimMax/2)+1);
for r = 0:floor(dimMax/2)
    Pf(1, r + 1) = nanmean( im( i{r+1} ) );
%     im3 = im; im3(i{r+1})=1000; 
%     imshow(im3)
end
