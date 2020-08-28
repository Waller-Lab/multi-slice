%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts 3D .TIF grayscale stack into a 3D matlab array
%
% July 06, 2020 by Shwetadwip Chowdhury
% 
% inputs:   
%           fname: path variable to .TIF stack containing sample phantom
%
%           (Optional inputs)
%           rows:       1D array of row indices to accept within overall 
%                       .TIF 3D phantom
%           cols:       1D array of col indices to accept within overall 
%                       .TIF 3D phantom
%           idx:        1D array of z-indices to accept within overall .TIF
%                       3D phantom
%
% outputs:
%           obj:        3D object from .TIF stack, scaled to be between 
%                       n_media and n_max
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj = tif2matrix(fname,rows,cols,idx)
cnt = 1;
info = imfinfo(fname);
num_images = numel(info);
im = imread(fname, 1);

if nargin < 4
    idx = 1:num_images;
end
if nargin==1
    rows = 1:size(im,1);
    cols = 1:size(im,2);
end


for k = idx
    if cnt == 1
        obj = zeros([length(rows),length(cols),length(idx)]);
    end
    im = imread(fname, k);
    obj(:,:,cnt) = im(rows,cols);
    cnt = cnt+1;
end
