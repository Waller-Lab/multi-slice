%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearly scales object from .TIF file so that refractive-index values
% fall between n_media and n_max
%
% July 06, 2020 by Shwetadwip Chowdhury
% 
% inputs:   
%           fname:      path variable to .TIF stack containing sample phantom
%           n_min:      refractive index of media (corresponding to '0'
%                       value in .TIF file)
%           n_max:      refractive index of feature within 3D phantom with
%                       maximum optical density

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

function obj = makePhantom(fname, n_min, n_max, rows, cols, idx)
if nargin < 4
    obj     = tif2matrix(fname);
elseif nargin < 6
    obj     = tif2matrix(fname,rows,cols);
else
    obj     = tif2matrix(fname,rows,cols,idx);
end

obj_max     = max(obj(:));
obj         = (n_max-n_min)/obj_max*(obj-obj_max)+n_max;