%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge two patches inside a common FOV using weighted averaging within 
% the overlap.
%
% July 24, 2020 by Shwetadwip Chowdhury
%
% inputs:
%           a_3D:           3D matrix of NaN, except for patch within
%                           volume containing object information
%           b_3D:           3D matrix of NaN, except for patch within
%                           volume containing object information
%
% outputs:
%           tot:            3D matrix of NaN, except for merged patch of
%                           a_3D and b_3D, where overlap is calculated via
%                           weighted average.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tot] = mergePatch(a_3D,b_3D)

% extracting just middle slice to identify overlap region between two
% patches
a = a_3D(:,:,round(end/2));
b = b_3D(:,:,round(end/2));

% finding overlap region
a_msk = ~isnan(a);
b_msk = ~isnan(b);
ovrlp = single((a_msk+b_msk)==2);

% setting coordinate axes
x = 1:size(a_msk,2);
y = 1:size(a_msk,1);
[xx,yy] = meshgrid(x,y);

% finding boundary of overlap regions, for both 'a' and 'b'. These will be
% used as the edges to pull overlap weighting from.
bndry_a = bndry(a_msk).*bndry(ovrlp);
bndry_b = bndry(b_msk).*bndry(ovrlp);
bndry_total = bndry(a_msk+b_msk>0);
bndry_a(bndry_total==1)=0;
bndry_b(bndry_total==1)=0;

if (max(bndry_a(:)) == 0 && max(bndry_b(:)) > 0)
    tot = a;
elseif (max(bndry_b(:)) == 0 && max(bndry_a(:)) > 0)
    tot = b;
else
    % identifying x and y coordinates of different boundary components
    ovrlp_x = xx(ovrlp==1);
    ovrlp_y = yy(ovrlp==1);
    bndry_a_x = xx(bndry_a==1);
    bndry_a_y = yy(bndry_a==1);
    bndry_b_x = xx(bndry_b==1);
    bndry_b_y = yy(bndry_b==1);
    
    ovrlp_xx    = repmat(ovrlp_x,[1,length(bndry_a_x)]);
    ovrlp_yy    = repmat(ovrlp_y,[1,length(bndry_a_x)]);
    bndry_a_xx  = repmat(bndry_a_x',[length(ovrlp_x),1]);
    bndry_a_yy  = repmat(bndry_a_y',[length(ovrlp_x),1]);
    
    % normalized distance from overlap boundary to nearest edge of mask 1
    dist_a_lin = min(sqrt((ovrlp_xx-bndry_a_xx).^2+(ovrlp_yy-bndry_a_yy).^2),[],2);
    dist_a = zeros(size(a))*NaN;
    dist_a(ovrlp==1) = dist_a_lin./max(dist_a_lin);
    
    ovrlp_xx    = repmat(ovrlp_x,[1,length(bndry_b_x)]);
    ovrlp_yy    = repmat(ovrlp_y,[1,length(bndry_b_x)]);
    bndry_b_xx  = repmat(bndry_b_x',[length(ovrlp_x),1]);
    bndry_b_yy  = repmat(bndry_b_y',[length(ovrlp_x),1]);
    
    % normalized distance from overlap boundary to nearest edge of mask 2
    dist_b_lin = min(sqrt((ovrlp_xx-bndry_b_xx).^2+(ovrlp_yy-bndry_b_yy).^2),[],2);
    dist_b = zeros(size(a))*NaN;
    dist_b(ovrlp==1) = dist_b_lin./max(dist_b_lin);
    
    % normalizing distance based on min distance between mask 1 & mask 2
    dist_a_norm = dist_a;
    dist_b_norm = dist_b;
    dist_a_norm(ovrlp==1) = dist_a(ovrlp==1)./(dist_a(ovrlp==1)+dist_b(ovrlp==1));
    dist_b_norm(ovrlp==1) = dist_b(ovrlp==1)./(dist_a(ovrlp==1)+dist_b(ovrlp==1));
    
    % extending this weighted overlap mask into 3D
    a_tot = a_3D;
    b_tot = b_3D;
    ovrlp3D = repmat(ovrlp==1,[1,1,size(a_3D,3)]);
    a_msk3D = repmat(a_msk,   [1,1,size(a_3D,3)]);
    b_msk3D = repmat(b_msk,   [1,1,size(a_3D,3)]);
    dist_a_norm_3D   = repmat(dist_a_norm,[1,1,size(a_3D,3)]);
    dist_b_norm_3D   = repmat(dist_b_norm,[1,1,size(a_3D,3)]);
    
    a_tot(ovrlp3D) = a_3D(ovrlp3D).*dist_a_norm_3D(ovrlp3D);
    b_tot(ovrlp3D) = b_3D(ovrlp3D).*dist_b_norm_3D(ovrlp3D);
    
    tot = zeros(size(a_3D))*NaN;
    tot(a_msk3D+b_msk3D>0) = 0;
    tot(a_msk3D) = a_tot(a_msk3D)+tot(a_msk3D);
    tot(b_msk3D) = b_tot(b_msk3D)+tot(b_msk3D);
end
end

% Finding boundary of a binary image
function bn = bndry(bin_im)
bin_im_xp = circshift(bin_im,[1,0]);
bin_im_xn = circshift(bin_im,[-1,0]);
bin_im_yp = circshift(bin_im,[0,1]);
bin_im_yn = circshift(bin_im,[0,-1]);
bn    =               bin_im~=bin_im_xp ...
    | bin_im~=bin_im_yp ...
    | bin_im~=bin_im_xn ...
    | bin_im~=bin_im_yn;
bn = bn.*bin_im;
end
