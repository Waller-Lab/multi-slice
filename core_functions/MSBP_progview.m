%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays a triframe cross-sectional view of a 3D volume, and may also
% plot a figure associated with a cost curve (optional)
%
% July 06, 2020 by Shwetadwip Chowdhury
% 
% inputs:   
%           obj:            3D matrix 
%           fignum:         number assigned to this figure
%           plot_range:     color axis with which to display values
%           pltCurve:       curve that may get plot (optional)
%           upto:           boundary up to which curve will be plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MSBP_progview(obj,fignum,plot_range,pltCurve, upto)

figure(fignum);
subplot(2,2,1);
imagesc(real(squeeze(obj(:,:,end/2)))); axis equal; axis tight;
caxis(plot_range); colormap gray; colorbar; title('x,y');
subplot(2,2,2);
imagesc(real(squeeze(obj(:,end/2,:)))); axis equal; axis tight;
caxis(plot_range); colormap gray; colorbar; title('x,z');
subplot(2,2,3);
imagesc(real(squeeze(obj(end/2,:,:)))); axis equal; axis tight;
caxis(plot_range); colormap gray; colorbar; title('y,z');

if nargin == 5
    subplot(2,2,4);
    plot(log10(1+pltCurve(1:upto)),':o'); title('cost function'); axis tight;
end