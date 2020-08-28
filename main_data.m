%%
% This file demonstrates the multi-slice beam propagation (MSBP) method on
% experimental C. elegans data. C. elegans are well known to be optically
% transparent, but multiple scattering organisms. This file will run the
% MSBP reconstruction protocol to reconstruct 3D refractive-index (RI)
% distributions of the C. elegans. This file includes: 1) downloading the
% experimental data and optical system parameters from saved .TIF and .MAT
% files; 2) decomposing the .TIF stack into patches, since MSBP is a
% computational intensive process and often cannot be applied on the full
% .TIF field-of-view due to practical hardware constraints; 3) A user-input
% field ('patchFOV_index') to choose which patch to reconstruct with MSBP;
% 4) iteratively reconstructing 3D refractive index within the patch C.
% elegans volume via gradient-based optimization (NOTE: this takes a long
% time!); and 5) saving the patch and important parameter variables into a
% .MAT file.
%
% Modification date: July 25, 2020
%
% Reference for MSBP:
% S. Chowdhury, M. Chen, R. Eckert, D. Ren, F. Wu, N. Repina, and L. 
% Waller, "High-resolution 3D refractive index microscopy of multiple-
% scattering samples from intensity images," Optica 6, 1211-1219 (2019) 
%
% Acknowledgements:
% This MSBP demonstration utilizes regularization using the 3D total 
% variation (TV) proximal operator. We gratefully acknowledge the work done 
% by the UNLocBoX team, which has provided an extensive open-source 
% compilation of several optimizers utilizing proximity operators. From
% the UNLocBoX toolbox, we use the 'prox_tv3d' function for 3D TV prox
% regularization. Full documentation of UNLocBoX can be found at:
%
% UNLocBoX toolbox:     https://epfl-lts2.github.io/unlocbox-html/
% download this toolbox and put into the same path as this .M file 
%
% Reference for 'prox_tv3d':  
% A. Beck and M. Teboulle. Fast gradient-based algorithms for constrained 
% total variation image denoising and deblurring problems. Image 
% Processing, IEEE Transactions on, 18(11):2419--2434, 2009


clear; close all;
clc;
%% adding relevant paths

folderPath1 = 'core_functions\';
folderPath2 = 'FOV_01_reconstruction_CElegan_resources\';
folderPath3 = 'unlocbox\';
folderPath4 = 'FOV_01_reconstruction_CElegan_outputs_V1\' ;

addpath(folderPath1);
addpath(folderPath2);
addpath(folderPath3);
init_unlocbox();

%% Downloading data from .TIF file and imaging parameters from .MAT file

totFOV_amp_acqs     = single(tif2matrix([folderPath2,'FOV_01_rawdata_CElegan.tif'])).^0.5;
load([folderPath2,'FOV_01_reconstruction_CElegan_params.mat']);
disp('finished loading data and data_parameters');

%% Setting parameters relevant to physical object volume (some are drawn directly from parameter .MAT file)

use_gpu     = true;                 % TRUE: use GPU device; FALSE: use computer CPU
ps          = dx;                   % pixel size (x,y,z) in object space (micron)
lambda      = lambda;               % central wavelength (micron)
NA          = NA;                   % numerical aperture of imaging and detection lens
n_imm       = 1.51;                 % refractive index of immersion media
n_m         = 1.33;                 % refractive index of media (for reconstruction purposes)
pdar        = 30;                   % padding size to avoid edge artifacts
z_plane     = 0;                    % center plane of reconstruction volume, where 0 um is object volume center

%% Setting up the patches of field-of-views (FOV) to divide up complete field-of-view into manageable chunks

% Setting parameters of patch FOVs
patchFOV_size   = 360;                          % lateral pixel dimension of a patch within the total FOV
grid_num        = 4;                            % number of patches per side of the total FOV
totalFOV_size   = size(totFOV_amp_acqs,1);      % size of the total FOV, which will basically be divided into a grid of (GRID_NUM)x(GRID_NUM) patches

% Calculating number of pixels in overlap region between two adjacent
% patches. PATCHFOV_SIZE and GRID_NUM must be selected so that FOV_OVERLAP
% is an integer!
FOV_overlap     = (grid_num*patchFOV_size-totalFOV_size)/(grid_num-1);
if FOV_overlap < 0
    error('patch size too small to span total FOV with user-inputed grid number');
elseif FOV_overlap ~= round(FOV_overlap)
    error(sprintf('patch size does not allow non-integer dimension length of overlap patch. Tune parameters ''patchFOV_size'' and ''grid_num'' to achieve integer value of ''FOV_overlap'', which is currently %2.2f pixels',FOV_overlap));
end
fprintf('%2.0f pixels overlap between adjacent patches in the FOV\n',FOV_overlap);

% calculating indices for each patch FOV in the object's total FOV
coordSt         = zeros(grid_num,patchFOV_size);
for k = 1:grid_num
  coordSt(k,:)  = (k-1)*(patchFOV_size-FOV_overlap)+(1:patchFOV_size);
end


%% Specifying and visualizing which patch in the total FOV to reconstruct.

patchFOV_index = 9;                    % user-input on which patch to reconstruct

% calculating row and column indices for user-selected patch
if patchFOV_index > grid_num^2
    error('patch index exceeds number of patches ');
end
rowInd = ceil(patchFOV_index/grid_num); 
colInd = mod(patchFOV_index,grid_num); if(colInd ==0); colInd = grid_num; end
rows   = coordSt(rowInd,:);
cols   = coordSt(colInd,:);

% visualizing the patch FOV within the total FOV, so that user can visually
% inspect patch's location and confirm that the patch is the one to select
figure(1)
imagesc(totFOV_amp_acqs(:,:,1)); axis equal; colormap gray; axis tight; caxis(1+[-1.0,1.4]); title(sprintf('patchFOV %2.0f',patchFOV_index));
hold on;
rectangle('Position',[cols(1),rows(1),patchFOV_size,patchFOV_size],...
         'LineWidth',3, 'EdgeColor','r')
hold off;

% cropping the total FOV to only consider the user-selected patch for
% subsequent MSBP reconstruction
amp_acqs    = totFOV_amp_acqs(rows,cols,:);

%% Setting spatial and frequency axes and propagation kernels

N           = size(amp_acqs,1)+...
                2*pdar;             % lateral pixel dimension of a padded patch within the object
x           = ps*[-N/2:N/2-1];      % 1D padded axis in x
[xx,yy]     = meshgrid(x,x);        % 2D padded grid in x/y

dfx         = 1/(N*ps);             % Fourier spacing of padded axis
fx          = dfx*[-N/2:N/2-1];     % 1D padded axis in fx
[fxx,fyy]   = meshgrid(fx,fx);      % 2D padded grid in fx/fy

fx          = ifftshift(fx);        % FFT shifting Fourier axes
fxx         = ifftshift(fxx);       % FFT shifting Fourier axes
fyy         = ifftshift(fyy);       % FFT shifting Fourier axes

% setting propagation kernels and pupil support
prop_phs            = 1i*2*pi*sqrt((n_imm/lambda)^2-(fxx.^2+fyy.^2));
NA_crop             = (fxx.^2 + fyy.^2 > (NA/lambda)^2);

% converting into GPU arrays if user targets gpu-enabling
if use_gpu
    xx              = gpuArray(xx);
    yy              = gpuArray(yy);
    fyy             = gpuArray(fyy);
    fyy             = gpuArray(fyy);
    prop_phs        = gpuArray(prop_phs);
    amp_acqs        = gpuArray(amp_acqs);
end
%% Downloading illumination k-vectors from .MAT file and accounting for system scan-angle orientation

fx_in           = -kx;                               
fy_in           = -ky;  

%% initializing forward model measurements and initial guess of reconstructed object

O               = 060;                  % axial dimension size of reconstruction space                                        
psz             = 0.5;                  % pixel size (z) in reconstructed object space(micron)
                                        % lateral pixel size is assumed to be same as variable 'ps'

reconObj        = single(0*randn([N,... % initialization of guess of reconstructed object (deltaRI, not RI), to be updated iteratively     
                              N,...
                              O,]));   	

if use_gpu
    reconObj     = gpuArray(reconObj);
end
%% optimization params for iterative reconstruction

maxiter         = 200;                  % number of iterations to run optimization protocol for
step_size       = 5e-4;                 % step size for gradient-based optimization protocol
plot_range      = [n_m-0.04,n_m+0.08];  % contrast to be used to show the reconstruction at each iteration
cost            = zeros(maxiter,1);     % cost function to evaluate convergence
      
reconObj_prox   = reconObj;             % used for Nesterov acceleration protocol for faster convergence
t_k             = 1;                    % parameter used for Nesterov acceleration

regParam        = 1e-3; %1e-3;                 % regularization parameter for 3D proxTV

%% initializing Figure windows to observe iterative process
close all;
% triframe cross-sectional views of the reconstructed object, as it undergoes
% iterative updates
figure('Name','Reconstruction result');
figNum = 1;
MSBP_progview(real(reconObj)+n_m,figNum,plot_range,cost, 0)
    
pause(0.01);

%% Running iterative optimization of object volume. Variable 'reconObj' is the final 3D refractive-index reconstruction!
gpu_fail_num = 0;

tic;
for iter = 1:maxiter
    
    pause(0.01);
    
    % randomly scramble angles and choose without replacement
    seq = randperm(length(fx_in));
    for illum_angle = 1:length(fx_in)
        
        % compute estimated exit field on the camera plane
        [efield,efield_vol]     = MultiSlice_Forward(reconObj, psz, xx, yy, dfx, prop_phs, NA_crop, lambda, fx_in(seq(illum_angle)), fy_in(seq(illum_angle)), z_plane, pdar, use_gpu);
        
        % compute gradient (and update refractive index at each layer)
        [reconObj,funcVal]      = BPM_update(reconObj, psz, efield, efield_vol, amp_acqs(:,:,seq(illum_angle)), prop_phs, NA_crop, lambda, z_plane, step_size, pdar);
        
        % compute accumulated error for current iteration
        cost(iter)         = cost(iter) + gather(funcVal);
        fprintf('illum_angle: %1.0d  iteration: %1.0d\n',illum_angle,iter)
    end
    
    % Prox operator is a memory-intensiver operator. If GPU crashes due to
    % memory requirements, use CPU instead. It will be slower but the
    % program won't crash.
    try
        reconObj_prox1 = prox_tv3d(real(reconObj), regParam);
    catch
        disp('running regularizer on CPU because GPU ran out of memory for this memory-intensive procedure');
        reconObj_prox1  = prox_tv3d(gather(real(reconObj)), regParam);
        reconObj_prox1  = gpuArray(reconObj_prox1);
        gpu_fail_num    = gpu_fail_num+1;   % counter to keep track of how many times GPU failed to regularize
    end
    
    if iter>1
        if cost(end) > cost(end-1)
            t_k   = 1;
            reconObj = reconObj_prox;
            continue;
        end
    end
    
    % Nesterov's update
    t_k1       = 0.5 * (1 + sqrt(1 + 4 * t_k^2));
    beta       = (t_k - 1)/t_k1;
    reconObj   = reconObj_prox1 + beta*(reconObj_prox1 - reconObj_prox);
    t_k        = t_k1;
    reconObj_prox = reconObj_prox1;
    fprintf('iteration: %d, error: %5.5e, elapsed time: %5.2f seconds\n',iter, cost(iter));
    
    MSBP_progview(real(reconObj)+n_m, figNum, plot_range, cost, iter)
    pause(0.01);

end
toc;
close_unlocbox();
%% Running Forward model on reconstructed object and comparing to raw measurements. Should be useful for troubleshooting purposes

efield_acqs_fwd = zeros([N-2*pdar,N-2*pdar,length(fx_in)]);
if use_gpu
    efield_acqs_fwd = gpuArray(efield_acqs_fwd);
end
for idx = 1:length(fx_in)
    [efield_fwd,~]              = MultiSlice_Forward(reconObj, psz, xx, yy, dfx, prop_phs, NA_crop, lambda, fx_in(idx), fy_in(idx), 0, pdar, use_gpu);     % Multi-slice forward model
    efield_acqs_fwd(:,:,idx)    = efield_fwd;                                                                 % storing efield output of forward model
    disp(['simulate data: ', num2str(idx)]);
end
amp_acqs_fwd = abs(efield_acqs_fwd); 
comp_acqs = cat(2,amp_acqs(:,:,1:length(fx_in)),amp_acqs_fwd);
sliderDisplayImVC2(comp_acqs, {'colormap gray', 'title(''(LEFT: measurement data)    (RIGHT: forward model confirmation)'')'});

%% In case you want to save reconstruted data and relevant parameters
% save([folderPath4 'FOV_01_reconstruction_CElegan_output_patch_' sprintf('%0.2d',patchFOV_index) '.mat'], ...
%     'reconObj', 'amp_acqs','cost', 'regParam','t_k','NA', 'lambda', 'ps', 'psz', ...
%     'step_size', 'n_imm', 'n_m', 'pdar','z_plane','N','O','rowInd','colInd','patchFOV_index', ...
%     'rows', 'cols', 'coordSt','totalFOV_size');
