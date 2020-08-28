%%
% This file demonstrates the multi-slice beam propagation (MSBP) method on
% a simulated phantom object. This file includes: 1) initializing a phantom
% object as a 3D distribution of refractive index; 2) initializing a spiral
% scan pattern of illumination angles; 3)Forward modeling the optical
% scattering process as light propagates through the phantom object for
% each illumination angle, according to the multi-slice framework. The
% output will be a set of simulated AMPLITUDE measurements; and 4)
% iteratively reconstructing the 3D refractive index within the volume via
% gradient-based optimization. 
%
% Note: no regularization is used
%
% Modification date: July 25, 2020
%
% reference:
% S. Chowdhury, M. Chen, R. Eckert, D. Ren, F. Wu, N. Repina, and L. 
% Waller, "High-resolution 3D refractive index microscopy of multiple-
% scattering samples from intensity images," Optica 6, 1211-1219 (2019) 


addpath('core_functions/');
%%
clear; 
close all;

%% Downloading phantom from .TIF file

n_media     = 1.33;         % refractive index of media
n_max       = 1.5;          % refractive index of max density feature 
pdar        = 120;            % padding size before running forward model (to avoid edge artifacts)
obj         = makePhantom('phantom.tif',n_media,n_max);     % making phantom

%% Setting parameters relevant to physical object volume

use_gpu     = true;                % TRUE: use GPU device; FALSE: use computer CPU
ps          = 0.1;                  % pixel size (x,y,z) in object space (micron)
lambda      = 0.532;                % central wavelength (micron)
NA          = 1.0;                  % numerical aperture of imaging and detection lens
n_imm       = 1.33;                 % refractive index of immersion media
n_m         = n_media;              % refractive index of media (for reconstruction purposes)
z_plane     = 0;                    % center plane of reconstruction volume, where 0 um is object volume center


%% Setting spatial and frequency axes and propagation kernels

N           = size(obj,1)+2*pdar;   % lateral pixel dimension of padded object
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
prop_crop           = (fxx.^2 + fyy.^2 > (n_imm/lambda)^2);

% converting into GPU arrays if user targets gpu-enabling
if use_gpu
    obj             = gpuArray(obj);
    xx              = gpuArray(xx);
    yy              = gpuArray(yy);
    fyy             = gpuArray(fyy);
    fyy             = gpuArray(fyy);
    prop_phs        = gpuArray(prop_phs);
end

%% Setting parameters relevant to illumination angles

N_k             = 200;          % number of illumination-angle acquisitions
revs            = 8;            % number of revolutions the spiral takes
outerCirc       = true;         % boolean on whether to include outer circle points
N_o             = 30;           % if outerCirc == true, number of points in outer circle

% freq initialization of illumination angles
[kx_in,ky_in]   = generateSpiralPath(N_k, revs, outerCirc, N_o);    % generating normalized spiral coordinates
kx_in           = NA/lambda*kx_in;                                  % scaling fx by NA/lambda to span pupil function
ky_in           = NA/lambda*ky_in;                                  % scaling fy by NA/lambda to span pupil function

fx_in = ky_in;
fy_in = kx_in;

plot(fx_in,fy_in,':o'); axis equal; axis tight;
%% Running forward model on object phantom to simulate instrument measurements

obj_pad         = padarray(obj,[pdar,pdar,0],n_media);            % padding 3D phantom array
illumAngles     = 1:length(fx_in);
efield_acqs     = zeros([size(obj,1),size(obj,2),length(illumAngles)]);

if use_gpu
    efield_acqs     = gpuArray(efield_acqs);
end

tic;
for idx = illumAngles
    [efield,~]              = MultiSlice_Forward(obj_pad, ps, xx, yy, dfx, ...
                                                 prop_phs, NA_crop, lambda, ...
                                                 fx_in(idx), fy_in(idx), 0, ...
                                                 pdar,use_gpu);     % Multi-slice forward model
    efield_acqs(:,:,idx)    = efield;                               % storing efield output of forward model
    disp(['simulate data: ', num2str(idx)]);
end
disp('Simulated acquisitions via MSBP forward model is complete, and ready for reconstruction');
toc;

amp_acqs        = abs(efield_acqs);     % taking only the amplitude of the 
                                        % stored efield simulated
                                        % measurements, since our
                                        % reconstruction uses
                                        % INTENSITY-ONLY measurements.
                                        
%% initializing forward model measurements and initial guess of reconstructed object

O               = 100;                   % axial dimension size of reconstruction space                                        
psz             = 0.1;                  % pixel size (z) in reconstructed object space(micron)
                                        % lateral pixel size is assumed to be same as variable 'ps'

reconObj        = 0*randn([N,... % initialization of guess of reconstructed object (deltaRI, not RI), to be updated iteratively     
                              N,...
                              O,]);   

if use_gpu
    reconObj     = gpuArray(reconObj);
end
%% optimization params for iterative reconstruction

maxiter         = 100;                  % number of iterations to run optimization protocol for
step_size       = 1e-3;                 % step size for gradient-based optimization protocol
plot_range      = [n_m-0.06,n_m+0.17];  % contrast to be used to show the reconstruction at each iteration
cost            = zeros(maxiter,1);     % cost function to evaluate convergence
      
reconObj_prox   = reconObj;             % used for Nesterov acceleration protocol for faster convergence
t_k             = 1;                    % parameter used for Nesterov acceleration

%% initializing Figure windows to observe iterative process
close all;

% triframe cross-sectional views of the true phantom (to be used as a visual
% benchmark to evaluate convergence accuracy)
figure('Name','True Phantom (padded)');
MSBP_progview(real(obj_pad),1,plot_range); 

% triframe cross-sectional views of the reconstructed object, as it goes
% iterative updates
figure('Name','Reconstruction result');
MSBP_progview(real(reconObj)+n_m,2,plot_range,cost, 0)
    
pause(0.01);

%% Running iterative optimization of object volume. Variable 'reconObj' is the final 3D refractive-index reconstruction!
tic;
for iter = (1:maxiter)
    
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
    
    reconObj_prox1 = reconObj;
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
    reconObj      = reconObj_prox1 + beta*(reconObj_prox1 - reconObj_prox);
    t_k        = t_k1;
    reconObj_prox = reconObj_prox1;
    fprintf('iteration: %d, error: %5.5e, elapsed time: %5.2f seconds\n',iter, cost(iter), toc());
    
    MSBP_progview(real(reconObj)+n_m, 2, plot_range, cost, iter)
    pause(0.01);

end
toc;
%% Running Forward model on reconstructed object and comparing to raw measurements. Should be useful for troubleshooting purposes

efield_acqs_fwd = zeros([size(obj,1),size(obj,2),length(illumAngles)]);
for idx = illumAngles
    [efield_fwd,~]              = MultiSlice_Forward(reconObj, psz, xx, yy, dfx, n_m, prop_phs, NA_crop, lambda, fx_in(idx), fy_in(idx), -z_plane, pdar);     % Multi-slice forward model
    efield_acqs_fwd(:,:,idx)    = efield_fwd;                                                                 % storing efield output of forward model
    disp(['simulate data: ', num2str(idx)]);
end
amp_acqs_fwd = abs(efield_acqs_fwd);
comp_acqs = cat(2,amp_acqs,amp_acqs_fwd);
sliderDisplayImVC2(comp_acqs, {'colormap gray', 'title(''(LEFT: measurement data)    (RIGHT: forward model confirmation)'')'});
