%%
% This function implements the forward model, as defined by the multi-slice beam propagation method.
%
% inputs:   obj:            3D matrix of \delta RI values
%           psz:            pixel size (z) in reconstructed object space(micron)
%           xx:             2D grid of x-coordinates
%           yy:             2D grid of y-coordinates
%           dfx:            Fourier spacing
%           prop_phs:       phase term of propagation kernel
%           NA_crop:        pupil support
%           lambda:         imaging wavelength
%           fx_in:          horizontal component of illumination wave-vector
%           fy_in:          vertical component of illumination wave-vector
%           z_plane:        center plane of reconstruction volume, where 0 um is object volume center
%           pdar:           padding used on the object
%           use_gpu:        % TRUE: use GPU device; FALSE: use computer CPU
%
% outputs:
%           efield:         electric-field output
%           efield_vol:     layer-specific electric-fields
%
% Author: Shwetadwip Chowdhury; July 25, 2020
% Thank you to Michael Chen and David Ren, for preliminary 
% versions of this code
%
% reference:
% S. Chowdhury, M. Chen, R. Eckert, D. Ren, F. Wu, N. Repina, and L. 
% Waller, "High-resolution 3D refractive index microscopy of multiple-
% scattering samples from intensity images," Optica 6, 1211-1219 (2019) 

function [efield,efield_vol] = MultiSlice_Forward(obj, psz, xx, yy, dfx, prop_phs, NA_crop, lambda, fx_in, fy_in, z_plane, pdar, use_gpu)

    % initalizing outpts
    efield_vol      = zeros(size(xx,1), size(xx,2), size(obj,3));  
    efield          = zeros(size(obj, 1), size(obj, 2), length(z_plane));
    
    % converting into GPU arrays if user targets gpu-enabling
    if use_gpu
        efield_vol  = gpuArray(efield_vol);
        efield      = gpuArray(efield);
    end
    
    % create an incident planewave
    fx_in_interp = fix(fx_in/dfx)*dfx;
    fy_in_interp = fix(fy_in/dfx)*dfx;
    U_in         = exp(1i * 2 * pi * (fx_in_interp * xx + fy_in_interp * yy));
    fU_current   = fft2(U_in);
    
    % multi-slice propagation
    prop_kernel  = exp(prop_phs * psz);
    for layerIdx = 1:size(obj,3)
        efield_vol(:, :, layerIdx)  = ifft2(fU_current.* prop_kernel);
        fU_current                  = fft2(efield_vol(:,:,layerIdx).* exp(1i*2*pi*(obj(:,:,layerIdx))*psz/lambda));
    end
    
    % back-propagation to volume center
    prop_kernel = conj(exp( (size(obj, 3)/2 - 1) .* prop_phs .* psz) );
    fU_current  = fU_current.* prop_kernel;
    
    % propagate the field to multiple focusing planes and crop by the NA  
    for zIdx = 1:length(z_plane)
        prop_kernel             = exp(z_plane(zIdx) .* prop_phs);
        prop_kernel(NA_crop)    = 0;
        efield(:,:,zIdx)        = ifft2(fU_current.* prop_kernel);
        efield                  = efield((pdar+1):(end-pdar),(pdar+1):(end-pdar),:);
    end   
    
end
