%%
% This function implements the back-projection gradient update component
% of the multi-slice beam propagation method.
%
% inputs:   obj:            3D matrix of \delta RI values
%           psz:            pixel size (z) in reconstructed object space(micron)
%           efield:         electric-field measurement simulated output from a Forward model
%           efield_vol:     layer-specific electric-fields, also output from a Forward model
%           amp_acqs:       amplitude measurements detected by imaging system
%           prop_phs:       phase term of propagation kernel
%           NA_crop:        pupil support
%           lambda:         imaging wavelength
%           z_plane:        center plane of reconstruction volume, where 0 um is object volume center
%           step_size:      step size for gradient-based optimization protocol
%           pdar:           padding used on the object
%           use_field:      TRUE: uses field-component for reconstruction 
%                           FALSE: uses amplitude-only component for reconstruction
%
% outputs:
%           obj:            object after gradient-update iteration
%           cost_value:     cost value
%
% Author: Shwetadwip Chowdhury; July 25, 2020
% Thank you to Michael Chen and David Ren, for preliminary 
% versions of this code
%
% reference:
% S. Chowdhury, M. Chen, R. Eckert, D. Ren, F. Wu, N. Repina, and L. 
% Waller, "High-resolution 3D refractive index microscopy of multiple-
% scattering samples from intensity images," Optica 6, 1211-1219 (2019) 

function [obj,cost_value] = BPM_update(obj, psz, efield, efield_vol, acq, prop_phs, NA_crop, lambda, z_plane, step_size, pdar, use_field)
    % compute residual 
    if use_field
        back_prop           = efield - acq;
        res                 = back_prop;
    else
        back_prop           = efield - abs(acq).* exp(1i*angle(efield));
        res                 = abs(efield) - abs(acq);
    end

    % compute cost
    cost_value              = norm(res(:))^2;
    
    %% compute gradient 
    back_prop(efield==0)    = 0;
    back_prop               = padarray(back_prop,[pdar,pdar,0]);
    
    % back-propagation to volume center
    for zIdx = 1:length(z_plane)
        prop_kernel             = conj(exp(1.* z_plane(zIdx) .* prop_phs));
        prop_kernel(NA_crop)    = 0;
        back_prop(:,:,zIdx)     = prop_kernel.* fft2(back_prop(:,:,zIdx));
    end
    
    % propagation to the last slice
    prop_kernel             = exp(1*(size(obj,3)/2 - 1)*prop_phs*psz);
    back_prop               = ifft2(back_prop.* prop_kernel);
    prop_kernel             = conj(exp(1*prop_phs*psz));
    
    % compute gradient at each slice and update the object
    for layerIdx = size(obj,3):-1:1
        grad                = conj(efield_vol(:,:,layerIdx)).* back_prop;
        transmission_layer  = conj(exp(1i*2*pi*obj(:,:,layerIdx)*psz/lambda));
        grad                = transmission_layer.* grad;
        grad                = (-1i*2*pi*psz/lambda) * grad;
        obj(:,:,layerIdx)   = obj(:,:,layerIdx) - step_size * grad;
        
        back_prop           = transmission_layer.* back_prop;
        if layerIdx>1
            back_prop       = ifft2(prop_kernel.* fft2(back_prop));
        end
    end
end