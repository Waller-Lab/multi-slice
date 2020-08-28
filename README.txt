This folder contains data collected via an angle-scanning optical imaging system on 11/2/2018 of a C. elegans sample. Details of this optical imaging system are contained in the paper:
S. Chowdhury, M. Chen, R. Eckert, D. Ren, F. Wu, N. Repina, and L. Waller, "High-resolution 3D refractive index microscopy of multiple-scattering samples from intensity images," Optica 6, 1211-1219 (2019) 



Folder 'FOV_01_reconstruction_CElegan_resources' contains helpful resources pertaining to the C. elegans data.   The folder contains:

'FOV_01_rawdata_CElegan.TIF':                       raw image stack collected by the imaging system. You can open it in ImageJ and scroll through the stack to get a feel for experimental data.   
'FOV_01_reconstruction_CElegan_params.mat':         contains parameters pertaining to the imaging system acquisition protocol and imaging parameters.

'FOV_01_reconstruction_CElegan.mat':                .MAT file containing workspace variables 'totalFOV', 'totalFOV_hc', and 'cmap'
                         'totalFOV':                raw reconstruction of the C. elegans head, and is basically the data corresponding to the gray-scale C. elegan figure in the Optica paper.
                      'totalFOV_hc':                halo-corrected version of 'totalFOV', corresponding to the colored C. elegans figure in the Optica paper
                             'cmap':                colormap used for the colored figures in the Optica paper

These figures were reconstructed with old multi-slice beam-propagation (MSBP) code from a while ago, when we were writing the Optica paper. Updated MSBP code outputs better reconstructions, which are in the following folders:




Folder 'FOV_01_reconstruction_CElegan_outputs_V1' contains reconstructed data (using updated MSBP code) on 16 specific patches of the C. elegans data from the .TIF file, dividing the whole field-of-view into a 4x4 grid:

'FOV_01_reconstruction_CElegan_output_patch_'+index+'.mat':     .MAT files containing reconstruction outputs after running MSBP algorithm on a patch of the C. elegans data indexed via variable 'index'. In this folder, 'index' ranges from 1 to 16 
                                                'reconObj':     Important workspace variable in the .MAT file. It is the 3D refractive index distribution of the indexed patch, as reconstructed via MSBP

'makeTotalFOV_V1.m':                                            MATLAB script that patches together the 'reconObj' 3D RI reconstruction outputs in each .MAT file, based on the 'index'
'FOV_01_reconstruction_CElegan_totalFOV_notuning.mat':          .MAT file that contains the output ('totalFOV') after running the patch script above on all 16 patches              
'FOV_01_reconstruction_CElegan_totalFOV.mat':                   .MAT file that contains the output ('totalFOV') after running the patch script above on all 16 patches, with some global tuning on the specific patches of homegenous RI




Folder 'FOV_01_reconstruction_CElegan_outputs_V2' contains reconstructed data on 2 larger patches that encompass the region of the FOV containing C. elegans data. Patches are NOT aligned on a grid in this case.

'backgroundSubtraction.m':              MATLAB script to conduct a low-order polynomial background subtraction on reconstruction patches. The code in 
'makeTotalFOV_V2.m':                    MATLAB script that combines reconstruction patches together

'FOV_01_reconstruction_CElegan_output_patch_'+index+'.mat':                 .MAT files containing reconstruction outputs after running MSBP algorithm on a patch of the C. elegans data indexed via variable 'index'. In this folder, 'index' is either 1 or 2. 
'FOV_01_reconstruction_CElegan_output_patch_'+index+'bkgdSubtrct.mat':      Same as .MAT file above, but after undergoing background subtraction
                                                           'reconObj':      Important workspace variable in the .MAT file. It is the 3D refractive index distribution of the indexed patch, as reconstructed via MSBP   
                                                           
'FOV_01_reconstruction_CElegan_totalFOV_bkgdSubtrct.mat':                   .MAT file that contains the output ('totalFOV') after running the patch script above on the two background-subtracted patches.