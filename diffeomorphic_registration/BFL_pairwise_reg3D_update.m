function [log_def_x, log_def_y, log_def_z, stats, warped_mov_im, backwarped_fix_im] = ...
    BFL_pairwise_reg3D(fix_im, mov_im, options)
% function [log_def_x, log_def_y, log_def_z, stats, warped_mov_im, backwarped_fix_im] = ...
%    BFL_pairwise_reg3D(fix_im, mov_im, options = [])

%input: 
% fix_im: a 3D double matrix (e.g. MRI volume of size 256x256x256)
% mov_im: a 3D double matrix (e.g. MRI volume of size 256x256x256)
%
% NOTE: This algorithm assumes fix_im and mov_im are of the same size!!!
% It implicitly assumes that the voxels are isometric (e.g. 1 mm^3). It
% will work if this is not the case, but the result may not be optimal.
% All testings of this implementation have been done on norm.mgz files
% created by FS. These are 1 mm^3 resolution, 256x256x256 size,
% skull-stripped, intensity normalized brain MRI scans.
%
%options: a structure with following fields:
%   *num_multires = <scalar, integer> num of levels in multi-resolution pyramid
%   *numiter = <scalar, integer> the number of (max.) "gradient descent" iterations the
%   algorithm will take at each resolution level. (default: 50)
%   *numiter_cell = <cell of scalars of size num_multires> alternative to numiter you can directly input the number of iterations
%   for each resolution level seperately. numiter_cell{1} will be the max.
%   number of iterations taken at the highest-resolution (e.g. the original
%   resolution of the images).
%   *sigma_diff = <scalar> determines the smoothness of the warp -- the higher the
%   smoother. this is the std. of the gaussian kernel used to smooth the velocity field.
%   the smallest it should be is probably 1. for many
%   applications we have found 2 to be optimal. (default: 2). 
%   *reg_weight = <scalar> this effectively determines the step size in the
%   optimization (gradient descent). a higher value will yield smaller
%   steps. (default: 150)
%   *log_def_x ( and log_def_y, log_def_z) = <a 3D matrix (size of
%   fix_im)>
%   all three of these should be input together. these
%   are initializations for the velocity fields (i.e. the parametrization
%   of the warp). 
%   *verbose = <scalar, 0 or nonzero> if nonzero, will display some
%   results. (default: 0)


% output:
% * [log_def_x, log_def_y, log_def_z] = 3 x <a 3D double matrix (size of
%   fix_im)> the velocity field of the output warp
% * stats = a structure with following fields:
%   MSE (mean square error between fixed image and warped moving image), backMSE (MSE between moving image and backward warp fixed image), 
%   harmoEner (harmonic energy of "forward" warp), backharmoEner (harmonic
%   energy of "backward" warp
% * warped_mov_im = the moving image warped (resampled onto the grid of the
% fixed image).
% * backwarped_fix_im = the fixed image resampled onto the moving image
% grid

%=========================================================================
%  Brain Fuse Lab 
%  Copyright (c) 2010-2011 Mert R. Sabuncu
%  A.A. Martinos Center for Biomedical Imaging, 
%  Mass. General Hospital, Harvard Medical School
%  CSAIL, Mass. Institute of Technology 
%  All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%    * Neither the names of the copyright holders nor the names of future
%      contributors may be used to endorse or promote products derived from this
%      software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
%ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
%
%=========================================================================


if (nargin<3)
    options = [];
end

if (~isfield(options,'num_multires'))
    options.num_multires = 4;
end
if (~isfield(options, 'numiter'))
    options.numiter = 150;
end

if (~isfield(options, 'rigidFlag'))
    options.rigidFlag = 1;
end


if (~isfield(options, 'sigma_diff'))
    options.sigma_diff = 2;
end

if ~isfield(options, 'verbose')
    options.verbose = 0;
end

if (~isfield(options, 'numiter_cell'))
    options.numiter_cell = cell(options.num_multires,1);
    for i = 1:options.num_multires
        options.numiter_cell{i,1} = options.numiter;
    end
end

if (~isfield(options, 'rigidParams'))
    options.rigidParams = [];
end



if ~isfield(options, 'reg_weight')
    options.reg_weight = 50; %% this determines the step size in the optimization/registration -- higher yields smaller step size
end


%%%% this is the minimum level at which registration is performed.
% 1 corresponds to the highest resolution -- i.e. the full resolution of
% images
% if your images are too big and your computational resources are limited,
% you might want to perform registration only upto a lower resoltion
% representation (ie based on sub-sampled images only).

if (~isfield(options, 'min_level'))
    options.min_level = 1;
end




stats = cell(numOfLevels,1);

for level = numOfLevels: -1: 1

    if (level < options.min_level)
        options.numiter = 0;
    else
        options.numiter = options.numiter_cell{level};
    end
    
    [log_def_x, log_def_y, log_def_z, stats{level,1}, warped_mov_im, backwarped_fix_im] = ...
        invconstdemonsreg3d_aux(double(pyramid1{level,1}), double(pyramid2{level,1}), options);
    
    if (level ~= 1)
        options.log_def_x = double(upscale(log_def_x, size(pyramid1{level - 1,1})));
        options.log_def_y = double(upscale(log_def_y, size(pyramid1{level - 1,1})));
        options.log_def_z = double(upscale(log_def_z, size(pyramid1{level - 1,1})));
    end
    display(['Finished registration at pyramid level: ' num2str(level)]);
    
end

return;




