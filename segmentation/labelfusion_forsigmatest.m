%function BFL_labelfusion_3D(SBJ_CELL, DATA_DIR, Output_dir, options)
function labelfusion_forsigmatest(varargin)

DATA_DIR = '/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFuseLab_original/data/';
% OutpRut_dir = '/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFuseLab_original/results/';

% addpath('../diffeomorphic_registration/');
% addpath('../misc/');
% addpath('../misc/toolbox_fast_marching/');
% addpath('../misc/toolbox_fast_marching/toolbox/');
% addpath('../misc/fminsearchbnd/');
% addpath('../segmentation/')

if (nargin ~=0)

    SBJ_CELL = varargin(1);
    SBJ_SEG_CELL = varargin(2);
    training_image_cell = varargin(3);
    training_image_cell = strsplit(training_image_cell{1},',');
    training_image_cell = setdiff(training_image_cell,SBJ_CELL);
    training_label_cell = varargin(4);
    training_label_cell = strsplit(training_label_cell{1},',');
    training_label_cell = setdiff(training_label_cell,SBJ_SEG_CELL);
    Output_dir = varargin(5);
    temp = varargin(6);
    sigma = str2num(temp{1})
else
    SBJ_CELL = {};
    SBJ_CELL{1} = 'orignal01.nii.gz';
    SBJ_SEG_CELL = {};
    SBJ_SEG_CELL = 'label01.nii.gz';
    training_image_cell = 'orignal01.nii.gz,orignal02.nii.gz,orignal04.nii.gz,orignal03.nii.gz';
    training_image_cell = strsplit(training_image_cell,',');
    training_image_cell = setdiff(training_image_cell,SBJ_CELL);
    
    training_label_cell = 'lable01.nii.gz,lable16.nii.gz,lable24.nii.gz,lable36.nii.gz,lable40.nii.gz';
    training_label_cell = strsplit(training_label_cell,',');
    training_label_cell = setdiff(training_label_cell,SBJ_SEG_CELL);
   sigma = 2; 
end

options.verbose = 1;

options.num_multires = 4;
options.numiter_cell = cell(options.num_multires,1);
options.reg_weight = 20;
options.TRAINING_DATA_DIR = DATA_DIR;
for i = 3:options.num_multires
    options.numiter_cell{i,1} = 150;
end
options.numiter_cell{2,1} = 30;
options.numiter_cell{1,1} = 150;
options.rigidFlag = 1; %%% This flag (if nonzero) tells the registration algorithm to pre-register w/ an affine transformation model

if (~isfield(options, 'min_level'))
    options.min_level = 1;
end
%%% this parameter determines the "nonlinearity/flexibility" of the registration
%typically it should be an integer greater than zero
%the smaller it is, the more nonlinear, deformable the registration model
%is
options.sigma_diff = sigma;

%the following parameter determines the weighting in the weighted local
%label fusion -- if it is very high (eg approaches infinitiy) then all
%training subjects get the same weight. smaller values make the training
%images that have a local similar appearance to the test image get weighed
%more.
options.sigma_labelfusion = 5;

% input:
% *SBJ_CELL : a cell with filenames of subject image data
% * DATA_DIR: the directory where the subject image data resides
% * Output_dir: the directory where the results will be written.
% * options: a structure with following potential fields:
%
%   (These are for segmentation)
%   *TRAINING_DATA_DIR: the directory where the training data data resides
%   *TRAINING_IMAGE_CELL: filenames for training images (eg)
%   *TRAINING_SEG_CELL: filenames for training label images (ie manual
%   segmentations)
%   * sigma_labelfusion : the sigma parameter that determines the weighting
%   for locally weighted label fusion (see [1]). If this is high, label
%   fusion is closer to majority voting where all training subjects have an
%   equal weight. if this is low, then training images that have an
%   intensity closer to the test image (the image to be segmented) have a
%   higher weight. (DEFAULT VAL: 10)

%   * dt_weight: this is the weight used for the distance transform --
%   which determines the "sharpness" of the propagated labels in label
%   fusion. (See [1], where this parameter is referred to as rho) (default:
%   1)
%   background_label = <int> the index of the background label (default: 0)

%   (These are for registration options)
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
%
%
%  This code is an implementation of "local label fusion" as described in
%  the following paper:
%  [1] "A Generative Model for Image Segmentation Based on Label Fusion," by MR
%  Sabuncu, BTT Yeo, K Van Leemput, P Golland and B Fischl. IEEE Transactions on Medical Imaging, vol. 29, no. 10, 2010.
%
%  There is one result per input subject and that is the automatic
%  segmentation. The file will be called: [Output_dir '/' SBJ_CELL{1,i} '_BFL_seg.mgz']
%




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


if (~isfield(options,'TRAINING_IMAGE_CELL'))
    options.TRAINING_IMAGE_CELL = training_image_cell;
end

if (~isfield(options,'TRAINING_SEG_CELL'))
    options.TRAINING_SEG_CELL = training_label_cell;
end

if (~isfield(options,'sigma_labelfusion'))
    options.sigma_labelfusion = 5;
end

if (~isfield(options,'dt_weight'))
    options.dt_weight = 1;
end

if (~isfield(options,'background_label'))
    options.background_label = 0;
end

if (~isfield(options,'outfile_postfix'))
    options.outfile_postfix = '';
end


NumTrainingSubjects = length(options.TRAINING_IMAGE_CELL);

if (length(options.TRAINING_IMAGE_CELL) ~= NumTrainingSubjects)
    error('The number of training images and training segmentations should be the same!');
end

NumTestSubjects = length(SBJ_CELL);

if (~isfield(options,'labels'))
    
    for s = 1:NumTrainingSubjects
        seg2 = MRIread([options.TRAINING_DATA_DIR '/' options.TRAINING_SEG_CELL{s}]);
        if (s == 1)
            labels = unique(seg2.vol(:));
        else
            labels = intersect(labels, unique(seg2.vol(:)));
        end
    end
else
    labels = options.labels;   
end

labels = union(labels, options.background_label);

numLabels = length(labels);


for s = 1:NumTestSubjects
    
    display(['Label fusion for : ' SBJ_CELL{s}]);
    vol1 = MRIread([DATA_DIR '/' SBJ_CELL{s}]);
    labelprob = zeros([size(vol1.vol), numLabels], 'single');
    
    for i = 1:NumTrainingSubjects
        
        display(['Processing training subject: ' options.TRAINING_IMAGE_CELL{i}]);
        
        vol2 = MRIread([options.TRAINING_DATA_DIR '/' options.TRAINING_IMAGE_CELL{i}]);
        seg2 = MRIread([options.TRAINING_DATA_DIR '/' options.TRAINING_SEG_CELL{i}]);
        
        if (isfield(options, 'rigidParamsCell'))
            options.rigidParams = options.rigidParamsCell{s};
        end

        [log_def_x, log_def_y, log_def_z, stats, wm] = ...
            BFL_pairwise_reg3D((vol1.vol), (vol2.vol), options);
        
        [def_x, def_y, def_z] = velocityfieldexp(single(log_def_x), single(log_def_y), single(log_def_z));
        
        mri_weight = exp(-(wm - vol1.vol).^2/2/options.sigma_labelfusion^2);
        
        dt_bk = zeros(size(vol1.vol), 'single');
        curlabelprob = zeros([size(vol1.vol), numLabels], 'single');
    
        for l = 1:numLabels
            cur_label = labels(l);
            if (cur_label ~= options.background_label)

                dt = fast_compute_distance_transform(seg2.vol, cur_label);
                
                dt = warpimage(single(dt), single(def_x), single(def_y),single(def_z));
                
                curlabelprob(:,:,:, l) = (dt.*options.dt_weight);
            else
                curlabelprob(:,:,:, l) = (dt_bk.*options.dt_weight);
            end
            
        end
        curlabelprob = curlabelprob - repmat(max(curlabelprob, [], 4), [1 1 1 numLabels]);
        curlabelprob = exp(curlabelprob);
        curlabelprob = curlabelprob./repmat(sum(curlabelprob, 4), [1 1 1 numLabels]);
        labelprob = labelprob + repmat(mri_weight, [1 1 1 numLabels]).*curlabelprob;

    end

    [Y, I] = max(labelprob, [], 4);
    snip_seg = (labels(I));
    clear Y I
    seg1 = seg2;
    seg1.vol = snip_seg;
    MRIwrite(seg1, [Output_dir{1} '/' SBJ_CELL{s} '_BFL_seg' options.outfile_postfix '.mgz']);
end


