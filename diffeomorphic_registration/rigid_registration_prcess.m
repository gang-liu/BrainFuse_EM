function options = rigid_registration_prcess(fix_im,mov_im,options)

numOfLevels = options.num_multires;
if (~isfield(options,'rigidFlag'))
    options.rigidFlag = 1;
end

if (~isfield(options, 'rigidParams'))
    options.rigidParams = [];
end

if ~isfield(options,'log_def_x')
    
   if (options.rigidFlag)
    if options.verbose
        display('Performing initial affine alignment...');
        tic
    end
    
    if (isempty(options.rigidParams))
        options_affine = options;
        options_affine.num_multires = 3;
        [affineMat, affineParams] = ...
            BFL_pairwise_affine_reg3D_multires(fix_im, mov_im, options_affine);
    else
        affineParams = options.rigidParams;
        affineMat = AffineParams2Mat_aux(affineParams);
    end

    [log_def_x, log_def_y, log_def_z] = AffineMat2VelocityField3D_aux(affineMat, size(fix_im));
    options.log_def_x = log_def_x;
    options.log_def_y = log_def_y;
    options.log_def_z = log_def_z;
    if options.verbose
        display('Initial affine alignment... done');
        display(['Translations (voxels): ' num2str(affineParams(7)) ' ' num2str(affineParams(8)) ' ' num2str(affineParams(9))]);
        display(['Rotations (rad): ' num2str(affineParams(4)) ' ' num2str(affineParams(5)) ' ' num2str(affineParams(6))]);
        display(['Log-Scale: ' num2str(affineParams(1)) ' ' num2str(affineParams(2)) ' ' num2str(affineParams(3))]);
        
        toc
    end
   else
    options.log_def_x = zeros(size(pyramid1{1,1}),class(pyramid1{numOfLevels,1}));
    options.log_def_y = zeros(size(pyramid1{1,1}),class(pyramid1{numOfLevels,1}));
    options.log_def_z = zeros(size(pyramid1{1,1}),class(pyramid1{numOfLevels,1}));
   end
end


pyramid_log_def_x = cell(numOfLevels, 1);
pyramid_log_def_y = cell(numOfLevels, 1);
pyramid_log_def_z = cell(numOfLevels, 1);

pyramid_log_def_x{1} = options.log_def_x;
pyramid_log_def_y{1} = options.log_def_y;
pyramid_log_def_z{1} = options.log_def_z;


for level = 2:1:numOfLevels



    size1 = size(pyramid_log_def_x{level-1,1});

    size1 = size1 - mod(size1,2);

    pyramid_log_def_x{level,1} = pyramid_log_def_x{level-1,1}(1:2:size1(1),1:2:size1(2), 1:2:size1(3));

    pyramid_log_def_x{level,1} = pyramid_log_def_x{level,1} + pyramid_log_def_x{level-1,1}(2:2:size1(1),1:2:size1(2),1:2:size1(3));

    pyramid_log_def_x{level,1} = pyramid_log_def_x{level,1} + pyramid_log_def_x{level-1,1}(1:2:size1(1),2:2:size1(2),1:2:size1(3));

    pyramid_log_def_x{level,1} = pyramid_log_def_x{level,1} + pyramid_log_def_x{level-1,1}(1:2:size1(1),1:2:size1(2),2:2:size1(3));

    pyramid_log_def_x{level,1} = pyramid_log_def_x{level,1} + pyramid_log_def_x{level-1,1}(1:2:size1(1),2:2:size1(2),2:2:size1(3));

    pyramid_log_def_x{level,1} = pyramid_log_def_x{level,1} + pyramid_log_def_x{level-1,1}(2:2:size1(1),1:2:size1(2),2:2:size1(3));

    pyramid_log_def_x{level,1} = pyramid_log_def_x{level,1} + pyramid_log_def_x{level-1,1}(2:2:size1(1),2:2:size1(2),1:2:size1(3));

    pyramid_log_def_x{level,1} = pyramid_log_def_x{level,1} + pyramid_log_def_x{level-1,1}(2:2:size1(1),2:2:size1(2),2:2:size1(3));



    pyramid_log_def_x{level,1} = pyramid_log_def_x{level,1}/8/2;


    pyramid_log_def_y{level,1} = pyramid_log_def_y{level-1,1}(1:2:size1(1),1:2:size1(2), 1:2:size1(3));

    pyramid_log_def_y{level,1} = pyramid_log_def_y{level,1} + pyramid_log_def_y{level-1,1}(2:2:size1(1),1:2:size1(2),1:2:size1(3));

    pyramid_log_def_y{level,1} = pyramid_log_def_y{level,1} + pyramid_log_def_y{level-1,1}(1:2:size1(1),2:2:size1(2),1:2:size1(3));

    pyramid_log_def_y{level,1} = pyramid_log_def_y{level,1} + pyramid_log_def_y{level-1,1}(1:2:size1(1),1:2:size1(2),2:2:size1(3));

    pyramid_log_def_y{level,1} = pyramid_log_def_y{level,1} + pyramid_log_def_y{level-1,1}(1:2:size1(1),2:2:size1(2),2:2:size1(3));

    pyramid_log_def_y{level,1} = pyramid_log_def_y{level,1} + pyramid_log_def_y{level-1,1}(2:2:size1(1),1:2:size1(2),2:2:size1(3));

    pyramid_log_def_y{level,1} = pyramid_log_def_y{level,1} + pyramid_log_def_y{level-1,1}(2:2:size1(1),2:2:size1(2),1:2:size1(3));

    pyramid_log_def_y{level,1} = pyramid_log_def_y{level,1} + pyramid_log_def_y{level-1,1}(2:2:size1(1),2:2:size1(2),2:2:size1(3));



    pyramid_log_def_y{level,1} = pyramid_log_def_y{level,1}/8/2;

    
    pyramid_log_def_z{level,1} = pyramid_log_def_z{level-1,1}(1:2:size1(1),1:2:size1(2), 1:2:size1(3));

    pyramid_log_def_z{level,1} = pyramid_log_def_z{level,1} + pyramid_log_def_z{level-1,1}(2:2:size1(1),1:2:size1(2),1:2:size1(3));

    pyramid_log_def_z{level,1} = pyramid_log_def_z{level,1} + pyramid_log_def_z{level-1,1}(1:2:size1(1),2:2:size1(2),1:2:size1(3));

    pyramid_log_def_z{level,1} = pyramid_log_def_z{level,1} + pyramid_log_def_z{level-1,1}(1:2:size1(1),1:2:size1(2),2:2:size1(3));

    pyramid_log_def_z{level,1} = pyramid_log_def_z{level,1} + pyramid_log_def_z{level-1,1}(1:2:size1(1),2:2:size1(2),2:2:size1(3));

    pyramid_log_def_z{level,1} = pyramid_log_def_z{level,1} + pyramid_log_def_z{level-1,1}(2:2:size1(1),1:2:size1(2),2:2:size1(3));

    pyramid_log_def_z{level,1} = pyramid_log_def_z{level,1} + pyramid_log_def_z{level-1,1}(2:2:size1(1),2:2:size1(2),1:2:size1(3));

    pyramid_log_def_z{level,1} = pyramid_log_def_z{level,1} + pyramid_log_def_z{level-1,1}(2:2:size1(1),2:2:size1(2),2:2:size1(3));



    pyramid_log_def_z{level,1} = pyramid_log_def_z{level,1}/8/2;

   
end

options.log_def_x = double(pyramid_log_def_x{numOfLevels,1});
options.log_def_y = double(pyramid_log_def_y{numOfLevels,1});
options.log_def_z = double(pyramid_log_def_z{numOfLevels,1});

clear pyramid_log_def_x pyramid_log_def_y pyramid_log_def_z

return;