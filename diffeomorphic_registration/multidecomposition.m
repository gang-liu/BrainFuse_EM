function [pyramid1,pyramid2] = multidecomposition(fix_im,mov_im,options)
if (~isfield(options,'num_multires'))
    options.num_multires = 4;
end

numOfLevels = options.num_multires;

pyramid1 = cell(numOfLevels, 1);

pyramid2 = cell(numOfLevels, 1);

pyramid1{1,1} = fix_im;

pyramid2{1,1} = mov_im;



for level = 2:1:numOfLevels



    size1 = size(pyramid1{level-1,1});

    size1 = size1 - mod(size1,2);

    pyramid1{level,1} = pyramid1{level-1,1}(1:2:size1(1),1:2:size1(2), 1:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(2:2:size1(1),1:2:size1(2),1:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(1:2:size1(1),2:2:size1(2),1:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(1:2:size1(1),1:2:size1(2),2:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(1:2:size1(1),2:2:size1(2),2:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(2:2:size1(1),1:2:size1(2),2:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(2:2:size1(1),2:2:size1(2),1:2:size1(3));

    pyramid1{level,1} = pyramid1{level,1} + pyramid1{level-1,1}(2:2:size1(1),2:2:size1(2),2:2:size1(3));



    pyramid1{level,1} = pyramid1{level,1}/8;



    size2 = size(pyramid2{level-1,1});

    size2 = size2 - mod(size2,2);

    pyramid2{level,1} = pyramid2{level-1,1}(1:2:size2(1),1:2:size2(2), 1:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(2:2:size2(1),1:2:size2(2),1:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(1:2:size2(1),2:2:size2(2),1:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(1:2:size2(1),1:2:size2(2),2:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(1:2:size2(1),2:2:size2(2),2:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(2:2:size2(1),1:2:size2(2),2:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(2:2:size2(1),2:2:size2(2),1:2:size2(3));

    pyramid2{level,1} = pyramid2{level,1} + pyramid2{level-1,1}(2:2:size2(1),2:2:size2(2),2:2:size2(3));



    pyramid2{level,1} = pyramid2{level,1}/8;



end
