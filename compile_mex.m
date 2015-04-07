%ITK_DIR = '/Users/local_msabuncu/Documents/work/code/ITK/';

ITK_DIR = '/autofs/cluster/con_003/users/mert/pubsw/InsightToolkit-3.20.0/';

if (~strcmp(comp_id(1:4), 'GLNX'))
    cd ./diffeomorphic_registration/
    compile_mex_example(ITK_DIR);
    cd ..
end
cd ./misc/toolbox_fast_marching/
compile_mex
cd ../..