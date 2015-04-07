#!/bin/bash
export MCR="/usr/pubsw/common/matlab/8.0"
echo $MCR
exedirorg=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_release_1015/segmentation/BFL_TEST/src


outputdir2=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_release_1015/r5/results17

train_images=orignal23.nii.gz,orignal25.nii.gz,orignal19.nii.gz,orignal17.nii.gz,orignal01.nii.gz
train_labels=lable23.nii.gz,lable25.nii.gz,lable19.nii.gz,lable17.nii.gz,lable01.nii.gz

numofsub=4
# outputdir2=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_release_1015/r"$atlas_num"/results"$randomtimes"/
declare -a subject_image=("orignal14.nii.gz" "orignal22.nii.gz" "orignal20.nii.gz" "orignal06")
declare -a subject_label=("lable14.nii.gz" "lable22.nii.gz" "lable20.nii.gz" "lable06.nii.gz")
n=0
until [ $n == $numofsub ]; do

     $exedirorg/run_BFL_TEST.sh $MCR ${subject_image[n]} ${subject_label[n]} "$train_images"  "$train_labels" "$outputdir2";

    (( n++ ));
    sleep 1;
done
exit




