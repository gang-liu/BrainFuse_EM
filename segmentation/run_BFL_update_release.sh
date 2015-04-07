#!/bin/bash
export MCR="/usr/pubsw/common/matlab/8.0"
echo $MCR
exedir=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFuseLab_release_temp/segmentation/BFL_update_release/src

declare -a subject_image=("orignal01.nii.gz" "orignal16.nii.gz" "orignal24.nii.gz" "orignal36.nii.gz" "orignal40.nii.gz")
declare -a subject_label=("lable01.nii.gz" "lable16.nii.gz" "lable24.nii.gz" "lable36.nii.gz" "lable40.nii.gz")
training_labels="lable01.nii.gz,lable16.nii.gz,lable24.nii.gz,lable36.nii.gz,lable40.nii.gz" 
training_images="orignal01.nii.gz,orignal16.nii.gz,orignal24.nii.gz,orignal36.nii.gz,orignal40.nii.gz" 
n=0;
## now loop through the above array
until [ $n == 5 ]; do 

#    echo subject_image(n)

    echo ${subject_image[n]}
    $exedir/run_BFL_update_release.sh $MCR ${subject_image[n]} ${subject_label[n]} $training_images $training_labels
#    pbsubmit -l nodes=1:ppn=4,vmem=28gb -m lukeliu -c "$exedir/run_BFL_fuison_update_sigma.sh  $MCR \"${subject_image[n]}\" \"${subject_label[n]}\" \"$training_images\" \"$training_labels\"";

((n++))

done



#  "exedir"/run_testbash.sh   $MCR "'1234','asdf'"




