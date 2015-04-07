#!/bin/bash
export MCR="/usr/pubsw/common/matlab/8.0"
echo $MCR
exedir=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_release_1015/segmentation/BFL_TEST/src

declare -a subject_image=("orignal01.nii.gz" "orignal02.nii.gz" "orignal03.nii.gz" "orignal04.nii.gz" "orignal05.nii.gz" "orignal06.nii.gz" "orignal07.nii.gz" "orignal08.nii.gz" "orignal09.nii.gz" "orignal10.nii.gz" "orignal11.nii.gz" "orignal12.nii.gz" "orignal13.nii.gz" "orignal14.nii.gz" "orignal15.nii.gz" "orignal16.nii.gz" "orignal17.nii.gz" "orignal18.nii.gz" "orignal19.nii.gz" "orignal20.nii.gz" "orignal21.nii.gz" "orignal22.nii.gz" "orignal23.nii.gz" "orignal24.nii.gz" "orignal25.nii.gz" "orignal26.nii.gz" "orignal27.nii.gz" "orignal28.nii.gz" "orignal29.nii.gz" "orignal30.nii.gz" "orignal31.nii.gz" "orignal32.nii.gz" "orignal33.nii.gz" "orignal34.nii.gz" "orignal35.nii.gz" "orignal36.nii.gz" "orignal37.nii.gz" "orignal38.nii.gz" "orignal39.nii.gz" "orignal40.nii.gz")
declare -a subject_label=("lable01.nii.gz" "lable02.nii.gz" "lable03.nii.gz" "lable04.nii.gz" "lable05.nii.gz" "lable06.nii.gz" "lable07.nii.gz" "lable08.nii.gz" "lable09.nii.gz" "lable10.nii.gz" "lable11.nii.gz" "lable12.nii.gz" "lable13.nii.gz" "lable14.nii.gz" "lable15.nii.gz" "lable16.nii.gz" "lable17.nii.gz" "lable18.nii.gz" "lable19.nii.gz" "lable20.nii.gz" "lable21.nii.gz" "lable22.nii.gz" "lable23.nii.gz" "lable24.nii.gz" "lable25.nii.gz" "lable26.nii.gz" "lable27.nii.gz" "lable28.nii.gz" "lable29.nii.gz" "lable30.nii.gz" "lable31.nii.gz" "lable32.nii.gz" "lable33.nii.gz" "lable34.nii.gz" "lable35.nii.gz" "lable36.nii.gz" "lable37.nii.gz" "lable38.nii.gz" "lable39.nii.gz" "lable40.nii.gz")
training_labels="lable01.nii.gz,lable02.nii.gz,lable03.nii.gz,lable04.nii.gz,lable05.nii.gz,lable06.nii.gz,lable07.nii.gz,lable08.nii.gz,lable09.nii.gz,lable10.nii.gz,lable11.nii.gz,lable12.nii.gz,lable13.nii.gz,lable14.nii.gz,lable15.nii.gz,lable16.nii.gz,lable17.nii.gz,lable18.nii.gz,lable19.nii.gz,lable20.nii.gz,lable21.nii.gz,lable22.nii.gz,lable23.nii.gz,lable24.nii.gz,lable25.nii.gz,lable26.nii.gz,lable27.nii.gz,lable28.nii.gz,lable29.nii.gz,lable30.nii.gz,lable31.nii.gz,lable32.nii.gz,lable33.nii.gz,lable34.nii.gz,lable35.nii.gz,lable36.nii.gz,lable37.nii.gz,lable38.nii.gz,lable39.nii.gz,lable40.nii.gz" 
training_images="orignal01.nii.gz,orignal02.nii.gz,orignal03.nii.gz,orignal04.nii.gz,orignal05.nii.gz,orignal06.nii.gz,orignal07.nii.gz,orignal08.nii.gz,orignal09.nii.gz,orignal10.nii.gz,orignal11.nii.gz,orignal12.nii.gz,orignal13.nii.gz,orignal14.nii.gz,orignal15.nii.gz,orignal16.nii.gz,orignal17.nii.gz,orignal18.nii.gz,orignal19.nii.gz,orignal20.nii.gz,orignal21.nii.gz,orignal22.nii.gz,orignal23.nii.gz,orignal24.nii.gz,orignal25.nii.gz,orignal26.nii.gz,orignal27.nii.gz,orignal28.nii.gz,orignal29.nii.gz,orignal30.nii.gz,orignal31.nii.gz,orignal32.nii.gz,orignal33.nii.gz,orignal34.nii.gz,orignal35.nii.gz,orignal36.nii.gz,orignal37.nii.gz,orignal38.nii.gz,orignal39.nii.gz,orignal40.nii.gz" 
n=30;
## now loop through the above array
until [ $n == 40 ]; do 

#    echo subject_image(n)

    echo ${subject_image[n]}
    $exedir/run_BFL_TEST.sh $MCR ${subject_image[n]} ${subject_label[n]} $training_images $training_labels
#    pbsubmit -l nodes=1:ppn=4,vmem=28gb -m lukeliu -c "$exedir/run_BFL_TEST.sh $MCR \"${subject_image[n]}\" \"${subject_label[n]}\" \"$training_images\" \"$training_labels\"";

((n++))

done



#  "exedir"/run_testbash.sh   $MCR "'1234','asdf'"




