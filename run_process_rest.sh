#!/bin/bash
export MCR="/usr/pubsw/common/matlab/8.0"
echo $MCR
# exedirorg=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_release_1015/segmentation/BFL_TEST/src
# exedirorg2=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFuseLab_update_sigma/segmentation/BFL_fuison_update_sigma/src
weighted_multisolution_dir=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_release_1015/segmentation/BFL_labelfusion_resolution_weighted/src
noweighted_multisolution_dir=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_release_1015/segmentation/BFL_labelfusion_resolution_noweighted/src
opsweighted_multisolution_dir=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_release_1015/segmentation/BFL_labelfusion_resolution_opsweighted/src


NUMPARAMS=$#
randomtimes=$1

atlas_num=$2;
numofsub=4


# outputdir2=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_release_1015/r"$atlas_num"/results"$randomtimes"/
# outputdir3=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFuseLab_update_sigma/r"$atlas_num"/results"$randomtimes"/

outputdir_weighted=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/results20141210/weighted_multireso/r"$atlas_num"/results"$randomtimes"/
outputdir_noweighted=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/results20141210/noweighted_multireso/r"$atlas_num"/results"$randomtimes"/
outputdir_opsweighted=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/results20141210/opsweighted_multireso/r"$atlas_num"/results"$randomtimes"/
mkdir "$outputdir_weighted" "$outputdir_noweighted" "$outputdir_opsweighted"

n=0;


filecontent=( `cat "/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_release_1015/r"$atlas_num"/results"$randomtimes"/atlasfilename.txt" `)

train_images=${filecontent[6]}
train_labels=${filecontent[7]}

subject_image=${filecontent[@]:15:4}
subject_label=${filecontent[@]:19:4}

subject_image=($subject_image)
subject_label=($subject_label)


until [ $n == $numofsub ]; do
     echo $exedirorg/run_BFL_TEST.sh $MCR ${subject_image[n]} ${subject_label[n]} "$train_images"  "$train_labels"  "$outputdir2"


#    $exedirorg/run_BFL_TEST.sh $MCR ${subject_image[n]} ${subject_label[n]} "$train_images"  "$train_labels"  "$outputdir2" 

#    pbsubmit -l nodes=1:ppn=3,vmem=21gb -m lukeliu -c "$exedirorg/run_BFL_TEST.sh   $MCR \"${subject_image[n]}\" \"${subject_label[n]}\" \"${train_images}\" \"${train_labels}\" \"${outputdir2}\"";
    pbsubmit -l nodes=1:ppn=3,vmem=21gb -m lukeliu -c "$weighted_multisolution_dir/run_BFL_labelfusion_resolution_weighted.sh   $MCR \"${subject_image[n]}\" \"${subject_label[n]}\" \"${train_images}\" \"${train_labels}\" \"${outputdir_weighted}\"";
    pbsubmit -l nodes=1:ppn=3,vmem=21gb -m lukeliu -c "$noweighted_multisolution_dir/run_BFL_labelfusion_resolution_noweighted.sh   $MCR \"${subject_image[n]}\" \"${subject_label[n]}\" \"${train_images}\" \"${train_labels}\" \"${outputdir_noweighted}\"";
    pbsubmit -l nodes=1:ppn=3,vmem=21gb -m lukeliu -c "$opsweighted_multisolution_dir/run_BFL_labelfusion_resolution_opsweighted.sh   $MCR \"${subject_image[n]}\" \"${subject_label[n]}\" \"${train_images}\" \"${train_labels}\" \"${outputdir_opsweighted}\"";



#    pbsubmit -l nodes=1:ppn=3,vmem=21gb -m lukeliu -c "$exedirorg2/run_BFL_fuison_update_sigma.sh   $MCR \"${subject_image[n]}\" \"${subject_label[n]}\" \"${train_images}\" \"${train_labels}\" \"${outputdir3}\""; 

(( n++ ));
    sleep 1;
done
exit










