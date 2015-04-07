#!/bin/bash

#  testbash.sh  randomtime_number atlas_number

exe_dir=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/BrainFusionLab_release_1015
# exe_dir2=/autofs/cluster/con_003/users/mert/lukeliu/ARRANGE/pubsw/codeupdate/testbash
n=1
until [ $n == 21 ]; do
# "$exe_dir2"/originalsigma_difftest.sh $n 5 9 
  "$exe_dir"/run_process_rest.sh $n 10 


echo $n
# the first parameter is the times it will be repeatly performed, and the second parameter is atlas number.
# the third paramter is sigma_diff
(( n++ ));
done


