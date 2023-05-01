#!/bin/bash

#for i in {1..100}
#do
#	bash /home/malosree/projects/def-gturecki/malosree/Permutation_analysis/bash_files/Step3.sh $i
#	bash /home/malosree/projects/def-gturecki/malosree/Permutation_analysis/bash_files/Step4.sh $i
#done

#cluster=(100 10 27 28 31 32 35 37 38 40 41 46 49 51 64 68 71 73 75 78 89 8 90 92 95 96 98 99)
#broad=(38 41 42 48 64 71 73 74 79 85 88 96 98)

#cluster=(46 49 64 71 75 78 89 96)
#cluster=(75 89)
cluster=(89)
for i in ${cluster[@]}
do
	 bash /home/malosree/projects/def-gturecki/malosree/Permutation_analysis/bash_files/Step3.sh $i
done
#for i in ${broad[@]}
#do
#	 bash /home/malosree/projects/def-gturecki/malosree/Permutation_analysis/bash_files/Step4.sh $i
#done
