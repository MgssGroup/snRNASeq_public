#!/bin/bash

#for i in {1..100}
#cluster=(100 17 19 26 38 3 40 41 47 49 52 55 57 58 61 62 63 65 68 70 78 79 96 97 99 9)
#cluster=(100 49)
cluster=(49)
for i in ${cluster[@]}
do
	bash /home/malosree/projects/def-gturecki/malosree/Permutation_analysis/bash_files/Step1.sh $i
	#bash /home/malosree/projects/def-gturecki/malosree/Permutation_analysis/bash_files/Step2.sh $i
done
