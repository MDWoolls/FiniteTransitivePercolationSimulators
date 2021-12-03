#sbatch --output="out3.txt" --job-name="perc3" run3.sh
#sbatch --output="outRAN.txt" --job-name="RAN" singlerun.sh
#sbatch --output="out2.txt" --job-name="perc2" run2.sh
#sbatch --output="out44.txt" --job-name="torrus" createtorrus.sh
#sbatch --output="out1.txt" --job-name="perc1" run1.sh
#sbatch --output="out443.txt" --job-name="torrus3" runtorrus.sh

sbatch --output="out_6_1.txt" --job-name="6_1_hypercube" singlerun2.sh 6 1 8 11
sbatch --output="out_6_2.txt" --job-name="6_2_hypercube" singlerun2.sh 6 2 9 11
sbatch --output="out_6_3.txt" --job-name="6_3_hypercube" singlerun2.sh 6 3 4 11

sbatch --output="outfrd_6_1.txt" --job-name="6_1_frdhypercube" singlerun.sh 6 1 8 11
sbatch --output="outfrd_6_2.txt" --job-name="6_2_frdhypercube" singlerun.sh 6 2 9 11
sbatch --output="outfrd_6_3.txt" --job-name="6_3_frdhypercube" singlerun.sh 6 3 4 11
