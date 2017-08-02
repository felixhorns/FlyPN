for i in $(seq 0 901)
do
    sbatch validate.sh $i validate_GH146_CombinatorialCode.out
done
	
