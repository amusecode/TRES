#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=astera      # Request nodes from the specific partition
#SBATCH --mem 250mb
#SBATCH --time 1:00:00
#SBATCH --job-name=TRES
#SBATCH --output=logs/array_%A_%a.out
#SBATCH --array=0-9
#SBATCH --output=TRES_job.%J.out
#SBATCH --error=TRES_job.%J.err

function clean_up {
    echo "### Running Clean_up ###"
    # all files should be deleted from the compute-node before job end
    # if there are more than one output files in the directory after the completion of one of the simulations, only delete the file
    # corresponding to that simulation. Otherwise, remove the whole directory    
    files=$(ls *)
    echo $files   
    num_files=$(ls *.hdf | wc -l)
    echo "$num_files"
    #if [ "$num_files" -eq "1" ]
    #then
    #    rm -rf "/hddstore/$USER"
    #    echo "Removed all files"
    #else
    #    #rm -rf "/hddstore/$USER/$SLURM_ARRAY_TASK_ID"
    #    rm "/hddstore/$USER/$FILE_NAME" 
    #    echo 'Removed'"$FILE_NAME"''
    #fi

    # copy the data from the node to the storage directory
    cp $FILE_NAME /zfs/helios/filer0/$USER/$STORAGE_FOLDER/
    rm "/hddstore/$USER/$FILE_NAME" 
    echo "Finished"; date
    
    # exit the script
    exit
}

# call "clean_up" function when this script exits, it is run even if
# SLURM cancels the job
trap 'clean_up' EXIT

##### pipeline below #####


# enter virtual environment
# this environment contains the necessary python packages 
module purge 
module load anaconda3/python3.7.3
source /home/$USER/Amuse-env/bin/activate

# create a directory on the node to store your data
mkdir -p /hddstore/$USER
OUTPUT_FOLDER=$(mktemp -d -p /hddstore/$USER)

# for each sbatch array, define your output filename
FILE_NAME='TRES_'"$SLURM_ARRAY_TASK_ID"'.hdf'
echo   $OUTPUT_FOLDER $FILE_NAME $SLURMD_NODENAME.
cd     $OUTPUT_FOLDER 

# create a directory where the data should be stored (for helios this is the /zfs/helios/filer0/ directory)       
STORAGE_FOLDER="TRES_data/TPS"
mkdir -p /zfs/helios/filer0/$USER/$STORAGE_FOLDER       

#perform simulation
TRES_DIR=/home/stoonen1/TRES
NUM_PER_JOB=3
echo python $TRES_DIR/TPS.py -n $NUM_PER_JOB -N $((NUM_PER_JOB*SLURM_ARRAY_TASK_ID)) -f $FILE_NAME -T 1
python $TRES_DIR/TPS.py -n $NUM_PER_JOB -N $((NUM_PER_JOB*SLURM_ARRAY_TASK_ID)) -f $FILE_NAME -T 1 
