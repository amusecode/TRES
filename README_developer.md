# TRES as developer
This is a short additional document that will help guide you through the installation and compilation process of TRES as a developer. For a more complete guide on how to use TRES, please have a look at the [README.md](https://github.com/amusecode/TRES/blob/main/README.md) file.

For the installation, it is important to know that TRES combines two different codes: SeBa, which is a stellar evolution code, and a dynamical code. The dynamical code is already integrated into TRES and thus doesn't need to be installed separately. SeBa however, needs to be installed through the AMUSE framework. We will explain how to do this shortly, but just make sure to satisfy the [AMUSE pre-requisites](https://amuse.readthedocs.io/en/latest/install/howto-install-AMUSE.html) and have installed the necessary python modules beforehand.

In case the you have access to a computer cluster, [Run TRES on cluster](#Run-TRES-on-cluster) gives some additional information on running TRES with slurm.
 
### Installation TRES
If you are interested in applying changes to the stellar evolution code SeBa, the code should be locally installed. To do this, you will have to install AMUSE directly from the source code:

```
git clone https://github.com/amusecode/amuse.git
```

It is important to install the amuse packages using the developer mode. We can do this by typing the following command in the root of the cloned repository:

```
pip install -e .
```

Then, SeBa can be build with:

```
make seba.code
```

If instead you wish to build all the community codes that amuse provides, we can write:

```
python setup.py develop_build
```

Now, you should be able to navigate to a directory called "amuse/src/amuse/community/seba", which contains the complete evolution code. To compile the code, write (exluding the comments preceded by a "#"):

```
make clean
make               # Here we compile all the C files
```

It is very important to know that anytime the code is compiled, the source code will be downloaded again, meaning any changes to SeBa will be overwritten. To prevent this from happening, we need to comment out a line in the Makefile. In line 47, put a "#" before "download src/SeBa".

TRES can simply be installed by cloning the github repository in the terminal:

```
git clone https://github.com/amusecode/TRES.git
```

## TRES development team
Members of the TRES development team are recommended to work on their own fork on github. To set this up:

0) After installing AMUSE and downloading TRES
1) First create a fork. Can be done easily on the github webinterface and creates the repository "username/TRES", where the username refers to your git account. 

2) Now we have to set up links to the official TRES repository and the forked one. Clone the forked TRES repository, which will be know as ‘origin’
```
git clone git@github.com:(username)/TRES.git
```
And add a link to the official TRES repository, which will be known as upstream
```
git remote add upstream git@github.com:amusecode/TRES.git
```
Also to list the current configured remote repositories for your fork:
```
git remote -v
```


Example of workflow:
1) Good practise to work on a new branch
```
git checkout -b (branch)
```
2) pull all updates from ‘origin’ to current branch
```
git pull origin main
```
3) pull all updates from upstream to current branch
```
git pull upstream main
```
4) make the desired changes to TRES
```
git add / git commit
git checkout main / git merge (branch)/ git branch -d (branch)
```
5) push changes to your own fork and follow that branch
```
git push origin main
```
6) If so desired, push changes to the official TRES repository
```
git push --set-upstream upstream (branch) 
```


Manage who has access to your fork via github:  "settings -> manage access -> invite a collaborator"

### Development tips & tricks
1) TRES makes use of object oriented programming. The structure of the object is as follows:
```
- self contains simulations parameters such as the initial and final time
- self.triple contains general parameters such as the time, the ID number 
- self.triple can have 1 or 2 children. Each child can have 1 or 2 children. 
- a child is either a star with parameters such as mass or radius, or a binary system with parameters such as semimajoraxis and eccentricity.
- By default TRES starts the simulations with the structure below. 


           Default triple 
==============================================
|                  self                      |                 
|                   |                        | 
|               self.triple                  |                 
|                |      |                    | 
|                |   self.triple.child1      |         
|         self.triple.child2                 | 
|              |      |                      | 
|              |  self.triple.child2.child1  |  
| self.triple.child2.child2                  | 
==============================================

where self.triple.child1 is the tertiary star, and self.triple.child2.child1 & self.triple.child2.child2 the inner two stars. self.triple.child 2 itself represents the inner binary, while self.triple the outer binary. 
If the triple experiences a merger, and reduces to a binary system, the structure becomes: 


               Binary
==============================================
|                self                        |              
|                  |                         | 
|              self.triple                   |                
|               |      |                     | 
|               |   self.triple.child1       |        
|        self.triple.child2                  | 
==============================================

where self.triple.child1 and self.triple.child2 are two stars and self.triple represents the binary. 

```




2) To receive more output, there are a number of 'REPORT' statements (REPORT_DEBUG, REPORT_DT, REPORT_SN_EVOLUTION, REPORT_TRIPLE_EVOLUTION) that can be set to True on the top of TRES.py. For the most extensive and generic option choose REPORT_DEBUG. This option will also create pdfs a txt file with output at every global (TRES) timestep, as well as create pdf of the time evolution of many parameters. 

To receive more output from SeBa, do the following: in setup_stellar_code(), comment out self.stellar_code = SeBa() and uncomment self.stellar_code = SeBa(redirection='none'). 

To receive more output from the secular code, there are two options. 1) in setup_secular_code(), comment out self.stellar_code = SecularTriple() and uncomment self.stellar_code = SecularTriple(redirection='none'). 2) in setup_secular_code(), set self.secular_code.parameters.verbose to True. 

3) To reduce the global (TRES) timestep (aka get more timestamps in the TRES terminal output, see above) set the maximum_time_step on the top of TRES.py

4) To use the detailed gyration radius and apsidal motion constant from SeBa, set GET_GYRATION_RADIUS_FROM_STELLAR_CODE and/or GET_AMC_FROM_STELLAR_CODE to True on the top of TRES.py

5) To start back up your simulation of a given triple at a specific time, use the input parameter tinit. For example to start the simulation of the default triple at 2.5Myr: 

```
amuse TRES.py --initial_time 2.5  
```
If you use this option, the stellar masses should reflect the masses on the zero-age main-sequence, while the orbital parameters should reflect values at the specified time. 
For the moment this only works for pre-mass transfer systems. 

6) To change the prescription of the Roche lobe, you can set self.secular_code.parameters.roche_radius_specification in setup_secular_code to 0, 1, or 2: 
```
0  Roche radius based on Eggleton's formula using the pericenter distance (i.e. including eccentricity factor)  [default]
1  Roche radius based on Sepinsky. Function of not only eccentricity but also stellar spins
2  Roche radius based on Eggleton's classical formula. Watch out this is only valid for circularized & synchronised systems
```

### Run TRES on cluster
Running computationally expensive simulations on a computer cluster can save a lot of time. However, clusters work somewhat different than your personal computer. There are two points we'd like to draw your attention to.

For starters, make sure if the pre-required packages for AMUSE are already installed. The easiest way to do this is try installing AMUSE and check where eventual errors occur. Unfortunately, on a cluster you will most likely not have sudo rights, so you'll have to figure out a way to install the missing packages. 

Second, clusters work with slurm. To run a simulation you must create a bash file that can be sumbitted as a slurm job. Down below is an example bash script for the helios cluster with additional comments for clarification. If you are using a different cluster, change the script accordingly.

```
#!/bin/bash

# first, we need to ask for resources with sbatch. The documentation can be found here: https://slurm.schedmd.com/sbatch.html
# this example uses the array functionality, which will submit an array of similar jobs that can be run simultaneously

#SBATCH --nodes 1
#SBATCH --tasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --time 1:00:00
#SBATCH -w helios-cn004
#SBATCH --job-name=TRES
#SBATCH --output=logs/array_%A_%a.out
#SBATCH --array=0-9
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

function clean_up {
    echo "### Running Clean_up ###"
    # all files should be deleted from the compute-node before job end
    # if there are more than one output files in the directory after the completion of one of the simulations, only delete the file
    # corresponding to that simulation. Otherwise, remove the whole directory
    num_files=$(ls *.hdf | wc -l)
    echo "$num_files"
    if [ "$num_files" -eq "1" ]
    then
        rm -rf "/hddstore/$USER"
        echo "Removed all files"
    else
        #rm -rf "/hddstore/$USER/$SLURM_ARRAY_TASK_ID"
        rm "/hddstore/$USER/$FILE_NAME" 
        echo 'Removed'"$FILE_NAME"''
    fi
    echo "Finished"; date
#    # - exit the script
    exit
}

# call "clean_up" function when this script exits, it is run even if
# SLURM cancels the job
trap 'clean_up' EXIT

##### pipeline below #####

# enter virtual environment
# this environment contains the necessary python packages 
. /home/fkummer/Amuse-env/bin/activate

# openmpi was not installed on helios and needed to be imported through a module 
module purge
module load openmpi/3.1.6


mkdir -p /hddstore/$USER                              # create a directory on the node to store your data
export OUTPUT_FOLDER="/hddstore/$USER"                # this simply creates an alias for the file/directory
cd $OUTPUT_FOLDER

export FILE_NAME='TRES_'"$SLURM_ARRAY_TASK_ID"'.hdf'  # for each sbatch array, define your output filename
export TRES="/home/fkummer/TRES"
python $TRES/TPS.py -n 10 --M_max 100 --M_min 15 --Q_min 0.1 --A_distr 5 --q_min 0.1 --E_max 0.9 --E_distr 3 --e_max 0.9 --e_distr 3 -z 0.0001 -f $FILE_NAME

export FOLDER_NAME="test_TRES/cpu_check"
mkdir -p /zfs/helios/filer0/$USER/$FOLDER_NAME        # create a directory where the data should be stored (for helios this is the /zfs/helios/filer0/ directory)                             
cp $FILE_NAME /zfs/helios/filer0/$USER/$FOLDER_NAME/  # copy the data from the node to the storage directory
```

