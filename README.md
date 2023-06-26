# README for GENERA

Software for the paper [doi]

## REQUIREMENTS
python 3 <br>
pytorch <br>
rdkit

## GENERAL SETUP

Some environment variables need to be exported: <br>
export GACONF=path/to/GACONF/download <br>
export EVOSRC=path/to/EVOSRC/download <br>
export SCRIPTS=path/to/SCRIPTS/download <br>

## WORKDIR SETUP
GENERA runs in a working directory where all temporary and output files will be saved. The working directory should start with a file called done_so_far that contains the starting population.
GENERA will update this file at each iteration with new entries. At the end of the execution, the done_so_far file will contain the generated molecules and with all the required objectives.


## THE DONE_SO_FAR FILE
The done_so_far file is the starting input file and will be the final output file. This file will contain each entry as a line. 
For each entry, there should be:
\<SMILES\> = <objective 1> <objective 2> … <objective n>

Where: <br>
\<SMILES\> is the SMILES string of the generated molecules devoid of stereochemical notation. <br>
\<objective n\> is the nth objective function for the Pareto multiobjective optimization.  The content of these lines depends on the chosen mode. See the paper for examples of objective functions.

 
## START COMMAND
The software should be called using the command:
$GACONF/pilot_local.csh workdir=<path_to_working_directory> recdir=<argument_for_specific_mode> minpop=1 mode=<mode_name> cont=yes nc=<number_of_CPUs> noprog=<stop_criterion_1> maxconfigs=<stop_criterion_2>

Where:
<ul>
<li> workdir should be set to the path (either absolute or relative) to the chosen working directory </li>
<li> mode is the name for the chosen mode. The software will look for a <mode_name>.csh file in the $GACONF folder. We provide the denovoDockpareto mode that reproduces Experiment 1 from the paper. </li>

<li> nc is the number of CPUs to be used. By default, nc equals the number of CPUs available on the current machine. </li>

<li>noprog is a stopping criterion on the GA iterations. The algorithm will stop running when fewer new molecules than noprog are produced.</li>

<li>maxconfings is a stopping criterion on GA iterations. The algorithm will stop when a total number of maxconfig molecules is found in the done_so_far file.</li>

</ul>
The user may rely on the currently provided $GACONF/local_denovoDockpareto.csh to change the score criteria defining the Pareto front. The resulting script may be saved under a new name local_<new_mode_name>.csh, which will run when calling pilot_local.csh with mode=new_mode_name as an argument. However, for consistency reasons, before the execution, two parameter files <new_mode_name>.rng and <new_mode_name>.pars need to be formally created in $GACONF. These are independent of the scoring details in <new_mode_name>.csh and therefore are best created by copying or linking to provided denovoDockpareto.{pars,rng} files.

## EXIT CONDITIONS

-noprog. When the number of novel valid molecules is below a given number, the pilot script will stop execution.

-reaching maxconfig. When the done_so_far file contains a number of rows greater or equal to maxconfigs, the pilot script will stop.
 
-Manual stopping. To stop the execution correctly and cleanly, the user should create a file called stop_now in the working directory. When the stop_now file is found in the working directory, the pilot will end the current iteration and stop.
NB: killing the pilot script execution is not guaranteed to stop all the running subprocesses, so please use the "stop_now file" method 


## EXAMPLE
We provide the necessary files and data to replicate the Experiment 1 from the paper. Remember that this mode has specific dependencies, namely the S4MPLE software (available at https://infochim.u-strasbg.fr/) and the ChemAxon package.



## CONTACTS
dhorvath@unistra.fr <br>
giuseppe.mangiatordi@ic.cnr.it <br>
giuseppe.lamanna@ic.cnr.it <br>
