#Running the test for Gaussian
###1. Change the PATH_TO_SAVE_CHK_FILE in H2O.gjf.
```text
vi H2O.gjf
```
###2. Run the sub.gaussian script to submit the job to queue, e.g.,
```text
sbatch sub.gaussian
```
Note that the format of the sub.gaussian script and running command, which depend on the job scheduling system 
(e.g., Slurm) on Linux clusters, maybe different from the given example.