This directory contains code to setup and run GRNBoost2 (via arboreto) on each
of the microglia datasets. I must emphasize that this process was rather janky.
There are known issues installing/running arboreto that required me to roll
back installations: see ```setup_arboerto.sh```

I had two further unresolved issues: some datasets would reproducibly stall
in their iterations, but not at the same iteration. So I had to repeatedly 
re-run these datasets to achieve the 100 iterations.

Next, I could not get dask slurm running, and so used the local cluster 
implementation and ran in slurm. It was unclear if this was an Alex problem or 
part of the aforementioned install issues. If you're going to continueusing 
GRNBoost2 it would be worth exploring this further, given how intense this 
process is.

Final note: the python scripts to run GRNBoost2 are loading the entire list of
microglia count matrices/metadata via rds2py and then looping over each (i.e.,
the process is serialized). It would be much smarter to have the datasets saved
individually as matrices to skip the bulky rds loading, and to better allow
splitting of tasks.

I ran everything in slurm (tmux + srun)

```
srun --ntasks=1 --cpus-per-task=40 --mem=100G --time=7-00:00:00 --pty bash -i
conda activate agg
python3 grnboost2_mcg_hg.py
```
