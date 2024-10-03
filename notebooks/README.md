# Process Mining results

Once you have run the mining across the entire PDB, run first the notebook called `process_results.ipynb` followed by the notebooks called `plot_results.ipynb` and `plot_sequence_conservation.ipynb` and `cross_reference_5_8.ipynb`:

- `process_results.ipynb`: reads the output for the mining and processes it but also adds other information, such as domain information. It also completes the degron, as sometimes the sequence in the structure files is cutoff. 

- `plot_results.ipynb`: plot the results for the figures used in the paper. 

- `plot_sequence_conservation.ipynb`: plot the sequence logo plots and alignments for the paper. 

- `cross_reference_5_8.ipynb`: Compare the overlaps of the 5 residue mining with that of the 8 residue mining. 

