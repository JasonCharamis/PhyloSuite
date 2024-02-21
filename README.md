**Scalable RAxML-NG-based phylogenetic analysis using Snakemake**

This is a Snakemake workflow for running scalable maximum likelihood (ML) phylogenetic analysis using RAxML-NG and associated tools (Pythia, ModelTest-NG). This workflow is considerably slower than the IQ-TREE-based one, but is much better in terms of accuracy, especially in difficult-to-analyze datasets.

This workflow performs all steps sequentially, from MSA to model selection and phylogeny inference. 
To include custom options, edit the config/config.yaml file.

Usage:
snakemake --cores 20 --snakefile Snakefile

All dependencies are downloaded and installed through conda.

Mafft
https://github.com/GSLBiotech/mafft

Trimal
https://github.com/inab/trimal

Pythia
https://github.com/tschuelia/PyPythia

ModelTest-NG
https://github.com/ddarriba/modeltest

RAxML-NG
https://github.com/amkozlov/raxml-ng

ETE3
http://etetoolkit.org/

**To run this fully automated pipeline using Docker containers use the following commands (Docker requires sudo rights)**:
1. git clone https://github.com/JasonCharamis/Snakemake-workflow-for-RAxML-NG-based-phylogenetic-analysis.git
2. cd Snakemake-workflow-for-RAxML-NG-based-phylogenetic-analysis/workflow/ && **sudo** docker build -t automated_phylogenetic_analysis:latest . 
3. Travel to the directory where you have the desired fasta files (better only those)
4. **sudo** docker run -it -v $(pwd):/workflow -w /workflow automated_phylogenetic_analysis:latest snakemake --cores 20 --snakefile Snakemake-workflow-for-RAxML-NG-based-phylogenetic-analysis/workflow/Snakefile --use-conda --conda-frontend mamba

Of course, to customize the run edit the config/config.yaml file. 

That's it! The pipeline will run automatically.



