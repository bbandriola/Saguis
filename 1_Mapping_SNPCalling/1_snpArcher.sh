# 1 -> download ref genome
wget --recursive --no-host-directories --cut-dirs=6 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/400/415/GCA_013400415.1_ASM1340041v1/ -P ./

# 2 -> run snparcher in the folder project
conda activate snparcher
# dryrun
snakemake --executor slurm -p --snakefile /home/bruna.andriola/Software/snpArcher/workflow/Snakefile --directory ./ --cores 48 --use-conda --workflow-profile ./profiles/genoma/ --dryrun
# actual run
snakemake -p -s /snpArcher/workflow/Snakefile -d ./ --cores 10 --use-conda --workflow-profile project/profiles/genoma
