# HostG
Graph convolutional neural network for host prediction


HostG is a GCN based model, which can learn the species masking feature via deep learning classifier, for new phage host prediction. To use HostG, you only need to input your contigs to the program.


# Required Dependencies
* Python 3.x
* Numpy
* Pytorch
* Networkx
* Pandas
* [Diamond](https://github.com/bbuchfink/diamond)
* BLAST
* MCL
* [Prodigal](https://github.com/hyattpd/Prodigal)

All these packages can be installed using Anaconda.

If you want to use the gpu to accelerate the program:
* cuda 10.1 
* Pytorch-gpu

## An easiler way to install
We recommend you to install all the package with [Anaconda](https://anaconda.org/)

After cloning this respository, you can use anaconda to install the **HostG.yaml**. This will install all packages you need with gpu mode (make sure you have installed cuda on your system to use the gpu version. Othervise, it will run with cpu version). The command is: `conda env create -f HostG.yaml`

# Prepare the dataset
```
cd dataset
bzip2 -d protein.fasta.bz2
bzip2 -d nucl.fasta.bz2
cd ..
```

# Usage (example)
Here we present an example to show how to run HostG. We support a file named "contigs.fa" in the Github folder and it contain contigs simulated from E. coli phage. The only command that you need to run is `python run_Speed_up.py --contigs test_contigs.fa --len 8000 --t [confidence(SoftMax value)]`. 

There are three parameters for the program: 1. `--contigs` is the path of your contigs file. 2. `--len` is the minimum length of the viral contigs that will be processed. The default length is 8000bp. 3. `--t `. This will output predictions only the confidence (SoftMax value) is larger than the threshold (from 0 to 1, default 0).  Both `--len` and '--t' are cutoffs decided by users. As shown in the paper, if you prefer higher accuracy, you can specify bigger values for len and t so that only the contigs longer than len and have confidence > t will be output.

The output file is **final_prediction.csv**. There are several column in this csv file: "contig_name, median_file_name, [taxa]".

# Database Extension
Since the limitation of storage on GitHub, we only upload part of bacteria genomes in *bacteria* folder. Information about these bacteria can be found in **database.csv** in the *dataset* folder. Thus, the model can only predict the reported labels for phages. 


All the bacteria genomes can be downloaded using the `datasets` binary. The guideline is on [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/quickstarts/command-line-tools/). Once you installed it on your systems, you can directly download the genomes by the given **All prokaryote.csv** file in *dataset* folder. remember to add the labels (Given in **All prokaryote.csv**) to **label.csv** to ensure the model can gain new labels from you downloaded dataset.

If you perfer website-based method to download all the genomes, you can try to upload all the accession in **All prokaryote.csv** to [NCBI batchentrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez)

**Note** If you do not know the taxonomy, you use can use this [link](https://github.com/KennthShang/PYlogeny) to generate a csv file by giving the accesion of the genomes. But be careful that you may not find all the taxonomy by using the script. Thus, you need to remove the unkown genomes to prevent some potential error.


We also describe and evaluate the extension ability in the paper [See references section]. If you want to extend the database, please follow the same format and upload your genomes into to *bacteria* folder. Also you need to update the **label.csv** in the *dataset* folder according to your new genomes. 


# Training on your new dataset
If you have new phage-host interactions and you want to train HostG on them. There are several steps you need to follow:
1. Add your phages into both **nucl.fasta** and **protein.fasta** file in *dataset* folder.
2. You need to generate the `gene_to_genome.csv` file about your new phages. The format is shown as below. You need to specify which phage (*contig_id*) the protein (*protein_id*) belongs to. If you do not know the keywords of the protein, you can use "hypothetical protein" instead. This will **NOT** influent the results of HostG. Then you should add all the information this new `gene_to_genomes.csv` into our `dataset_gene_to_genomes.csv`.
![image](https://user-images.githubusercontent.com/22445402/131283164-6f67621c-3d40-4647-a0d1-4d2fd820ab15.png)
4. Add your host genomes in *bacteria*.
5. Add both labels for phages and hosts into **label.csv**

Then, we you predict for other contigs, HostG will integrate these new informations for prediction.




# Format of the file

1. The format of the input file should be a fasta file. Since the script will pass contigs with non-ACGT characters, which means those non-ACGT contigs will be remained unpredict.
2. The format of the output file is a csv file which contain the prediction of each phage. *contig_name* is the accession from the input. *idx* is the temporary name used in HostG 

![image](https://user-images.githubusercontent.com/22445402/131282066-e8c9743f-2b56-431e-84d3-cecca893aea1.png)

** Note ** Because HostG run on each taxonomy iteratively, there might exist some consistency problem for the prediction at each taxa level. We supply a script to make sure the prediction is consistent to NCBI taxonomy tree. Thus, there will be two output files: final_prediction_consistency.csv and final_prediction.csv (raw output). This is because some of the user might want to use the taxonomy from GTDB. User can replace the `label.csv` in dataset folder with GTDB taxonomies and use the raw output as the prediction.


# References
The pre-print version is avaliable on https://arxiv.org/abs/2105.13570


## Contact
If you have any questions, please email us: jyshang2-c@my.cityu.edu.hk


## Notes
if the program output an error (which is caused by your machine):
`Error: mkl-service + Intel(R) MKL: MKL_THREADING_LAYER=INTEL is incompatible with libgomp.so.1 library.`
You can type in the command `export MKL_SERVICE_FORCE_INTEL=1` before runing *run_Speed_up.py*
