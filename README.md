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

All these packages can be installed using Anaconda.

If you want to use the gpu to accelerate the program:
* cuda 10.1 
* Pytorch-gpu

## An easiler way to install
We recommend you to install all the package with [Anaconda](https://anaconda.org/)

After cloning this respository, you can use anaconda to install the **HostG.yaml**. This will install all packages you need with gpu mode (make sure you have installed cuda on your system).

# Usage (example)
Here we present an example to show how to run PhaGCN. We support a file named "contigs.fa" in the Github folder and it contain contigs simulated from E. coli phage. The only command that you need to run is `python run_Speed_up.py --contigs test_contigs.fa --len 8000 --t [confidence]`. 

There are three parameters for the program: 1. `--contigs` is the path of your contigs file. 2. `--len` is the length of the contigs you want to predict. As shown in our paper, with the length of contigs increases, the recall and precision increase. The default length is 8000bp. `--t [confidence]`. This will output predictions only the confidence is larger than the threshold. We recommend you to choose a proper length according to your needs.

The output file is **final_prediction.csv**. There are several column in this csv file: "contig_name, median_file_name, [taxa]".

# Database Extension
Since the limitation of storage on GitHub, we only upload part of bacteria genomes in *bacteria* folder. Information about these bacteria can be found in **bacteria.csv** in the *dataset* folder.

If you want to extend the database, please follow the same format and upload your genomes into to *bacteria* folder. Also you need to update the **label.csv** in the *dataset* folder according to your new genomes. 

# Notice
If you want to use HostG, you need to take care of:
1. The script will pass contigs with non-ACGT characters, which means those non-ACGT contigs will be remained unpredict.


# References
The paper was submitted to *Bioinformatics* and the pre-print version is avaliable on https://arxiv.org/abs/2105.13570


## Contact
If you have any questions, please email us: jyshang2-c@my.cityu.edu.hk
