# MicroDynamics
A bioinformatics pipeline for studying microbial community landscape along chronic human disease progression path

## Setup
The setup process for the pipeline requires the following steps:
### Download pipeline
```bash
git clone https://github.com/liluacrobat/MicroDynamics.git
```

### Requirement
The following software is required
* MATLAB

Pleae add all the subforders of the MicroDynamics working directory to the path list of Matlab

## Input
There are two input files for the MicroDynamics pipeline: an OTU table and a meta file of disease behavior. In the following example, a gut microbiom data set of Crohn's disease is used. 
#### OTU table of 16S rRNA sequences
```
OTU_ID Sample_1 Sample_2 Sample_3
OTU1  15  30  12
OTU2  152 116 130
OTU3  27  208 74
```
#### Meta file of disease behavior
```
Sample_ID Disease_behavior
Sample_1  HC
Sample_2  Penetrating (B3)
Sample_3  Inflammatory (B1)
```

## Usage
### 1. Preprocesssing
Load the OTU table and meta data. Exclude samples without enough sequencing depth (default:10,000). 
```
[Table_otu, Table_clinic] = script_data_processing(filen_otu, file_meta, params)
```
#### Optional arguments  
```
params       
    -- min_count
         The minimum sequence counts used as a threshold to retain a microbiome sample [default: 10,000]
    -- pseudo_count
         A small number added to the relative abundance before 10-base log
         transformation [default: 10^-6]
    -- last_tax
         Flag indicating whether the last column of the OTU table is taxonomy
         or not. If the last column of the table is the taxonomy, you
         specify 1 [default: 0]
    -- col_label
         The column of the clinical information used for feature
         selection
    -- mapping
         Mapping from class categories to numerical labels
```

### 2. Identifying disease associated OTUs
Perform feature selection within the LOGO framework [1].
```
Feature_Table = script_feature_LOGO(Table_otu, Table_clinic, params)
```
#### Optional arguments 
```
params  
    -- sigma
         Kernel width (k nearest neighbor) [default: 10]
    -- lambda
         Regularization parameter [default: 10^-4]
    -- threshold
         Threshold of feature weight
```
The regularization parameter lambda can be optimized using 10-fold cross-validation.
```
[opt_lambda, ACC_LOGO] = script_param_LOGO(Table_otu, Table_clinic, params)
```
#### Optional arguments
```
params       
    -- lam_ls
         Range of the regularization parameter lambda [default: 10^-5~100]
    -- sigma
         Kernel width (k nearest neighbor) [default: 10]
    -- folds
         Number of folds for cross-validation [default: 10]
```

### 3. Random sampling-based consensus clustering
Perform random sampling-based consensus clustering to group samples with similar microbial community composition.
```
cidx = script_consensus_clustering(Feature_Table, params)
```
#### Optional arguments
```
params       
     -- med          : Method used for clustering 
                         options 
                            km: k-means
                            dmm: dirichlet multinomial mixtures
                            cst: community state typing with partition around medoids
     -- alpha        : Parameter for Dirichlet process prior
     -- cluster_num  : Number of clusters (available when med=km)
     -- iters        : Number of iterations for consensus clustering [default: 1000]
```
When k-means or community state typing (CST) with partition around medoids is used for clustering, the number of clusters can be automaticly determined using
```
[cluster_num, eva_Gap] = script_gap_statistic(Feature_Table,med)
 ```
#### Optional argument med
Specify the method used for clustering
```
'km': k-means
'cst': community state typing with partition around medoids
```

### 4. Embedded structuring learning
Learn a principal tree using the DDRTree method [2].
```
params    
    -- sigma
         Bandwidth parameter
    -- lambda
         Regularization parameter for inverse graph embedding
    -- col_label
         The column of the clinical information used for feature
         selection
```
The parameters employed by DDRTree is tuned by using the elbow method.
```
[var_opt, CurveLength, Error_mse] = script_elbow_DDRTree(Feature_Table, params)
```
#### Optional arguments
```
params      
    -- sigma
         Bandwidth parameter
    -- lambda
         Regularization parameter for inverse graph embedding
    -- f_sig
         Optimize the bandwidth [0: optimize lambda; 1: optimize sigma]
    -- var_ls
         Range of the variable to be tuned
```

### 5. Visualization
Visualize the principal tree learned from data.
```
script_visualization(PrincipalTree, annotations, levels, params)
```
#### Optional arguments
```
params        
    -- FaceColor
         Facecolor for annotations
    -- order
         Order used to sort labels
    -- post_flag
         Flag of post processing [default: 0]
           0: principal tree before post-processing
           1: principal tree after post-processing
    -- prog_flag
         Plot samples of the extracted progression paths [default: 0]
```

## Example
We provide two tab-delimited files of a human gut microbiome data set [3] in the example directory: a OTU table file CD_16S_OTU.txt, and a meta data file CD_meta_clinic.txt. 

Due to the storage limitation of GitHub, please download the CD_16S_OTU.tsv file from https://drive.google.com/file/d/1OhjjGS5Kw5G4ImOlzy8G8HluHM1aNjUN/view?usp=sharing and put the file in the 'example/data' directory. 

We can learn a microbial progression model by running Demo_MicroDynamics.m. The precalculated results are stored in the 'example/precalculated/' directory. 

## Reference
[1] Yijun Sun, Sinisa Todorovic, and Steve Goodison. "Local-learning-based feature selection for high-dimensional data analysis." *IEEE Transactions on Pattern Analysis and Machine Intelligence*, 32(9), 1610-1626, 2010.  
[2] Qi Mao, Li Wang, Steve Goodison, and Yijun Sun. "Dimensionality reduction via graph structure learning." In *Proceedings of the 21th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*, pp. 765-774, 2015.  
[3] Jonas Halfvarson, Colin J. Brislawn, Regina Lamendella, Yoshiki Vázquez-Baeza, William A. Walters, Lisa M. Bramer, Mauro D'amato et al. "Dynamics of the human gut microbiome in inflammatory bowel disease." *Nature Microbiology*, 2(5), 1-7, 2017.  
