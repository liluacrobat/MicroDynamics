# MicroDynamics
A pipeline for delineating community dynamics of human microbiome associated with disease progression

## Setup
The setup process for the proposed pipeline requires the following steps:
### Download the pipeline
```bash
git clone https://github.com/liluacrobat/MicroDynamics.git
```

### Requirement
The following software is required
* MATLAB

## Input
There are two input files for the MicroDynamics pipeline including an OTU table and a meta file of disease behavior. The following example is a human gut microbiom data set of Crohn's disease. 
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
         Number of observation (sequence) count to apply as the minimum
         total observation count of a sample for that sample to be retained.
         If you want to include samples with sequencing depth higher than
         or equal to 10,000, you specify 10,000 [default: 10,000]
    -- pseudo_count
         A small number added to the relative abundance before 10-base log
         transformation [default: 10^-6]
    -- last_tax
         Flag of whether the last column of the OTU table is taxonomy
         or not. If the last column of the table is the taxonomy, you
         specify 1 [default: 0]
    -- col_label
         The Column of the clinical information used for feature
         selection
    -- mapping
         Mapping from class categories to numerical labels
```

### 2. Identifying disease associated OTUs
Feature selection within LOGO [1] framework.
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

### 3. Random sampling based consensus clustering
Perform random sampling based consensus clustering to group samples with similar microbial community composition.
```
cidx = script_consensus_clustering(Feature_Table, params)
```
#### Optional arguments
```
params        
    -- cluster_num  
         Number of clusters
    -- iters        
         Number of iterations
```
The number of clusters can be estimated based on gap statistics.
```
[cluster_num, eva_Gap] = script_kmeans_gap(Feature_Table)
```

### 4. Embedded structuring learning
Learn a principal tree using DDRTree [2] method.
```
params    
    -- sigma
         Bandwidth parameter
    -- lambda
         Regularization parameter for inverse graph embedding
    -- col_label
         The Column of the clinical information used for feature
         selection
```
The parameters employed by DDRTree can be tuned using elbow method.
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
           0 : principal tree before post-processing
           1 : principal tree after post-processing
    -- prog_flag
         Plot samples of the extracted progression paths [default: 0]
```

## Example
We provide two tab delimited files of a human gut microbiome data set [3] in the example directory: a OTU table file CD_16S_OTU.txt, and a meta data file CD_meta_clinic.txt. 

Due to the storage limitation of GitHub, please download the CD_16S_OTU.tsv file from https://drive.google.com/file/d/1OhjjGS5Kw5G4ImOlzy8G8HluHM1aNjUN/view?usp=sharing and put the file under the 'example/data' directory. 

Then, we can learn a microbial progression model by running Demo_MicroDynamics.m. The precalculated results are stored in the 'example/precalculated/' directory. 

## Reference
[1] Sun, Yijun, Sinisa Todorovic, and Steve Goodison. "Local-learning-based feature selection for high-dimensional data analysis." *IEEE Transactions on Pattern Analysis and Machine Intelligence* 32, no. 9 (2009): 1610-1626.  
[2] Mao, Qi, Li Wang, Steve Goodison, and Yijun Sun. "Dimensionality reduction via graph structure learning." In *Proceedings of the 21th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*, pp. 765-774. 2015.  
[3] Halfvarson, Jonas, Colin J. Brislawn, Regina Lamendella, Yoshiki Vázquez-Baeza, William A. Walters, Lisa M. Bramer, Mauro D'amato et al. "Dynamics of the human gut microbiome in inflammatory bowel disease." *Nature Microbiology* 2, no. 5 (2017): 1-7.  
