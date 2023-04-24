# SRD129_Bordetella_bronchiseptica

**Last updated 24Jan2022**

This is the repository for scripts and data files pertaining to the research paper "Shifts in the swine nasal microbiota following Bordetella bronchiseptica challenge in a longitudinal study."

Fastq files are located in Bioproject PRJNA525911.

R version 4.1.2 (2020-06-22) and accompanying packages were used to run the scripts for this research.

### **Nomenclature**
BB = *Bordetella bronchiseptica*

## **Table of contents**
| Chapter | Description |
| -- | -- |
| [data](https://github.com/nielsend/SRD129_Bbronchiseptica/tree/master/data) | Includes data files needed to carry out R analysis |
| [scripts](https://github.com/nielsend/SRD129_Bbronchiseptica/tree/master/scripts) | Text file of commands used for mothur, R scripts for 16S rRNA analysis|

## **Scripts description and the order to run them**
| Order | Description |
| -- | -- |
| 1a | Process 16S sequence data and generate output for R scripts. |
| 1b | Create OTU table from mothur output to use for creating phyloseq objects. |
| 2  | Generate phyloseq object to use for 3_alpha_beta_diversity_BB.R. Run [adonis](https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/adonis) function with distance matrices to assess how variation is attributed to different experimental treatments or uncontrolled covariates. |
| 3  | Run alpha (Shannon, Inverse Simpson) and beta diversity (generating NMDS, pairwise comparisons) analyses, data visualization. |
| 4  | Identify differentially abundant bacterial taxa (genus level) between groups within each day, data visualization. Also, analyze Bordetella counts |
| 5  | Generate a list of percent total genera found in each treatment group per day for each tissue, data visualization.
| 6  | Plot the F-statistic from PERMANOVA pairwise comparisons of Control and BB groups over time. This displays the magnitude of change in the nasal bacterial community structure of the BB group relative to Control. Also includes scripts for data visualization. |

