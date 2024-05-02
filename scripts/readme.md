### Github repository for 

## Dataset jki_seq2_seq4 (16S rRNA gene amplicon sequencing) pipeline analyses
Performed in R v.4.1.3 [R Core Team](https://www.r-project.org)

#### Color code
- Site
  - ![#1d5b65](https://placehold.co/15x15/1d5b65/1d5b65.png) `#1d5b65` (Silty)
  - ![#a34828](https://placehold.co/15x15/a34828/a34828.png) `#a34828` (Sandy)
- Climate
  - ![#f8b195](https://placehold.co/15x15/f8b195/f8b195.png) `#f8b195` (Dry)
  - ![#a8e6c3](https://placehold.co/15x15/a8e6c3/a8e6c3.png) `#a8e6c3` (Wet)
- Season
  - ![#377881](https://placehold.co/15x15/377881/377881.png) `#377881` (Silty dry)
  - ![#6baea5](https://placehold.co/15x15/6baea5/6baea5.png) `#6baea5` (Silty wet)
  - ![#d1711f](https://placehold.co/15x15/d1711f/d1711f.png) `#d1711f` (Sandy dry)
  - ![#efaf7a](https://placehold.co/15x15/efaf7a/efaf7a.png) `#efaf7a` (Sandy wet)
- Microhabitats
  - ![#b05644](https://placehold.co/15x15/b05644/b05644.png) `#b05644` (BS)
  - ![#d9b967](https://placehold.co/15x15/d9b967/d9b967.png) `#d9b967` (RH)
  - ![#57896A](https://placehold.co/15x15/57896A/57896A.png) `#57896A` (RP)
- Rotations
  - ![#e7e5cc](https://placehold.co/15x15/e7e5cc/e7e5cc.png) `#e7e5cc` (W1)
  - ![#9cc184](https://placehold.co/15x15/9cc184/9cc184.png) `#9cc184` (W2)
  - ![#447243](https://placehold.co/15x15/447243/447243.png) `#447243` (W3)
  - ![#1e3d14](https://placehold.co/15x15/1e3d14/1e3d14.png) `#1e3d14` (WM)
- Phylum
  - ![#628dbd](https://placehold.co/15x15/1e3d14/1e3d14.png) `#1e3d14` (Actinobacteriota)
  - ![#a8e6c3](https://placehold.co/15x15/a8e6c3/a8e6c3.png) `#a8e6c3` (Bacteroidota)
  - ![#f8b195](https://placehold.co/15x15/f8b195/f8b195.png) `#f8b195` (Patescibacteriota)
  - ![#80c698](https://placehold.co/15x15/80c698/80c698.png) `#80c698` (Proteobacteria)
- Greenhouse
  - ![#549ccb](https://placehold.co/15x15/549ccb/549ccb.png) `#549ccb` (Disease)
  - ![#c6ba9c](https://placehold.co/15x15/c6ba9c/c6ba9c.png) `#c6ba9c` (Root biomass)
  - ![#80c698](https://placehold.co/15x15/80c698/80c698.png) `#80c698` (Shoot biomass)

### Packages
library("dada2")\
library("ShortRead")\
library("Biostrings")\
library("phyloseq")\
library("microbiome")\
library("vegan")\
library("DESeq2") \
library("dplyr")\
library("stringr")\
library("ggpubr")\
library("tidyr")\
library("ggplot2")\
library("readxl")\
library("RColorBrewer")
  
### 1. Dada2
This script uses the raw data obtained from the BioProject SRA.\
A R server is required. \
Database used: [SILVA 138 SSU](https://www.arb-silva.de/documentation/release-138/) 

### 2. Creating phyloseq object
This script creates phyloseq object based on these files:

- jki_seq1_metadata.csv
- jki_seq1_otu.xlsx
- jki_seq1_seqs.xlsx
- jki_seq1_taxa.xlsx

### 3. Rename NA
This script replace NA for the latest taxonomy found for an ASV.

### 4. Clean dataset
This script removes unwanted taxonomic groups from the dataset.
- Root NAs (Domain, Phylum)
- Eukaryotes (Domain)
- Cloroplasts (Order)
- Mitochondria (Family)

### 5. Data selection
This script selects group of samples to be analyzed separetely.

### 6. Rarefaction
This script performs rarefaction based on the minimum sequences.

### 7. Alpha diversity
This script calculates alpha diversity based on the rarefied data.

### 8. Ordination 
This scripts creates MDS plots and calculates PERMANOVA and ANOSIM.

### 9. DESeq2
This script performs differential abundance (DA) between two groups.
