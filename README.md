**********************************
# Copyright and License Information
**********************************
[![](https://img.shields.io/badge/Python-3.5.2-brightgreen)](https://www.python.org/downloads/release/python-352/)
![](https://img.shields.io/badge/matplotlib-3.3.3-blue)
[![](https://img.shields.io/badge/R-4.1.0-orange)](https://cloud.r-project.org/src/base/R-4/R-4.1.0.tar.gz)
![](https://img.shields.io/badge/ggplot2-3.3.5-red)

Copyright (C) 2021
Northwest A&F University,
Mingzhi Liao, Xianyong Lan
  
Authors: Hongfei Liu

Contact: lhf563@nwafu.edu.cn

This project is avaliable for the comparison of different circRNA software packages predicted from short-read illumina sequencing datasets.
All of data and source code are free and you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

******************
# Table of contents
******************
* [Background](#background)
* [Usage](#usage)
  * [config](#config)
  * [results](#results)
  * [scripts](#scripts)
    * [circRNA detection](#circrna-detection)
    * [Downstream analysis](#downstream-analysis)
  * [shell\_scripts](#shell_scripts)
* [Citation](#citation)


# Background
  Circular RNA is generally formed by the "back-splicing" process between the upstream splice acceptor and the downstream donor in/not in the regulation of the corresponding RNA-binding proteins or cis-elements. Therefore, more and more software packages that have been developed based on the identification of the back-spliced junction (BSJ) reads. However, recent studies have developed two software tools that can detect circRNA candidates by constructing k-mer table or/and de bruijn graph rather than reads mapping.
  
  Here, we compared the precision, sensitivity and detection efficiency between software tools based on different algorithms. Eleven representative detection tools with two types of algorithm were selected for the overall pipeline of analysis of RNA-seq datasets with/without RNase R treatment in two cell lines. Precision, sensitivity, AUC, F1 score and detection efficiency metrics were assessed to compare prediction tools. Meanwhile, the sensitivity and distribution of highly expressed circRNAs before and after RNase R treatment were also revealed by their enrichment, unaffected and depleted candidate frequencies. Eventually, we found that compared to the k-mer based tools, [CIRI2][CIRI2] and [KNIFE][KNIFE] with reads mapping based had relatively superior and balanced detection performance regardless of the cell line or RNase R (-/+) datasets. In summary, the novel k-mer based software show dominant performance on sensitivity and computational efficiency in circRNA discovery. This study may provide new insights into development and application in circRNA detection tools.

  The real datasets are available from the NCBI Gene Expression Omnibus (GEO) database (BioProject: [PRJNA231724](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA231724), GEO: [GSE53327](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53327)) followed by different RNase R treatment in two cell lines.
# Usage
- bash.sh: one-step shell script of general circRNA-seq analysis pipeline for all software tools. As for the specific and detailed usage of this shell script and/or other scripts, please read the description in these scripts or try to run them.
- SRR_list.txt: the file contains accesson_id of all fastq files
## config
The required config files for specific software
- CDBG_config.ini: config file of [CircDBG](https://github.com/lxwgcool/CircDBG), which contains reference fasta file, annotation file, reads1/2, and other required or optional parameters. Of which, Reference, GTF, Reads1/2, and options in Parameter section is important and required.
- CM_config.ini: config file of [CircMarker](https://github.com/lxwgcool/CircMarker) like [CircDBG](https://github.com/lxwgcool/CircDBG)
- paired_sample: config file of [segemehl][segemehl], which maily contains the absolute or relative path of results of BSJ reads given by STAR alignment
## results  
The raw, filtered or/and annotated predicted candidates by each software under different dataset (can be downloaded from figshare: https://doi.org/10.6084/m9.figshare.19090640.v1).
- circ_candidates.bed ([CIRCexplorer2][CIRCexplorer2]) or other different files exclude followings: raw identified circRNA bed or other format files for each software package
- circ_candidates_convert.bed: coonverted circRNA bed format file (genome coordinate converted to uniform 0-based format)
- circularRNA_known.txt: annotated circRNA information file generated from [CIRCexplorer2][CIRCexplorer2] annotate moudle
- circRNA_known_annotated.txt: it contains annotated circRNA list (circularRNA_known.txt) merged by the known circRNA retrieved from [circBase][circBase] and [circAtlas][circAtlas] database
- BackgroundPlotdf.csv: The dataset of plotting upset plot of circRNA candidates on background datasets with different depths
- BackgroundRaw.csv: Raw datasets of expression matrix of circRNA predicted by all software on background datasets with different depths
- performance.csv: recorded prediction performance indices of all software tools
## scripts
### circRNA detection
- get_reads_length.py: get the average read length of all datasets
- summarize.py: summarize the raw results (*_filt.sngl.bed) of [segemehl][segemehl] into summary bed file (*.sum.bed)
### Downstream analysis
The following scripts are used for data clean, analysis and visualization with python or R
- conversion.r: convert raw genome coordinate (0-based or 1-based) to 0-based.
- annotation.r: merge annotated circRNA list (circularRNA_known.txt) by the known circRNA retrieved from [circBase][circBase] and [circAtlas][circAtlas] database
- pr-RnaseR.py: data clean and downstream analysis python script for positive, mixed, and real datasets
- background_analysis.py: data clean and downstream analysis python script for background datasets
## shell_scripts  
The general prediction pipeline for each software from short-read RNA-seq
# Citation
If you find this code useful in your research, please cite:  
  
  Liu H, Akhatayeva Z, Pan C, Liao M, Lan X. Comprehensive comparison of two types of algorithm for circRNA detection from short-read RNA-Seq. Bioinformatics. 2022 Apr 28:btac302. doi: [10.1093/bioinformatics/btac302](https://doi.org/10.1093/bioinformatics/btac302).

[CIRCexplorer2]: https://circexplorer2.readthedocs.io/en/latest/
[CIRI2]: https://sourceforge.net/projects/ciri/files/CIRI2/
[KNIFE]: https://github.com/blawney/knife_circ_rna
[segemehl]: https://www.bioinf.uni-leipzig.de/Software/segemehl/
[circBase]: http://www.circbase.org/
[circAtlas]: http://circatlas.biols.ac.cn/