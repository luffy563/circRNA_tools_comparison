*********************************
Copyright and License Information
**********************************
Copyright (C) 2021
Northwest A&F University,
Mingzhi Liao, Xianyong Lan
  
Authors: Hongfei Liu

This project is avaliable for the comparison of different circRNA software packages predicted from short-read illumina sequencing datasets.
All of data and source code are free and you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

******************
Table of contents
******************
Background
Circular RNA is generally formed by the “back-splicing” process between the upstream splice acceptor and the downstream donor in/not in the regulation of the corresponding RNA-binding proteins or cis-elements. Therefore, more and more software packages that have been developed based on the identification of the back-spliced junction reads. However, recent studies have developed two software tools that can detect circRNA candidates by constructing k-mer table or/and de bruijn graph rather than reads mapping.
Here, we compared the precision, sensitivity and detection efficiency between software tools based on different algorithms. Eleven representative detection tools with two types of algorithm were selected for the overall pipeline of analysis of RNA-seq datasets with/without RNase R treatment in two cell lines. Precision, sensitivity, AUC, F1 score and detection efficiency metrics were assessed to compare prediction tools. Meanwhile, the sensitivity and distribution of highly expressed circRNAs before and after RNase R treatment were also revealed by their enrichment, unaffected and depleted candidate frequencies. Eventually, we found that compared to the k-mer based tools, CIRI2 and KNIFE with reads mapping based had relatively superior and balanced detection performance regardless of the cell line or RNase R (-/+) datasets. In summary, the novel k-mer based software show dominant performance on sensitivity and computational efficiency in circRNA discovery. This study may provide new insights into development and application in circRNA detection tools.
The real datasets are available from the NCBI Gene Expression Omnibus (GEO) database (BioProject: PRJNA231724, GEO: GSE53327) followed by different RNase R treatment in two cell lines.
Useage
bash.sh: one-step shell script of general circRNA-seq analysis pipeline for all software tools
results
the raw, filtered or/and annotated predicted candidates by each software under different dataset (can be downloaded from figshare).
circ_candidates.bed: raw identified circRNA bed format file for each software package
circ_candidates_convert.bed: coonverted circRNA bed format file (genome coordinate converted to uniform 0-based format)
circularRNA_known.txt: annotated circRNA information file generated from CIRCexplorer2 annotate moudle
circRNA_known_annotated.txt: it contains annotated circRNA list (circularRNA_known.txt) merged by the known circRNA retrieved from circBase and circAtlas database
performance.csv: recorded prediction performance indices of all software tools

scripts: downstream analysis scripts for data clean, analysis and visualization with python or R
conversion.r: convert raw genome coordinate (0-based or 1-based) to 0-based.
annotation.r: merge annotated circRNA list (circularRNA_known.txt) by the known circRNA retrieved from circBase and circAtlas database
pr-RnaseR.py: data clean and downstream analysis python script

shell_scripts: general prediction pipeline for each software from short reads



