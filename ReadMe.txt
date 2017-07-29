Program name:	Con_Sea_Identification_and_PIC_Calculation.pl

Function:	1) Identification of the Con_Island and Con_Sea regions.
	2) Extraction of the sequences of each Con_Island and Con_Sea regions.
	3) Calculation of SNP No. , InDel No. and PIC (potentially informative characters) values of each Con_Sea regions in each sample pairs.

Definition:	 Based on alignment result, conservative sites among the multi-cp genomes were labeled. Regions containing over several (default: 50) continuous conserved sites in a cross-genus consensus genome sequence were identified and defined as "Con_Islands", while the regions between two adjacent Con_Islands were named "Con_Seas".

Citation:	Chao Feng, Meizhen Xu, Chen Feng, Eric J. B. von Wettberg, Ming Kang. The complete chloroplast genome of Primulina and two novel strategies for development of high polymorphic loci for phylogenetic studies. BMC Evolutionary Biology 2017;

Author:	Chao Feng (chaofeng@scbg.ac.cn)

Release date:	Jul. 18th, 2017

Version:	Version 1.0

Usage:	perl Con_Sea_Identification_and_PIC_Calculation.pl $inputfile $Con_Island_size
		$inputfile: The result multiple sequence alignment (in fq format)
		$Con_Island_size: Minimum length of Con_Island
	E.g.: perl Con_Sea_Identification_and_PIC_Calculation.pl File_S1.fas 50

Outputfile:	Con_Island.txt: The info of Con_Island, including No., start site, end site and length.
	Con_Island.seq.txt: Consensus sequences of each Con_Island regions.
	Con_Sea.txt: The info of Con_Sea, including No., start site, end site and length.
	Con_Sea.seq.txt: The sequences of each Con_Sea regions for each sample.
	Con_Sea.PIC.txt: SNP No. , InDel No. and PIC values of each Con_Sea regions in each sample pair.
