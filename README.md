# Assignment 1: Part 1 
### Introduction: Biological and Analytical Background

The primary goal of genome assembly and reference-based alignment is to reconstruct an accurate and contiguous representation of an organism’s genome that enables reliable identification of genetic variation. In bacterial pathogens such as *Salmonella enterica*, whole-genome sequencing (WGS) provides the resolution required to detect variation at the single-nucleotide, gene, and structural levels, supporting applications including comparative genomics, phylogenetic inference, and antimicrobial resistance surveillance (McDermott et al., 2016). A high-quality assembly is essential for these analyses, as fragmentation or sequence inaccuracies can obscure biologically meaningful differences when genomes are aligned to reference sequences. Consequently, genome assembly is not merely a preprocessing step, but a foundational analytical stage that directly influences the accuracy of reference-based comparisons and downstream interpretation.  

Multiple assembly strategies have been developed to address these goals, each involving distinct trade-offs between contiguity, accuracy, and computational complexity. Short-read sequencing approaches offer high per-base accuracy but often produce fragmented assemblies due to limited read length, particularly in repetitive regions and plasmids (Taylor et al., 2019). In contrast, long-read sequencing technologies such as Oxford Nanopore sequencing generate reads capable of spanning complex genomic regions, enabling assemblies that approach complete bacterial chromosomes (Xu et al., 2020). However, these advantages are offset by higher raw error rates, especially insertion–deletion errors in homopolymeric regions, which can propagate into assemblies and affect reference-aligned variant detection (Xu et al., 2020; McDermott et al., 2016). Benchmarking studies have further demonstrated that assembly outcomes are sensitive to the choice of assembly algorithm and parameter settings, with statistically significant differences observed among long-read assemblers in sequence accuracy and variant-level metrics, even when overall biological conclusions remain similar (Chen et al., 2020).   

In addition to sequencing technology, genome assembly approaches differ in whether genomes are reconstructed de novo or guided by an existing reference. De novo assembly avoids reference bias and enables discovery of novel genomic content but is highly sensitive to sequencing depth and coverage uniformity and may fail to recover complete genomes in low-abundance or heterogeneous datasets despite long-read data (Gauthier et al., 2025). In contrast, reference-guided assembly leverages prior genomic knowledge to improve sensitivity and completeness when a closely related reference is available, although this approach introduces the risk of reference bias. Together, these trade-offs emphasize that assembly strategy selection is context-dependent and must be aligned with the biological question and downstream analytical goal. Overall, these approaches differ in their ability to balance assembly contiguity, base-level accuracy, sensitivity to low-coverage regions, and susceptibility to reference bias, highlighting that no single method is universally optimal.  

Given the availability of Oxford Nanopore long-read data, a long-read–based assembly strategy is commonly used for genome reconstruction. (Xu et al., 2020; Taylor et al., 2019). Such an approach is expected to improve assembly contiguity and enable resolution of repetitive genomic regions compared to short-read–only methods (Taylor et al., 2019); however, it also introduces challenges related to higher sequencing error rates and parameter sensitivity, particularly for consensus polishing and variant detection (Xu et al., 2020; Chen et al., 2020). Parameters governing read filtering, assembly stringency, and polishing depth can influence the balance between contiguity and base-level accuracy, underscoring the need for careful evaluation of these choices in downstream analyses (Chen et al., 2020).  

### Proposed Methods  

The planned workflow follows a feasible genomics pipeline consisting of read quality assessment, reference assembly, consensus polishing, reference-based alignment, variant detection, and visualization.  

1.	*Sequencing Data and Quality Assessment:* Raw Oxford Nanopore sequencing reads (FASTQ format) generated using R10 chemistry (Q20+, N50 ~5–15 kb) for Salmonella enterica will be used as input for all downstream analyses. Read quality and length distributions will be assessed using NanoPlot v1.43.0 to confirm expected read length and quality profiles prior to assembly. Adapter trimming and filtering will not be applied initially, as modern long-read assemblers are robust to moderate sequencing error rates.

2.	*Genome Assembly:* De novo assembly of the Salmonella enterica genome will be performed using Flye v2.9, a long-read assembler optimized for Nanopore sequencing data. Flye will be executed using the --nano-hq parameter to reflect the high-accuracy R10 chemistry, with default genome size estimation and multithreading enabled. Flye was selected due to its strong performance on bacterial genomes and ability to resolve repetitive regions using long reads, as demonstrated in recent benchmarking studies.(Wick and Holt, 2021)

3.	*Assembly Polishing and Quality Evaluation:* Conducted using Medaka v1.11.3, which applies a neural network–based model trained on Oxford Nanopore data to improve base-level accuracy. Assembly quality will be evaluated using QUAST v5.2.0, assessing metrics including total assembly length, N50, number of contigs, and GC content. These metrics will be compared to expected values for Salmonella enterica to identify potential assembly fragmentation or structural inconsistencies. (Luan et al, 2024)

4.	*Reference Genome Alignment:* The assembled genome will be aligned to a publicly available Salmonella enterica reference genome downloaded from NCBI using Minimap2 v2.26. Alignment will be performed using the -ax asm5 preset, which is optimized for closely related bacterial genomes with low divergence, producing accurate whole-genome alignments. Resulting SAM files will be converted to sorted and indexed BAM files using Samtools v1.19 to facilitate downstream analysis and visualization. (Wick and Holt, 2021; Banovic et al, 2024)

5.	*Variant Calling and Visualization:* Identified using Sniffles2 v2.2, enabling detection of single-nucleotide variants, insertions, deletions, and structural variants from long-read alignments. Identified variants will be summarized and visualized using IGV v2.16, allowing manual inspection of alignment depth, variant support, and potential assembly errors. (Hall et al, 2024) 

6.	*Workflow Management and Reproducibility:* All analyses will be conducted in a Unix-based environment, with scripts version-controlled using Git v2.43 and hosted in a public GitHub repository. Software versions, parameters, and intermediate outputs will be documented to ensure reproducibility and facilitate transition to final methods and code implementation in Part 2 of the assignment






 ### References 

- Banović Đeri, B., Nešić, S., Vićić, I., Samardžić, J., & Nikolić, D. (2024). Benchmarking of five NGS mapping tools for the reference alignment of bacterial outer membrane vesicles–associated small RNAs. *Frontiers in Microbiology, 15*, Article 1401985.

- Chen, Z., Erickson, D. L., & Meng, J. (2025). Benchmarking long-read assemblers for genomic analyses of bacterial pathogens using Oxford Nanopore sequencing. *International Journal of Molecular Sciences, 26*(1).

- Gauthier, J., Mohammadi, S., Kukavica-Ibrulj, I., Boyle, B., Landgraff, C., Goodridge, L., White, K., Chapman, B., & Levesque, R. C. (2025). Leveraging artificial intelligence community analytics and nanopore metagenomic surveillance to monitor early enteropathogen outbreaks. *Frontiers in Public Health, 13*, Article 1675080.

- Hall, M. B., Wick, R. R., Judd, L. M., Nguyen, A. N., Steinig, E. J., Xie, O., Davies, M., Seemann, T., Stinear, T. P., & Coin, L. J. M. (2024). Benchmarking reveals superiority of deep learning variant callers on bacterial nanopore sequence data. *eLife, 13*, e98300.

- Hong, Y.-P., Chen, B.-H., Wang, Y.-W., Teng, R.-H., Wei, H.-L., & Chiou, C.-S. (2024). The usefulness of nanopore sequencing in whole-genome sequencing-based genotyping of *Listeria monocytogenes* and *Salmonella enterica* serovar Enteritidis. *Microbiology Spectrum, 12*(7), e00509-24.

- Luan, T., Commichaux, S., Hoffmann, M., Jayeola, V., Jang, J. H., Pop, M., Rand, H., & Luo, Y. (2024). Benchmarking short and long read polishing tools for nanopore assemblies: Achieving near-perfect genomes for outbreak isolates. *BMC Genomics, 25*, 679.

- McDermott, P. F., Tyson, G. H., Kabera, C., Chen, Y., Li, C., Folster, J. P., Ayers, S. L., Lam, C., Tate, H. P., & Zhao, S. (2016). Whole-genome sequencing for detecting antimicrobial resistance in nontyphoidal *Salmonella*. *Antimicrobial Agents and Chemotherapy, 60*(9), 5515–5520.

- Taylor, T. L., Volkening, J. D., DeJesus, E., Simmons, M., Dimitrov, K. M., Tillman, G. E., & Afonso, C. L. (2019). Rapid, multiplexed, whole genome and plasmid sequencing of foodborne pathogens using long-read nanopore technology. *Scientific Reports, 9*, 16350.

- Wick, R. R., Howden, B. P., & Stinear, T. P. (2025). Autocycler: Long-read consensus assembly for bacterial genomes. *Bioinformatics, 41*(9), btaf474.

- Wick, R. R., & Holt, K. E. (2021). Benchmarking of long-read assemblers for prokaryote whole genome sequencing (Version 4). *F1000Research, 8*, 2138.

- Xu, F., Ge, C., Luo, H., Li, S., Wiedmann, M., Deng, X., Zhang, G., Stevenson, A., Baker, R. C., & Tang, S. (2020). Evaluation of real-time nanopore sequencing for *Salmonella* serotype prediction. *Food Microbiology, 89*, 103452.

