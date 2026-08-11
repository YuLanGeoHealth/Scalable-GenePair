# Scalable-GenePair
This repo stores code to reproduce the analysis in the manuscript, "Comparing subsampling strategies for efficient pairwise analysis of large pathogen genomic and spatial datasets: an application to Mycobacterium tuberculosis transmission". 

-full_cc_gen.R: Use the case-control approach to subsample the full dataset for GenePair analysis. 

-full_dc_gen.R: Use the divide-and-conquer approach to subsample the full dataset for GenePair analysis.

-sim_cc_gen.R: Use the case-control approach to subsample the benchmark dataset for the simulation study. 

-sim_dc_gen.R: Use the divide-and-conquer approach to subsample the benchmark dataset for the simulation study. 

-sim_cc_results.R: Summary results using case-control approaches from the simulation.

-sim_dc_results.R: Combine and summary results using case-control approaches from the simulation.

Please refer to https://github.com/warrenjl/GenePair for GenePair codes on pairwise analysis on the subsampling dataset. 

The data used in this study contain confidential information and cannot be deposited in a public repository. Sequence data have been deposited in GenBank and are publicly available in the NCBI BioProject database (accession number PRJNA1311115). Aggregate-level demographic, clinical, and lineage data are provided in the appendix of our previous publication via [https://www.thelancet.com/cms/10.1016/j.lanmic.2026.101369/attachment/46b2ca56-601f-48d8-9b47-a199e15a81a3/mmc1.pdf.](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(26)00024-8/fulltext)

Reference: Wu CY, Chen YA, Ioerger TR, Lan Y, Bai RY, Li MH, Lee CH, Lin JN, Lee SS, Chien ST, Warren JL. Lineage-specific transmission and spatial clustering of Mycobacterium tuberculosis in Kaohsiung, Taiwan, in 2019–23: a population-based genomic study. The Lancet Microbe. 2026 Jul 2.
