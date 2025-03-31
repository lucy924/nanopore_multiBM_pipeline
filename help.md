# Helpful information with this pipeline


## Run dag only
snakemake --dag get_input_for_minknow | dot -Tsvg > dag.get_input_for_minknow.svg
snakemake --dag get_report_for_sample | dot -Tsvg > dag.get_report_for_sample.svg


For when I'm ready to publish:  
https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility

panel_metadata layout

| ID  | Gene name | Variant Type (snv, sv, mod) | Is this record the whole gene? | Is variant in coding region? | SNP ID     | Result Options | Illumina EPIC ID | DNA methylation | BCG response characteristic | BCG non response characteristic | Notes                                                                 | WARNINGS | Ref                  | chrom | start pos  | end pos    | length | strand |
| --- | --------- | --------------------------- | ------------------------------ | ---------------------------- | --------- | ------- | ---------------- | --------------- | --------------------------- | ------------------------------ | -------------------------------------------------------------------- | -------- | -------------------- | ----- | ---------- | ---------- | ------- | ------ |
| 1   | ATG2B     | snv                         | No                             | Yes                          | rs3759601 | G/C     |                  |                 | G                           | CC                             | Allele not associated with cancer recurrence in Asian population but population did not have any CC genotype |          | Buffen, K. et al. A... | chr14 | 96,311,131 | 96,311,131 | 1       | -      |
| 2   | NOS2      | sv                          | No                             | No                           |           | (CCTTT)n|                  |                 | <13 repeats                 |                                | NOS2-promoter microsatellite (CCTTT)n, if <13 repeats likely to benefit from treatment. primers used for genomic location came from here: https://journals.sagepub.com/doi/10.1177/20587392211052948 |          | Ryk, C. et al. Outc... | chr17 | 27,803,068 | 27,803,263 | 196     | -      |
| 4   | NOS3      | snv                         | No                             | Yes                          | rs1799983 | T/A/G   |                  |                 | GG                          |                                | homozygous GG had decreased risk - note that GT/TT patients were grouped together |          | Ryk, C. et al. Outc... | chr7  | 150,999,023 | 150,999,023 | 1       | +      |
| 5   | NOS3      | snv                         | No                             | No                           | rs2070744 | C/T     |                  |                 | TT                          |                                | note that CT/CC patients were grouped together. SNP in promoter region |          | Ryk, C. et al. Outc... | chr7  | 150,992,991 | 150,992,991 | 1       | +      |
| 7   | FLT1      | mod                         | No                             | No                           |           |         | cg26544530       | promoter        | reduced methylation         | increased methylation         | PLACEHOLDER - check it actually corresponds to the gene                |          | Ilijazi, D. et al. ... | chr13 | 28,896,826 | 28,896,827 | 2       | -      |
| 8   | ARID1A    |                             | Yes                            |                              |           |         |                  | wildtype        | mutations                   | any/all mutations associated with BCG non-response                    |          |          | Pietzak, E. J. et a... | chr1  | 26,728,912 | 26,780,756 | 51,845  | +      |
| 9   | CHST11    | mod                         | Yes                            |                              |           |         |                  | promoter        |                            | increased methylation         | has 5 CpG sites: Validated as predictive biomarker in study            |          | Ilijazi, D. et al. ... | chr12 | 104,456,948 | 104,762,014 | 305,067 | +      |


"Whole genes" in panel are from gene start to gene end. 2000bp upstream and 1000bp downstream to account for flanking functional regions are added as part of the pipeline.

hg38 is used.

https://nanoporetech.com/document/adaptive-sampling  
As a rule of thumb, for a library with a normal distribution with an N50 of ~8 kb, aim for adding 20 kb of buffer.

# To Do
- [ ] Before more sequencing happens need to ensure references probes are all covered. Code in make_panel_bed misses out some genomic locations as they are not in my current epic probe genomic locations file
- [ ] write script to combine multiple sample results for biomarker classification section

In setup need to give sh scripts execution permission (i.e. make_adaptive_ref.sh)

To get input for minknow:
`snakemake get_input_for_minknow`  
Lowest flanking should be 1000bp

Basecalling - this should already be done by minknow, basecalled at presumably HAC during sequencing and alignment output already.

MinKNOW running basecalling for testing:  
- pod5 files here: `/external/data/lucy/bm_minknow_output/test1`
- set up instructions here: `/home/dejlu879/ProjectProtocol/bm_pipeline_dev/minknow_instructions`
- had to change a bunch of permissions to 777 so minknow could access them
- settings:
    - HAC
    - mCG_hmCG
    - 5Hz
    - FloPRO114M
    - alignment to hg38
- outputs tree: 
    - /external/data/lucy/bm_minknow_output/test1
        - PAK79113_95da30ad_e3be95f9_6.pod5
        - PAK79113_95da30ad_e3be95f9_7.pod5
        - basecalling
            - fail
                - bam_runid_*.bam
                - fastq_runid_*.fastq.gz
            - ont_basecall_client_log-2025-02-10_16-26-02.log
            - pass
                - bam_runid_*.bam
                - fastq_runid_*.fastq.gz
            - sequencing_summary.txt
            - sequencing_telemetry.js


Test1 was with 2 pod5 files /external/data/lucy/bm_minknow_output/test1/basecalling/pass
Test2 is with early_pod5s: /external/data/lucy/20241107-BiomarkerProject/early_pod5s/basecalling/pass  
nextflow output:  
    ```
    Completed at: 13-Feb-2025 15:16:10
    Duration    : 56m 51s
    CPU hours   : 9.6
    Succeeded   : 894
    ```

Methatlas - will need to think on a better algorithm for this kind of data, and be mindful of licenses/permissions/acknowledgments.  
Methatlas notes:  
- Platelet isn't available as they don't have gDNA. Future investigations may reveal the ability to use mDNA signatures for this
- https://www.sciencedirect.com/science/article/pii/S0167527323016960
- Th1 isn't available with this reference
- Dendritic cells aren't available with this reference
- M2 macrophages specifically aren't available with this reference

### MethylCIBERSORT set up  
```
cd workflow/scripts/MethylCIBERSORT
mamba create -n mCS_venv -c conda-forge r-base=3.6.3 r-quadprog r-matrixstats r-nmf r-doparallel r-glmnet r-caret bioconda::bioconductor-minfi bioconda::bioconductor-multtest
mamba activate mCS_venv
R CMD INSTALL MethylCIBERSORT_0.2.0.tar.gz
Rscript -e 'options(repos = c(CRAN = "https://cran.stat.auckland.ac.nz")); if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")'
```

### 1.4 Immune infiltrate ratios of interest

- LMR, lymphocyte-to-monocyte ratio
  - LMR > 3.25 predictive of progression after BCG (more so than NLR)
- NLR, Neutrophil to Lymphocyte Ration
- PLR, platelet-to lymphocyte ratio
- Th1 infiltrate
  - https://pubmed.ncbi.nlm.nih.gov/37304240/
  - entry 70 in lit review
- tumour infiltrating dendritic cells and M2 macrophage infiltration are associated with poor response to BCG therapy
- increasing CD4T and CD8T cell proportions were associated with a statistically significant decreased hazard of tumor recurrence or death
  - peripheral blood
- higher levels of transcripts reflecting the abundance of B cells, tumor-associated M1/M2 macrophages, CD8+ T cells, and Tregs. - in BCG unresponsive tumours
  - High TAM infiltration
- CD8 cells in BCG unresponsive vs responsive tumours
- PDL1+ subset of Treg cells (FOXP3+) associated with recurrence after BCG therapy
- DC CD83+ associated with poor response to BCG
- high CD4T cell count in BCG pretreatment microenvironment have a significantly prolonged recurrence-free survival compared to patients with a low CD4T cell count

Other ratios of interest
- IL-6/IL-10 ratio > 0.10 associated with lower recurrence rate after BCG (I think this is just cytokine testing)
  - Can maybe get with methylation of promoters of these two genes

Methatlas does not work. methylcibersort does.  

Advise running with snakemake --use-conda --conda-frontend mamba

Installing methylcibersort  
had to open an R console to install  
install.packages("quadprog")  

Cibersortx docker token: 
New token: 62c6ce7d617349ca44a7ee8efb7df65c, valid until 2025-08-18 11:23:15


ignoring hmC for now (not integrated into mC, just removed from analysis altogether).  
If we want to add it back in, look at convert_DSS.py

Immune infiltrate stuff isn't working at the moment.  
Focusing on scoring parts.

Modification classification of "reduced methylation" and "increased methylation" may need to be specifically thresholded during the biomarker scoring in pre-clinical trial depending on the scoring method

Threshold for saying region is methylated: 80%
Error's +/- 10% methylation  
with a >10% difference between the two groups  

you might get odd errors if there are too many bam files.
If you get something like this:  
```
ERROR ~ Error executing process > 'ingress:catSortBams (1)'

Caused by:
  Process `ingress:catSortBams (1)` terminated with an error exit status (1)

Command executed:

  samtools cat -b <(find input_bams -name 'reads*.bam' | sort)     | samtools sort - -@ 2 --write-index -o reads.bam##idx##reads.bam.bai

Command exit status:
  1

Command output:
  (empty)

Work dir:
  /home/dejlu879/ProjectProtocol/bm_pipeline_dev/results/test2/wf-humvar/work/e2/0c8a013e3255204a2112f7277ecac6
  ```
This might be the problem.

## Modification thresholds
80% methylated is considered methylated
80% methylated across a region is considered a methylated region


## Exp Ratio
Input has to be exactly as I've done it, in "Expression Ratio Components" referring to the gene names in the panel.
These numbers are representative, actual correlation of methylation to expression will need to be determined for each component of the expression ratio.

## Running the file
snakemake --dag run_sample_panel_for_BM_classifier | dot -Tsvg > dag.run_sample_panel_for_BM_classifier.svg
snakemake run_sample_panel_for_BM_classifier --use-conda

snakemake run_sample_panel_for_BM_classifier --use-conda --report --benchmark-extended

*Figuring out SLC11A1 microsatellite alleles*  
TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGTGTGTG|TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGTGTGTG/TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGTGTGTG|TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGT/TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGTGTGTG|TGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGTG/TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGT|TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGT/TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGT|TGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGTG/TGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGTG|TGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGTG/

TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGTGTGTG|TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGT/TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGT|TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGT/TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGT|TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGT/TGTGTGTGTGTACGTGTGTGTGTACGTGTGTGTGTGTGTGTGT|TGTGTGTGTGTACGTGTGTGTGTGTGTGTGTGTG

## Panel development
For FLT1's single location, I have based the response/non-response characteristic on the paper Fig2A for estimating highest and lowest methylation in this gene (2.0 == 100% methylated, -0.6 == 35% methylated respectively), and then used these numbers to figure out whereabouts the "differences" might encompass on Fig4, FLT1.  
Based on attempting to do something along the lines of +/- 10% for error with a 10% difference between responder and non-responder, I am setting a cap for the methylation in the scoring panel for 60% or below methylation belonging to responders and 70% or higher belonging to non-responders. This is calculated by 50 + 10 and 80 - 10. Any other error being added uses the associated HR/SD.  
This method will need reviewing.  
This is a ballpark set of numbers for testing the pipeline, these numbers should be calculated based on individual clinical trials for a whole panel together.  

## Advice:
https://community.nanoporetech.com/posts/adaptive-sampling-reads
  Note other reply: pad your target regions by your expected read length
  
  Hi Daria,
  Just like Chris said, you should be considering your padding/buffer regions according to your fragment length distribution. We are actually looking into this in detail now trying to compile a range of best conditions and settings to use with adaptive sampling to optimize enrichment. What I would advise at the moment is to buffer your target region possibly up to N90 or even N99 of your fragment length distribution if you are targeting 1% of the genome/total DNA sequence in the sample. This will also be dependent on the size of your target region, so we expect that as you target larger regions it wouldn't be as beneficial to use such large buffer regions.
  
  To answer your questions:
  If the first bases are mapped within you bed file in the adaptive sampling settings, it will sequence the remaining of the strand. For the bed file in the adaptive sampling, add your regions of interest + buffer size. However, make sure in the alignment settings you load a bed file carrying exclusively the target region, without the buffer.
  I advise extending your buffer regions up to N90 of your fragment length distribution to try and maximize the enrichment from adaptive sampling.
  
  More detailed information on this feature will come out soon.
  
Miguel

snakemake
8.27.1

To run:  
mamba create -n bm_min_venv conda=24.7.1 peppy=0.40.7 snakemake=8.27.1 
conda config --set channel_priority strict
