# Rumex meiotic drive
This repository contains scripts used in the project of identifying  meiotic drive in _Rumex_ using pooled sequencing.

There're three folders including scripts for different analysis, in each script there're helpful comments and usage.
* transcriptome (analyze transcriptome alignments)
    * ```get_freq.py``` : allele frequency quantification without mapping bias correction
    * ```get_freq_nobias.py``` : allele frequency quantification with mapping bias correction
* genome (analyze genome alignments)
    * ```filt_vcf_pool1.py``` : filter and extract allele depths in pool 1 from VCF
    * ```freq_window.py``` : allele frequency quantification in sliding windows without mapping bias correction
    * ```freq_window_nobias.py``` : allele frequency quantification in sliding windows with mapping bias correction (still revising the bias correction formula)

* AF_DP (plot the average allele frequency at different depths)
    * ```get_freq_depth.py``` : get allele frequency and the depth at each informative site
    * ```get_depth_data.py``` : calculate the average allele frequency at different depths


