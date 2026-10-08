# Systematic benchmark of ancient DNA read mapping

Code and pipelines behind:

> Oliva A, Tobler R, Cooper A, Llamas B, Souilmi Y. **Systematic benchmark of ancient DNA read mapping.** *Briefings in Bioinformatics* (2021). [doi:10.1093/bib/bbab076](https://doi.org/10.1093/bib/bbab076)

Follow-up comparing the BWA-mem settings proposed by Xu et al. (2021) with the best settings from this benchmark:

> Oliva A, Tobler R, Llamas B, Souilmi Y. **Additional evaluations show that specific BWA-aln settings still outperform BWA-mem for ancient DNA data alignment.** *Ecology and Evolution* (2021). [doi:10.1002/ece3.8297](https://doi.org/10.1002/ece3.8297)

## What the study did

Ancient DNA (aDNA) reads are short (~30–80 bp) and damaged, which makes them hard to map and pushes them towards the reference allele (reference bias). The usual aDNA mapping settings had barely changed in ten years, so we tested them against newer options.

- **30 mapping strategies** across four mappers: BWA-aln, BWA-mem, Bowtie2 and NovoAlign, including an IUPAC-augmented reference for NovoAlign.
- **Simulated reads with a known true position**, generated from human chromosome 22 for individuals from three populations (1000 Genomes), cut to aDNA-like lengths and given post-mortem damage with gargammel.
- **Each strategy scored on** mapping precision and accuracy, reference bias at known SNPs, the downstream effect on population-genetic statistics (D-statistics, PCA), and run time.

**Main findings**
- Well-tuned BWA-aln remained one of the most reliable choices for aDNA.
- Specific NovoAlign and BWA-mem settings also reached high precision with low reference bias.
- Filtering out reads with low mapping quality reduced reference bias for every mapper.
- Unbiased NovoAlign results needed the IUPAC reference. This shows the value of putting population variation into the reference, as pangenome graphs do.

## Repository layout

```
Simulations/
  CreateReads/      simulate aDNA reads
    [Pipeline]-Create_aDNA_reads(Mitty).py   full simulation pipeline (exported from Jupyter)
    gen_reads.sh                             Mitty wrapper: filter variants, generate and corrupt reads
  SubSampleReads/   downsample mapped reads to an empirical aDNA read-length distribution
    Pipeline.sh     entry point (SLURM job script)
    GetMapReads.sh, CreateTSV.sh, TSVlength.py, main.py, PicardRun.sh
    Motala1_read_length_dist_chr22.txt       empirical length distribution (Motala, doi:10.1038/nature13673)
Analysis/
  CalculateMappingDistance.R   distance between true and mapped position for every read,
                               corrected for strand and soft/hard clipping
```

## Running it

**1. Simulate reads** (`Simulations/CreateReads`)
- Requires [Mitty](https://github.com/sbg/Mitty), [gargammel](https://github.com/grenaud/gargammel), Python 3 with Biopython, bgzip and tabix.
- Inputs: the GRCh37 reference (`human_g1k_v37`) and the 1000 Genomes phase 3 chr22 VCF, filtered to the individual you want to simulate. The script defaults to NA19471.
- Edit the variables at the top of the pipeline script, then run it. Reads are generated from 20 to 169 bp in 15 equal-sized length bins, then damaged with gargammel.

**2. Map** the simulated FASTQ files with each mapper and parameter set you want to compare. The parameter sets are listed in the paper.

**3. Downsample** (`Simulations/SubSampleReads`)
- Requires Picard, SAMtools and Python 3.7, with all BAM/SAM files in one directory.
- Run `bash Pipeline.sh`, or submit it with `sbatch` on a SLURM cluster.

**4. Score** (`Analysis`)
- Requires R with `data.table` and `stringi`.
- Put the SAM files in `sam/`, then run `CalculateMappingDistance.R`. It writes one `MapDist_<sample>.txt` per input, with true position, mapped position, MAPQ, CIGAR and mapping distance for each read.

## Citation

If you use this code, please cite the *Briefings in Bioinformatics* paper above.

## License

GPL-3.0. See [LICENSE](LICENSE).
