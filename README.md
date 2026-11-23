# nf-aligning

nf-aligning is a Nextflow pipeline that aligns DNase-seq and ATAC-seq short reads to a reference genome and calls peaks with [hotspot3](https://github.com/vierstralab/hotspot3). Reads can come from local FASTQ files or be downloaded from the Sequence Read Archive (SRA) by run accession.

The main steps are:
- **Download** *(SRA mode only)*: fetch each run with `prefetch` and `fasterq-dump`, and its metadata with `ffq`.
- **Adapter trimming** with `fastp`, with quality and length filtering disabled.
- **Alignment** with `bwa aln`, followed by `bwa sampe` (paired-end) or `bwa samse` (single-end).
- **Read filtering**: `bin/filter_reads.py` **flags** (does not remove) reads failing filters: QC-fail (`0x200`) for MAPQ < 10, more than 2 mismatches, unpaired or improperly paired reads, or insert size > 750 bp; non-nuclear (`0x1000`) for reads off the nuclear chromosomes.
- **Duplicate marking**: runs are sorted and merged per `sample_id`, and duplicates are **marked** (not removed) with Picard `MarkDuplicatesWithMateCigar`.
- **Filtered CRAM**: unmapped and QC-fail reads are dropped (`samtools view -F 516`), and with `save_cram_mode = "nuclear"` only reads on nuclear chromosomes are kept (default `all`). The result is saved as `<sample_id>.filtered.cram`. Sequencing stats are computed on the duplicate-marked reads before this filter, and insert size metrics (paired-end only) after it.
- **Peak calling** with [hotspot3](https://github.com/vierstralab/hotspot3) on the filtered CRAM, at each FDR in `params.fdrs`.

The per-sample outputs (filtered CRAM, peaks, normalized density, fit statistics, p-values) are the input to [nf-index](https://github.com/vierstralab/nf-index), which builds the consensus DHS index.

⚠️ Check the duplication rate of each library. hotspot3 counts duplicates, and above ~50% duplication its background model breaks down and peak calling can produce unexpected results. Downsample such libraries first, see [Libraries with a high duplication rate](#libraries-with-a-high-duplication-rate).

Coding agents: read [AGENTS.md](AGENTS.md) first.

## Table of contents

- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [1. Clone the repository](#1-clone-the-repository)
  - [2. Build the container](#2-build-the-container)
  - [3. Set the parameters](#3-set-the-parameters)
  - [4. Configure the compute environment](#4-configure-the-compute-environment)
- [Parameters](#parameters)
- [Usage](#usage)
  - [Align local FASTQ files](#align-local-fastq-files)
  - [Align reads downloaded from SRA](#align-reads-downloaded-from-sra)
  - [QC workflows](#qc-workflows)
  - [Subsampling workflows](#subsampling-workflows)
  - [Libraries with a high duplication rate](#libraries-with-a-high-duplication-rate)
  - [Run pipeline steps separately](#run-pipeline-steps-separately)
- [Output files](#output-files)
  - [Main outputs](#main-outputs)
  - [Peak and hotspot calls](#peak-and-hotspot-calls)
- [Authors](#authors)

# Installation

## Prerequisites
- [Nextflow](https://www.nextflow.io/)
- Apptainer/Singularity, to build and run the pipeline's container
- A compute cluster or other Nextflow executor (see [4. Configure the compute environment](#4-configure-the-compute-environment))

## 1. Clone the repository
```bash
git clone https://github.com/vierstralab/nf-aligning.git
cd nf-aligning
```

## 2. Build the container

All tools run inside a single Apptainer/Singularity image, so no conda environment is needed on the host. The image contains bwa, samtools, Picard, fastp, sra-tools, ffq, preseq, [hotspot3](https://github.com/vierstralab/hotspot3), hotspot1 (used only for the SPOT score) and kentUtils.

```bash
cd container
singularity build --fakeroot ../aligning-pipeline.sif aligning-pipeline.def
```

The build takes about 15 minutes. It writes the image to the repository root, where `params.container` (`${projectDir}/aligning-pipeline.sif`) expects it; if you put it elsewhere, point `params.container` at it. Rebuild the image whenever a file in `container/` changes.

`--fakeroot` builds without root access, but your cluster account needs subuid/subgid ranges configured (ask your cluster admin if the build fails with a fakeroot permission error). With root, use `sudo singularity build ../aligning-pipeline.sif aligning-pipeline.def` instead.

The build files in `container/`:

| File | Purpose |
|---|---|
| `aligning-pipeline.def` | Apptainer/Singularity recipe |
| `environment.yml` | Conda env with the alignment, QC and hotspot3 dependencies; hotspot3 itself is installed with `pip` during the build |

## 3. Set the parameters

Fill in the reference files in `params.config` (see [Parameters](#parameters)). `params.config` is a symlink to the human GRCh38 config, `params.human.config`; `params.mouse.config` is the mm10 version. To switch to mouse, run `ln -sf params.mouse.config params.config`. The paths in both point to the Altius cluster and need to be replaced with your own.

## 4. Configure the compute environment

The pipeline doesn't depend on a specific cluster. `-profile Altius` in the examples below is the SLURM profile used at Altius; add a profile for your own executor to `nextflow.config`, for example:

```groovy
profiles {
    slurm {
        process.executor = 'slurm'
        process.queue = 'your-queue-name'
        singularity.enabled = true
        singularity.autoMounts = true
    }
}
```

and run with `-profile slurm` in place of `-profile Altius`. See the Nextflow docs on [executors](https://www.nextflow.io/docs/latest/executor.html) and [config profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles).

# Parameters

Set in `params.config`; the paths there are examples from the human GRCh38 configuration.

| Parameter | Description |
|---|---|
| `outdir` | Output directory (relative path recommended) |
| `samples_file` | Tab-separated samples file (see [Usage](#usage)) |
| `container` | Path to the pipeline's container image |
| `threads` | CPU threads for download, trimming, alignment, sorting and CRAM conversion |
| `genome_fasta_file` | Reference genome FASTA (with `.fai`). The BWA index must sit next to it, named after the FASTA without its extension (`bwa index -p GRCh38_no_alts GRCh38_no_alts.fa`) |
| `nuclear_chroms` | Text file with nuclear chromosome names, one per line |
| `save_cram_mode` | `all` (default) keeps all reads that pass the filters in the output CRAMs; `nuclear` keeps only reads on nuclear chromosomes |
| `fdrs` | Comma-separated FDR thresholds for peak calling |
| `save_debug` | Keep the hotspot3 `debug/` folder |
| `nuclear_chrom_sizes` | Nuclear chromosome sizes, passed to hotspot3 |
| `mappable_bases` | Mappable bases (BED, gzipped), passed to hotspot3 |

Used only by `subsampling.nf -entry spot1score` (`chrominfo` and `spot1_mappable_file` must still exist for any `subsampling.nf` entry):

| Parameter | Description |
|---|---|
| `chrominfo` | Chromosome coordinates (BED) for hotspot1 |
| `spot1_mappable_file` | Uniquely mappable positions (BED) for reads of length `readlength` |
| `readlength` | Read length (k-mer size) of `spot1_mappable_file` |
| `subsampling_spot1_depth` | Number of reads to subsample to before scoring (default 5,000,000; keep the default) |

# Usage

## Align local FASTQ files

```bash
nextflow run main.nf -profile Altius --samples_file samples.tsv
```

Runs trimming, alignment, filtering and hotspot3 peak calling, writing per-sample outputs to `<outdir>/<sample_id>/` (see [Output files](#output-files)).

`samples.tsv` is tab-separated with a header row:

| Column | Description |
|---|---|
| `sample_id` | Sample identifier; runs with the same `sample_id` are merged after alignment |
| `align_id` | Unique identifier of the run, used in file names. Rows with a repeated `align_id` are dropped |
| `reads1` | Absolute path to the read 1 FASTQ |
| `reads2` | Absolute path to the read 2 FASTQ (paired-end only) |
| `type` | `paired` or `single` |
| `adapterP7` | (Optional) Adapter sequence for read 1; detected by `fastp` if omitted |
| `adapterP5` | (Optional) Adapter sequence for read 2, paired-end only; detected by `fastp` if omitted |

Adapter sequences are cut at the first character that is not A, C, G or T, with a warning.

```tsv
sample_id	align_id	reads1	reads2	type	adapterP7	adapterP5
LIB001	AG001	/data/reads/AG001_R1.fastq.gz	/data/reads/AG001_R2.fastq.gz	paired
LIB002	AG002	/data/reads/AG002.fastq.gz		single
```

## Align reads downloaded from SRA

```bash
nextflow run main.nf -profile Altius -entry alignFromSRA --samples_file sra_runs.tsv
```

Downloads each run with `prefetch` and `fasterq-dump`, detects whether it is paired- or single-end from the downloaded files, then trims, aligns and calls peaks as above. `sra_runs.tsv` has two columns:

| Column | Description |
|---|---|
| `sample_id` | Sample identifier; runs with the same `sample_id` are merged after alignment |
| `align_id` | SRA run accession (e.g. `SRR26815981`) |

```tsv
sample_id	align_id
SRX22512992	SRR26815981
SRX22512992	SRR26815982
```

Downloaded FASTQs are saved to `<outdir>/<sample_id>/<align_id>/`, and the `ffq` run metadata to `<outdir>/<sample_id>/stats/<align_id>/<align_id>_info.json`. Runs that fail to download are skipped and the pipeline continues. Empty placeholder FASTQs are saved for them, so check that every run in `sra_runs.tsv` has non-empty FASTQs.

## QC workflows

`stats.nf` takes existing CRAMs. In the tables below, the nuclear filter is `samtools view -F 516` (drops unmapped and QC-fail reads) on `nuclear_chroms`. `nuclearStats` writes `<sample_id>.total_sequencing_stats.txt`, the same name as the main pipeline, so give it its own `--outdir`.

| Entry | Input columns | Steps | Output in `<outdir>/<sample_id>/` |
|---|---|---|---|
| default | `sample_id`, `cram_file`, `cram_index` | Read counts and Picard duplication metrics | `<sample_id>.sequencing_stats.txt` (filtered aligned reads; duplicates, i.e. flag `1024`), `<sample_id>.spotdups.txt` and `<sample_id>.percent_dup.txt` (Picard duplication metrics recomputed after clearing the existing duplicate flags) |
| `-entry nuclearStats` | `sample_id`, `cram_file`, `cram_index` | Nuclear filter, then `bin/bamcounts.py` | `<sample_id>.total_sequencing_stats.txt` |
| `-entry preseq` | `sample_id`, `cram_file`, `cram_index`, `read_type` (`paired` or `single`) | `preseq c_curve` on reads on `nuclear_chroms` | `<outdir>/<sample_id>.preseq_hist.txt` (directly in `<outdir>`) |

## Subsampling workflows

`subsampling.nf` takes existing CRAMs and calls peaks or computes the SPOT score at a matched depth. `subsampleToFrac` writes files with the same names as the main pipeline (`<sample_id>.total_sequencing_stats.txt` and the hotspot3 outputs), so give it its own `--outdir`.

| Entry | Input columns | Steps | Output in `<outdir>/<sample_id>/` |
|---|---|---|---|
| `-entry subsampleToFrac` | `sample_id`, `cram_file`, `cram_index`, `frac` | Nuclear filter → keep fraction `frac` of read pairs (`samtools view --subsample <frac> --subsample-seed 42`) → re-mark duplicates → hotspot3, `bamcounts.py` stats and CRAM | `<sample_id>.subsampled_pairs.bam(.bai)`, `<sample_id>.marked.cram(.crai)`, `<sample_id>.total_sequencing_stats.txt` and the hotspot3 outputs |
| `-entry spot1score` | `sample_id`, `cram_file` (index at `<cram_file>.crai`) | Read 1 only → nuclear filter → subsample to `subsampling_spot1_depth` reads → hotspot1 SPOT score | `r1.spot.out` (not prefixed with `sample_id`) |

## Libraries with a high duplication rate

hotspot3 fits the background assuming that nonspecific cleavages follow a negative binomial distribution. Duplicates are marked but kept in the filtered CRAM and counted by hotspot3, so in libraries with a high duplication rate this assumption breaks down and peak calling can produce unexpected results. We downsampled such libraries to a duplication rate below 50%, where hotspot3 still works reliably:

1. Compute the library complexity curve with `nextflow run stats.nf -entry preseq` on the filtered CRAMs. `<outdir>/<sample_id>.preseq_hist.txt` lists the number of distinct reads (pairs, for paired-end data) at increasing sequencing depth. The duplication rate at a given depth is `1 - distinct_reads / total_reads`; take the largest depth at which it is below 50% and divide it by the library size to get the fraction to keep:
```python
import pandas as pd

curve = pd.read_table('<sample_id>.preseq_hist.txt')  # total_reads, distinct_reads
total, distinct = curve.iloc[:, 0], curve.iloc[:, 1]
dup_rate = 1 - distinct / total  # the last row is the full library
frac = total[dup_rate < 0.5].max() / total.max()
```
2. For libraries with `frac < 1`, run `nextflow run subsampling.nf -entry subsampleToFrac` with that `frac`. It keeps the fraction of read pairs, re-marks duplicates and calls peaks on the downsampled reads.

## Run pipeline steps separately

The steps can also be run on their own, which is mostly useful for debugging or for finishing samples that stopped partway (if `-resume` can't pick them up).

Align already trimmed FASTQs:
```bash
nextflow run aligning.nf -profile Altius --samples_file trimmed.tsv
```
`trimmed.tsv` has columns `sample_id`, `reads1`, `reads2` and `type`, as in [Align local FASTQ files](#align-local-fastq-files). Outputs are the alignment files listed in [Main outputs](#main-outputs).

Call peaks on existing CRAM (or BAM) files:
```bash
nextflow run hotspots_calling.nf -profile Altius --samples_file crams.tsv
```
`crams.tsv` has columns `sample_id`, `cram_file` and optionally `cram_index` (defaults to `<cram_file>.crai`, created if missing; set it to the `.bai` for BAM files). Outputs are the hotspot3 files listed in [Main outputs](#main-outputs).

# Output files

Results are written to `<outdir>/<sample_id>/`.

## Main outputs

```
<outdir>/<sample_id>/
├── <align_id>/                                     # downloaded FASTQs (SRA mode only)
├── <sample_id>.filtered.cram(.crai)                # filtered, duplicate-marked alignments
├── <sample_id>.total_sequencing_stats.txt          # sequencing stats (bin/bamcounts.py)
├── stats/
│   ├── <align_id>/
│   │   ├── fastp.json, fastp.html                  # adapter trimming reports
│   │   └── <align_id>_info.json                    # SRA run metadata from ffq (SRA mode only)
│   ├── <sample_id>.MarkDuplicates.picard           # Picard duplicate marking metrics
│   └── <sample_id>.CollectInsertSizeMetrics.picard(.pdf)  # Picard insert size metrics (paired-end only)
│
│   # hotspot3
├── <sample_id>.cutcounts.bed.gz(.tbi)              # per-base cut counts
├── <sample_id>.total_cutcounts                     # total cut count
├── <sample_id>.fit_stats.bed.gz(.tbi)              # background model fit statistics per segment
├── <sample_id>.fit_stats.bed.gz.iter<N>(.tbi)      # fit statistics of intermediate refits (only for samples that needed them)
├── <sample_id>.pvals.parquet                       # per-base p-values
├── <sample_id>.normalized_density.bw               # normalized signal density
├── <sample_id>.background.bw                       # estimated background level
├── <sample_id>.peak_calling.log                    # hotspot3 log
├── fdr<X>/                                         # peak and hotspot calls, one folder per FDR in `fdrs`
└── debug/                                          # hotspot3 intermediate files (emptied unless `save_debug = true`)
```

See the [hotspot3 README](https://github.com/vierstralab/hotspot3#output-files) for the format of each hotspot3 file and how to interpret the fit statistics.

## Peak and hotspot calls

Each `fdr<X>/` folder contains:

| File | Content |
|---|---|
| `<sample_id>.hotspots.fdr<X>.bed.gz(.tbi)` (+ `.bb`) | Hotspots: contiguous regions of significant bases |
| `<sample_id>.peaks.fdr<X>.bed.gz(.tbi)` (+ `.bb`) | Peaks: summits within hotspots, trimmed to the width of the peak |

Use peaks for downstream analyses such as building the DHS index.

# Authors

Developed by:

- Sergey Abramov
- Alexandr Boytsov
- Athanasios Teodosiadis
- Jemma Nelson

Altius Institute for Biomedical Sciences
