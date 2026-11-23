# nf-aligning
 
A Nextflow pipeline for aligning DNase-seq/ATAC-seq style short-read data to a reference genome with BWA. Reads can come from local FASTQ files already on the server, or be pulled directly from the Sequence Read Archive (SRA) by accession.
 
The pipeline handles adapter trimming (`fastp`), alignment (`bwa`), and downstream hotspot calling, and can optionally bootstrap its input by downloading raw reads from SRA before running the same trim/align/hotspot-calling steps.
 
## Pipeline Steps
 
1. **Download** *(SRA mode only, `-entry alignFromSRA`)* — fetch raw reads for each `align_id` from the SRA using `prefetch` and `fasterq-dump`.
2. **trimReads** — trim adapter sequences from the reads with `fastp`, quality and length filtering disabled.
3. **alignReads** — align each run with `bwa aln` followed by `bwa sampe` (paired-end) or `bwa samse` (single-end). `bin/filter_reads.py` then **flags** (does not remove) reads failing filters: QC-fail (`0x200`) for MAPQ < 10, more than 2 mismatches, unpaired/improperly paired reads, or insert size > 750 bp; non-nuclear (`0x1000`) for reads off the nuclear chromosomes. Runs are sorted, merged per `sample_id`, and duplicates are **marked** (not removed) with Picard `MarkDuplicatesWithMateCigar`. From the marked BAM:
   - per-sample sequencing stats are written to `<sample_id>.total_sequencing_stats.txt` (`bin/bamcounts.py`);
   - unmapped and QC-fail reads are dropped (`samtools view -F 516`), and when `save_cram_mode = "nuclear"` only reads on nuclear chromosomes are kept;
   - insert-size metrics (paired-end only) are computed on these filtered reads;
   - the filtered reads are saved as `<sample_id>.filtered.cram`.
4. **call_hotspots** — call hotspots with [hotspot3](https://github.com/vierstralab/hotspot3) on the filtered CRAM. hotspot3 is a peak caller for control-free DNase-seq/ATAC-seq data: it fits a local negative-binomial background, tests each base, and reports *hotspots* (contiguous significant regions) and *peaks* (summits within them) at each FDR in `params.fdrs`. It has mainly been tested on datasets with ≥ 10 million tags. See the [hotspot3 README](https://github.com/vierstralab/hotspot3) for the method, options and output formats. This step runs in both the default and `alignFromSRA` workflows. If you don't need it, remove the `| call_hotspots` line from the workflow in `main.nf`.
Additional QC counts (filtered aligned reads, duplicate reads and percent duplication) can be computed from the CRAMs with `nextflow run stats.nf`.
 
### Method summary
 
Reads were trimmed of adapter sequences (`fastp --disable_quality_filtering --disable_length_filtering --detect_adapter_for_pe`) and aligned to either the GRCh38/hg38 or mm10 genome assembly using `bwa aln`.
 
## Requirements
 
- [Nextflow](https://www.nextflow.io/)
- A container runtime (Docker/Singularity/Apptainer) able to run the pipeline's container image, referenced via `params.container`
- A compute cluster/executor configured in `nextflow.config` — see [Compute environment](#compute-environment) below
## Compute environment
 
This pipeline does not require access to any specific cluster. `-profile Altius` is just the example profile used at Altius (SLURM-based); adjust or add a profile in `nextflow.config` to match your own cluster/executor.
 
See the [Nextflow executor documentation](https://www.nextflow.io/docs/latest/executor.html) and [config profiles documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) for how to define a profile for your environment.
 
A minimal example for a SLURM cluster:
 
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
 
Run with `-profile slurm` (or whatever you name it) in place of `-profile Altius`.
 
## Container
 
The pipeline runs entirely inside a single Apptainer/Singularity container — no conda environment needs to be installed on the host. The image bundles bwa, samtools, Picard, fastp, sra-tools/ffq, hotspot1 (used only by the optional SPOT score), [hotspot3](https://github.com/vierstralab/hotspot3) and kentUtils, so alignment, subsampling and hotspot calling all run from the same image.
 
The build files live in `container/`:
 
| File | Purpose |
|---|---|
| `container/aligning-pipeline.def` | Apptainer/Singularity recipe for the image |
| `container/environment.yml` | Conda env for alignment and QC (bwa, samtools, Picard, fastp, sra-tools, …), pinned to Python 3.8 |
| `container/environment.hotspot3.yml` | Separate conda env for hotspot3, which needs Python 3.9+ and newer libraries than the pinned alignment env |
 
They are kept separate so that installing hotspot3 doesn't change the pinned versions of the alignment tools. Both are installed side by side in the image and are on `PATH`.
 
### Building the `.sif` image
 
```bash
cd container
singularity build --fakeroot ../aligning-pipeline.sif aligning-pipeline.def
```
 
Run the build from inside `container/` so the recipe can find the `environment*.yml` files. The image is written to the repository root, which is where `params.container` (`${projectDir}/aligning-pipeline.sif`) expects it.
 
This takes roughly **~15 minutes** (it compiles hotspot1/kentUtils from source and resolves both conda environments).
 
`--fakeroot` lets you build without root access; your cluster account needs subuid/subgid ranges configured for this to work (ask your cluster admin if the build fails with a fakeroot-related permission error). If you have actual root/sudo instead, use `sudo singularity build ../aligning-pipeline.sif aligning-pipeline.def`.
 
If you put the image somewhere else, point `params.container` at it. Rebuild the image any time a file in `container/` changes.
 
## Setup
 
1. Clone the repository and build the container image referenced by `params.container` (see [Container](#container)).
2. Fill in the required file paths in `params.config` (reference genome, chrom sizes, hotspot annotation files, etc. — see [Parameters](#parameters) below).
3. Configure a profile for your compute environment in `nextflow.config` (see [Compute environment](#compute-environment) above).
4. Prepare a tab-separated samples file and point `params.samples_file` at it (see [Input files](#input-files) below).
## Usage
 
### Align local FASTQ files
 
```bash
nextflow run main.nf -profile Altius
```
 
This runs the default workflow:
 
1. Reads `params.samples_file` and groups records by `align_id`.
2. Trims adapters with `fastp` (`trimReads`).
3. Aligns, filters, marks duplicates and writes `<sample_id>.filtered.cram` (`alignReads`).
4. Calls hotspots and peaks on the filtered CRAM with [hotspot3](https://github.com/vierstralab/hotspot3) (`call_hotspots`). These are the outputs `nf-index` needs.
 
### Align reads downloaded from SRA
 
Use this mode to fetch raw reads for one or more SRA run accessions and run them through the same trimming/alignment/hotspot-calling pipeline, without needing to stage FASTQs on disk first.
 
```bash
nextflow run main.nf -profile Altius -entry alignFromSRA
```
 
This runs the `alignFromSRA` workflow:
 
1. Reads `params.samples_file` and, for each unique `align_id`, downloads the run with `prefetch` and converts it to FASTQ with `fasterq-dump` (paired- and single-end output are both attempted; whichever exists is used).
2. Also pulls run metadata via `ffq` and publishes it as `<outdir>/<sample_id>/stats/<align_id>/<align_id>_info.json`, next to that run's `fastp` reports.
3. Detects whether the run is paired- or single-end based on which FASTQ files are non-empty, and drops any run for which no reads were downloaded (e.g. failed/invalid accessions — download errors are ignored rather than failing the pipeline).
4. Trims and aligns exactly as in the default workflow, then calls hotspots with [hotspot3](https://github.com/vierstralab/hotspot3) on the filtered CRAM.
Downloaded FASTQs are published to `<outdir>/<sample_id>/`, so re-running the pipeline can reuse previously downloaded data instead of re-fetching it from SRA.
 
### Call hotspots on existing CRAMs
 
Runs [hotspot3](https://github.com/vierstralab/hotspot3) directly on existing CRAM (or BAM) files, for example to re-call hotspots with different settings without re-aligning:
 
```bash
nextflow run hotspots_calling.nf -profile Altius --samples_file crams.tsv
```
 
`crams.tsv` is tab-separated with columns `sample_id`, `cram_file` and optionally `cram_index` (defaults to `<cram_file>.crai`). If the index is missing it is created on the fly. Outputs are the [hotspot3](https://github.com/vierstralab/hotspot3) files listed under [Output](#output), written to `<outdir>/<sample_id>/`.
 
### QC statistics from existing CRAMs
 
```bash
nextflow run stats.nf -profile Altius --samples_file crams.tsv
```
 
Needs columns `sample_id`, `cram_file`, `cram_index`. Writes `<sample_id>.sequencing_stats.txt` (filtered aligned read count and duplicate count, i.e. reads with SAM flag `1024`) and `<sample_id>.spotdups.txt` / `<sample_id>.percent_dup.txt` (Picard duplication metrics recomputed from scratch).
 
## Input files
 
`params.samples_file` is a **tab-separated** file with a header row. The columns required depend on the entry point:
 
### Default workflow (local FASTQs)
 
| Column | Description |
|---|---|
| `sample_id` | Sample/group identifier used for output paths |
| `align_id` | Unique identifier for this alignment (used for grouping and file naming) |
| `reads1` | Path to read 1 FASTQ |
| `reads2` | Path to read 2 FASTQ (required only when `type` is `paired`) |
| `type` | `paired` or `single` |
| `adapterP7` | (Optional) Adapter sequence for read 1; auto-detected by `fastp` if omitted |
| `adapterP5` | (Optional) Adapter sequence for read 2, paired-end only; auto-detected if omitted |
 
Adapter sequences are sanitized automatically — any trailing ambiguous (non-ACGT) bases are stripped, with a warning printed to the console.
 
Example:
 
```tsv
sample_id align_id reads1 reads2 type adapterP7 adapterP5
LIB001 AG001 /data/reads/AG001_R1.fastq.gz /data/reads/AG001_R2.fastq.gz paired
LIB002 AG002 /data/reads/AG002.fastq.gz single
```
 
### `alignFromSRA` workflow
 
| Column | Description |
|---|---|
| `sample_id` | Sample/group identifier used for output paths |
| `align_id` | SRA run accession (e.g. `SRR12345678`) to download |
 
Example:
 
```tsv
sample_id align_id
SRX22512992 SRR26815981
SRX22512992 SRR26815982
```
 
## Parameters
 
Key parameters set in `params.config` (paths below are examples from a human GRCh38 configuration and should be adjusted for your genome/environment):
 
| Parameter | Description |
|---|---|
| `outdir` | Output directory (relative path recommended) |
| `samples_file` | Path to the tab-separated samples file |
| `container` | Path to the pipeline's container image |
| `threads` | CPU threads for read download/trimming steps |
| `genome_fasta_file` | Reference genome FASTA (indexed for BWA) |
| `nuclear_chroms` / `nuclear_chrom_sizes` | Nuclear chromosome name/size lists, used to restrict analysis to nuclear DNA |
| `save_cram_mode` | Which reads to retain in output CRAMs (`all` or `nuclear`) |
| `fdrs` | Comma-separated FDR thresholds for hotspot calling |
| `save_debug` | Whether to keep intermediate debug output |
| `mappable_bases` | Mappable bases (BED, gzipped) passed to hotspot3 |
 
Used only by the optional SPOT score (`subsampling.nf -entry spot1score`, see [Additional workflows](#additional-workflows)):
 
| Parameter | Description |
|---|---|
| `chrominfo` | Chromosome info BED for hotspot1 |
| `hotspots_dir` | Path to the hotspot1 distribution (the copy inside the container is at `/hotspot1`) |
| `readlength` | Read length (k-mer size) of the hotspot1 mappability file |
| `subsampling_spot1_depth` | Number of reads to subsample to before scoring |
 
## Output
 
Results are published under `<outdir>/<sample_id>/`. Raw FASTQs are nested under a subdirectory named by `align_id` (the SRA run accession in SRA mode); all downstream output files are prefixed with `sample_id` instead:
 
```
<outdir>/<sample_id>/
├── <align_id>/                                          # raw run FASTQs (SRA mode only), e.g. SRR13520417/
├── <sample_id>.filtered.cram(.crai)                      # filtered, duplicate-marked alignments (main alignment output)
├── <sample_id>.background.bw                            # local background signal (bigWig)
├── <sample_id>.cutcounts.bed.gz(.tbi)                    # per-base cut counts, tabix-indexed
├── <sample_id>.fit_stats.tsv.gz(.tbi, .iter1, .iter1.tbi)  # background model fit statistics
├── <sample_id>.normalized_density.bw                     # normalized signal density (bigWig)
├── <sample_id>.peak_calling.log                          # peak calling log
├── <sample_id>.per_segment.background.bw                 # per-segment background signal
├── <sample_id>.pvals.parquet/                            # per-base p-values (partitioned parquet)
├── <sample_id>.thresholds.bw                             # significance thresholds (bigWig)
├── <sample_id>.total_cutcounts                           # total cut count
├── <sample_id>.total_sequencing_stats.txt                # sequencing depth/stats summary
├── debug/                                                # intermediate/debug output (only when `save_debug = true`)
├── fdr0.001/                                             # hotspot calls thresholded at FDR 0.001
├── fdr0.01/                                              # hotspot calls thresholded at FDR 0.01
├── fdr0.05/                                              # hotspot calls thresholded at FDR 0.05
└── stats/
    ├── <align_id>/                                       # one subdirectory per run
    │   ├── fastp.json, fastp.html                        # adapter-trimming reports
    │   └── <align_id>_info.json                          # SRA run metadata from `ffq` (SRA mode only)
    ├── <sample_id>.CollectInsertSizeMetrics.picard(.pdf)  # Picard insert size metrics (paired-end only)
    └── <sample_id>.MarkDuplicates.picard                  # Picard duplicate-marking metrics
```
 
### Which step produces what
 
- **Alignment** (always): `<sample_id>.filtered.cram`, `<sample_id>.total_sequencing_stats.txt` and `stats/`.
- **Hotspot calling** ([hotspot3](https://github.com/vierstralab/hotspot3)): everything else.
### Peaks and hotspots
 
Each `fdr<X>/` directory contains:
 
| File | Content |
|---|---|
| `<sample_id>.hotspots.fdr<X>.bed.gz` (+ `.bb`) | Hotspots: contiguous regions of significant bases |
| `<sample_id>.peaks.fdr<X>.bed.gz` (+ `.bb`) | Peaks: summits within hotspots, trimmed to peak width. Column 4 is the `sample_id` |
 
- Use **peaks** for downstream analyses such as building the DHS index.
- Which `fdr*` directories are created is set by `params.fdrs` (comma-separated FDR thresholds).
- The fit-stats file name depends on the hotspot3 version. Current releases write `<sample_id>.fit_stats.bed.gz`.
 
 
## Notes
 
- SRA downloads use `errorStrategy 'ignore'`: if `prefetch`/`fasterq-dump` fails for a given accession, the pipeline continues rather than aborting, and that run is silently dropped from downstream processing.
- `align_id` must be unique across the samples file — records are deduplicated on this column before processing.
- With the `Altius` profile, tasks that fail for reasons other than running out of memory (exit codes 134–143, retried up to 3 times) are ignored, so a finished run can be missing samples. Check `.nextflow.log` for ignored tasks and confirm every `sample_id` has a `.filtered.cram` (and `fdr*/` directories, if hotspots were expected).
- To re-call hotspots without re-aligning, use `hotspots_calling.nf` on the CRAMs (see [Call hotspots on existing CRAMs](#call-hotspots-on-existing-crams)). The `fromCutcounts` entry in that file is defunct: it reads from a hard-coded directory.
## Additional workflows
 
These entry points are not part of the main run. They take existing CRAMs/BAMs and are mostly used for QC and for depth-matched comparisons. All read a tab-separated `--samples_file` with a header row. Several of them write files with the same names as the main pipeline (for example `<sample_id>.total_sequencing_stats.txt` and hotspot3 outputs), so give them their own `--outdir` to avoid overwriting a main run.
 
### `stats.nf`: QC counts
 
| Entry | Input columns | Steps | Output in `<outdir>/<sample_id>/` |
|---|---|---|---|
| default | `sample_id`, `cram_file`, `cram_index` | `collect_basic_stats` and `percent_dup` | `<sample_id>.sequencing_stats.txt` (filtered aligned reads; duplicates, i.e. flag `1024`), `<sample_id>.spotdups.txt` and `<sample_id>.percent_dup.txt` (Picard duplication metrics recomputed after clearing existing duplicate flags) |
| `-entry nuclearStats` | `sample_id`, `cram_file`, `cram_index` | nuclear filter (`samtools view -F 516` on `nuclear_chroms`), then `bin/bamcounts.py` | `<sample_id>.total_sequencing_stats.txt` |
| `-entry preseq` | `sample_id`, `cram_file`, `cram_index`, `read_type` (`paired` or `single`) | `preseq c_curve` on nuclear reads | `<outdir>/<sample_id>.preseq_hist.txt` (directly in `<outdir>`, not per sample) |
 
`preseq` is not in the container: `run_preseq` uses a hard-coded conda environment (`/home/sabramov/miniconda3/envs/super-index`). Change the `conda` line in `stats.nf` to an environment with `preseq` and `samtools` before using it.
 
### `subsampling.nf`: depth-matched hotspots and SPOT score
 
| Entry | Input columns | Steps | Output in `<outdir>/<sample_id>/` |
|---|---|---|---|
| `-entry subsampleToFrac` | `sample_id`, `cram_file`, `cram_index`, `frac` | nuclear filter → keep fraction `frac` of read pairs (`samtools view --subsample <frac> --subsample-seed 42`, pairs kept together) → re-mark duplicates → hotspot3, `bamcounts.py` stats and CRAM conversion | `<sample_id>.subsampled_pairs.bam(.bai)`, `<sample_id>.marked.cram(.crai)`, `<sample_id>.total_sequencing_stats.txt` and the hotspot3 outputs |
| `-entry filterAndCallHotspots` | `sample_id`, `bam_file`, `bam_index` | nuclear filter → hotspot3 | hotspot3 outputs. Use this for BAMs that have not been through `filter_nuclear`, e.g. marked BAMs from another pipeline |
| `-entry spot1score` | `sample_id`, `cram_file` (index must be `<cram_file>.crai`) | `preprocessBams` (below) → hotspot1 → SPOT score | `r1.spot.out` (hotspot1 SPOT score; the file name is not prefixed with `sample_id`) |
 
`subsampleToFrac` needs a conda environment for its subsampling step: `subsample_with_pairs_frac` uses `conda params.conda`, which isn't set in the params files. Pass `--conda <env with samtools>` or change that line to `container "${params.container}"`.
 
`preprocessBams` is the sub-workflow used by `spot1score`, not an entry point. It keeps read 1 of each pair (`samtools view -F 140`, which also drops unmapped reads and reads with an unmapped mate), applies the nuclear filter, and randomly subsamples to `subsampling_spot1_depth` reads (`bin/random_sample.py`, seed 42; samples with fewer reads are kept whole).
 
`spot1score` runs hotspot1 through `bin/runhotspot.bash`, which also expects a mappability file at `<genome prefix>.K<readlength>.mappable_only.bed` next to `genome_fasta_file`. `hotspots_dir` in the params files points to an Altius path; set it to `/hotspot1` to use the copy inside the container.
 
