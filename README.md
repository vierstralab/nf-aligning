# nf-aligning
 
A Nextflow pipeline for aligning DNase-seq/ATAC-seq style short-read data to a reference genome with BWA. Reads can come from local FASTQ files already on the server, or be pulled directly from the Sequence Read Archive (SRA) by accession.
 
The pipeline handles adapter trimming (`fastp`), alignment (`bwa`), and downstream hotspot calling, and can optionally bootstrap its input by downloading raw reads from SRA before running the same trim/align/hotspot-calling steps.
 
## Pipeline Steps
 
1. **Download** *(SRA mode only, `-entry alignFromSRA`)* — fetch raw reads for each `align_id` from the SRA using `prefetch` and `fasterq-dump`.
2. **trimReads** — trim adapter sequences from the reads with `fastp`, quality and length filtering disabled.
3. **alignReads** — align each run with `bwa aln` followed by `bwa sampe` (paired-end) or `bwa samse` (single-end). `bin/filter_reads.py` then **flags** (does not remove) reads failing filters: QC-fail (`0x200`) for MAPQ < 10, more than 2 mismatches, unpaired/improperly paired reads, or insert size > 750 bp; non-nuclear (`0x1000`) for reads off the nuclear chromosomes. Runs are sorted, merged per `sample_id`, and duplicates are **marked** (not removed) with Picard `MarkDuplicatesWithMateCigar`. From the marked BAM:
   - per-sample sequencing stats are written to `<sample_id>.total_sequencing_stats.txt` (`bin/bamcounts.py`);
   - unmapped and QC-fail reads are dropped (`samtools view -F 516`), and when `save_cram_mode = "nuclear"` only reads on nuclear chromosomes are kept;
   - insert-size metrics (paired-end only) and MACS2 peak calling (if `do_macs`) are run on these filtered reads;
   - the filtered reads are saved as `<sample_id>.filtered.cram`.
4. **call_hotspots** — call hotspots with [hotspot3](https://github.com/vierstralab/hotspot3) on the filtered CRAM. hotspot3 is a peak caller for control-free DNase-seq/ATAC-seq data: it fits a local negative-binomial background, tests each base, and reports *hotspots* (contiguous significant regions) and *peaks* (summits within them) at each FDR in `params.fdrs`. It has mainly been tested on datasets with ≥ 10 million tags. See the [hotspot3 README](https://github.com/vierstralab/hotspot3) for the method, options and output formats. **This step runs automatically only in the `alignFromSRA` workflow.** For local FASTQs it is commented out in `main.nf`; run it afterwards with `hotspots_calling.nf` (see [Call hotspots on existing CRAMs](#call-hotspots-on-existing-crams)).
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
 
The pipeline runs entirely inside a single Apptainer/Singularity container built from `aligning-pipeline.def` — no conda environment needs to be installed on the host. The image bundles hotspot1, [hotspot3](https://github.com/vierstralab/hotspot3), kentUtils, and all of their Python/R dependencies, so alignment, subsampling, and hotspot calling all run from the same container.
 
### Building the `.sif` image
 
From the repository root:
 
```bash
singularity build --fakeroot aligning-pipeline.sif aligning-pipeline.def
```
 
This takes roughly **~15 minutes** to build (it compiles hotspot1/kentUtils from source and resolves two separate conda environments along the way).
 
`--fakeroot` lets you build without root access; your cluster account needs subuid/subgid ranges configured for this to work (ask your cluster admin if the build fails with a fakeroot-related permission error). If you have actual root/sudo instead, use `sudo singularity build aligning-pipeline.sif aligning-pipeline.def`.
 
Once built, point `params.container` at the resulting `aligning-pipeline.sif` path. Every process that references `params.container` runs from that single image — no separate conda setup or module loads needed for subsampling, hotspot calling, or anything else downstream.
 
Rebuild the image any time `aligning-pipeline.def` or one of the `environment*.yml` files it references changes.
 
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
> **Note:** hotspot calling is **not** run in this mode (`call_hotspots` is commented out at the end of the default workflow in `main.nf`). To get peaks and the other [hotspot3](https://github.com/vierstralab/hotspot3) outputs needed by `nf-index`, run [`hotspots_calling.nf`](#call-hotspots-on-existing-crams) on the resulting CRAMs.
 
### Align reads downloaded from SRA
 
Use this mode to fetch raw reads for one or more SRA run accessions and run them through the same trimming/alignment/hotspot-calling pipeline, without needing to stage FASTQs on disk first.
 
```bash
nextflow run main.nf -profile Altius -entry alignFromSRA
```
 
This runs the `alignFromSRA` workflow:
 
1. Reads `params.samples_file` and, for each unique `align_id`, downloads the run with `prefetch` and converts it to FASTQ with `fasterq-dump` (paired- and single-end output are both attempted; whichever exists is used).
2. Also pulls run metadata via `ffq` and publishes it alongside the FASTQs under `<outdir>/<sample_id>/stats`.
3. Detects whether the run is paired- or single-end based on which FASTQ files are non-empty, and drops any run for which no reads were downloaded (e.g. failed/invalid accessions — download errors are ignored rather than failing the pipeline).
4. Trims and aligns exactly as in the default workflow, then calls hotspots with [hotspot3](https://github.com/vierstralab/hotspot3) on the filtered CRAM.
Downloaded FASTQs are published to `<outdir>/<sample_id>/`, so re-running the pipeline can reuse previously downloaded data instead of re-fetching it from SRA.
 
### Call hotspots on existing CRAMs
 
Runs [hotspot3](https://github.com/vierstralab/hotspot3) directly on CRAM (or BAM) files, for example after the default local-FASTQ workflow or to re-call hotspots without re-aligning:
 
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
| `chrom_sizes` / `chrom_sizes_bed` | Genome-wide chromosome sizes |
| `do_macs` | Whether to run MACS2 peak calling |
| `save_cram_mode` | Which reads to retain in output CRAMs (`all` or `nuclear`) |
| `fdrs` | Comma-separated FDR thresholds for hotspot calling |
| `save_debug` | Whether to keep intermediate debug output |
| `mappable_bases` / `mappable` | Mappability annotation files used in hotspot calling |
| `centers` | Center-site annotation file for hotspot calling |
| `density_buckets` | Chromatin density bucket file |
| `chrominfo` | Chromosome info file for the SPOT1 score step |
| `hotspots_dir` | Path to the hotspot1 distribution |
| `readlength` | Read length used by hotspot1 |
| `subsampling_spot1_depth` | Read-depth target used when subsampling for the hotspot1 SPOT score |
 
## Output
 
Results are published under `<outdir>/<sample_id>/`. Raw FASTQs are nested under a subdirectory named by `align_id` (the SRA run accession in SRA mode); all downstream output files are prefixed with `sample_id` instead:
 
```
<outdir>/<sample_id>/
├── <align_id>/                                          # raw run FASTQs (SRA mode), named by SRA run accession, e.g. SRR13520417/
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
    ├── <align_id>/                                       # per-run stats subdirectory
    ├── <align_id>_info.json                              # SRA run metadata (from `ffq`)
    ├── <sample_id>.CollectInsertSizeMetrics.picard(.pdf)  # Picard insert size metrics
    ├── <sample_id>.MarkDuplicates.picard                  # Picard duplicate-marking metrics
    └── macs2.<sample_id>_*                                # MACS2 peak calling output (narrowPeak, summits, pileup/lambda bedgraphs, xls)
```
 
### Which step produces what

- **Alignment** (always): `<sample_id>.filtered.cram`, `<sample_id>.total_sequencing_stats.txt` and `stats/`.
- **Hotspot calling** ([hotspot3](https://github.com/vierstralab/hotspot3)): everything else. These files only appear after hotspot calling has run, either in SRA mode or via `hotspots_calling.nf`.

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
## TODO
 
- [ ] Document the stats and subsampling scripts (`spot1score`, `subsampleToFrac`, `filterAndCallHotspots`, `preprocessBams`) once cleaned up.
 
