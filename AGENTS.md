# AGENTS.md — nf-aligning

Guidance for AI coding agents (Claude Code and others) working in this repository. Human-facing docs are in `README.md`; this file is the operational summary an agent needs to run, debug, or modify the pipeline safely.

## What this repo is

A Nextflow (DSL2) pipeline that takes DNase-seq / ATAC-seq short reads to aligned, filtered CRAMs and [hotspot3](https://github.com/vierstralab/hotspot3) peak calls. It is the **first stage** of the processing used for the paper [TODO: add link]; its per-sample outputs (peaks, normalized density, fit stats, p-value parquet, CRAM) are the inputs to the second stage, [`nf-index`](https://github.com/vierstralab/nf-index), which builds the consensus DHS index and the DHS × sample matrices.

```
FASTQ (local or SRA) ─► fastp trim ─► bwa aln/sampe|samse ─► filter_reads.py (flag) + sort
   ─► merge runs per sample_id ─► Picard MarkDuplicatesWithMateCigar
   ─► samtools view -F 516 [nuclear chroms] ─► <sample_id>.filtered.cram
   ─► hotspot3 ─► peaks / hotspots / density / p-values per FDR
```

## Entry points

Run everything from the repo root (or point `nextflow run` at it). Always pass a profile; add `-resume` when re-running.

| Command | Input (`--samples_file`, TSV) | What it does |
|---|---|---|
| `nextflow run main.nf -profile <p>` | `sample_id, align_id, reads1, reads2, type[, adapterP7, adapterP5]` | Trim → align → stats → CRAM → hotspot3. |
| `nextflow run main.nf -profile <p> -entry alignFromSRA` | `sample_id, align_id` (align_id = SRR accession) | Download from SRA → trim → align → CRAM → hotspot3. |
| `nextflow run aligning.nf -profile <p>` | `sample_id, reads1, reads2, type` (already trimmed) | Align → stats → CRAM. For debugging or finishing a stopped run. |
| `nextflow run hotspots_calling.nf -profile <p>` | `sample_id, cram_file[, cram_index]` | hotspot3 only, on existing CRAM/BAM. Re-calls peaks without re-aligning. |
| `nextflow run stats.nf -profile <p>` | `sample_id, cram_file, cram_index` | `collect_basic_stats` (filtered_aligned and duplicate counts) + `percent_dup`. |
| `nextflow run stats.nf -profile <p> -entry nuclearStats` | same | `filter_nuclear` + `bamcounts.py` stats. |
| `nextflow run stats.nf -profile <p> -entry preseq` | `sample_id, cram_file, cram_index, read_type` | `preseq c_curve` on nuclear reads. |
| `nextflow run subsampling.nf -profile <p> -entry spot1score` | `sample_id, cram_file` | R1 only → `filter_nuclear` → subsample to `subsampling_spot1_depth` → hotspot1 SPOT score. |
| `nextflow run subsampling.nf -profile <p> -entry subsampleToFrac` | `sample_id, cram_file, cram_index, frac` | Subsample pairs to a fraction → re-mark dups → hotspot3 + stats + `<sample_id>.marked.cram`. |

Some `stats.nf` and `subsampling.nf` entries (`nuclearStats`, `subsampleToFrac`) reuse output file names from the main run, so run them with a separate `--outdir`.

Code below a `// defunc` comment (and files starting with one) is unmaintained: don't use it, fix it, or document it unless asked.

## Configuration

- `nextflow.config` includes `params.config` into `params {}`. `params.config` is a symlink to `params.human.config` (GRCh38). For mouse (mm10), `ln -sf params.mouse.config params.config`. `-c params.mouse.config` does not work (bare assignments outside `params {}` don't become params); override single values with `--<param>`. `params.mouse.config` has `samples_file = ""`.
- All reference paths in the params files are absolute paths on the Altius `/net/seq/...` filesystem. On any other system every path must be replaced.
- `params.container` defaults to `${projectDir}/aligning-pipeline.sif`. Build files are in `container/`: `aligning-pipeline.def`, `aligning-pipeline.def` and `environment.yml` (one conda env, `aligning`, for alignment, QC and hotspot3; hotspot3 is `pip install`ed on top during the build). Build from inside that folder: `cd container && singularity build --fakeroot ../aligning-pipeline.sif aligning-pipeline.def` (~15 min). Every process runs in this image; hotspot1 is at `/hotspot1`.
- `helpers.nf:get_container()` adds `--bind` options for the directory holding a reference file (and its symlink target). It is evaluated when a `.nf` file loads, so every file it is called on must exist. If a process can't see a reference file inside the container, check this first.
- Profiles: `Altius` (SLURM, apptainer via `module load apptainer`, singularity on). `new_cluster` and `old_cluster` only set queues and are meant to be combined, e.g. `-profile Altius,new_cluster`. Outside Altius, add a profile instead of editing `Altius`.
- Key params: `outdir` (relative path expected), `samples_file`, `threads`, `genome_fasta_file` (BWA index prefix is the FASTA path without its extension, `.fai` also needed), `nuclear_chroms` (file, one name per line), `nuclear_chrom_sizes`, `mappable_bases`, `fdrs` (comma-separated, default `0.05,0.01,0.001`), `save_cram_mode` (`all` (default) | `nuclear`), `save_debug`. `chrominfo`, `spot1_mappable_file`, `readlength` and `subsampling_spot1_depth` (default 5M, keep it) are used only by `subsampling.nf -entry spot1score`, but `chrominfo` and `spot1_mappable_file` go through `get_container()` when `subsampling.nf` loads, so they must exist for every entry in it.

## Processing details (easy to get wrong)

- `bin/filter_reads.py` **does not remove reads**. It sets flags: 0x200 (QC fail) for MAPQ < 10, NM > 2, unpaired or improper pairs, or insert > 750; 0x1000 for non-nuclear. Removal happens later with `samtools view -F 516` (unmapped + QC-fail): `filter_passing_reads` (all chromosomes) when `save_cram_mode = "all"` (default), `filter_nuclear` (only `nuclear_chroms`) when `"nuclear"`. The `stats.nf` and `subsampling.nf` entries always use `filter_nuclear`.
- Duplicates are **marked, not removed** (MarkDuplicatesWithMateCigar, `MINIMUM_DISTANCE=300`). The CRAM keeps them with flag 1024, and hotspot3 counts them. Above 50% duplication the hotspot3 background fit breaks down, so such libraries are downsampled below 50% with `subsampling.nf -entry subsampleToFrac`, with `frac` estimated from the `stats.nf -entry preseq` complexity curve (README, "Libraries with a high duplication rate").
- Runs are merged by `sample_id` after alignment. `align_id` must be unique; rows are de-duplicated on it.
- `insert_size` runs only for paired-end data, on the nuclear-filtered BAM.
- hotspot3 runs on the filtered CRAM with `--save_density --reference <fasta> --debug`; `debug/` is emptied afterwards unless `save_debug = true`. `--debug` also makes hotspot3 write `<sample_id>.fit_stats.bed.gz.iter<N>` into the sample folder for samples that needed refits; these are always published.
- For hotspot3 options, output formats and fit-stats flags (`refit_with_constraint`, `success_fit`, `max_bg_reached`), read the [hotspot3 README and AGENTS.md](https://github.com/vierstralab/hotspot3) rather than guessing.

## Outputs (per sample, under `<outdir>/<sample_id>/`)

- `<sample_id>.filtered.cram` + `.crai`: the main alignment deliverable.
- `<sample_id>.total_sequencing_stats.txt` (from `bamcounts.py`, on the duplicate-marked BAM before filtering).
- `stats/`: `<sample_id>.MarkDuplicates.picard`, `<sample_id>.CollectInsertSizeMetrics.picard(.pdf)`, and per run `<align_id>/fastp.{json,html}` plus `<align_id>/<align_id>_info.json` in SRA mode.
- hotspot3 outputs: `fdr<X>/<sample_id>.{peaks,hotspots}.fdr<X>.{bed.gz,bed.gz.tbi,bb}`, `<sample_id>.normalized_density.bw`, `<sample_id>.background.bw`, `<sample_id>.fit_stats.bed.gz`, `<sample_id>.pvals.parquet`, `<sample_id>.cutcounts.bed.gz`, `<sample_id>.total_cutcounts`, `<sample_id>.peak_calling.log`.

The nf-index `samples_file` columns map to these outputs: `cram_file` → `.filtered.cram`; `peaks_file_0.001fdr` / `peaks_file_0.01fdr` → `fdr0.001/…peaks…` / `fdr0.01/…peaks…`; `normalized_density_bw`; `hotspot3_fit_stats_file` → `fit_stats`; `hotspot3_pvals_parquet` → `pvals.parquet`.

## Gotchas

- **Silent failures.** The `Altius` profile retries exit codes 134–143 up to 3 times and otherwise uses `errorStrategy 'ignore'`. `merge_bam` and SRA download are also `ignore`. A finished run can therefore be missing samples. After every run, compare the sample IDs in `samples_file` with the directories that contain a `.filtered.cram` (and `fdr*/` if hotspots were expected), and check `.nextflow.log` for `FAILED`/`IGNORED` tasks.

- Picard's `READ_NAME_REGEX` assumes Illumina-style read names. SRA-renamed reads won't parse for optical duplicates (a warning, not an error).

## Working rules for agents

- Do not launch full runs without being asked. They submit SLURM jobs at scale (`queueSize = 1200`). For testing, use a 1–2 row `samples_file`, a scratch `outdir`, and `-resume`.
- Never delete `work/` or `outdir` contents unless told to. `-resume` depends on `work/`.
- Before changing a process, grep for which `.nf` files `include` it. `aligning.nf` processes are reused by `stats.nf` and `subsampling.nf`.
- Keep `README.md` in sync when you change entry points, samples-file columns, params, or output names.
- Syntax check without running: `nextflow run <file>.nf -preview` or `nextflow config` to resolve params.
