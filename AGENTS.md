# AGENTS.md — nf-aligning
 
Guidance for AI coding agents (Claude Code and others) working in this repository. Human-facing docs are in `README.md`; this file is the operational summary an agent needs to run, debug, or modify the pipeline safely.
 
## What this repo is
 
A Nextflow (DSL2) pipeline that takes DNase-seq / ATAC-seq short reads to aligned, filtered CRAMs and hotspot3 peak calls. It is the **first stage** of the processing used for the paper; its per-sample outputs (peaks, normalized density, fit stats, p-value parquet, CRAM) are the inputs to the second stage, [`nf-index`](https://github.com/vierstralab/nf-index), which builds the consensus DHS index and the DHS × sample matrices.
 
```
FASTQ (local or SRA) ─► fastp trim ─► bwa aln/sampe|samse ─► filter_reads.py (flag) + sort
   ─► merge runs per sample_id ─► Picard MarkDuplicatesWithMateCigar
   ─► samtools view -F 516 [nuclear chroms] ─► <sample_id>.filtered.cram
   ─► hotspot3 ─► peaks / hotspots / density / p-values per FDR
```
 
## Entry points (what actually runs)
 
Run everything from the repo root (or point `nextflow run` at it). Always pass a profile; add `-resume` when re-running.
 
| Command | Input (`--samples_file`, TSV) | What it does |
|---|---|---|
| `nextflow run main.nf -profile <p>` | `sample_id, align_id, reads1, reads2, type[, adapterP7, adapterP5]` | Trim → align → stats → CRAM. **Does not call hotspots**: `call_hotspots` is commented out at the end of the default workflow in `main.nf`. |
| `nextflow run main.nf -profile <p> -entry alignFromSRA` | `sample_id, align_id` (align_id = SRR accession) | Download from SRA → trim → align → CRAM → **hotspot3**. |
| `nextflow run hotspots_calling.nf -profile <p>` | `sample_id, cram_file[, cram_index]` | hotspot3 only, on existing CRAM/BAM. Use this to add peak calls after the default (local FASTQ) workflow. |
| `nextflow run stats.nf -profile <p>` | `sample_id, cram_file, cram_index` | `collect_basic_stats` (filtered_aligned and duplicate counts) + `percent_dup`. |
| `nextflow run stats.nf -entry nuclearStats` | same | Nuclear filter + `bamcounts.py` stats. |
| `nextflow run subsampling.nf -entry spot1score` | `sample_id, cram_file` | R1 only → nuclear → subsample to `subsampling_spot1_depth` → HotSpot1 SPOT score. |
| `nextflow run subsampling.nf -entry subsampleToFrac` | `sample_id, cram_file, cram_index, frac` | Subsample pairs to a fraction → re-mark dups → hotspot3. Needs `--conda` (see Gotchas). |
| `nextflow run subsampling.nf -entry filterAndCallHotspots` | `sample_id, bam_file, bam_index` | Nuclear filter → hotspot3. |
 
Treat these as legacy and avoid them unless asked: `encode_download.nf` (hard-coded metadata path); the `fromCutcounts` / `tmp_debug` entries in `hotspots_calling.nf` (marked defunc, hard-coded `prev_run_dir`). `hotspots_calling.hotspot3.nf` is a byte-identical copy of `hotspots_calling.nf`. `main.nf` imports `hotspots_calling.nf`, so edit that one.
 
## Configuration
 
- `nextflow.config` includes `params.config` into `params {}`. **`params.config` is the human GRCh38 config** (identical to `params.human.config`). For mouse (mm10), copy `params.mouse.config` over `params.config` or pass `-c`/`--param` overrides. Note that `params.mouse.config` has `samples_file = ""`.
- All reference paths in the params files are absolute paths on the Altius `/net/seq/...` filesystem. On any other system every path must be replaced.
- `params.container` defaults to `${projectDir}/aligning-pipeline.sif`. Build it with `singularity build --fakeroot aligning-pipeline.sif aligning-pipeline.def` (~15 min). Almost every process runs in this image, which contains bwa, samtools, picard, fastp, sra-tools, ffq, macs2, HotSpot1, HotSpot2, kentUtils and hotspot3.
- `helpers.nf:get_container()` adds `--bind` options for the directories holding the reference FASTA, chrom sizes and nuclear chrom list (including symlink targets). If a process can't see a reference file inside the container, check this first.
- Profiles: `Altius` (SLURM, apptainer via `module load apptainer`, singularity on). `new_cluster` and `old_cluster` only set queues and are meant to be combined, e.g. `-profile Altius,new_cluster`. Outside Altius, add a profile instead of editing `Altius`.
Key params: `outdir` (relative path expected), `samples_file`, `threads`, `genome_fasta_file` (the BWA index must sit next to it with the same basename), `nuclear_chroms`, `chrom_sizes`, `nuclear_chrom_sizes`, `mappable_bases`, `fdrs` (comma-separated, default `0.05,0.01,0.001`), `save_cram_mode` (`nuclear` | `all`), `do_macs`, `save_debug`.
 
## Processing details agents get wrong
 
- `bin/filter_reads.py` **does not remove reads**. It sets flags: 0x200 (QC fail) for MAPQ < 10, NM > 2, unpaired or improper pairs, or insert > 750; 0x1000 for non-nuclear. Removal happens later in `filter_nuclear` with `samtools view -F 516` (unmapped + QC-fail), restricted to `nuclear_chroms` when `save_cram_mode == "nuclear"`.
- Duplicates are **marked, not removed** (MarkDuplicatesWithMateCigar, `MINIMUM_DISTANCE=300`). The CRAM keeps them with flag 1024.
- Runs are merged by `sample_id` after alignment. `align_id` must be unique; rows are de-duplicated on it.
- `insert_size` runs only for paired-end data. `macs2` runs only if `do_macs`. Both use the nuclear-filtered BAM.
- hotspot3 runs on the filtered CRAM with `--save_density --reference <fasta>`. `debug/` is emptied unless `save_debug = true`.
## Outputs (per sample, under `<outdir>/<sample_id>/`)
 
- `<sample_id>.filtered.cram` + `.crai`: the main alignment deliverable.
- `<sample_id>.total_sequencing_stats.txt` (from `bamcounts.py`).
- `stats/`: `<sample_id>.MarkDuplicates.picard`, `<sample_id>.CollectInsertSizeMetrics.picard(.pdf)`, `macs2.<sample_id>*`, `<align_id>/fastp.{json,html}`, and `<align_id>_info.json` in SRA mode.
- hotspot3 outputs (SRA mode or `hotspots_calling.nf`): `fdr<X>/<sample_id>.peaks.fdr<X>.bed.gz`, `fdr<X>/<sample_id>.hotspots.fdr<X>.bed.gz`, `<sample_id>.normalized_density.bw`, `<sample_id>.fit_stats.*.gz`, `<sample_id>.pvals.parquet`, `<sample_id>.cutcounts.bed.gz`, `<sample_id>.total_cutcounts`, `<sample_id>.peak_calling.log`.
The nf-index `samples_file` columns map to these outputs: `cram_file` → `.filtered.cram`; `peaks_file_0.001fdr` / `peaks_file_0.01fdr` → `fdr0.001/…peaks…` / `fdr0.01/…peaks…`; `normalized_density_bw`; `hotspot3_fit_stats_file` → `fit_stats`; `hotspot3_pvals_parquet` → `pvals.parquet`.
 
## Gotchas
 
- **Silent failures.** The `Altius` profile retries exit codes 134–143 up to 3 times and otherwise uses `errorStrategy 'ignore'`. `merge_bam` and SRA download are also `ignore`. A finished run can therefore be missing samples. After every run, compare the sample IDs in `samples_file` with the directories that contain a `.filtered.cram` (and `fdr*/` if hotspots were expected), and check `.nextflow.log` for `FAILED`/`IGNORED` tasks.
- **No hotspots in the default entry.** Run `hotspots_calling.nf` on the CRAMs afterwards, or uncomment `call_hotspots` in `main.nf`.
- `subsampling.nf`'s `subsample_with_pairs_frac` uses `conda params.conda`, which is not defined in `params.config`. Pass `--conda <env>` or switch it to `container`.
- `stats.nf:run_preseq` hard-codes a personal conda path (`/home/sabramov/miniconda3/envs/super-index`).
- `mark_duplicates` in `aligning.nf` uses `${ag_id}` in `tag`/`publishDir` while its input is named `sample_id`. Check that `stats/<sample_id>.MarkDuplicates.picard` is actually published before relying on it.
- Picard's `READ_NAME_REGEX` assumes Illumina-style read names. SRA-renamed reads won't parse for optical duplicates (a warning, not an error).
## Working rules for agents
 
- Do not launch full runs without being asked. They submit SLURM jobs at scale (`queueSize = 1200`). For testing, use a 1–2 row `samples_file`, a scratch `outdir`, and `-resume`.
- Never delete `work/` or `outdir` contents unless told to. `-resume` depends on `work/`.
- Before changing a process, grep for which `.nf` files `include` it. `aligning.nf` processes are reused by `stats.nf`, `subsampling.nf` and `encode_download.nf`.
- Keep `README.md` in sync when you change entry points, samples-file columns, params, or output names.
- Syntax check without running: `nextflow run <file>.nf -preview` or `nextflow config` to resolve params.
 
