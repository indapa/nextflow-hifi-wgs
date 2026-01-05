# GitHub Copilot Instructions for nextflow-hifi-wgs

## Project Overview

This is a Nextflow DSL2 pipeline for secondary analysis of PacBio HiFi sequencing reads. The pipeline performs read alignment, variant calling, phasing, and annotation.

## Technology Stack

- **Workflow Engine**: Nextflow DSL2
- **Containerization**: Docker containers (with Wave and Fusion support)
- **Primary Tools**: pbmm2 (alignment), DeepVariant (variant calling), HiPhase/WhatsHap (phasing), Ensembl VEP (annotation)

## Code Structure and Conventions

### Pipeline Organization

- **Main workflow file**: `main.nf` contains workflow definitions and orchestration
- **Modules directory**: `modules/` contains process definitions organized by tool:
  - `modules/pbtools/` - pbmm2 alignment and HiPhase phasing processes
  - `modules/deepvariant/` - DeepVariant and DeepTrio variant calling processes
  - `modules/samtools/` - BAM file statistics processes
  - `modules/ensemblvep/` - Ensembl VEP annotation processes
  - `modules/whatshap/` - WhatsHap trio phasing processes
- **Configuration**: `nextflow.config` contains process execution settings and parameters

### Nextflow DSL2 Conventions

1. **Process definitions**: Each process should be defined in its module file with clear inputs, outputs, and emit labels
2. **Channel handling**: Use tuple channels with explicit element naming, e.g., `tuple val(sample_id), path(bam), path(bai)`
3. **Process directives**:
   - Use `tag` for sample identification in logs: `tag "$sample_id"`
   - Use `publishDir` for output management with mode 'copy'
   - Use `label` for resource management: `label 'high_memory'`
   - Specify `container` for each process with version-pinned images

4. **Workflow entry points**: Define multiple workflows with descriptive names:
   - Default workflow runs `ALIGN_DEEP_VARIANT_HIPHASE_VEP_SYT1`
   - Use `-entry <workflow_name>` to run alternative workflows

### Coding Style

1. **Indentation**: Use 4 spaces for indentation in Nextflow scripts
2. **Comments**: Add clear comments explaining complex channel transformations and pipeline steps
3. **Naming conventions**:
   - Use snake_case for process names: `pbmm2_align_syt1_region`
   - Use descriptive channel names with `_ch` suffix: `input_bams_ch`, `hiphase_input_ch`
   - Use clear variable names for tuples and maps

4. **Process scripts**: Use bash script blocks with clear step numbering and explanations
5. **Error handling**: Add safety checks at workflow start for required parameters

### Channel Transformations

Common patterns used in this pipeline:
- `.join()` - Combine channels by key
- `.map{}` - Transform channel elements
- `.splitCsv(header: true)` - Parse CSV samplesheets
- `.groupTuple(by: 0)` - Group elements by first tuple element

### Docker Container Usage

- All processes use containerized tools with version-pinned images
- Enable Docker, Wave, and Fusion in `nextflow.config`
- Container images are primarily from quay.io and community repositories

## Input/Output Conventions

### Samplesheet Format

1. **Alignment workflow** (`samplesheet` parameter):
   ```csv
   sample_id,bam_file
   sample1,/path/to/sample1.bam
   ```
   BAM files should contain unaligned HiFi reads.

2. **HiPhase workflow** (`hiphase_samplesheet` parameter):
   ```csv
   sample_id,vcf,vcf_tbi,bam,bai
   sample1,/path/to/sample1.vcf,/path/to/sample1.vcf.tbi,/path/to/sample1.bam,/path/to/sample1.bai
   ```
   Requires aligned BAM files and variant calls from DeepVariant.

3. **DeepTrio workflow** (`trio_samplesheet` parameter):
   ```csv
   family_id,sample_id,bam,bai,role
   family1,child1,/path/to/child.bam,/path/to/child.bai,child
   family1,father1,/path/to/father.bam,/path/to/father.bai,father
   family1,mother1,/path/to/mother.bam,/path/to/mother.bai,mother
   ```

### Output Files

Outputs are published to directories specified in `nextflow.config`:
- `aligned_output_dir` - Aligned BAM files
- `deepvariant_output_dir` - Variant call files (VCF)
- `hiphase_output_dir` - Phased VCF files
- `bam_stats_output_dir` - BAM statistics

## Running and Testing

### Basic Execution

```bash
nextflow run main.nf -c nextflow.config
```

### Alternative Workflows

```bash
# For HiPhase and VEP only
nextflow run main.nf -entry HIPHASE_VEP_ONLY --hiphase_samplesheet path/to/samplesheet.csv

# For DeepTrio variant calling
nextflow run main.nf -entry RUN_DEEPTRIO --trio_samplesheet path/to/trio_samplesheet.csv
```

### Resume Capabilities

Nextflow supports resume functionality with `-resume` flag to restart from last successful step.

## When Making Changes

1. **Adding new processes**: Create in appropriate module directory with clear inputs/outputs
2. **Modifying workflows**: Ensure channel compatibility between connected processes
3. **Container updates**: Pin version numbers to ensure reproducibility
4. **Parameter changes**: Update `nextflow.config` with sensible defaults
5. **Documentation**: Update README.md for user-facing changes

## Common Patterns to Follow

1. **Safety checks**: Validate required parameters at workflow start
2. **Channel debugging**: Use `.view()` during development (commented out in production)
3. **Stub sections**: Provide stub implementations for testing workflow logic
4. **Resource management**: Use labels in `nextflow.config` for different resource requirements

## Key Files to Review

- `main.nf` - Workflow orchestration and entry points
- `nextflow.config` - Configuration and parameters
- `README.md` - User documentation
- `modules/*/main.nf` - Process implementations
