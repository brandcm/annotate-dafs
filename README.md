# annotate-dafs

[![Tests](https://github.com/brandcm/annotate-dafs/actions/workflows/tests.yml/badge.svg)](https://github.com/brandcm/annotate-dafs/actions/workflows/tests.yml)

`annotate-dafs` is a Python command-line tool for annotating VCF files with derived allele frequencies (DAFs) using ancestral allele information and alternate allele frequencies.

Ancestral alleles can be provided either:

- in an existing VCF INFO field using `--aa-field`, or
- from an ancestral FASTA sequence using `--fasta`.

Alternate allele frequencies can be obtained from:

1. an INFO field such as `AF` (default), using the VCF ALT alleles;
2. genotypes for explicitly selected samples using `--samples` or `--sample-file`; or
3. an INFO allele-frequency field paired with a corresponding INFO field containing the alternate alleles using `--af-field` and `--alts-field`.

The package can also download and prepare the **human ancestral sequence FASTAs** provided by Ensembl for GRCh37/hg19 and GRCh38/hg38. These ancestral sequences are inferred by Ensembl Compara from primate EPO multiple-sequence alignments using Ortheus.

The downloaded FASTAs are therefore specifically the Ensembl ancestral sequences provided for the human reference assemblies. `annotate-dafs annotate` itself is not restricted to human data and can use a user-provided ancestral FASTA for other systems, provided that its sequence names and coordinates correspond to the input VCF.

## Installation

Create a Python environment and install the package:

```bash
conda create -n annotate-dafs python=3.12
conda activate annotate-dafs
pip install .
```

For development, install the package in editable mode:

```bash
pip install -e .
```

## Download ancestral sequences

### Ensembl ancestral sequence provenance

Ensembl ancestral sequences are inferred from Enredo-Pecan-Ortheus (EPO) multiple-sequence alignments. Ortheus uses the phylogenetic relationships among aligned genomes to infer ancestral sequences at internal nodes of the species tree.

The human ancestral sequence distributed for GRCh37 is based on Ensembl's 6-primate EPO alignment, which includes human, chimpanzee, gorilla, orangutan, macaque, and marmoset. The downloaded FASTA represents ancestral sequence calls mapped to the human reference coordinate system; it is not a reconstructed genome from a single ancestral species.

The `download-ancestral` command currently downloads the human ancestral sequence resources distributed by Ensembl. The `annotate` command is more general and can instead use an appropriate user-supplied ancestral FASTA for other organisms or reference systems.

For additional information, see the Ensembl documentation on ancestral sequences and multiple-genome alignments.

### Download and prepare the FASTA

Ensembl ancestral sequences can be downloaded and prepared directly with the `download-ancestral` command.

For GRCh37/hg19:

```bash
annotate-dafs download-ancestral \
    --assembly GRCh37 \
    --output-dir ~/ancestral_sequences
```

The GRCh37 ancestral sequence is retrieved from Ensembl release 75, the last Ensembl release providing the GRCh37 ancestral sequence.

For GRCh38/hg38, specify the Ensembl release:

```bash
annotate-dafs download-ancestral \
    --assembly GRCh38 \
    --release <release> \
    --output-dir ~/ancestral_sequences
```

The output directory is created automatically if it does not already exist. The command downloads the Ensembl ancestral sequence archive, combines the individual FASTA files, normalizes Ensembl sequence headers, BGZF-compresses the resulting FASTA, and creates the indexes required by `pysam`.

For GRCh37, the resulting files are:

```text
homo_sapiens_ancestor_GRCh37.fa.gz
homo_sapiens_ancestor_GRCh37.fa.gz.fai
homo_sapiens_ancestor_GRCh37.fa.gz.gzi
```

Existing output files are not overwritten by default. Use `--force` to replace them.

Ensembl ancestral FASTA headers such as:

```text
ANCESTOR_for_chromosome:GRCh37:22:1:51304566:1
```

are normalized to their sequence identifier:

```text
22
```

Supercontig identifiers are similarly retained in normalized form, for example `GL000191.1`.

## Ancestral allele sequence convention

Ancestral allele calls from Ensembl use the following conventions:

| Character | Meaning |
|---|---|
| `A`, `C`, `T`, `G` | High-confidence ancestral allele call |
| `a`, `c`, `t`, `g` | Low-confidence ancestral allele call |
| `N` | Ancestral state is unresolved |
| `-` | The extant species contains an insertion at this position |
| `.` | No coverage in the alignment |

By default, both high- and low-confidence ancestral allele calls are accepted when calculating DAF. Lowercase calls are converted to uppercase in the output.

Use `--exclude-low-confidence-ancestral` to treat lowercase ancestral alleles as unresolved:

```bash
annotate-dafs annotate \
    --fasta ~/ancestral_sequences/homo_sapiens_ancestor_GRCh37.fa.gz \
    --exclude-low-confidence-ancestral \
    --vcf input.vcf \
    --output output.vcf
```

Missing or unresolved ancestral states (`N`, `-`, or `.`) are not annotated.

FASTA-based ancestral allele retrieval currently supports SNVs only.

## Basic usage

### Using an ancestral allele INFO field

```bash
annotate-dafs annotate \
    --aa-field AA \
    --vcf input.vcf \
    --output output.vcf
```

### Using an ancestral FASTA

```bash
annotate-dafs annotate \
    --fasta ~/ancestral_sequences/homo_sapiens_ancestor_GRCh37.fa.gz \
    --vcf input.vcf \
    --output output.vcf
```

Chromosome names with and without a leading `chr` prefix are supported when matching VCF chromosomes to FASTA sequences. For example, VCF chromosome `chr22` can be matched to FASTA sequence `22`.

When a resolved ancestral allele is retrieved from a FASTA and DAF is successfully calculated, the ancestral allele used for the calculation is written to the output VCF `AA` INFO field. An existing `AA` value is replaced by the FASTA-derived allele for that variant.

The corresponding INFO definitions are:

```text
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele used to calculate DAF">
##INFO=<ID=DAF,Number=.,Type=Float,Description="Derived allele frequency">
```

## Allele frequency sources

By default, `annotate-dafs annotate` uses the `AF` INFO field and the VCF ALT alleles.

A different INFO field can be selected with `--af-field`:

```bash
annotate-dafs annotate \
    --aa-field AA \
    --af-field EAS_AF \
    --vcf input.vcf \
    --output output.vcf
```

### Calculate allele frequencies from selected samples

When `--samples` or `--sample-file` is provided, allele frequencies are calculated directly from the selected samples' genotypes rather than using an INFO AF field.

For example:

```bash
annotate-dafs annotate \
    --aa-field AA \
    --samples sample1 sample2 sample3 \
    --vcf input.vcf \
    --output output.vcf
```

Or provide one sample ID per line in a file:

```bash
annotate-dafs annotate \
    --aa-field AA \
    --sample-file samples.txt \
    --vcf input.vcf \
    --output output.vcf
```

Phased and unphased genotypes are supported, and missing alleles are excluded from the frequency calculation.

### Use alternate alleles from an INFO field

Some VCFs contain allele frequencies corresponding to a set of alternate alleles that differs from the VCF ALT field. In this case, `--alts-field` can be used together with the corresponding `--af-field`.

For example:

```bash
annotate-dafs annotate \
    --aa-field AA \
    --af-field CUSTOM_AF \
    --alts-field CUSTOM_ALTS \
    --vcf input.vcf \
    --output output.vcf
```

The alleles in `CUSTOM_ALTS` and frequencies in `CUSTOM_AF` are paired by position and must contain the same number of values. The VCF REF allele is assumed to be the reference allele for both representations.

When samples are explicitly selected, genotype allele indices are interpreted against the VCF ALT field, so `--alts-field` is not used.

A record with no value for `--af-field`, a missing value within it, or no value for `--alts-field` when that option is used, is written to the output without a DAF annotation. A record whose `--alts-field` and `--af-field` contain different numbers of values is an error.

## DAF calculation

DAF is calculated by comparing the ancestral allele with the reference and alternate alleles at each variant.

- When the **ancestral allele is the REF allele**, the ALT allele frequencies are reported as DAF.
- When the **ancestral allele is an ALT allele**, that allele is assigned a DAF of 0, while the frequencies of the other ALT alleles are retained.
- When the **ancestral allele does not match the REF or any ALT allele**, all ALT alleles are assigned a DAF of 1.0.
- When the **ancestral allele is missing or unresolved**, no DAF annotation is added.
- When the **record has no alternate allele**, it is monomorphic for the reference allele, which is reported with a DAF of 0 when it is ancestral and 1 when it is derived.

The DAF INFO field is declared with `Number=.` because the number of values follows the allele frequency source: it is the number of VCF ALT alleles by default, the number of alleles in `--alts-field` when that is used, and one for a record without an alternate allele.

For example, a variant with `REF=A`, `ALT=G`, and `AF=0.10` produces:

| Ancestral allele | DAF |
|---|---:|
| A | 0.10 |
| G | 0.00 |

For a multiallelic variant with `REF=A`, `ALT=G,T`, and `AF=0.10,0.20`:

| Ancestral allele | DAF |
|---|---|
| A | 0.10, 0.20 |
| G | 0.00, 0.20 |
| T | 0.10, 0.00 |

## Command-line options

`annotate-dafs` provides two subcommands:

```text
annotate
download-ancestral
```

### `annotate`

| Option | Description |
|---|---|
| `--aa-field` | INFO field containing the ancestral allele. Required when `--fasta` is not used. |
| `--fasta` | FASTA file containing the ancestral sequence. Required when `--aa-field` is not used. |
| `--exclude-low-confidence-ancestral` | Treat lowercase ancestral alleles from a FASTA as unresolved and skip DAF annotation. |
| `--vcf` | Input VCF file to annotate. Required. |
| `--output` | Output VCF file. Required and must differ from the input VCF. |
| `--af-field` | INFO field containing alternate allele frequencies. Default: `AF`. |
| `--alts-field` | INFO field containing alternate alleles corresponding to `--af-field`. |
| `--samples` | Sample IDs whose genotypes should be used to calculate allele frequencies. |
| `--sample-file` | File containing one sample ID per line whose genotypes should be used to calculate allele frequencies. |
| `--daf-field` | INFO field name for the derived allele frequency annotation. Default: `DAF`. |
| `--daf-field-description` | Description for the DAF INFO field. Default: `Derived allele frequency`. |
| `--log-level` | Logging level: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. Default: `INFO`. |

### `download-ancestral`

| Option | Description |
|---|---|
| `--assembly` | Genome assembly: `GRCh37` or `GRCh38`. Required. |
| `--release` | Ensembl release. Required for GRCh38 and not permitted for GRCh37. |
| `--output-dir` | Directory in which the prepared ancestral FASTA and indexes are written. Required. |
| `--force` | Overwrite existing ancestral FASTA output files. |
| `--log-level` | Logging level: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. Default: `INFO`. |

Use:

```bash
annotate-dafs --help
annotate-dafs annotate --help
annotate-dafs download-ancestral --help
```

to display the complete command-line interface.

## Requirements and testing

`annotate-dafs` requires Python 3.10 or later and uses `pysam` and `vcfpy`.

Tests can be run from the repository root with:

```bash
pytest -q
```

## References

- [Ensembl ancestral sequences](https://grch37.ensembl.org/info/genome/compara/ancestral_sequences.html)
- [Ensembl release 75](https://www.ebi.ac.uk/about/news/updates-from-data-resources/ensembl-75/)
- [Ensembl multiple genome alignments](https://www.ensembl.org/info/genome/compara/multiple_genome_alignments.html)

## License

`annotate-dafs` is distributed under the MIT License. See [LICENSE](LICENSE) for details.
