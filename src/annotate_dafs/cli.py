# src/annotate_dafs/cli.py

"""Command-line interface for annotate-dafs."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pysam
import vcfpy

from .ancestral import (
    download_ancestral_fasta,
    retrieve_ancestral_allele_from_fasta,
    retrieve_ancestral_allele_from_info,
)
from .frequencies import (
    calculate_daf,
    get_alt_alleles_and_frequencies,
)
from .io import (
    add_info_field,
    load_samples_from_file,
    validate_samples,
)


LOGGER = logging.getLogger(__name__)


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments shared by annotate-dafs commands."""
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
        help="Logging level [default: INFO].",
    )


def add_annotate_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the annotate command."""
    ancestral_group = parser.add_mutually_exclusive_group(required=True)

    ancestral_group.add_argument(
        "--aa-field",
        help="INFO field containing the ancestral allele.",
    )

    ancestral_group.add_argument(
        "--fasta",
        type=Path,
        help="FASTA file containing the ancestral sequence.",
    )

    parser.add_argument(
        "--exclude-low-confidence-ancestral",
        action="store_true",
        help="Treat lowercase ancestral alleles as unresolved and skip DAF annotation.",
    )

    sample_group = parser.add_mutually_exclusive_group()

    sample_group.add_argument(
        "--samples",
        nargs="+",
        help="Sample IDs to use when calculating allele frequencies.",
    )

    sample_group.add_argument(
        "--sample-file",
        type=Path,
        help="File containing one sample ID per line.",
    )

    parser.add_argument(
        "--vcf",
        type=Path,
        required=True,
        help="Input VCF file.",
    )

    parser.add_argument(
        "--daf-field",
        default="DAF",
        help="INFO field name for DAF [default: DAF].",
    )

    parser.add_argument(
        "--daf-field-description",
        default="Derived allele frequency",
        help="Description for the DAF INFO field.",
    )

    parser.add_argument(
        "--af-field",
        default="AF",
        help="INFO field containing alternate allele frequencies.",
    )

    parser.add_argument(
        "--alts-field",
        help="Optional INFO field containing alternate alleles.",
    )

    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output VCF file.",
    )

    add_common_arguments(parser)


def add_download_ancestral_arguments(
    parser: argparse.ArgumentParser,
) -> None:
    """Add arguments for the download-ancestral command."""
    parser.add_argument(
        "--assembly",
        choices=("GRCh37", "GRCh38"),
        required=True,
        help="Human reference genome assembly.",
    )

    parser.add_argument(
        "--release",
        type=int,
        help="Ensembl release for GRCh38.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for the downloaded ancestral FASTA.",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files.",
    )

    add_common_arguments(parser)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Annotate VCFs with derived allele frequencies and manage "
            "ancestral sequence data."
        )
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Annotate a VCF with derived allele frequencies.",
        description=(
            "Annotate a VCF with derived allele frequency (DAF) using "
            "ancestral allele information."
        ),
    )
    add_annotate_arguments(annotate_parser)

    download_parser = subparsers.add_parser(
        "download-ancestral",
        help="Download and prepare an Ensembl ancestral sequence FASTA.",
        description=(
            "Download, combine, BGZF-compress, and index an Ensembl "
            "ancestral sequence FASTA."
        ),
    )
    add_download_ancestral_arguments(download_parser)

    args = parser.parse_args(argv)

    if args.command == "annotate":
        if args.output.resolve() == args.vcf.resolve():
            parser.error("--output must be different from --vcf.")

    if args.command == "download-ancestral":
        if args.assembly == "GRCh37" and args.release is not None:
            parser.error(
                "--release cannot be specified for GRCh37; "
                "the GRCh37 ancestral sequence is retrieved from "
                "Ensembl release 75."
            )

        if args.assembly == "GRCh38" and args.release is None:
            parser.error("--release is required for GRCh38.")

    return args


def configure_logging(log_level: str) -> None:
    """Configure application logging."""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format="%(levelname)s: %(message)s",
    )


def process_variant(
    variant: vcfpy.Record,
    aa_field: str | None,
    fasta: pysam.FastaFile | None,
    exclude_low_confidence_ancestral: bool,
    daf_field: str,
    output_aa_field: str,
    af_field: str,
    alts_field: str | None,
    samples: set[str] | None,
) -> None:
    """Calculate and add DAF and the ancestral allele to a single VCF variant."""
    if fasta is not None:
        ancestral_allele = retrieve_ancestral_allele_from_fasta(
            fasta=fasta,
            chromosome=variant.CHROM,
            position=variant.POS,
            reference=variant.REF,
            exclude_low_confidence=exclude_low_confidence_ancestral,
        )
    else:
        if aa_field is None:
            raise ValueError("An ancestral allele source was not provided.")

        ancestral_allele = retrieve_ancestral_allele_from_info(
            variant,
            aa_field,
            exclude_low_confidence=exclude_low_confidence_ancestral,
        )

    if ancestral_allele is None:
        return

    allele_frequencies = get_alt_alleles_and_frequencies(
        variant=variant,
        af_field=af_field,
        alts_field=alts_field,
        samples=samples,
    )

    if allele_frequencies is None:
        return

    alternate_alleles, alternate_frequencies = allele_frequencies

    daf = calculate_daf(
        ancestral_allele=ancestral_allele,
        reference_allele=variant.REF,
        alternate_alleles=alternate_alleles,
        alternate_frequencies=alternate_frequencies,
    )

    if daf:
        variant.INFO[output_aa_field] = ancestral_allele
        variant.INFO[daf_field] = daf


def run_annotate(args: argparse.Namespace) -> int:
    """Run the DAF annotation workflow."""
    samples: set[str] | None = None

    if args.samples is not None:
        samples = set(args.samples)
    elif args.sample_file is not None:
        samples = set(load_samples_from_file(args.sample_file))

    fasta: pysam.FastaFile | None = None

    if args.fasta is not None:
        fasta = pysam.FastaFile(str(args.fasta))

    output_aa_field = args.aa_field if args.aa_field is not None else "AA"

    try:
        reader = vcfpy.Reader.from_path(args.vcf)

        if samples is not None:
            validate_samples(
                samples=samples,
                available_samples=reader.header.samples.names,
            )

        add_info_field(
            header=reader.header,
            field_id=output_aa_field,
            description="Ancestral allele used to calculate DAF",
            number="1",
            field_type="String",
        )

        # the number of DAF values follows the allele frequency source, which
        # for --alts-field and for monomorphic records is not the ALT count
        add_info_field(
            header=reader.header,
            field_id=args.daf_field,
            description=args.daf_field_description,
            number=".",
        )

        writer = vcfpy.Writer.from_path(
            args.output,
            reader.header,
        )

        processed = 0
        annotated = 0

        try:
            for variant in reader:
                process_variant(
                    variant=variant,
                    aa_field=args.aa_field,
                    fasta=fasta,
                    exclude_low_confidence_ancestral=(
                        args.exclude_low_confidence_ancestral
                    ),
                    daf_field=args.daf_field,
                    output_aa_field=output_aa_field,
                    af_field=args.af_field,
                    alts_field=args.alts_field,
                    samples=samples,
                )

                if args.daf_field in variant.INFO:
                    annotated += 1

                writer.write_record(variant)
                processed += 1

        finally:
            writer.close()

    finally:
        if fasta is not None:
            fasta.close()

    LOGGER.info(
        "Processed %d variants; annotated %d with DAF.",
        processed,
        annotated,
    )

    return 0


def run_download_ancestral(args: argparse.Namespace) -> int:
    """Run the ancestral FASTA download workflow."""
    download_ancestral_fasta(
        assembly=args.assembly,
        output_directory=args.output_dir,
        release=args.release,
        force=args.force,
    )

    return 0


def main(argv: list[str] | None = None) -> int:
    """Run annotate-dafs."""
    args = parse_args(argv)
    configure_logging(args.log_level)

    if args.command == "annotate":
        return run_annotate(args)

    if args.command == "download-ancestral":
        return run_download_ancestral(args)

    raise ValueError(f"Unknown command: {args.command}")
