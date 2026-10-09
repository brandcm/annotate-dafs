# src/annotate_dafs/frequencies.py

"""Allele-frequency and derived-allele-frequency calculations."""

from __future__ import annotations

import re
from collections import Counter
from typing import Sequence

import vcfpy

from .io import parse_alt_alleles


def calculate_allele_frequencies(
    genotype_calls: Sequence[Sequence[int | None]],
    number_of_alts: int,
) -> list[float] | None:
    """Calculate alternate allele frequencies from genotype allele calls."""
    allele_counts = Counter()
    total_alleles = 0

    for genotype in genotype_calls:
        for allele in genotype:
            if allele is None:
                continue

            if allele < 0 or allele > number_of_alts:
                raise ValueError(
                    f"Genotype contains invalid allele index {allele}; "
                    f"expected an index from 0 to {number_of_alts}."
                )

            allele_counts[allele] += 1
            total_alleles += 1

    if total_alleles == 0:
        return None

    return [
        allele_counts[allele_index] / total_alleles
        for allele_index in range(1, number_of_alts + 1)
    ]


def get_genotype_calls(
    variant: vcfpy.Record,
    samples: set[str],
) -> list[list[int | None]]:
    """Retrieve genotype allele calls for the requested samples."""
    parsed_calls = []

    for call in variant.calls:
        if call.sample not in samples:
            continue

        genotype = call.data.get("GT", "")
        alleles: list[int | None] = []

        for allele in re.split(r"[|/]", genotype):
            if allele in {"", "."}:
                alleles.append(None)
            else:
                alleles.append(int(allele))

        parsed_calls.append(alleles)

    return parsed_calls


def parse_af_values(
    variant: vcfpy.Record,
    af_field: str,
) -> list[float] | None:
    """Retrieve and validate alternate allele frequencies from INFO."""
    af_values = variant.INFO.get(af_field)

    if af_values is None:
        return None

    if not isinstance(af_values, list):
        af_values = [af_values]

    frequencies = []

    for value in af_values:
        if value is None or value == ".":
            return None

        frequency = float(value)

        if not 0 <= frequency <= 1:
            raise ValueError(
                f"Invalid allele frequency {frequency} at "
                f"{variant.CHROM}:{variant.POS}; expected a value in [0, 1]."
            )

        frequencies.append(frequency)

    return frequencies


def get_alt_alleles_and_frequencies(
    variant: vcfpy.Record,
    af_field: str,
    samples: set[str] | None,
    alts_field: str | None = None,
) -> tuple[list[str], list[float]] | None:
    """Get alternate alleles and their frequencies."""

    if samples is not None:
        genotype_calls = get_genotype_calls(variant, samples)

        frequencies = calculate_allele_frequencies(
            genotype_calls,
            number_of_alts=len(variant.ALT),
        )

        if frequencies is None:
            return None

        alleles = [alt.value.upper() for alt in variant.ALT]

        return alleles, frequencies

    # a monomorphic record has no alternate alleles, and so no frequencies to read
    if alts_field is None and not variant.ALT:
        return [], []

    frequencies = parse_af_values(variant, af_field)

    if frequencies is None:
        return None

    if alts_field is None:
        alleles = [alt.value.upper() for alt in variant.ALT]
    elif variant.INFO.get(alts_field) is None:
        return None
    else:
        alleles = parse_alt_alleles(variant, alts_field)

    if len(alleles) != len(frequencies):
        raise ValueError(
            f"Number of alternate alleles ({len(alleles)}) does not match "
            f"number of allele frequencies ({len(frequencies)}) at "
            f"{variant.CHROM}:{variant.POS}."
        )

    return alleles, frequencies


def calculate_daf(
    ancestral_allele: str,
    reference_allele: str,
    alternate_alleles: Sequence[str],
    alternate_frequencies: Sequence[float] | None,
) -> list[float] | None:
    """Calculate derived allele frequency for each alternate allele."""
    if alternate_frequencies is None:
        return None

    if len(alternate_frequencies) != len(alternate_alleles):
        raise ValueError(
            "Number of allele frequencies does not match number of "
            "alternate alleles."
        )

    ancestral_allele = ancestral_allele.upper()
    reference_allele = reference_allele.upper()
    alternate_alleles = [allele.upper() for allele in alternate_alleles]

    # a record without alternate alleles is monomorphic for the reference allele,
    # which is derived when it is not the ancestral allele
    if not alternate_alleles:
        return [0.0] if ancestral_allele == reference_allele else [1.0]

    if ancestral_allele == reference_allele:
        return list(alternate_frequencies)

    if ancestral_allele in alternate_alleles:
        ancestral_index = alternate_alleles.index(ancestral_allele)

        return [
            0.0 if index == ancestral_index else frequency
            for index, frequency in enumerate(alternate_frequencies)
        ]

    return [1.0] * len(alternate_alleles)