from pathlib import Path

import pytest
import vcfpy

from annotate_dafs.frequencies import (
    calculate_allele_frequencies,
    calculate_daf,
    get_alt_alleles_and_frequencies,
    get_genotype_calls,
)


DATA_DIR = Path(__file__).parent / "data"


def read_records(filename: str) -> list[vcfpy.Record]:
    """Read all records from a test VCF."""
    reader = vcfpy.Reader.from_path(DATA_DIR / filename)
    return list(reader)


@pytest.fixture
def genotype_records() -> list[vcfpy.Record]:
    """Records from the genotype test VCF."""
    return read_records("genotypes.vcf")


@pytest.fixture
def af_records() -> list[vcfpy.Record]:
    """Records from the allele-frequency test VCF."""
    return read_records("allele_frequencies.vcf")


def test_get_genotype_calls_unphased(genotype_records):
    """Unphased genotypes are parsed into individual allele calls."""
    variant = genotype_records[2]

    calls = get_genotype_calls(
        variant,
        {"sample1", "sample2", "sample3", "sample4"},
    )

    assert calls == [[0, 1], [0, 1], [0, 0], [None, None]]


def test_get_genotype_calls_phased(genotype_records):
    """Phased genotypes are parsed into individual allele calls."""
    variant = genotype_records[3]

    calls = get_genotype_calls(
        variant,
        {"sample1", "sample2", "sample3", "sample4"},
    )

    assert calls == [[0, 0], [0, 1], [1, 1], [0, 1]]


def test_get_genotype_calls_selected_samples(genotype_records):
    """Only requested samples are returned."""
    variant = genotype_records[0]

    calls = get_genotype_calls(
        variant,
        {"sample4"},
    )

    assert calls == [[0, 1]]


def test_calculate_allele_frequencies_biallelic(genotype_records):
    """Allele frequencies are calculated correctly for a biallelic variant."""
    variant = genotype_records[0]

    calls = get_genotype_calls(
        variant,
        {"sample1", "sample2", "sample3", "sample4"},
    )

    frequencies = calculate_allele_frequencies(calls, len(variant.ALT))

    assert frequencies == pytest.approx([0.125])


def test_calculate_allele_frequencies_missing(genotype_records):
    """Missing genotypes are excluded from allele-frequency calculations."""
    variant = genotype_records[2]

    calls = get_genotype_calls(
        variant,
        {"sample1", "sample2", "sample3", "sample4"},
    )

    frequencies = calculate_allele_frequencies(calls, len(variant.ALT))

    assert frequencies == pytest.approx([2 / 6])


def test_calculate_allele_frequencies_partial_missing():
    """Missing alleles are excluded while called alleles are retained."""
    frequencies = calculate_allele_frequencies(
        [[0, None], [None, 1]],
        number_of_alts=1,
    )

    assert frequencies == pytest.approx([0.5])


def test_calculate_allele_frequencies_multiallelic(genotype_records):
    """Allele frequencies are calculated independently for multiple ALT alleles."""
    variant = genotype_records[4]

    calls = get_genotype_calls(
        variant,
        {"sample1", "sample2", "sample3", "sample4"},
    )

    frequencies = calculate_allele_frequencies(calls, len(variant.ALT))

    assert frequencies == pytest.approx([1 / 6, 1 / 6])


def test_info_af_is_used_when_samples_not_selected(af_records):
    """INFO/AF is used when no samples are explicitly selected."""
    variant = af_records[0]

    alleles, frequencies = get_alt_alleles_and_frequencies(
        variant,
        samples=None,
        af_field="AF",
    )

    assert alleles == ["G"]
    assert frequencies == pytest.approx([0.01])


def test_selected_samples_override_info_af(genotype_records, af_records):
    """Explicit sample selection takes precedence over INFO/AF."""
    genotype_variant = genotype_records[0]
    af_variant = af_records[0]

    # The INFO/AF value is 0.01, whereas sample4 has one ALT allele
    # out of two called alleles, giving an allele frequency of 0.5.
    assert af_variant.INFO["AF"] == [0.01]

    alleles, frequencies = get_alt_alleles_and_frequencies(
        genotype_variant,
        samples={"sample4"},
        af_field="AF",
    )

    assert alleles == ["G"]
    assert frequencies == pytest.approx([0.5])


def test_calculate_daf_ancestral_is_reference():
    """ALT frequencies are derived frequencies when REF is ancestral."""
    daf = calculate_daf(
        ancestral_allele="A",
        reference_allele="A",
        alternate_alleles=["G"],
        alternate_frequencies=[0.125],
    )

    assert daf == pytest.approx([0.125])


def test_calculate_daf_ancestral_is_alternate():
    """An ancestral ALT has zero derived allele frequency."""
    daf = calculate_daf(
        ancestral_allele="G",
        reference_allele="A",
        alternate_alleles=["G"],
        alternate_frequencies=[0.125],
    )

    assert daf == pytest.approx([0.0])


def test_calculate_daf_ancestral_is_multiallelic_alternate():
    """Ancestral ALT is zero while other ALT alleles retain their frequencies."""
    daf = calculate_daf(
        ancestral_allele="A",
        reference_allele="C",
        alternate_alleles=["A", "T"],
        alternate_frequencies=[0.01, 0.03],
    )

    assert daf == pytest.approx([0.0, 0.03])


def test_calculate_daf_ancestral_not_observed():
    """All ALT alleles are treated as derived when the ancestral allele is absent."""
    daf = calculate_daf(
        ancestral_allele="G",
        reference_allele="C",
        alternate_alleles=["A", "T"],
        alternate_frequencies=[0.01, 0.03],
    )

    assert daf == pytest.approx([1.0, 1.0])


def test_calculate_daf_frequency_length_mismatch():
    """A mismatch between ALT alleles and frequencies raises an error."""
    with pytest.raises(
        ValueError,
        match="Number of allele frequencies does not match",
    ):
        calculate_daf(
            ancestral_allele="A",
            reference_allele="A",
            alternate_alleles=["G", "T"],
            alternate_frequencies=[0.1],
        )


def test_calculate_daf_missing_frequencies():
    """Missing allele frequencies produce missing DAF."""
    daf = calculate_daf(
        ancestral_allele="A",
        reference_allele="A",
        alternate_alleles=["G"],
        alternate_frequencies=None,
    )

    assert daf is None


def test_custom_alts_field_is_used_with_info_af(af_records):
    """A custom ALT field is paired with the corresponding INFO AF field."""
    variant = af_records[0]

    variant.INFO["CUSTOM_AF"] = [0.25]
    variant.INFO["CUSTOM_ALTS"] = ["T"]

    alleles, frequencies = get_alt_alleles_and_frequencies(
        variant,
        samples=None,
        af_field="CUSTOM_AF",
        alts_field="CUSTOM_ALTS",
    )

    assert alleles == ["T"]
    assert frequencies == pytest.approx([0.25])

@pytest.fixture
def monomorphic_records() -> list[vcfpy.Record]:
    """Records from the monomorphic test VCF."""
    return read_records("monomorphic.vcf")


def test_calculate_daf_monomorphic_ancestral_is_reference():
    """A monomorphic reference allele that is ancestral has zero derived frequency."""
    daf = calculate_daf(
        ancestral_allele="A",
        reference_allele="A",
        alternate_alleles=[],
        alternate_frequencies=[],
    )

    assert daf == pytest.approx([0.0])


def test_calculate_daf_monomorphic_ancestral_is_not_reference():
    """A monomorphic reference allele that is derived has a derived frequency of one."""
    daf = calculate_daf(
        ancestral_allele="G",
        reference_allele="A",
        alternate_alleles=[],
        alternate_frequencies=[],
    )

    assert daf == pytest.approx([1.0])


def test_monomorphic_record_has_no_alleles_or_frequencies(monomorphic_records):
    """A record without an ALT allele yields empty alleles and frequencies."""
    variant = monomorphic_records[0]

    alleles, frequencies = get_alt_alleles_and_frequencies(
        variant,
        samples=None,
        af_field="AF",
    )

    assert alleles == []
    assert frequencies == []


def test_monomorphic_record_uses_custom_alts_field(monomorphic_records):
    """A custom ALT field is used even when the record has no ALT allele."""
    variant = monomorphic_records[2]

    alleles, frequencies = get_alt_alleles_and_frequencies(
        variant,
        samples=None,
        af_field="CUSTOM_AF",
        alts_field="CUSTOM_ALTS",
    )

    assert alleles == ["T"]
    assert frequencies == pytest.approx([0.4])


def test_missing_allele_frequency_is_treated_as_missing(af_records):
    """A missing allele frequency value produces no frequencies."""
    variant = af_records[0]

    variant.INFO["AF"] = [None]

    assert get_alt_alleles_and_frequencies(
        variant,
        samples=None,
        af_field="AF",
    ) is None


def test_missing_alts_field_is_treated_as_missing(af_records):
    """A record without the requested ALT field produces no frequencies."""
    variant = af_records[0]

    variant.INFO["CUSTOM_AF"] = [0.25]

    assert get_alt_alleles_and_frequencies(
        variant,
        samples=None,
        af_field="CUSTOM_AF",
        alts_field="CUSTOM_ALTS",
    ) is None


def test_allele_frequency_out_of_range(af_records):
    """An allele frequency outside [0, 1] raises an error."""
    variant = af_records[0]

    variant.INFO["AF"] = [1.5]

    with pytest.raises(ValueError, match="Invalid allele frequency"):
        get_alt_alleles_and_frequencies(
            variant,
            samples=None,
            af_field="AF",
        )


def test_allele_and_frequency_length_mismatch(af_records):
    """A custom ALT field of a different length to its AF field raises an error."""
    variant = af_records[0]

    variant.INFO["CUSTOM_AF"] = [0.25]
    variant.INFO["CUSTOM_ALTS"] = ["T", "C"]

    with pytest.raises(
        ValueError,
        match="does not match",
    ):
        get_alt_alleles_and_frequencies(
            variant,
            samples=None,
            af_field="CUSTOM_AF",
            alts_field="CUSTOM_ALTS",
        )


def test_calculate_allele_frequencies_all_missing():
    """All-missing genotypes produce no allele frequencies."""
    frequencies = calculate_allele_frequencies(
        [[None, None], [None, None]],
        number_of_alts=1,
    )

    assert frequencies is None


def test_calculate_allele_frequencies_invalid_allele_index():
    """A genotype allele index beyond the ALT count raises an error."""
    with pytest.raises(ValueError, match="invalid allele index"):
        calculate_allele_frequencies(
            [[0, 2]],
            number_of_alts=1,
        )

def test_omitted_alts_field_with_zero_frequency_is_monomorphic(monomorphic_records):
    """An omitted ALT field with an allele frequency of zero means no alternate allele."""
    variant = monomorphic_records[4]

    alleles, frequencies = get_alt_alleles_and_frequencies(
        variant,
        samples=None,
        af_field="CUSTOM_AF",
        alts_field="CUSTOM_ALTS",
    )

    assert alleles == []
    assert frequencies == []


def test_omitted_alts_field_with_a_frequency_is_treated_as_missing(monomorphic_records):
    """An omitted ALT field alongside a non-zero allele frequency produces no frequencies."""
    variant = monomorphic_records[6]

    assert get_alt_alleles_and_frequencies(
        variant,
        samples=None,
        af_field="CUSTOM_AF",
        alts_field="CUSTOM_ALTS",
    ) is None
