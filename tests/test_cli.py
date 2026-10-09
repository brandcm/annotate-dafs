from pathlib import Path

import pysam
import pytest
import vcfpy

from annotate_dafs.cli import main


DATA_DIR = Path(__file__).parent / "data"


def read_records(path: Path) -> list[vcfpy.Record]:
    """Read all records from a VCF file."""
    reader = vcfpy.Reader.from_path(path)
    return list(reader)


def test_cli_annotates_daf_from_info(tmp_path):
    """CLI annotates DAF using ancestral alleles and INFO/AF."""
    input_vcf = DATA_DIR / "allele_frequencies.vcf"
    output_vcf = tmp_path / "output.vcf"

    exit_code = main(
        [
            "annotate",
            "--aa-field",
            "AA",
            "--vcf",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    reader = vcfpy.Reader.from_path(output_vcf)

    assert "DAF" in reader.header.info_ids()

    records = list(reader)

    assert records[0].INFO["DAF"] == [0.01]
    assert records[1].INFO["DAF"] == [0.1]
    assert records[2].INFO["DAF"] == [0.0]
    assert records[3].INFO["DAF"] == [0.0]
    assert records[4].INFO["DAF"] == [0.0, 0.03]
    assert records[5].INFO["DAF"] == [0.08, 0.09]


def test_cli_uses_selected_samples(tmp_path):
    """CLI uses genotype-derived AF when samples are explicitly selected."""
    input_vcf = DATA_DIR / "genotypes.vcf"
    output_vcf = tmp_path / "output.vcf"

    exit_code = main(
        [
            "annotate",
            "--aa-field",
            "AA",
            "--samples",
            "sample4",
            "--vcf",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    records = read_records(output_vcf)

    assert records[0].INFO["DAF"] == [0.5]


def test_cli_uses_fasta(tmp_path):
    """CLI annotates DAF and AA using ancestral alleles from a FASTA."""
    input_vcf = DATA_DIR / "cli_input.vcf"
    ancestral_fasta = DATA_DIR / "ancestral.fa"
    output_vcf = tmp_path / "output.vcf"

    exit_code = main(
        [
            "annotate",
            "--fasta",
            str(ancestral_fasta),
            "--vcf",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    records = read_records(output_vcf)

    assert records[0].INFO["AA"] == "C"
    assert records[0].INFO["DAF"] == [1.0]


def test_cli_fasta_replaces_existing_ancestral_allele(tmp_path):
    """FASTA ancestry replaces an existing AA value used by the input VCF."""
    input_vcf = tmp_path / "input.vcf"
    ancestral_fasta = tmp_path / "ancestral.fa"
    output_vcf = tmp_path / "output.vcf"

    input_vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr22,length=10>\n"
        '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">\n'
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr22\t10\t.\tG\tA\t30\tPASS\tAA=G;AF=0.2\n"
    )

    ancestral_fasta.write_text(
        ">chr22\n"
        "AAAAAAAAAA\n"
    )
    pysam.faidx(str(ancestral_fasta))

    exit_code = main(
        [
            "annotate",
            "--fasta",
            str(ancestral_fasta),
            "--vcf",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    records = read_records(output_vcf)

    assert records[0].INFO["AA"] == "A"
    assert records[0].INFO["DAF"] == [0.0]


def test_cli_fasta_adds_ancestral_allele_header(tmp_path):
    """FASTA mode adds a String ancestral-allele INFO definition."""
    input_vcf = tmp_path / "input.vcf"
    ancestral_fasta = tmp_path / "ancestral.fa"
    output_vcf = tmp_path / "output.vcf"

    input_vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr22,length=10>\n"
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr22\t10\t.\tG\tA\t30\tPASS\tAF=0.2\n"
    )

    ancestral_fasta.write_text(
        ">chr22\n"
        "AAAAAAAAAA\n"
    )
    pysam.faidx(str(ancestral_fasta))

    exit_code = main(
        [
            "annotate",
            "--fasta",
            str(ancestral_fasta),
            "--vcf",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    reader = vcfpy.Reader.from_path(output_vcf)

    aa_lines = [
        line
        for line in reader.header.get_lines("INFO")
        if line.mapping["ID"] == "AA"
    ]

    assert len(aa_lines) == 1
    assert aa_lines[0].mapping["Number"] == 1
    assert aa_lines[0].mapping["Type"] == "String"


def test_cli_excludes_low_confidence_ancestral(tmp_path):
    """CLI skips DAF annotation for lowercase ancestral alleles when requested."""
    input_vcf = DATA_DIR / "cli_input.vcf"
    ancestral_fasta = tmp_path / "ancestral.fa"
    output_vcf = tmp_path / "output.vcf"

    ancestral_fasta.write_text(">chr22\nAAAAAAAAAa\n")
    pysam.faidx(str(ancestral_fasta))

    exit_code = main(
        [
            "annotate",
            "--fasta",
            str(ancestral_fasta),
            "--exclude-low-confidence-ancestral",
            "--vcf",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    records = read_records(output_vcf)

    assert "DAF" not in records[0].INFO


def test_cli_adds_daf_header(tmp_path):
    """CLI adds the DAF INFO definition to the output VCF header."""
    input_vcf = DATA_DIR / "cli_input.vcf"
    output_vcf = tmp_path / "output.vcf"

    exit_code = main(
        [
            "annotate",
            "--aa-field",
            "AA",
            "--vcf",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    reader = vcfpy.Reader.from_path(output_vcf)

    assert "DAF" in reader.header.info_ids()


def test_cli_uses_custom_alts_field(tmp_path):
    """CLI uses a custom ALT field with the corresponding INFO AF field."""
    input_vcf = DATA_DIR / "cli_input.vcf"
    output_vcf = tmp_path / "custom_alts_output.vcf"

    exit_code = main(
        [
            "annotate",
            "--aa-field",
            "AA",
            "--vcf",
            str(input_vcf),
            "--af-field",
            "CUSTOM_AF",
            "--alts-field",
            "CUSTOM_ALTS",
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    records = read_records(output_vcf)

    # VCF ALT is G, but CUSTOM_ALTS is T with CUSTOM_AF=0.4.
    # AA=A, so T is derived and should retain its frequency.
    assert records[0].INFO["DAF"] == pytest.approx([0.4])


def test_cli_annotates_monomorphic_records(tmp_path):
    """CLI annotates DAF for records without an ALT allele."""
    input_vcf = DATA_DIR / "monomorphic.vcf"
    output_vcf = tmp_path / "monomorphic_output.vcf"

    exit_code = main(
        [
            "annotate",
            "--aa-field",
            "AA",
            "--vcf",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    records = read_records(output_vcf)

    # The reference allele is ancestral, so no derived allele is present.
    assert records[0].INFO["DAF"] == pytest.approx([0.0])

    # The reference allele is derived and fixed in every genotype.
    assert records[1].INFO["DAF"] == pytest.approx([1.0])

    assert records[2].INFO["DAF"] == pytest.approx([0.0])
    assert records[3].INFO["DAF"] == pytest.approx([0.25])


def test_cli_annotates_monomorphic_records_with_custom_alts_field(tmp_path):
    """A custom ALT field is annotated on records without an ALT allele."""
    input_vcf = DATA_DIR / "monomorphic.vcf"
    output_vcf = tmp_path / "monomorphic_custom_alts_output.vcf"

    exit_code = main(
        [
            "annotate",
            "--aa-field",
            "AA",
            "--vcf",
            str(input_vcf),
            "--af-field",
            "CUSTOM_AF",
            "--alts-field",
            "CUSTOM_ALTS",
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    records = read_records(output_vcf)

    # The record has no custom allele frequency, so it is not annotated.
    assert "DAF" not in records[0].INFO

    # The custom ALT field is annotated although the record has no ALT allele.
    assert records[2].INFO["DAF"] == pytest.approx([0.4])

    # The custom ALT field holds more alleles than the ALT column.
    assert records[3].INFO["DAF"] == pytest.approx([0.1, 0.2])


def test_cli_annotates_monomorphic_records_from_fasta(tmp_path):
    """CLI annotates DAF for records without an ALT allele using an ancestral FASTA."""
    input_vcf = tmp_path / "monomorphic_fasta_input.vcf"
    ancestral_fasta = tmp_path / "ancestral.fa"
    output_vcf = tmp_path / "monomorphic_fasta_output.vcf"

    input_vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr22,length=20>\n"
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
        "chr22\t10\t.\tA\t.\t30\tPASS\tAN=2\tGT\t0/0\n"
        "chr22\t20\t.\tG\t.\t30\tPASS\tAN=2\tGT\t0/0\n"
    )

    ancestral_fasta.write_text(">chr22\n" + "A" * 20 + "\n")
    pysam.faidx(str(ancestral_fasta))

    exit_code = main(
        [
            "annotate",
            "--fasta",
            str(ancestral_fasta),
            "--vcf",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    records = read_records(output_vcf)

    # The ancestral allele matches REF, so the reference allele is not derived.
    assert records[0].INFO["AA"] == "A"
    assert records[0].INFO["DAF"] == pytest.approx([0.0])

    # The ancestral allele differs from REF, so the reference allele is derived.
    assert records[1].INFO["AA"] == "A"
    assert records[1].INFO["DAF"] == pytest.approx([1.0])

def test_cli_daf_header_permits_any_number_of_values(tmp_path):
    """The DAF INFO definition permits a variable number of values."""
    input_vcf = DATA_DIR / "cli_input.vcf"
    output_vcf = tmp_path / "daf_header_output.vcf"

    exit_code = main(
        [
            "annotate",
            "--aa-field",
            "AA",
            "--vcf",
            str(input_vcf),
            "--output",
            str(output_vcf),
        ]
    )

    assert exit_code == 0

    reader = vcfpy.Reader.from_path(output_vcf)

    daf_lines = [
        line
        for line in reader.header.get_lines("INFO")
        if line.mapping["ID"] == "DAF"
    ]

    assert len(daf_lines) == 1
    assert daf_lines[0].mapping["Number"] == "."
    assert daf_lines[0].mapping["Type"] == "Float"


def test_cli_download_ancestral_grch37(tmp_path, monkeypatch):
    """Download the GRCh37 ancestral FASTA with the expected arguments."""
    calls = []

    def mock_download_ancestral_fasta(
        assembly,
        output_directory,
        release=None,
        force=False,
    ):
        calls.append(
            {
                "assembly": assembly,
                "output_directory": output_directory,
                "release": release,
                "force": force,
            }
        )

    monkeypatch.setattr(
        "annotate_dafs.cli.download_ancestral_fasta",
        mock_download_ancestral_fasta,
    )

    output_directory = tmp_path / "ancestral_sequences"

    exit_code = main(
        [
            "download-ancestral",
            "--assembly",
            "GRCh37",
            "--output-dir",
            str(output_directory),
        ]
    )

    assert exit_code == 0
    assert calls == [
        {
            "assembly": "GRCh37",
            "output_directory": output_directory,
            "release": None,
            "force": False,
        }
    ]


def test_cli_download_ancestral_grch38(tmp_path, monkeypatch):
    """Download a GRCh38 ancestral FASTA for the requested Ensembl release."""
    calls = []

    def mock_download_ancestral_fasta(
        assembly,
        output_directory,
        release=None,
        force=False,
    ):
        calls.append(
            {
                "assembly": assembly,
                "output_directory": output_directory,
                "release": release,
                "force": force,
            }
        )

    monkeypatch.setattr(
        "annotate_dafs.cli.download_ancestral_fasta",
        mock_download_ancestral_fasta,
    )

    output_directory = tmp_path / "ancestral_sequences"

    exit_code = main(
        [
            "download-ancestral",
            "--assembly",
            "GRCh38",
            "--release",
            "115",
            "--output-dir",
            str(output_directory),
        ]
    )

    assert exit_code == 0
    assert calls == [
        {
            "assembly": "GRCh38",
            "output_directory": output_directory,
            "release": 115,
            "force": False,
        }
    ]


def test_cli_download_ancestral_force(tmp_path, monkeypatch):
    """Pass the force option to the ancestral FASTA downloader."""
    calls = []

    def mock_download_ancestral_fasta(
        assembly,
        output_directory,
        release=None,
        force=False,
    ):
        calls.append(
            {
                "assembly": assembly,
                "output_directory": output_directory,
                "release": release,
                "force": force,
            }
        )

    monkeypatch.setattr(
        "annotate_dafs.cli.download_ancestral_fasta",
        mock_download_ancestral_fasta,
    )

    output_directory = tmp_path / "ancestral_sequences"

    exit_code = main(
        [
            "download-ancestral",
            "--assembly",
            "GRCh37",
            "--output-dir",
            str(output_directory),
            "--force",
        ]
    )

    assert exit_code == 0
    assert calls == [
        {
            "assembly": "GRCh37",
            "output_directory": output_directory,
            "release": None,
            "force": True,
        }
    ]


def test_cli_download_ancestral_grch37_rejects_release(tmp_path):
    """Reject an explicit Ensembl release for GRCh37."""
    with pytest.raises(SystemExit) as error:
        main(
            [
                "download-ancestral",
                "--assembly",
                "GRCh37",
                "--release",
                "75",
                "--output-dir",
                str(tmp_path),
            ]
        )

    assert error.value.code == 2


def test_cli_download_ancestral_grch38_requires_release(tmp_path):
    """Require an Ensembl release when downloading GRCh38."""
    with pytest.raises(SystemExit) as error:
        main(
            [
                "download-ancestral",
                "--assembly",
                "GRCh38",
                "--output-dir",
                str(tmp_path),
            ]
        )

    assert error.value.code == 2
