from pathlib import Path
import gzip


def open_variant_file(path):
    """Open plain or gzipped VCF files as text."""

    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"Variant file not found: {path}"
        )

    if path.name.lower().endswith(".vcf.gz"):
        return gzip.open(path, "rt")

    if path.name.lower().endswith(".vcf"):
        return open(path, "r")

    raise ValueError(
        "Unsupported file type. Use .vcf or .vcf.gz."
    )


def load_vcf(path):
    """
    Stream a VCF file record by record.

    Supports .vcf and .vcf.gz files.
    """

    with open_variant_file(path) as file:

        header = None

        for line_number, line in enumerate(file, start=1):

            line = line.rstrip("\r\n")

            if not line:
                continue

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.split("\t")
                continue

            if header is None:
                raise ValueError(
                    "VCF header (#CHROM) was not found."
                )

            fields = line.split("\t")

            if len(fields) != len(header):

                # Some simple test files may contain
                # whitespace instead of tabs.
                fields = line.split()

            if len(fields) != len(header):
                raise ValueError(
                    f"Malformed VCF record on line "
                    f"{line_number}: expected "
                    f"{len(header)} columns but found "
                    f"{len(fields)}."
                )

            yield dict(
                zip(header, fields)
            )