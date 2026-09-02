from src.variants.variant_evidence import (
    parse_variant_evidence,
    calculate_variant_evidence,
)


VCF_COLUMNS = {
    "#CHROM",
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT",
}


def get_sample_names(row):
    """Return sample columns from a VCF record."""

    return [
        key
        for key in row
        if key not in VCF_COLUMNS
    ]


def debug_variants(
    vcf_records,
    max_reports=None,
):
    """
    Debug VCF records.

    Supports single- and multi-sample VCFs.

    The loader streams records, so the complete VCF
    does not need to be held in memory.

    max_reports limits stored reports while all
    variants are still analyzed.
    """

    reports = []

    variants_analyzed = 0
    suspicious_records = 0

    sample_names_seen = set()

    for row in vcf_records:

        variants_analyzed += 1

        sample_names = get_sample_names(
            row
        )

        if not sample_names:
            sample_names = [None]

        for sample_name in sample_names:

            if sample_name is not None:

                sample_names_seen.add(
                    sample_name
                )

            evidence = parse_variant_evidence(
                row,
                sample_name=sample_name,
            )

            reasons = calculate_variant_evidence(
                evidence
            )

            if reasons:
                suspicious_records += 1

            chrom = row.get(
                "#CHROM",
                row.get("CHROM")
            )

            position = row.get(
                "POS"
            )

            reference = row.get(
                "REF"
            )

            alternate = row.get(
                "ALT"
            )

            report = {
                "variant": (
                    f"{chrom}:"
                    f"{position} "
                    f"{reference}>"
                    f"{alternate}"
                ),

                "sample": (
                    sample_name
                    or "No sample"
                ),

                "chromosome": chrom,

                "position": position,

                "reference": reference,

                "alternate": alternate,

                "genotype": evidence.get(
                    "genotype"
                ),

                "depth": evidence.get(
                    "depth"
                ),

                "ref_reads": evidence.get(
                    "ref_reads"
                ),

                "alt_reads": evidence.get(
                    "alt_reads"
                ),

                "allele_fraction": evidence.get(
                    "allele_fraction"
                ),

                "mapping_quality": evidence.get(
                    "mapping_quality"
                ),

                "quality": evidence.get(
                    "quality"
                ),

                "filter": evidence.get(
                    "filter"
                ),

                "artifact_evidence": len(
                    reasons
                ),

                "reasons": reasons,

                "requires_investigation": bool(
                    reasons
                ),
            }

            if (
                max_reports is None
                or len(reports) < max_reports
            ):

                reports.append(
                    report
                )

    summary = {
        "variants_analyzed":
            variants_analyzed,

        "reports_returned":
            len(reports),

        "samples":
            sorted(sample_names_seen),

        "suspicious_records":
            suspicious_records,
    }

    return reports, summary