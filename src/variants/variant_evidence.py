def parse_format_evidence(row):
    """
    Extract common variant evidence from a VCF record.

    Supports:
    DP = read depth
    AD = reference/alternate read counts
    AF = alternate allele frequency
    MQ = mapping quality
    """

    evidence = {
        "depth": None,
        "ref_reads": None,
        "alt_reads": None,
        "allele_fraction": None,
        "mapping_quality": None,
    }

    format_field = row.get("FORMAT")

    if not isinstance(format_field, str):
        return evidence

    excluded_columns = {
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

    sample_columns = [
        key
        for key in row.index
        if key not in excluded_columns
    ]

    if not sample_columns:
        return evidence

    sample = row[sample_columns[0]]

    # Pandas can represent missing values as floats.
    if not isinstance(sample, str):
        return evidence

    format_keys = format_field.split(":")
    sample_values = sample.split(":")

    values = dict(
        zip(format_keys, sample_values)
    )

    # -------------------------
    # DEPTH
    # -------------------------

    if "DP" in values:

        try:
            evidence["depth"] = int(
                values["DP"]
            )
        except (ValueError, TypeError):
            pass

    # -------------------------
    # ALLELE DEPTH
    # -------------------------

    if "AD" in values:

        try:

            ad = values["AD"].split(",")

            if len(ad) >= 2:

                evidence["ref_reads"] = int(
                    ad[0]
                )

                evidence["alt_reads"] = int(
                    ad[1]
                )

        except (ValueError, TypeError):
            pass

    # -------------------------
    # ALLELE FREQUENCY
    # -------------------------

    if "AF" in values:

        try:
            evidence["allele_fraction"] = float(
                values["AF"]
            )
        except (ValueError, TypeError):
            pass

    # -------------------------
    # MAPPING QUALITY
    # -------------------------

    if "MQ" in values:

        try:
            evidence["mapping_quality"] = float(
                values["MQ"]
            )
        except (ValueError, TypeError):
            pass

    # -------------------------
    # CALCULATE AF FROM AD
    # -------------------------

    if (
        evidence["allele_fraction"] is None
        and evidence["ref_reads"] is not None
        and evidence["alt_reads"] is not None
    ):

        total = (
            evidence["ref_reads"]
            + evidence["alt_reads"]
        )

        if total > 0:

            evidence["allele_fraction"] = (
                evidence["alt_reads"]
                / total
            )

    return evidence


def calculate_variant_evidence(
    depth,
    allele_fraction,
    mapping_quality,
    min_depth=10,
    min_vaf=0.20,
    max_vaf=0.80,
    min_mapping_quality=20,
):
    """
    Identify suspicious evidence for a variant.
    """

    reasons = []

    if depth is not None:

        if depth < min_depth:
            reasons.append(
                "low sequencing depth"
            )

    if allele_fraction is not None:

        if allele_fraction < min_vaf:
            reasons.append(
                "weak alternate-allele support"
            )

        if allele_fraction > max_vaf:
            reasons.append(
                "unusual allele balance"
            )

    if mapping_quality is not None:

        if mapping_quality < min_mapping_quality:
            reasons.append(
                "low mapping quality"
            )

    return reasons