def to_int(value):
    try:
        if value is None or value == ".":
            return None

        return int(value)

    except (TypeError, ValueError):
        return None


def to_float(value):
    try:
        if value is None or value == ".":
            return None

        return float(value)

    except (TypeError, ValueError):
        return None


def parse_info(info):
    """Parse the VCF INFO field."""

    result = {}

    if not info or info == ".":
        return result

    for item in info.split(";"):

        if "=" in item:

            key, value = item.split(
                "=",
                1
            )

            result[key] = value

        else:

            result[item] = True

    return result


def parse_sample(
    format_field,
    sample_field,
):
    """Parse FORMAT and sample values."""

    if (
        not format_field
        or format_field == "."
        or not sample_field
        or sample_field == "."
    ):
        return {}

    keys = format_field.split(":")
    values = sample_field.split(":")

    result = {}

    for key, value in zip(
        keys,
        values
    ):
        result[key] = value

    return result


def parse_variant_evidence(
    row,
    sample_name=None,
):
    """
    Extract available evidence from a VCF record.

    Evidence sources:

    INFO:
        DP, AF, MQ

    FORMAT/sample:
        GT, DP, AD, AF, MQ
    """

    info = parse_info(
        row.get("INFO", ".")
    )

    standard_columns = {
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
        for key in row
        if key not in standard_columns
    ]

    # Select requested sample.
    if sample_name is not None:

        sample_field = row.get(
            sample_name
        )

    # Otherwise use first sample.
    elif sample_columns:

        sample_name = sample_columns[0]

        sample_field = row.get(
            sample_name
        )

    else:

        sample_field = None

    sample = parse_sample(
        row.get("FORMAT", "."),
        sample_field,
    )

    # ----------------------------
    # Depth
    # ----------------------------

    depth = to_int(
        sample.get("DP")
    )

    if depth is None:
        depth = to_int(
            info.get("DP")
        )

    # ----------------------------
    # Allele depth
    # ----------------------------

    ref_reads = None
    alt_reads = None

    ad = sample.get("AD")

    if ad and ad != ".":

        parts = ad.split(",")

        if len(parts) >= 2:

            ref_reads = to_int(
                parts[0]
            )

            alt_reads = to_int(
                parts[1]
            )

    # ----------------------------
    # Allele fraction
    # ----------------------------

    allele_fraction = None

    af = sample.get("AF")

    if af is None:
        af = info.get("AF")

    if af is not None:

        af = str(af)

        if "," in af:
            af = af.split(",")[0]

        allele_fraction = to_float(
            af
        )

    # Calculate from AD if AF absent.
    if (
        allele_fraction is None
        and ref_reads is not None
        and alt_reads is not None
    ):

        total = (
            ref_reads
            + alt_reads
        )

        if total > 0:

            allele_fraction = (
                alt_reads / total
            )

    # ----------------------------
    # Mapping quality
    # ----------------------------

    mapping_quality = to_float(
        sample.get("MQ")
    )

    if mapping_quality is None:

        mapping_quality = to_float(
            info.get("MQ")
        )

    # ----------------------------
    # Other evidence
    # ----------------------------

    quality = to_float(
        row.get("QUAL")
    )

    genotype = sample.get("GT")

    filter_status = row.get(
        "FILTER",
        "."
    )

    return {
        "sample": sample_name,
        "depth": depth,
        "ref_reads": ref_reads,
        "alt_reads": alt_reads,
        "allele_fraction": allele_fraction,
        "mapping_quality": mapping_quality,
        "quality": quality,
        "filter": filter_status,
        "genotype": genotype,
    }


def calculate_variant_evidence(
    evidence,
    min_depth=10,
    min_vaf=0.20,
    max_vaf=0.80,
    min_mapping_quality=20,
    min_quality=20,
):
    """
    Identify evidence requiring investigation.

    This does NOT declare a variant false.
    """

    reasons = []

    depth = evidence.get(
        "depth"
    )

    vaf = evidence.get(
        "allele_fraction"
    )

    mq = evidence.get(
        "mapping_quality"
    )

    quality = evidence.get(
        "quality"
    )

    filter_status = evidence.get(
        "filter"
    )

    # Low depth
    if (
        depth is not None
        and depth < min_depth
    ):
        reasons.append(
            "low sequencing depth"
        )

    # Allele support
    if vaf is not None:

        if vaf < min_vaf:

            reasons.append(
                "weak alternate-allele support"
            )

        elif vaf > max_vaf:

            reasons.append(
                "unusual allele balance"
            )

    # Mapping quality
    if (
        mq is not None
        and mq < min_mapping_quality
    ):
        reasons.append(
            "low mapping quality"
        )

    # Variant quality
    if (
        quality is not None
        and quality < min_quality
    ):
        reasons.append(
            "low variant quality"
        )

    # Caller filter
    if filter_status not in {
        None,
        "",
        ".",
        "PASS",
    }:

        reasons.append(
            f"variant caller filter: "
            f"{filter_status}"
        )

    return reasons