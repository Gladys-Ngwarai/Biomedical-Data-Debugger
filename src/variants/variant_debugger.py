from src.variants.variant_evidence import (
    parse_format_evidence,
    calculate_variant_evidence,
)


def debug_variants(vcf):

    reports = []

    for _, row in vcf.iterrows():

        evidence = parse_format_evidence(
            row
        )

        reasons = calculate_variant_evidence(
            depth=evidence["depth"],
            allele_fraction=evidence[
                "allele_fraction"
            ],
            mapping_quality=evidence[
                "mapping_quality"
            ],
        )

        chrom = row.get(
            "#CHROM",
            row.get("CHROM")
        )

        position = row.get("POS")
        ref = row.get("REF")
        alt = row.get("ALT")

        reports.append(
            {
                "variant": (
                    f"{chrom}:{position} "
                    f"{ref}>{alt}"
                ),
                "chromosome": chrom,
                "position": position,
                "reference": ref,
                "alternate": alt,
                "depth": evidence["depth"],
                "allele_fraction": evidence[
                    "allele_fraction"
                ],
                "mapping_quality": evidence[
                    "mapping_quality"
                ],
                "artifact_evidence": len(
                    reasons
                ),
                "reasons": reasons,
                "requires_investigation": (
                    len(reasons) >= 1
                ),
            }
        )

    return reports