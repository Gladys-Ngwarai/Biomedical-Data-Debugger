from src.variants.load_vcf import load_vcf
from src.variants.variant_debugger import debug_variants


def main():

    path = "data/raw/test_variants.vcf"

    records = load_vcf(path)

    reports, summary = debug_variants(records)

    print()
    print("DNA Variant Debugger")
    print("=" * 40)

    print(
        f"Variants analyzed: "
        f"{summary['variants_analyzed']}"
    )

    print(
        f"Variants requiring investigation: "
        f"{summary['suspicious_records']}"
    )

    if summary["samples"]:
        samples = ", ".join(summary["samples"])
    else:
        samples = "None"

    print(
        f"Samples detected: {samples}"
    )

    print()

    for report in reports:

        print(
            f"Variant: {report['variant']}"
        )

        print(
            f"Sample: {report['sample']}"
        )

        print(
            f"Genotype: {report['genotype']}"
        )

        print(
            f"Depth: {report['depth']}"
        )

        print(
            f"Reference reads: "
            f"{report['ref_reads']}"
        )

        print(
            f"Alternate reads: "
            f"{report['alt_reads']}"
        )

        print(
            f"Allele fraction: "
            f"{report['allele_fraction']}"
        )

        print(
            f"Mapping quality: "
            f"{report['mapping_quality']}"
        )

        print(
            f"Quality: "
            f"{report['quality']}"
        )

        print(
            f"Filter: "
            f"{report['filter']}"
        )

        print(
            f"Evidence count: "
            f"{report['artifact_evidence']}"
        )

        if report["requires_investigation"]:

            print(
                "Status: Requires investigation"
            )

            print(
                "Reasons: "
                + ", ".join(
                    report["reasons"]
                )
            )

        else:

            print(
                "Status: No major artifact evidence"
            )

        print("-" * 40)


if __name__ == "__main__":
    main()