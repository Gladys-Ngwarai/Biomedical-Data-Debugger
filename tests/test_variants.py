from src.variants.load_vcf import load_vcf
from src.variants.variant_debugger import debug_variants


def main():

    path = "data/raw/test_variants.vcf"

    vcf = load_vcf(path)

    print()
    print("DNA Variant Debugger")
    print("=" * 40)

    print(
        f"Variants loaded: {len(vcf)}"
    )

    reports = debug_variants(vcf)

    print()

    for report in reports:

        print(
            f"Variant: "
            f"{report['variant']}"
        )

        print(
            f"Depth: "
            f"{report['depth']}"
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
            f"Evidence count: "
            f"{report['artifact_evidence']}"
        )

        if report["reasons"]:

            print(
                "Reasons: "
                + ", ".join(
                    report["reasons"]
                )
            )

        else:

            print(
                "Status: "
                "No major artifact evidence"
            )

        print("-" * 40)


if __name__ == "__main__":
    main()