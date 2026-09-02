import pandas as pd


def load_vcf(path):
    """
    Load a VCF file into a pandas DataFrame.

    This lightweight loader focuses on standard variant
    fields and FORMAT/sample evidence.
    """

    records = []

    with open(path, "r") as file:

        for line in file:

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):

                header = line.strip().split("\t")

                continue

            fields = line.strip().split("\t")

            if len(fields) < 8:
                continue

            record = dict(
                zip(header, fields)
            )

            records.append(record)

    return pd.DataFrame(records)