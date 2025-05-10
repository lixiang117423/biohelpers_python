#!/usr/bin/env python3

import argparse

import pandas as pd


def process_gene(file_gtf, species, distance):
    df_gene = file_gtf[file_gtf["feature"] == "gene"].copy()
    df_gene["species"] = species
    df_gene["temp_chr"] = df_gene["chr"].str.replace("Chr", "")
    df_gene["temp_id"] = (
        df_gene["attribute"].str.replace("g", "").astype(int) * distance
    )

    df_gene = df_gene.sort_values(["temp_chr", "temp_id"])
    max_id = df_gene["temp_id"].max()
    len_all = len(str(max_id))
    df_gene["len_gene"] = df_gene["temp_id"].astype(str).str.len()
    df_gene["len_diff"] = len_all - df_gene["len_gene"]
    df_gene["new_id"] = (
        df_gene["species"]
        + df_gene["temp_chr"].str.zfill(2)
        + "g"
        + df_gene["len_diff"].apply(lambda x: "0" * x)
        + df_gene["temp_id"].astype(str)
    )
    df_gene["attribute_new"] = (
        "ID=" + df_gene["new_id"] + ";" + "Name=" + df_gene["new_id"]
    )
    df_gene_final = df_gene[
        [
            "chr",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
            "new_id",
            "attribute_new",
        ]
    ]
    return df_gene_final


def process_mrna(file_gtf, df_gene_final):
    df_mrna = file_gtf[file_gtf["feature"] == "transcript"].copy()
    df_mrna["gene_id"] = df_mrna["attribute"].str.split("\\.").str[0]
    df_gene_temp = df_gene_final[["chr", "attribute", "new_id"]].rename(
        columns={"attribute": "gene_id"}
    )
    df_mrna = pd.merge(df_mrna, df_gene_temp, on=["chr", "gene_id"], how="left")
    df_mrna["mrna_id"] = (
        df_mrna["attribute"].str.split("\\.").str[1].str.replace("t", "mRNA")
    )
    df_mrna["mrna_id_final"] = df_mrna["new_id"] + "." + df_mrna["mrna_id"]
    df_mrna["attribute_new"] = (
        "ID="
        + df_mrna["mrna_id_final"]
        + ";"
        + "Name="
        + df_mrna["mrna_id_final"]
        + ";"
        + "Parent="
        + df_mrna["new_id"]
    )
    return df_mrna, df_mrna[["chr", "attribute", "attribute_new"]].rename(
        columns={"attribute": "mrna_id"}
    )


def process_others(file_gtf, df_mrna_temp):
    df_other = file_gtf[~file_gtf["feature"].isin(["gene", "transcript"])].copy()
    df_other["mrna_id"] = df_other["attribute"].str.split('"').str[1]
    df_other = pd.merge(df_other, df_mrna_temp, on=["chr", "mrna_id"])
    df_other["attribute"] = df_other["attribute_new"]
    return df_other


def main():
    parser = argparse.ArgumentParser(
        description="Process Braker GTF to GFF3 with customized gene/mRNA/feature annotation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--inout", required=True, help="Input GTF file")
    parser.add_argument(
        "-s",
        "--species",
        required=True,
        default="Os",
        help="Species name string (for prefix, e.g. Os)",
    )
    parser.add_argument(
        "-d",
        "--distance",
        type=int,
        default=10,
        help="Gene id multiplier/distance, like Os01g000010 and Os01g000020.",
    )
    parser.add_argument("-o", "--output", required=True, help="Output GFF3 file")

    args = parser.parse_args()

    col_names = [
        "chr",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
    file_gtf = pd.read_table(
        args.inout, delimiter="\t", header=None, comment="#", names=col_names
    )

    df_gene_final = process_gene(file_gtf, args.species, args.distance)
    df_mrna, df_mrna_temp = process_mrna(file_gtf, df_gene_final)
    df_other = process_others(file_gtf, df_mrna_temp)

    df_gene_out = df_gene_final[
        [
            "chr",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute_new",
        ]
    ].rename(columns={"attribute_new": "attribute"})
    df_mrna_out = df_mrna[
        [
            "chr",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute_new",
        ]
    ].rename(columns={"attribute_new": "attribute"})
    df_other_out = df_other[
        [
            "chr",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]
    ]

    df_final_temp = pd.concat(
        [df_gene_out, df_mrna_out, df_other_out], axis=0, ignore_index=True
    )[["chr", "source", "feature", "start", "end", "attribute"]]
    df_final = (
        pd.merge(
            file_gtf,
            df_final_temp,
            on=["chr", "source", "feature", "start", "end"],
            how="left",
        )
        .drop(columns="attribute_x")
        .rename(columns={"attribute_y": "attribute"})
    )
    df_final_all = df_final[
        [
            "chr",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]
    ]

    df_final_all.to_csv(args.output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
