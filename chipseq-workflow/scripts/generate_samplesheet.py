from pathlib import Path
import pandas as pd

ROOT = "/users/j/a/jaileru/hepatic-nutrient-regulation"
metadata = pd.read_csv(f"{ROOT}/data/metadata/SraRunTable.csv")

chip = metadata.loc[
    metadata["Assay Type"].eq("ChIP-Seq")
].copy()

chip["diet"] = chip["source_name"].apply(
    lambda x: "CR" if "calorie restricted" in x.lower() else "HFD"
)

chip["antibody"] = chip["chip_antibody"].str.strip().replace({
    "RXRa": "RXRA",
    "PPARa": "PPARA",
})

chip["sample"] = chip["diet"] + "_" + chip["antibody"]

chip = chip.sort_values(["sample", "Run"])

chip["replicate"] = (
    chip.groupby("sample").cumcount() + 1
)

chip["fastq_1"] = chip["Run"].apply(
    lambda run: str(
        (Path(f"{ROOT}/chipseq-workflow/data/fastq") / f"{run}.fastq.gz").resolve()
    )
)

chip["fastq_2"] = ""

# Match each ChIP sample to a diet-specific IgG control.
chip["control"] = chip["diet"] + "_IgG"
chip["control_replicate"] = chip["replicate"]

# IgG samples themselves do not have controls.
is_igg = chip["antibody"].eq("IgG")
chip.loc[is_igg, ["antibody", "control", "control_replicate"]] = pd.NA


columns = [
    "sample", "fastq_1", "fastq_2", "replicate",
    "antibody", "control", "control_replicate"
]

samplesheet = chip[columns]


samplesheet.to_csv(
    f"{ROOT}/chipseq-workflow/samplesheet.csv",
    index=False
)
