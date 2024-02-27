import os
from typing import List, Optional, Union, Literal
import pandas as pd
import scanpy as sc
import anndata as ad

from ..tools._filter import sample_size_filter
from ..utils import get_logger


def load_tcga(
    file_path: Optional[str] = None,
    metadata: Optional[str] = None,
    genes: Optional[Union[List[str], str]] = None,
    sample_type: Optional[str] = None,
    tumor_type: Optional[str] = None,
    filter_samples: Optional[int] = None,
    drop_duplicates: bool = True,
    backed: Optional[Literal["r", "r+"]] = None,
) -> ad.AnnData:
    """Load TCGA dataset.
    Parameters:
        file_path: Path to the TCGA RNAseq file.
        metadata: Path to the TCGA metadata file, downloaded from cBioportal.
        genes: List of genes to be used for the analysis.
        sample_type: Type of samples to be used for the analysis, e.g. 'Normal' - Solid Tissue Normal, 'Tumor' - All non-normal tissue. By default, all samples are used.
        tumor_type: Type of tumor to be used for the analysis, e.g. 'AML' - Acute Myeloid Leukemia. By default, all samples are used.
        filter_samples: Minimal number of samples per project. By default, projects with more than 50 samples are used.
        drop_duplicates: If True, drop duplicate samples. By default, duplicates are dropped and only the first entry is kept.
        backed: If 'r', load AnnData in backed mode instead of fully loading it into memory (memory mode). If you want to modify backed attributes of the AnnData object, you need to choose 'r+'. By default, the AnnData object is fully loaded into memory.

    Returns:
        anndata.AnnData: TCGA dataset.
    """

    dataloader_logger = get_logger("dataloader")

    # check if the file_path is provided, if not, use the TCGA RNAseq file in .datasets folder
    if file_path is None:
        # Navigate to the directory containing the TCGA RNAseq file
        file_path = os.path.join(
            os.path.dirname(__file__),
            "datasets",
            "TCGA_log2TPMplus1_protein_coding_transcripts_20221113.h5ad",
        )
        dataloader_logger.info(
            "No file path provided. Using the TCGA RNAseq file in the .datasets folder."
        )

    # check if backed is provided and if so, raise an unimplemented error
    if backed is not None:
        raise NotImplementedError(
            "Backed mode is not yet implemented. Please use None -- fully loading into the memory."
        )

    adata = sc.read(file_path, cache=True, backed=backed)
    adata.obs["project_id"] = [
        pid.split("TCGA-")[1] for pid in adata.obs["project.project_id"]
    ]

    # drop duplicates
    if drop_duplicates:
        dataloader_logger.info("Dropping duplicates...")
        adata = adata[
            ~adata.obs.duplicated(subset=(["submitter_id", "samples.sample_type"]))
        ]

    # check if the TCGA sample index is in the correct format
    adata.obs_names = [ob.split(".")[0][:-1] for ob in adata.obs_names]
    assert (
        adata.obs_names.duplicated().sum() == 0
    ), "TCGA sample index is not unique. Drop duplicates and try again."

    # Load metadata, if not provided, use the cBioportal metadata in .datasets folder
    if metadata is None:
        # Navigate to the directory containing the cBioportal metadata file
        metadata = os.path.join(
            os.path.dirname(__file__),
            "datasets",
            "combined_study_clinical_data_cBioPortal.tsv",
        )

    # Load metadata, if provided, and filter samples based on sample_type and tumor_type
    metadata = pd.read_csv(metadata, sep="\t", index_col=2)
    metadata.columns = [
        "cbioportal." + "_".join(col.lower().split(" ")) for col in metadata.columns
    ]

    # check if all patient samples in the metadata are present in the adata object
    patient_id_metadata = metadata["cbioportal.patient_id"].unique()
    patient_id_RNAseq = adata.obs["submitter_id"].unique()

    missing_samples = len(set(patient_id_metadata).difference(patient_id_RNAseq))
    assert (
        missing_samples / len(patient_id_RNAseq) < 0.15
    ), "Too many patient samples are missing in the RNAseq data. Check the metadata whether the keys are patient_id."
    dataloader_logger.info(
        f"{missing_samples} patient samples are missing in the RNAseq data."
    )

    # check if sample_type is provided and if so, filter the adata.obs
    if sample_type is not None:
        dataloader_logger.info(f"Take {sample_type} samples only.")
        if sample_type.lower() == "tumor":
            adata = adata[adata.obs["samples.sample_type"] != "Solid Tissue Normal"]
        elif sample_type.lower() == "normal":
            adata = adata[adata.obs["samples.sample_type"] == "Solid Tissue Normal"]

    # check if tumor_type is provided and if so, filter the adata.obs
    if tumor_type is not None:
        dataloader_logger.info(f"Take patients with {tumor_type}")
        adata = adata[adata.obs["project_id"] == tumor_type]

    if genes is not None:
        dataloader_logger.info(f"Take gene expression data of {genes}")
        adata = adata[:, genes]

    # filter out patients with < cutoff samples
    if filter_samples is not None:
        dataloader_logger.info(
            f"Filtering out patients with < {filter_samples} samples..."
        )
        filtered_patients = sample_size_filter(
            adata.obs, groupby="project_id", cutoff=filter_samples
        )
        adata = adata[filtered_patients.index]

    adata.obs = adata.obs.join(metadata, how="left")

    # rename the subtype column: cbioportal.subtype -> subtype
    adata.obs["subtype"] = adata.obs["cbioportal.subtype"]
    # add normal key to the subtype column
    adata.obs.loc[
        adata.obs["samples.sample_type"] == "Solid Tissue Normal", "subtype"
    ] = "Normal"

    # rename the subtype column: subtype(SKCM) -> subtype, as cbioportal.subtype is not always available
    dataloader_logger.debug(
        "Adding a subtype column, SKCM are divided into SKCM-Primary and SKCM-Metastatic and all other subtypes are obtained from cBioportal"
    )

    for idx, row in adata.obs[["subtype(SKCM)"]].dropna().iterrows():
        adata.obs.loc[idx, "subtype"] = row["subtype(SKCM)"]

    return adata
