"""
Dagster pipeline for downloading RNA-seq FASTQ files from NCBI GEO/SRA.

GSE Experiments:
- GSE192721: scRNA-seq from HSCT patients (15 runs)
- GSE244263: scRNA-seq from ALS patients (40 runs)
- GSE227835: scRNA-seq from Myasthenia Gravis patients (40 runs)
- GSE218936: Data not available in SRA
- GSE200743: Data not available in SRA
"""

from dagster import (
    asset,
    Definitions,
    AssetExecutionContext,
    OpExecutionContext,
    Config,
    graph_asset,
    op,
    job,
    Out,
    In,
    graph,
    DynamicOut,
    DynamicOutput,
)
from dataclasses import dataclass
from typing import List, Dict, Optional
import subprocess
import os
import json


# ============================================================================
# Configuration
# ============================================================================

@dataclass
class GSEExperiment:
    """Represents a GEO Series experiment with SRA metadata."""
    gse_id: str
    bioproject: str
    sra_study: str
    description: str
    run_accessions: List[str]


# RNA-seq experiments from NCBI GEO
GSE_EXPERIMENTS: Dict[str, GSEExperiment] = {
    "GSE192721": GSEExperiment(
        gse_id="GSE192721",
        bioproject="PRJNA792834",
        sra_study="SRP352749",
        description="Immune landscape in HSCT patients - scRNA-seq",
        run_accessions=[
            "SRR17352796", "SRR17352797", "SRR17352798", "SRR17352799",
            "SRR17352800", "SRR17352801", "SRR17352802", "SRR17352803",
            "SRR17352804", "SRR17352805", "SRR17352806", "SRR17352807",
            "SRR17352808", "SRR17352809", "SRR17352810"
        ]
    ),
    "GSE244263": GSEExperiment(
        gse_id="GSE244263",
        bioproject="PRJNA1022022",
        sra_study="SRP463782",
        description="scRNA-seq from ALS patients PBMC",
        run_accessions=[
            "SRR26213475", "SRR26213476", "SRR26213477", "SRR26213478",
            "SRR26213479", "SRR26213480", "SRR26213481", "SRR26213482",
            "SRR26213483", "SRR26213484", "SRR26213485", "SRR26213486",
            "SRR26213487", "SRR26213488", "SRR26213489", "SRR26213490",
            "SRR26213491", "SRR26213492", "SRR26213493", "SRR26213494",
            "SRR26213495", "SRR26213496", "SRR26213497", "SRR26213498",
            "SRR26213499", "SRR26213500", "SRR26213501", "SRR26213502",
            "SRR26213503", "SRR26213504", "SRR26213505", "SRR26213506",
            "SRR26213507", "SRR26213508", "SRR26213509", "SRR26213510",
            "SRR26213511", "SRR26213512", "SRR26213513", "SRR26213514"
        ]
    ),
    "GSE227835": GSEExperiment(
        gse_id="GSE227835",
        bioproject="PRJNA964702",
        sra_study="SRP435413",
        description="scRNA-seq from Myasthenia Gravis patients PBMC",
        run_accessions=[
            "SRR24385872", "SRR24385873", "SRR24385874", "SRR24385875",
            "SRR24385876", "SRR24385877", "SRR24385878", "SRR24385879",
            "SRR24385880", "SRR24385881", "SRR24385882", "SRR24385883",
            "SRR24385884", "SRR24385885", "SRR24385886", "SRR24385887",
            "SRR24385888", "SRR24385889", "SRR24385890", "SRR24385891",
            "SRR24385892", "SRR24385893", "SRR24385894", "SRR24385895",
            "SRR24385896", "SRR24385897", "SRR24385898", "SRR24385899",
            "SRR24385900", "SRR24385901", "SRR24385902", "SRR24385903",
            "SRR24385904", "SRR24385905", "SRR24385906", "SRR24385907",
            "SRR24385908", "SRR24385909", "SRR24385910", "SRR24385911"
        ]
    ),
}

# Experiments without SRA data available
UNAVAILABLE_EXPERIMENTS = ["GSE218936", "GSE200743"]


# ============================================================================
# Ops (Operations)
# ============================================================================

@op(out=DynamicOut())
def get_srr_accessions(context: OpExecutionContext, gse_id: str):
    """Get all SRR accessions for a given GSE experiment."""
    if gse_id not in GSE_EXPERIMENTS:
        context.log.warning(f"{gse_id} not found or has no SRA data available")
        return
    
    experiment = GSE_EXPERIMENTS[gse_id]
    context.log.info(f"Processing {gse_id}: {experiment.description}")
    context.log.info(f"Found {len(experiment.run_accessions)} runs")
    
    for srr in experiment.run_accessions:
        yield DynamicOutput(value=srr, mapping_key=srr)


@op
def prefetch_sra(context: OpExecutionContext, srr_accession: str) -> str:
    """Download SRA file using prefetch from SRA Toolkit."""
    sra_bin = r"C:\Users\danie\sratoolkit.3.3.0-win64\bin"
    prefetch_exe = os.path.join(sra_bin, "prefetch.exe")
    
    sra_cache = os.path.join(os.path.expanduser("~"), "ncbi", "sra")
    os.makedirs(sra_cache, exist_ok=True)
    
    sra_file = os.path.join(sra_cache, f"{srr_accession}", f"{srr_accession}.sra")
    lock_file = sra_file + ".lock"
    if os.path.exists(lock_file):
        context.log.info(f"Removing stale lock file: {lock_file}")
        os.remove(lock_file)
    
    cmd = [prefetch_exe, "--progress", "--force", "all", "--max-size", "100g", srr_accession]
    context.log.info(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=7200)
        context.log.info(f"Successfully prefetched {srr_accession}")
        context.log.info(f"stdout: {result.stdout}")
        return srr_accession
    except subprocess.CalledProcessError as e:
        context.log.error(f"Failed to prefetch {srr_accession}: {e.stderr}")
        raise


@op
def convert_to_fastq(context: OpExecutionContext, srr_accession: str) -> Dict[str, str]:
    """Convert SRA file to FASTQ using fasterq-dump."""
    sra_bin = r"C:\Users\danie\sratoolkit.3.3.0-win64\bin"
    fasterq_exe = os.path.join(sra_bin, "fasterq-dump.exe")
    
    output_dir = os.environ.get("FASTQ_OUTPUT_DIR", r"C:\Users\danie\ncbi_rnaseq_data\fastq")
    os.makedirs(output_dir, exist_ok=True)
    
    cmd = [
        fasterq_exe,
        "--split-files",
        "--outdir", output_dir,
        "--threads", "4",
        "--progress",
        "--force",
        srr_accession
    ]
    context.log.info(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=14400)
        context.log.info(f"Successfully converted {srr_accession} to FASTQ")
        context.log.info(f"stdout: {result.stdout}")
        
        fastq_files = {
            "r1": os.path.join(output_dir, f"{srr_accession}_1.fastq"),
            "r2": os.path.join(output_dir, f"{srr_accession}_2.fastq"),
        }
        return fastq_files
    except subprocess.CalledProcessError as e:
        context.log.error(f"Failed to convert {srr_accession}: {e.stderr}")
        raise


@op
def compress_fastq(context: OpExecutionContext, fastq_files: Dict[str, str]) -> Dict[str, str]:
    """Compress FASTQ files with gzip."""
    compressed = {}
    for key, filepath in fastq_files.items():
        if os.path.exists(filepath):
            cmd = ["gzip", "-f", filepath]
            context.log.info(f"Compressing {filepath}")
            subprocess.run(cmd, check=True)
            compressed[key] = f"{filepath}.gz"
    return compressed


@op
def collect_results(context: OpExecutionContext, fastq_results: List[Dict[str, str]]) -> Dict:
    """Collect all FASTQ download results."""
    context.log.info(f"Collected {len(fastq_results)} FASTQ file pairs")
    return {"total_samples": len(fastq_results), "files": fastq_results}


# ============================================================================
# Graph Definition
# ============================================================================

@graph
def download_gse_fastq():
    """
    Graph for downloading FASTQ files from a GSE experiment.
    
    Flow:
    1. Get SRR accessions from GSE ID
    2. Prefetch each SRA file (parallel) - returns SRR accession
    3. Convert SRA to FASTQ using accession (parallel)
    4. Compress FASTQ files (parallel)
    5. Collect results
    """
    srr_accessions = get_srr_accessions()
    
    prefetched = srr_accessions.map(prefetch_sra)
    fastq_files = prefetched.map(convert_to_fastq)
    compressed_files = fastq_files.map(compress_fastq)
    
    return collect_results(compressed_files.collect())


# ============================================================================
# Assets
# ============================================================================

@asset
def gse192721_metadata() -> Dict:
    """Metadata for GSE192721 - HSCT immune landscape scRNA-seq."""
    exp = GSE_EXPERIMENTS["GSE192721"]
    return {
        "gse_id": exp.gse_id,
        "bioproject": exp.bioproject,
        "sra_study": exp.sra_study,
        "description": exp.description,
        "num_runs": len(exp.run_accessions),
        "run_accessions": exp.run_accessions,
    }


@asset
def gse244263_metadata() -> Dict:
    """Metadata for GSE244263 - ALS PBMC scRNA-seq."""
    exp = GSE_EXPERIMENTS["GSE244263"]
    return {
        "gse_id": exp.gse_id,
        "bioproject": exp.bioproject,
        "sra_study": exp.sra_study,
        "description": exp.description,
        "num_runs": len(exp.run_accessions),
        "run_accessions": exp.run_accessions,
    }


@asset
def gse227835_metadata() -> Dict:
    """Metadata for GSE227835 - Myasthenia Gravis PBMC scRNA-seq."""
    exp = GSE_EXPERIMENTS["GSE227835"]
    return {
        "gse_id": exp.gse_id,
        "bioproject": exp.bioproject,
        "sra_study": exp.sra_study,
        "description": exp.description,
        "num_runs": len(exp.run_accessions),
        "run_accessions": exp.run_accessions,
    }


@asset
def all_experiments_summary(
    gse192721_metadata: Dict,
    gse244263_metadata: Dict,
    gse227835_metadata: Dict,
) -> Dict:
    """Summary of all available GSE experiments."""
    return {
        "available_experiments": [
            gse192721_metadata,
            gse244263_metadata,
            gse227835_metadata,
        ],
        "unavailable_experiments": UNAVAILABLE_EXPERIMENTS,
        "total_runs": sum([
            gse192721_metadata["num_runs"],
            gse244263_metadata["num_runs"],
            gse227835_metadata["num_runs"],
        ])
    }


# ============================================================================
# Jobs
# ============================================================================

download_gse192721_job = download_gse_fastq.to_job(
    name="download_gse192721_fastq",
    description="Download FASTQ files for GSE192721 (HSCT immune landscape)",
    config={
        "ops": {
            "get_srr_accessions": {
                "inputs": {"gse_id": "GSE192721"}
            }
        }
    }
)

download_gse244263_job = download_gse_fastq.to_job(
    name="download_gse244263_fastq",
    description="Download FASTQ files for GSE244263 (ALS PBMC scRNA-seq)",
    config={
        "ops": {
            "get_srr_accessions": {
                "inputs": {"gse_id": "GSE244263"}
            }
        }
    }
)

download_gse227835_job = download_gse_fastq.to_job(
    name="download_gse227835_fastq",
    description="Download FASTQ files for GSE227835 (Myasthenia Gravis)",
    config={
        "ops": {
            "get_srr_accessions": {
                "inputs": {"gse_id": "GSE227835"}
            }
        }
    }
)


# ============================================================================
# Dagster Definitions
# ============================================================================

defs = Definitions(
    assets=[
        gse192721_metadata,
        gse244263_metadata,
        gse227835_metadata,
        all_experiments_summary,
    ],
    jobs=[
        download_gse192721_job,
        download_gse244263_job,
        download_gse227835_job,
    ],
)


if __name__ == "__main__":
    # Print summary of experiments
    print("=" * 60)
    print("NCBI RNA-seq FASTQ Download Pipeline")
    print("=" * 60)
    
    for gse_id, exp in GSE_EXPERIMENTS.items():
        print(f"\n{gse_id}:")
        print(f"  BioProject: {exp.bioproject}")
        print(f"  SRA Study:  {exp.sra_study}")
        print(f"  Runs:       {len(exp.run_accessions)}")
        print(f"  Desc:       {exp.description}")
    
    print(f"\nUnavailable: {', '.join(UNAVAILABLE_EXPERIMENTS)}")
    print("\nTo run: dagster dev -f ncbi_fastq_pipeline.py")
