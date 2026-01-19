"""
NCBI RNA-seq FASTQ Downloader - Qt6 GUI Application
Automates downloading RNA-seq FASTQ files from NCBI GEO/SRA with Dagster integration.
"""

import sys
import os
import json
import subprocess
import threading
from dataclasses import dataclass, asdict
from typing import List, Dict, Optional
from datetime import datetime
import urllib.request
import xml.etree.ElementTree as ET

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QTabWidget, QTableWidget, QTableWidgetItem, QPushButton, QLineEdit,
    QLabel, QProgressBar, QTextEdit, QGroupBox, QSplitter, QFrame,
    QHeaderView, QMessageBox, QFileDialog, QComboBox, QSpinBox,
    QCheckBox, QStatusBar, QToolBar, QMenu, QMenuBar, QDialog,
    QFormLayout, QDialogButtonBox, QScrollArea
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer, QSize
from PyQt6.QtGui import QAction, QIcon, QFont, QColor, QPalette, QPixmap

# Optional: matplotlib for visualizations
try:
    import matplotlib
    matplotlib.use('QtAgg')
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# Optional: DE analysis dependencies
try:
    import pandas as pd
    import numpy as np
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    import seaborn as sns
    HAS_DESEQ = True
except ImportError:
    HAS_DESEQ = False


# ============================================================================
# Data Models
# ============================================================================

@dataclass
class SRAExperiment:
    gse_id: str
    bioproject: str = ""
    sra_study: str = ""
    title: str = ""
    organism: str = ""
    experiment_type: str = ""
    platform: str = ""
    num_samples: int = 0
    total_runs: int = 0
    total_bases: int = 0
    run_accessions: List[str] = None
    status: str = "pending"
    
    def __post_init__(self):
        if self.run_accessions is None:
            self.run_accessions = []


# ============================================================================
# Worker Threads
# ============================================================================

class NCBIFetchWorker(QThread):
    """Worker thread for fetching experiment info from NCBI."""
    progress = pyqtSignal(str)
    finished_fetch = pyqtSignal(object)
    error = pyqtSignal(str)
    
    def __init__(self, gse_id: str):
        super().__init__()
        self.gse_id = gse_id
    
    def run(self):
        try:
            self.progress.emit(f"Fetching metadata for {self.gse_id}...")
            experiment = self.fetch_gse_info(self.gse_id)
            self.finished_fetch.emit(experiment)
        except Exception as e:
            self.error.emit(str(e))
    
    def fetch_gse_info(self, gse_id: str) -> SRAExperiment:
        """Fetch GSE experiment info from NCBI."""
        import re
        
        # Fetch from GEO with full view to get BioProject
        geo_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}&targ=self&form=xml&view=full"
        
        try:
            with urllib.request.urlopen(geo_url, timeout=30) as response:
                xml_data = response.read().decode('utf-8')
        except Exception as e:
            raise Exception(f"Failed to fetch GEO data: {e}")
        
        experiment = SRAExperiment(gse_id=gse_id)
        bioproject = None
        
        try:
            root = ET.fromstring(xml_data)
            ns = {'geo': 'http://www.ncbi.nlm.nih.gov/geo/info/MINiML'}
            
            # Parse Series info
            series = root.find('.//geo:Series', ns)
            if series is not None:
                title_elem = series.find('geo:Title', ns)
                if title_elem is not None:
                    experiment.title = title_elem.text or ""
                
                # Get sample count
                samples = series.findall('.//geo:Sample', ns)
                experiment.num_samples = len(samples)
                
                # Find BioProject in Relations
                for relation in series.findall('.//geo:Relation', ns):
                    rel_type = relation.get('type', '')
                    target = relation.get('target', '')
                    if 'BioProject' in rel_type:
                        bioproject = target.split('/')[-1] if '/' in target else target
                        experiment.bioproject = bioproject
            
            # Also search in raw XML for PRJNA (always check, some GSE don't have it in Relations)
            match = re.search(r'(PRJNA\d+)', xml_data)
            if match:
                bioproject = match.group(1)
                experiment.bioproject = bioproject
            
            # If no BioProject found, try fetching from HTML page
            if not bioproject:
                try:
                    html_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}"
                    with urllib.request.urlopen(html_url, timeout=30) as response:
                        html_data = response.read().decode('utf-8')
                    match = re.search(r'(PRJNA\d+)', html_data)
                    if match:
                        bioproject = match.group(1)
                        experiment.bioproject = bioproject
                except Exception:
                    pass
            
            # Try to get platform
            platform = root.find('.//geo:Platform', ns)
            if platform is not None:
                title_elem = platform.find('geo:Title', ns)
                if title_elem is not None:
                    experiment.platform = title_elem.text or ""
        except ET.ParseError:
            pass
        
        # Fetch SRA info using BioProject if available
        self.progress.emit(f"Fetching SRA data for {gse_id}...")
        
        if bioproject:
            sra_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={bioproject}[BioProject]&retmax=1000"
        else:
            sra_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={gse_id}&retmax=1000"
        
        try:
            with urllib.request.urlopen(sra_url, timeout=30) as response:
                sra_xml = response.read().decode('utf-8')
            
            root = ET.fromstring(sra_xml)
            id_list = root.find('IdList')
            if id_list is not None:
                ids = [id_elem.text for id_elem in id_list.findall('Id')]
                
                # Fetch run accessions
                if ids:
                    self.progress.emit(f"Fetching {len(ids)} run accessions...")
                    experiment.run_accessions = self.fetch_run_accessions(ids[:100])
                    experiment.total_runs = len(experiment.run_accessions)
        except Exception:
            pass
        
        experiment.experiment_type = "RNA-Seq"
        experiment.status = "ready" if experiment.run_accessions else "no_sra_data"
        
        return experiment
    
    def fetch_run_accessions(self, sra_ids: List[str]) -> List[str]:
        """Fetch SRR accessions from SRA IDs."""
        if not sra_ids:
            return []
        
        id_str = ",".join(sra_ids[:50])
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={id_str}&rettype=runinfo&retmode=text"
        
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read().decode('utf-8')
            
            runs = []
            lines = data.strip().split('\n')
            if len(lines) > 1:
                header = lines[0].split(',')
                run_idx = header.index('Run') if 'Run' in header else 0
                for line in lines[1:]:
                    parts = line.split(',')
                    if len(parts) > run_idx and parts[run_idx].startswith('SRR'):
                        runs.append(parts[run_idx])
            return runs
        except Exception:
            return []


class DownloadWorker(QThread):
    """Worker thread for downloading SRA/FASTQ files."""
    progress = pyqtSignal(str, int)  # message, percentage
    finished_download = pyqtSignal(str, bool)  # srr_id, success
    log = pyqtSignal(str)
    
    def __init__(self, srr_accessions: List[str], output_dir: str, sra_bin_path: str):
        super().__init__()
        self.srr_accessions = srr_accessions
        self.output_dir = output_dir
        self.sra_bin_path = sra_bin_path
        self.is_cancelled = False
    
    def run(self):
        total = len(self.srr_accessions)
        for i, srr in enumerate(self.srr_accessions):
            if self.is_cancelled:
                break
            
            pct = int((i / total) * 100)
            self.progress.emit(f"Downloading {srr} ({i+1}/{total})", pct)
            
            success = self.download_srr(srr)
            self.finished_download.emit(srr, success)
        
        self.progress.emit("Download complete", 100)
    
    def download_srr(self, srr: str) -> bool:
        """Download a single SRR accession."""
        prefetch = os.path.join(self.sra_bin_path, "prefetch.exe")
        fasterq = os.path.join(self.sra_bin_path, "fasterq-dump.exe")
        
        fastq_dir = os.path.join(self.output_dir, "fastq")
        os.makedirs(fastq_dir, exist_ok=True)
        
        # SRA files go to default NCBI location: ~/ncbi/sra/
        sra_cache = os.path.join(os.path.expanduser("~"), "ncbi", "sra")
        
        try:
            sra_cache = os.path.join(os.path.expanduser("~"), "ncbi", "sra")
            sra_file = os.path.join(sra_cache, f"{srr}.sra")
            
            sra_file_subdir = os.path.join(sra_cache, srr, f"{srr}.sra")
            
            if not os.path.exists(sra_file) and not os.path.exists(sra_file_subdir):
                self.log.emit(f"[{srr}] Running prefetch...")
                cmd = [prefetch, "--progress", "--max-size", "100g", srr]
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
                
                if result.returncode != 0:
                    self.log.emit(f"[{srr}] Prefetch failed: {result.stderr}")
                    return False
                
                if not os.path.exists(sra_file) and not os.path.exists(sra_file_subdir):
                    self.log.emit(f"[{srr}] Prefetch completed but SRA file not found!")
                    return False
            else:
                self.log.emit(f"[{srr}] SRA file already exists, skipping prefetch")
            
            self.log.emit(f"[{srr}] Converting to FASTQ...")
            
            # Convert to FASTQ - use the SRR accession, not the path
            # fasterq-dump works better with accession when file is in default location
            cmd = [fasterq, "--split-files", "--outdir", fastq_dir, "--progress", "--force", srr]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
            
            if result.returncode != 0:
                self.log.emit(f"[{srr}] FASTQ conversion failed: {result.stderr}")
                return False
            
            self.log.emit(f"[{srr}] Completed successfully")
            return True
            
        except subprocess.TimeoutExpired:
            self.log.emit(f"[{srr}] Timeout")
            return False
        except Exception as e:
            self.log.emit(f"[{srr}] Error: {e}")
            return False
    
    def cancel(self):
        self.is_cancelled = True


class DagsterWorker(QThread):
    """Worker thread for running Dagster."""
    output = pyqtSignal(str)
    finished_dagster = pyqtSignal(bool)
    
    SRA_BIN_PATH = r"C:\Users\danie\sratoolkit.3.3.0-win64\bin"
    
    def __init__(self, pipeline_path: str, job_name: str):
        super().__init__()
        self.pipeline_path = pipeline_path
        self.job_name = job_name
        self.process = None
    
    def run(self):
        try:
            cmd = [sys.executable, "-m", "dagster", "job", "execute", "-f", self.pipeline_path, "-j", self.job_name]
            self.output.emit(f"Running: {' '.join(cmd)}")
            
            env = os.environ.copy()
            env["PATH"] = self.SRA_BIN_PATH + os.pathsep + env.get("PATH", "")
            env["FASTQ_OUTPUT_DIR"] = os.path.join(os.path.expanduser("~"), "ncbi_rnaseq_data", "fastq")
            
            self.process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                env=env
            )
            
            for line in self.process.stdout:
                self.output.emit(line.strip())
            
            self.process.wait()
            self.finished_dagster.emit(self.process.returncode == 0)
            
        except Exception as e:
            self.output.emit(f"Error: {e}")
            self.finished_dagster.emit(False)
    
    def stop(self):
        if self.process:
            self.process.terminate()


class SalmonWorker(QThread):
    """Worker thread for running Salmon quantification via WSL on Windows."""
    progress = pyqtSignal(str, int)
    log = pyqtSignal(str)
    finished_quant = pyqtSignal(bool, str)
    
    TRANSCRIPTOME_URLS = {
        "Human (GRCh38)": "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "Mouse (GRCm39)": "https://ftp.ensembl.org/pub/release-110/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
    }
    
    USE_WSL = True  # Use WSL for Salmon on Windows
    
    def __init__(self, samples: dict, output_dir: str, organism: str):
        super().__init__()
        self.samples = samples
        self.output_dir = output_dir
        self.organism = organism
        self.salmon_dir = os.path.join(output_dir, "salmon")
        self.index_dir = os.path.join(self.salmon_dir, "index", organism.split()[0].lower())
    
    def win_to_wsl_path(self, win_path: str) -> str:
        """Convert Windows path to WSL path."""
        path = win_path.replace("\\", "/")
        if len(path) >= 2 and path[1] == ":":
            drive = path[0].lower()
            path = f"/mnt/{drive}{path[2:]}"
        return path
    
    def run_salmon_cmd(self, cmd: list, timeout: int = 3600) -> subprocess.CompletedProcess:
        """Run salmon command, using WSL if on Windows."""
        if self.USE_WSL and sys.platform == "win32":
            wsl_cmd = ["wsl", "-d", "Ubuntu-24.04", "--", "salmon"] + cmd[1:]
            return subprocess.run(wsl_cmd, capture_output=True, text=True, timeout=timeout)
        else:
            return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    
    def run(self):
        try:
            os.makedirs(self.salmon_dir, exist_ok=True)
            
            if not self.check_index():
                self.progress.emit("Downloading transcriptome...", 5)
                if not self.download_transcriptome():
                    self.finished_quant.emit(False, "")
                    return
                
                self.progress.emit("Building Salmon index...", 15)
                if not self.build_index():
                    self.finished_quant.emit(False, "")
                    return
            
            quant_dir = os.path.join(self.salmon_dir, "quant")
            os.makedirs(quant_dir, exist_ok=True)
            
            total = len(self.samples)
            for i, (sample, files) in enumerate(self.samples.items()):
                pct = 25 + int((i / total) * 65)
                self.progress.emit(f"Quantifying {sample} ({i+1}/{total})", pct)
                
                if not self.quantify_sample(sample, files, quant_dir):
                    self.log.emit(f"Warning: Failed to quantify {sample}")
            
            self.progress.emit("Combining results...", 95)
            count_matrix_path = self.combine_counts(quant_dir)
            
            self.progress.emit("Complete!", 100)
            self.finished_quant.emit(True, count_matrix_path)
            
        except Exception as e:
            self.log.emit(f"Error: {e}")
            import traceback
            self.log.emit(traceback.format_exc())
            self.finished_quant.emit(False, "")
    
    def check_index(self) -> bool:
        """Check if Salmon index exists."""
        return os.path.exists(os.path.join(self.index_dir, "complete_ref_lens.bin"))
    
    def download_transcriptome(self) -> bool:
        """Download reference transcriptome."""
        url = self.TRANSCRIPTOME_URLS.get(self.organism)
        if not url:
            self.log.emit(f"Unknown organism: {self.organism}")
            return False
        
        transcriptome_dir = os.path.join(self.salmon_dir, "transcriptomes")
        os.makedirs(transcriptome_dir, exist_ok=True)
        
        filename = os.path.basename(url)
        self.transcriptome_path = os.path.join(transcriptome_dir, filename)
        
        if os.path.exists(self.transcriptome_path):
            self.log.emit(f"Transcriptome already downloaded: {filename}")
            return True
        
        self.log.emit(f"Downloading {filename}...")
        try:
            urllib.request.urlretrieve(url, self.transcriptome_path)
            self.log.emit("Download complete")
            return True
        except Exception as e:
            self.log.emit(f"Download failed: {e}")
            return False
    
    def build_index(self) -> bool:
        """Build Salmon index from transcriptome."""
        os.makedirs(self.index_dir, exist_ok=True)
        
        if self.USE_WSL and sys.platform == "win32":
            transcriptome = self.win_to_wsl_path(self.transcriptome_path)
            index_dir = self.win_to_wsl_path(self.index_dir)
        else:
            transcriptome = self.transcriptome_path
            index_dir = self.index_dir
        
        cmd = [
            "salmon", "index",
            "-t", transcriptome,
            "-i", index_dir,
            "-k", "31"
        ]
        self.log.emit(f"Building index: {' '.join(cmd)}")
        
        try:
            result = self.run_salmon_cmd(cmd, timeout=3600)
            if result.returncode == 0:
                self.log.emit("Index built successfully")
                return True
            else:
                self.log.emit(f"Index build failed: {result.stderr}")
                return False
        except Exception as e:
            self.log.emit(f"Index build error: {e}")
            return False
    
    def quantify_sample(self, sample: str, files: dict, quant_dir: str) -> bool:
        """Run Salmon quant on a single sample."""
        if not files.get("r1"):
            self.log.emit(f"✗ {sample}: No R1 file found, skipping")
            return False
        
        sample_out = os.path.join(quant_dir, sample)
        
        if self.USE_WSL and sys.platform == "win32":
            index_dir = self.win_to_wsl_path(self.index_dir)
            r1 = self.win_to_wsl_path(files["r1"])
            r2 = self.win_to_wsl_path(files["r2"]) if files.get("r2") else None
            out_dir = self.win_to_wsl_path(sample_out)
        else:
            index_dir = self.index_dir
            r1 = files["r1"]
            r2 = files.get("r2")
            out_dir = sample_out
        
        if r2:
            cmd = [
                "salmon", "quant",
                "-i", index_dir,
                "-l", "A",
                "-1", r1,
                "-2", r2,
                "-o", out_dir,
                "-p", "4",
                "--validateMappings"
            ]
        else:
            cmd = [
                "salmon", "quant",
                "-i", index_dir,
                "-l", "A",
                "-r", r1,
                "-o", out_dir,
                "-p", "4",
                "--validateMappings"
            ]
        
        self.log.emit(f"Quantifying {sample}...")
        
        try:
            result = self.run_salmon_cmd(cmd, timeout=7200)
            if result.returncode == 0:
                self.log.emit(f"✓ {sample} complete")
                return True
            else:
                self.log.emit(f"✗ {sample} failed: {result.stderr[:200]}")
                return False
        except Exception as e:
            self.log.emit(f"✗ {sample} error: {e}")
            return False
    
    def combine_counts(self, quant_dir: str) -> str:
        """Combine Salmon quant.sf files into a count matrix."""
        import glob
        
        count_data = {}
        samples = []
        
        for sample_dir in glob.glob(os.path.join(quant_dir, "*")):
            if not os.path.isdir(sample_dir):
                continue
            
            quant_file = os.path.join(sample_dir, "quant.sf")
            if not os.path.exists(quant_file):
                continue
            
            sample = os.path.basename(sample_dir)
            samples.append(sample)
            
            with open(quant_file, 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        gene = parts[0].split('.')[0]
                        count = float(parts[4])
                        
                        if gene not in count_data:
                            count_data[gene] = {}
                        count_data[gene][sample] = int(count)
        
        output_path = os.path.join(self.salmon_dir, "count_matrix.csv")
        
        with open(output_path, 'w') as f:
            f.write("Gene," + ",".join(samples) + "\n")
            for gene in sorted(count_data.keys()):
                counts = [str(count_data[gene].get(s, 0)) for s in samples]
                f.write(f"{gene},{','.join(counts)}\n")
        
        self.log.emit(f"Count matrix saved: {len(count_data)} genes, {len(samples)} samples")
        return output_path


class GEOSupplementaryWorker(QThread):
    """Worker thread for fetching GEO supplementary files list."""
    finished_fetch = pyqtSignal(list)
    error = pyqtSignal(str)
    
    def __init__(self, gse_id: str):
        super().__init__()
        self.gse_id = gse_id
    
    def run(self):
        try:
            files = self.fetch_supplementary_files()
            self.finished_fetch.emit(files)
        except Exception as e:
            self.error.emit(str(e))
            self.finished_fetch.emit([])
    
    def fetch_supplementary_files(self) -> list:
        """Fetch list of supplementary files from GEO."""
        import re
        
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={self.gse_id}"
        
        with urllib.request.urlopen(url, timeout=30) as response:
            html = response.read().decode('utf-8')
        
        files = []
        
        ftp_pattern = r'href="(ftp://ftp\.ncbi\.nlm\.nih\.gov/geo/series/[^"]+/suppl/([^"]+))"'
        matches = re.findall(ftp_pattern, html)
        
        for url, filename in matches:
            http_url = url.replace("ftp://", "https://")
            files.append({
                "name": filename,
                "url": http_url,
                "size": "Unknown"
            })
        
        http_pattern = r'href="(https://[^"]*geo[^"]*suppl[^"]*)"[^>]*>([^<]+)<'
        http_matches = re.findall(http_pattern, html)
        
        for url, filename in http_matches:
            if filename not in [f["name"] for f in files]:
                files.append({
                    "name": filename.strip(),
                    "url": url,
                    "size": "Unknown"
                })
        
        return files


# ============================================================================
# Visualization Widget
# ============================================================================

class VisualizationWidget(QWidget):
    """Widget for displaying experiment visualizations."""
    
    def __init__(self):
        super().__init__()
        self.experiments: List[SRAExperiment] = []
        self.init_ui()
    
    def init_ui(self):
        layout = QVBoxLayout(self)
        
        if not HAS_MATPLOTLIB:
            label = QLabel("Install matplotlib for visualizations:\npip install matplotlib")
            label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(label)
            return
        
        # Chart selector
        control_layout = QHBoxLayout()
        control_layout.addWidget(QLabel("Chart Type:"))
        
        self.chart_combo = QComboBox()
        self.chart_combo.addItems([
            "Runs per Experiment",
            "Samples per Experiment", 
            "Download Status",
            "Data Summary"
        ])
        self.chart_combo.currentIndexChanged.connect(self.update_chart)
        control_layout.addWidget(self.chart_combo)
        
        self.refresh_btn = QPushButton("Refresh")
        self.refresh_btn.clicked.connect(self.update_chart)
        control_layout.addWidget(self.refresh_btn)
        control_layout.addStretch()
        
        layout.addLayout(control_layout)
        
        # Matplotlib canvas
        self.figure = Figure(figsize=(10, 6), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)
        
        # Stats summary
        self.stats_label = QLabel()
        self.stats_label.setStyleSheet("font-size: 12px; padding: 10px; background: #2d2d2d; border-radius: 5px;")
        layout.addWidget(self.stats_label)
    
    def set_experiments(self, experiments: List[SRAExperiment]):
        self.experiments = experiments
        self.update_chart()
    
    def update_chart(self):
        if not HAS_MATPLOTLIB or not self.experiments:
            return
        
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        
        chart_type = self.chart_combo.currentIndex()
        
        if chart_type == 0:  # Runs per Experiment
            self.plot_runs_per_experiment(ax)
        elif chart_type == 1:  # Samples per Experiment
            self.plot_samples_per_experiment(ax)
        elif chart_type == 2:  # Download Status
            self.plot_download_status(ax)
        elif chart_type == 3:  # Data Summary
            self.plot_data_summary(ax)
        
        self.figure.tight_layout()
        self.canvas.draw()
        self.update_stats()
    
    def plot_runs_per_experiment(self, ax):
        gse_ids = [exp.gse_id for exp in self.experiments]
        runs = [exp.total_runs for exp in self.experiments]
        
        colors = ['#4ade80' if exp.status == 'ready' else '#f87171' for exp in self.experiments]
        bars = ax.bar(gse_ids, runs, color=colors, edgecolor='white', linewidth=0.5)
        
        ax.set_xlabel('GSE Experiment', fontsize=10)
        ax.set_ylabel('Number of Runs', fontsize=10)
        ax.set_title('SRA Runs per GSE Experiment', fontsize=12, fontweight='bold')
        ax.tick_params(axis='x', rotation=45)
        
        for bar, run in zip(bars, runs):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                   str(run), ha='center', va='bottom', fontsize=9)
        
        ax.set_facecolor('#1a1a2e')
        self.figure.patch.set_facecolor('#16213e')
        ax.tick_params(colors='white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')
        for spine in ax.spines.values():
            spine.set_color('white')
    
    def plot_samples_per_experiment(self, ax):
        gse_ids = [exp.gse_id for exp in self.experiments]
        samples = [exp.num_samples for exp in self.experiments]
        
        colors = plt.cm.viridis([i/len(gse_ids) for i in range(len(gse_ids))])
        ax.barh(gse_ids, samples, color=colors, edgecolor='white', linewidth=0.5)
        
        ax.set_xlabel('Number of Samples', fontsize=10)
        ax.set_ylabel('GSE Experiment', fontsize=10)
        ax.set_title('Samples per GSE Experiment', fontsize=12, fontweight='bold')
        
        ax.set_facecolor('#1a1a2e')
        self.figure.patch.set_facecolor('#16213e')
        ax.tick_params(colors='white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')
        for spine in ax.spines.values():
            spine.set_color('white')
    
    def plot_download_status(self, ax):
        status_counts = {}
        for exp in self.experiments:
            status_counts[exp.status] = status_counts.get(exp.status, 0) + 1
        
        labels = list(status_counts.keys())
        sizes = list(status_counts.values())
        colors = {'ready': '#4ade80', 'pending': '#fbbf24', 'downloading': '#60a5fa', 
                  'completed': '#22c55e', 'error': '#ef4444', 'no_sra_data': '#6b7280'}
        pie_colors = [colors.get(s, '#888888') for s in labels]
        
        wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%',
                                          colors=pie_colors, explode=[0.02]*len(labels))
        ax.set_title('Experiment Status Distribution', fontsize=12, fontweight='bold', color='white')
        
        for text in texts + autotexts:
            text.set_color('white')
        
        self.figure.patch.set_facecolor('#16213e')
    
    def plot_data_summary(self, ax):
        total_runs = sum(exp.total_runs for exp in self.experiments)
        total_samples = sum(exp.num_samples for exp in self.experiments)
        ready = sum(1 for exp in self.experiments if exp.status == 'ready')
        
        categories = ['Experiments', 'Total Runs', 'Total Samples', 'Ready']
        values = [len(self.experiments), total_runs, total_samples, ready]
        
        colors = ['#60a5fa', '#4ade80', '#fbbf24', '#22c55e']
        bars = ax.bar(categories, values, color=colors, edgecolor='white', linewidth=0.5)
        
        ax.set_title('RNA-seq Data Summary', fontsize=12, fontweight='bold')
        
        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                   str(val), ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        ax.set_facecolor('#1a1a2e')
        self.figure.patch.set_facecolor('#16213e')
        ax.tick_params(colors='white')
        ax.title.set_color('white')
        for spine in ax.spines.values():
            spine.set_color('white')
    
    def update_stats(self):
        if not self.experiments:
            self.stats_label.setText("No experiments loaded")
            return
        
        total_runs = sum(exp.total_runs for exp in self.experiments)
        total_samples = sum(exp.num_samples for exp in self.experiments)
        ready = sum(1 for exp in self.experiments if exp.status == 'ready')
        
        stats = f"""
        <b>Summary Statistics:</b><br>
        📊 Total Experiments: {len(self.experiments)}<br>
        🧬 Total SRA Runs: {total_runs}<br>
        🔬 Total Samples: {total_samples}<br>
        ✅ Ready for Download: {ready}<br>
        ⏳ Pending: {len(self.experiments) - ready}
        """
        self.stats_label.setText(stats)


# ============================================================================
# Main Window
# ============================================================================

class NCBIRNASeqGUI(QMainWindow):
    """Main application window."""
    
    def __init__(self):
        super().__init__()
        self.experiments: List[SRAExperiment] = []
        self.workers = []
        self.sra_bin_path = r"C:\Users\danie\sratoolkit.3.3.0-win64\bin"
        self.output_dir = os.path.join(os.path.expanduser("~"), "ncbi_rnaseq_data")
        self.pipeline_path = os.path.join(os.path.dirname(__file__), "ncbi_fastq_pipeline.py")
        
        self.init_ui()
        self.apply_dark_theme()
        self.load_default_experiments()
    
    def init_ui(self):
        self.setWindowTitle("NCBI RNA-seq FASTQ Downloader")
        self.setMinimumSize(1200, 800)
        
        # Menu bar
        self.create_menu_bar()
        
        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        
        # Tabs
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)
        
        # Tab 1: Experiments
        self.create_experiments_tab()
        
        # Tab 2: Download
        self.create_download_tab()
        
        # Tab 3: Visualization
        self.create_visualization_tab()
        
        # Tab 4: Dagster
        self.create_dagster_tab()
        
        # Tab 5: Processing
        self.create_processing_tab()
        
        # Tab 6: DE Analysis
        self.create_de_analysis_tab()
        
        # Status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")
    
    def create_menu_bar(self):
        menubar = self.menuBar()
        
        # File menu
        file_menu = menubar.addMenu("File")
        
        export_action = QAction("Export Experiments", self)
        export_action.triggered.connect(self.export_experiments)
        file_menu.addAction(export_action)
        
        import_action = QAction("Import Experiments", self)
        import_action.triggered.connect(self.import_experiments)
        file_menu.addAction(import_action)
        
        file_menu.addSeparator()
        
        exit_action = QAction("Exit", self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)
        
        # Settings menu
        settings_menu = menubar.addMenu("Settings")
        
        sra_path_action = QAction("Set SRA Toolkit Path", self)
        sra_path_action.triggered.connect(self.set_sra_path)
        settings_menu.addAction(sra_path_action)
        
        output_dir_action = QAction("Set Output Directory", self)
        output_dir_action.triggered.connect(self.set_output_dir)
        settings_menu.addAction(output_dir_action)
        
        # Help menu
        help_menu = menubar.addMenu("Help")
        
        about_action = QAction("About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)
    
    def create_experiments_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Input section
        input_group = QGroupBox("Add GSE Experiment")
        input_layout = QHBoxLayout(input_group)
        
        input_layout.addWidget(QLabel("GSE ID:"))
        self.gse_input = QLineEdit()
        self.gse_input.setPlaceholderText("e.g., GSE192721")
        self.gse_input.returnPressed.connect(self.add_experiment)
        input_layout.addWidget(self.gse_input)
        
        self.add_btn = QPushButton("Add")
        self.add_btn.clicked.connect(self.add_experiment)
        input_layout.addWidget(self.add_btn)
        
        self.batch_btn = QPushButton("Add Multiple")
        self.batch_btn.clicked.connect(self.add_batch_experiments)
        input_layout.addWidget(self.batch_btn)
        
        layout.addWidget(input_group)
        
        # Experiments table
        self.exp_table = QTableWidget()
        self.exp_table.setColumnCount(8)
        self.exp_table.setHorizontalHeaderLabels([
            "GSE ID", "Title", "BioProject", "Platform", 
            "Samples", "Runs", "Status", "Actions"
        ])
        header = self.exp_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(5, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(6, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(7, QHeaderView.ResizeMode.Fixed)
        self.exp_table.setColumnWidth(7, 100)
        self.exp_table.setAlternatingRowColors(True)
        layout.addWidget(self.exp_table)
        
        # Buttons
        btn_layout = QHBoxLayout()
        
        self.fetch_all_btn = QPushButton("Fetch All Info")
        self.fetch_all_btn.clicked.connect(self.fetch_all_experiments)
        btn_layout.addWidget(self.fetch_all_btn)
        
        self.clear_btn = QPushButton("Clear All")
        self.clear_btn.clicked.connect(self.clear_experiments)
        btn_layout.addWidget(self.clear_btn)
        
        btn_layout.addStretch()
        
        self.exp_count_label = QLabel("Experiments: 0 | Runs: 0")
        btn_layout.addWidget(self.exp_count_label)
        
        layout.addLayout(btn_layout)
        
        self.tabs.addTab(tab, "📋 Experiments")
    
    def create_download_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Settings
        settings_group = QGroupBox("Download Settings")
        settings_layout = QFormLayout(settings_group)
        
        self.output_dir_edit = QLineEdit(self.output_dir)
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.browse_output_dir)
        output_layout = QHBoxLayout()
        output_layout.addWidget(self.output_dir_edit)
        output_layout.addWidget(browse_btn)
        settings_layout.addRow("Output Directory:", output_layout)
        
        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 8)
        self.threads_spin.setValue(2)
        settings_layout.addRow("Parallel Downloads:", self.threads_spin)
        
        self.max_runs_spin = QSpinBox()
        self.max_runs_spin.setRange(0, 1000)
        self.max_runs_spin.setValue(2)
        self.max_runs_spin.setSpecialValueText("All runs")
        settings_layout.addRow("Max Runs (0 = all):", self.max_runs_spin)
        
        self.skip_existing_check = QCheckBox("Skip already downloaded")
        self.skip_existing_check.setChecked(True)
        settings_layout.addRow("", self.skip_existing_check)
        
        layout.addWidget(settings_group)
        
        # Progress
        progress_group = QGroupBox("Download Progress")
        progress_layout = QVBoxLayout(progress_group)
        
        self.download_progress = QProgressBar()
        progress_layout.addWidget(self.download_progress)
        
        self.current_download_label = QLabel("Ready to download")
        progress_layout.addWidget(self.current_download_label)
        
        layout.addWidget(progress_group)
        
        # Log
        log_group = QGroupBox("Download Log")
        log_layout = QVBoxLayout(log_group)
        
        self.download_log = QTextEdit()
        self.download_log.setReadOnly(True)
        self.download_log.setFont(QFont("Consolas", 9))
        log_layout.addWidget(self.download_log)
        
        layout.addWidget(log_group)
        
        # Buttons
        btn_layout = QHBoxLayout()
        
        self.start_download_btn = QPushButton("▶ Start Download")
        self.start_download_btn.clicked.connect(self.start_download)
        btn_layout.addWidget(self.start_download_btn)
        
        self.stop_download_btn = QPushButton("⏹ Stop")
        self.stop_download_btn.clicked.connect(self.stop_download)
        self.stop_download_btn.setEnabled(False)
        btn_layout.addWidget(self.stop_download_btn)
        
        self.show_files_btn = QPushButton("📂 Show Downloaded Files")
        self.show_files_btn.clicked.connect(self.show_downloaded_files)
        btn_layout.addWidget(self.show_files_btn)
        
        btn_layout.addStretch()
        
        layout.addLayout(btn_layout)
        
        self.tabs.addTab(tab, "⬇️ Download")
    
    def create_visualization_tab(self):
        self.viz_widget = VisualizationWidget()
        self.tabs.addTab(self.viz_widget, "📊 Visualization")
    
    def create_dagster_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Info
        info_group = QGroupBox("Dagster Pipeline Integration")
        info_layout = QVBoxLayout(info_group)
        
        info_label = QLabel("""
        <b>Dagster Integration</b><br><br>
        This tab allows you to run the Dagster pipeline for orchestrated data processing.<br>
        The pipeline will download and convert FASTQ files with proper dependency management.
        """)
        info_layout.addWidget(info_label)
        
        self.pipeline_path_edit = QLineEdit(self.pipeline_path)
        pipeline_layout = QHBoxLayout()
        pipeline_layout.addWidget(QLabel("Pipeline:"))
        pipeline_layout.addWidget(self.pipeline_path_edit)
        info_layout.addLayout(pipeline_layout)
        
        layout.addWidget(info_group)
        
        # Job selection
        job_group = QGroupBox("Select Job")
        job_layout = QHBoxLayout(job_group)
        
        self.job_combo = QComboBox()
        self.job_combo.addItems([
            "download_gse192721_fastq",
            "download_gse244263_fastq", 
            "download_gse227835_fastq"
        ])
        job_layout.addWidget(self.job_combo)
        
        self.run_dagster_btn = QPushButton("▶ Run Dagster Job")
        self.run_dagster_btn.clicked.connect(self.run_dagster_job)
        job_layout.addWidget(self.run_dagster_btn)
        
        self.launch_ui_btn = QPushButton("🌐 Launch Dagster UI")
        self.launch_ui_btn.clicked.connect(self.launch_dagster_ui)
        job_layout.addWidget(self.launch_ui_btn)
        
        layout.addWidget(job_group)
        
        # Output
        output_group = QGroupBox("Dagster Output")
        output_layout = QVBoxLayout(output_group)
        
        self.dagster_output = QTextEdit()
        self.dagster_output.setReadOnly(True)
        self.dagster_output.setFont(QFont("Consolas", 9))
        output_layout.addWidget(self.dagster_output)
        
        layout.addWidget(output_group)
        
        self.tabs.addTab(tab, "⚙️ Dagster")
    
    def create_processing_tab(self):
        """Create the FASTQ Processing tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        # Scan for downloaded files section
        scan_group = QGroupBox("Downloaded FASTQ Files")
        scan_layout = QVBoxLayout(scan_group)
        
        scan_controls = QHBoxLayout()
        self.fastq_dir_edit = QLineEdit(os.path.join(self.output_dir, "fastq"))
        scan_controls.addWidget(QLabel("FASTQ Directory:"))
        scan_controls.addWidget(self.fastq_dir_edit)
        
        browse_fastq_btn = QPushButton("Browse")
        browse_fastq_btn.clicked.connect(self.browse_fastq_dir)
        scan_controls.addWidget(browse_fastq_btn)
        
        scan_btn = QPushButton("🔍 Scan for FASTQ Files")
        scan_btn.clicked.connect(self.scan_fastq_files)
        scan_controls.addWidget(scan_btn)
        scan_layout.addLayout(scan_controls)
        
        self.fastq_table = QTableWidget()
        self.fastq_table.setColumnCount(4)
        self.fastq_table.setHorizontalHeaderLabels(["Sample", "R1 File", "R2 File", "Size (MB)"])
        header = self.fastq_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        self.fastq_table.setMaximumHeight(200)
        scan_layout.addWidget(self.fastq_table)
        
        self.fastq_summary_label = QLabel("Click 'Scan' to find downloaded FASTQ files")
        scan_layout.addWidget(self.fastq_summary_label)
        
        layout.addWidget(scan_group)
        
        # Salmon Quantification section
        salmon_group = QGroupBox("Salmon Quantification (FASTQ → Count Matrix)")
        salmon_layout = QVBoxLayout(salmon_group)
        
        salmon_info = QLabel("""
        <b>Salmon pseudo-alignment</b> - Fast and accurate transcript quantification<br>
        • Downloads reference transcriptome (first run only)<br>
        • Builds Salmon index (first run only)<br>
        • Quantifies each sample and generates count matrix
        """)
        salmon_layout.addWidget(salmon_info)
        
        salmon_controls = QHBoxLayout()
        
        salmon_controls.addWidget(QLabel("Organism:"))
        self.organism_combo = QComboBox()
        self.organism_combo.addItems(["Human (GRCh38)", "Mouse (GRCm39)"])
        salmon_controls.addWidget(self.organism_combo)
        
        self.check_salmon_btn = QPushButton("Check Salmon")
        self.check_salmon_btn.clicked.connect(self.check_salmon_installation)
        salmon_controls.addWidget(self.check_salmon_btn)
        
        self.run_salmon_btn = QPushButton("▶ Run Salmon Quantification")
        self.run_salmon_btn.clicked.connect(self.run_salmon_quantification)
        salmon_controls.addWidget(self.run_salmon_btn)
        
        salmon_controls.addStretch()
        salmon_layout.addLayout(salmon_controls)
        
        self.salmon_progress = QProgressBar()
        salmon_layout.addWidget(self.salmon_progress)
        
        self.salmon_status_label = QLabel("Ready")
        salmon_layout.addWidget(self.salmon_status_label)
        
        layout.addWidget(salmon_group)
        
        # GEO Supplementary Files section
        geo_group = QGroupBox("GEO Supplementary Files (Pre-computed Count Matrices)")
        geo_layout = QVBoxLayout(geo_group)
        
        geo_info = QLabel("""
        <b>Download processed data directly from GEO</b><br>
        Many datasets include pre-computed count matrices as supplementary files.
        """)
        geo_layout.addWidget(geo_info)
        
        geo_controls = QHBoxLayout()
        geo_controls.addWidget(QLabel("GSE ID:"))
        self.geo_supp_input = QLineEdit()
        self.geo_supp_input.setPlaceholderText("e.g., GSE192721")
        geo_controls.addWidget(self.geo_supp_input)
        
        self.fetch_supp_btn = QPushButton("🔍 Fetch Supplementary Files")
        self.fetch_supp_btn.clicked.connect(self.fetch_geo_supplementary)
        geo_controls.addWidget(self.fetch_supp_btn)
        geo_controls.addStretch()
        geo_layout.addLayout(geo_controls)
        
        self.supp_files_table = QTableWidget()
        self.supp_files_table.setColumnCount(4)
        self.supp_files_table.setHorizontalHeaderLabels(["Filename", "Size", "Type", "Action"])
        header = self.supp_files_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.Fixed)
        self.supp_files_table.setColumnWidth(3, 120)
        self.supp_files_table.setMaximumHeight(150)
        geo_layout.addWidget(self.supp_files_table)
        
        layout.addWidget(geo_group)
        
        # Processing log
        log_group = QGroupBox("Processing Log")
        log_layout = QVBoxLayout(log_group)
        
        self.processing_log = QTextEdit()
        self.processing_log.setReadOnly(True)
        self.processing_log.setFont(QFont("Consolas", 9))
        self.processing_log.setMaximumHeight(150)
        log_layout.addWidget(self.processing_log)
        
        layout.addWidget(log_group)
        
        # Initialize
        self.discovered_samples = []
        self.geo_supp_files = []
        
        self.tabs.addTab(tab, "🔬 Processing")
    
    def browse_fastq_dir(self):
        """Browse for FASTQ directory."""
        dir_path = QFileDialog.getExistingDirectory(self, "Select FASTQ Directory", self.fastq_dir_edit.text())
        if dir_path:
            self.fastq_dir_edit.setText(dir_path)
    
    def scan_fastq_files(self):
        """Scan directory for downloaded FASTQ files."""
        import glob
        
        fastq_dir = self.fastq_dir_edit.text()
        if not os.path.exists(fastq_dir):
            QMessageBox.warning(self, "Not Found", f"Directory not found: {fastq_dir}")
            return
        
        fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq")) + \
                      glob.glob(os.path.join(fastq_dir, "*.fastq.gz")) + \
                      glob.glob(os.path.join(fastq_dir, "*.fq")) + \
                      glob.glob(os.path.join(fastq_dir, "*.fq.gz"))
        
        samples = {}
        import re
        for f in fastq_files:
            basename = os.path.basename(f)
            
            match = re.match(r'^(SRR\d+)_(\d+)\.', basename)
            if match:
                sample_base = match.group(1)
                read_num = match.group(2)
                if sample_base not in samples:
                    samples[sample_base] = {"r1": None, "r2": None, "r3": None}
                if read_num == "1":
                    samples[sample_base]["r1"] = f
                elif read_num == "2":
                    samples[sample_base]["r2"] = f
                elif read_num == "3":
                    samples[sample_base]["r3"] = f
            elif "_1.fastq" in basename or "_1.fq" in basename or "_R1" in basename:
                sample_name = basename.split("_1")[0].split("_R1")[0]
                if sample_name not in samples:
                    samples[sample_name] = {"r1": None, "r2": None}
                samples[sample_name]["r1"] = f
            elif "_2.fastq" in basename or "_2.fq" in basename or "_R2" in basename:
                sample_name = basename.split("_2")[0].split("_R2")[0]
                if sample_name not in samples:
                    samples[sample_name] = {"r1": None, "r2": None}
                samples[sample_name]["r2"] = f
            else:
                sample_name = basename.replace(".fastq.gz", "").replace(".fastq", "").replace(".fq.gz", "").replace(".fq", "")
                if sample_name not in samples:
                    samples[sample_name] = {"r1": None, "r2": None}
                samples[sample_name]["r1"] = f
        
        for sample in samples:
            files = samples[sample]
            if not files.get("r1") and files.get("r2"):
                files["r1"] = files["r2"]
                files["r2"] = files.get("r3")
        
        self.discovered_samples = samples
        self.fastq_table.setRowCount(len(samples))
        
        total_size = 0
        for i, (sample, files) in enumerate(samples.items()):
            self.fastq_table.setItem(i, 0, QTableWidgetItem(sample))
            
            r1_name = os.path.basename(files["r1"]) if files["r1"] else "-"
            r2_name = os.path.basename(files["r2"]) if files["r2"] else "-"
            self.fastq_table.setItem(i, 1, QTableWidgetItem(r1_name))
            self.fastq_table.setItem(i, 2, QTableWidgetItem(r2_name))
            
            size = 0
            if files["r1"] and os.path.exists(files["r1"]):
                size += os.path.getsize(files["r1"]) / (1024 * 1024)
            if files["r2"] and os.path.exists(files["r2"]):
                size += os.path.getsize(files["r2"]) / (1024 * 1024)
            total_size += size
            self.fastq_table.setItem(i, 3, QTableWidgetItem(f"{size:.1f}"))
        
        self.fastq_summary_label.setText(f"Found {len(samples)} samples, {total_size:.1f} MB total")
        self.processing_log.append(f"Scanned {fastq_dir}: found {len(samples)} samples")
    
    def check_salmon_installation(self):
        """Check if Salmon is installed (via WSL on Windows)."""
        try:
            if sys.platform == "win32":
                result = subprocess.run(["wsl", "-d", "Ubuntu-24.04", "--", "salmon", "--version"], 
                                       capture_output=True, text=True, timeout=15)
            else:
                result = subprocess.run(["salmon", "--version"], capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                version = result.stdout.strip() or result.stderr.strip()
                self.salmon_status_label.setText(f"Salmon found (WSL): {version}")
                self.processing_log.append(f"Salmon installed in WSL: {version}")
                QMessageBox.information(self, "Salmon Found", f"Salmon is installed in WSL:\n{version}")
            else:
                self.show_salmon_install_instructions()
        except FileNotFoundError:
            self.show_salmon_install_instructions()
        except Exception as e:
            self.processing_log.append(f"Error checking Salmon: {e}")
            self.show_salmon_install_instructions()
    
    def show_salmon_install_instructions(self):
        """Show Salmon installation instructions."""
        msg = """
        <b>Salmon not found in WSL!</b><br><br>
        <b>To install (run in WSL terminal):</b><br>
        <code>sudo apt-get update && sudo apt-get install -y salmon</code><br><br>
        <b>Or use conda in WSL:</b><br>
        <code>conda install -c bioconda salmon</code><br><br>
        <b>Alternative:</b> Use the GEO Supplementary Files section to download pre-computed count matrices.
        """
        self.salmon_status_label.setText("Salmon not found in WSL - see instructions")
        QMessageBox.warning(self, "Salmon Not Found", msg)
    
    def run_salmon_quantification(self):
        """Run Salmon quantification on discovered samples."""
        if not self.discovered_samples:
            QMessageBox.warning(self, "No Samples", "Please scan for FASTQ files first.")
            return
        
        try:
            if sys.platform == "win32":
                result = subprocess.run(["wsl", "-d", "Ubuntu-24.04", "--", "salmon", "--version"], 
                                       capture_output=True, text=True, timeout=15)
            else:
                result = subprocess.run(["salmon", "--version"], capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                self.show_salmon_install_instructions()
                return
        except:
            self.show_salmon_install_instructions()
            return
        
        organism = self.organism_combo.currentText()
        self.processing_log.append(f"Starting Salmon quantification for {len(self.discovered_samples)} samples...")
        self.processing_log.append(f"Organism: {organism}")
        
        self.salmon_worker = SalmonWorker(
            samples=self.discovered_samples,
            output_dir=self.output_dir,
            organism=organism
        )
        self.salmon_worker.progress.connect(self.on_salmon_progress)
        self.salmon_worker.log.connect(lambda msg: self.processing_log.append(msg))
        self.salmon_worker.finished_quant.connect(self.on_salmon_finished)
        self.salmon_worker.start()
        
        self.run_salmon_btn.setEnabled(False)
        self.salmon_status_label.setText("Running Salmon quantification...")
    
    def on_salmon_progress(self, message: str, percent: int):
        """Handle Salmon progress updates."""
        self.salmon_progress.setValue(percent)
        self.salmon_status_label.setText(message)
    
    def on_salmon_finished(self, success: bool, count_matrix_path: str):
        """Handle Salmon completion."""
        self.run_salmon_btn.setEnabled(True)
        if success:
            self.salmon_status_label.setText(f"✅ Complete! Count matrix: {count_matrix_path}")
            self.processing_log.append(f"Count matrix saved to: {count_matrix_path}")
            
            reply = QMessageBox.question(
                self, "Quantification Complete",
                f"Count matrix generated!\n\n{count_matrix_path}\n\nLoad into DE Analysis tab?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            if reply == QMessageBox.StandardButton.Yes:
                self.load_count_matrix_from_path(count_matrix_path)
        else:
            self.salmon_status_label.setText("❌ Quantification failed")
    
    def load_count_matrix_from_path(self, path: str):
        """Load count matrix and switch to DE tab."""
        try:
            if HAS_DESEQ:
                sep = '\t' if path.endswith('.tsv') else ','
                self.count_matrix = pd.read_csv(path, index_col=0, sep=sep)
                
                samples = self.count_matrix.columns.tolist()
                n_half = len(samples) // 2
                conditions = ["Control"] * n_half + ["Treatment"] * (len(samples) - n_half)
                
                self.sample_metadata = pd.DataFrame({
                    "sample": samples,
                    "condition": conditions
                }).set_index("sample")
                
                self.de_summary_label.setText(f"Loaded: {len(self.count_matrix)} genes, {len(samples)} samples")
                self.run_deseq_btn.setEnabled(True)
                
                self.tabs.setCurrentIndex(self.tabs.count() - 1)
                
                QMessageBox.information(self, "Data Loaded", 
                    f"Loaded count matrix with {len(self.count_matrix)} genes and {len(samples)} samples.\n"
                    f"Click 'Run DESeq2 Analysis' to analyze.")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to load count matrix: {e}")
    
    def fetch_geo_supplementary(self):
        """Fetch list of supplementary files from GEO."""
        gse_id = self.geo_supp_input.text().strip().upper()
        if not gse_id.startswith("GSE"):
            QMessageBox.warning(self, "Invalid ID", "Please enter a valid GSE ID (e.g., GSE192721)")
            return
        
        self.processing_log.append(f"Fetching supplementary files for {gse_id}...")
        self.fetch_supp_btn.setEnabled(False)
        
        self.geo_fetch_worker = GEOSupplementaryWorker(gse_id)
        self.geo_fetch_worker.finished_fetch.connect(self.on_geo_supp_fetched)
        self.geo_fetch_worker.error.connect(lambda e: self.processing_log.append(f"Error: {e}"))
        self.geo_fetch_worker.start()
    
    def on_geo_supp_fetched(self, files: list):
        """Handle fetched supplementary files."""
        self.fetch_supp_btn.setEnabled(True)
        self.geo_supp_files = files
        
        self.supp_files_table.setRowCount(len(files))
        
        for i, file_info in enumerate(files):
            self.supp_files_table.setItem(i, 0, QTableWidgetItem(file_info["name"]))
            self.supp_files_table.setItem(i, 1, QTableWidgetItem(file_info.get("size", "Unknown")))
            
            file_type = "Unknown"
            name_lower = file_info["name"].lower()
            if "count" in name_lower or "expression" in name_lower or "matrix" in name_lower:
                file_type = "Count Matrix"
            elif name_lower.endswith(".txt.gz") or name_lower.endswith(".csv.gz"):
                file_type = "Data Table"
            elif "raw" in name_lower:
                file_type = "Raw Data"
            self.supp_files_table.setItem(i, 2, QTableWidgetItem(file_type))
            
            download_btn = QPushButton("Download")
            download_btn.setStyleSheet("background-color: #0f3460; color: white; padding: 2px 8px;")
            download_btn.clicked.connect(lambda checked, f=file_info: self.download_geo_file(f))
            self.supp_files_table.setCellWidget(i, 3, download_btn)
        
        self.processing_log.append(f"Found {len(files)} supplementary files")
        
        if not files:
            QMessageBox.information(self, "No Files", "No supplementary files found for this GSE.")
    
    def download_geo_file(self, file_info: dict):
        """Download a GEO supplementary file."""
        url = file_info["url"]
        filename = file_info["name"]
        
        save_path = os.path.join(self.output_dir, "geo_data", filename)
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        
        self.processing_log.append(f"Downloading {filename}...")
        
        try:
            urllib.request.urlretrieve(url, save_path)
            self.processing_log.append(f"✅ Downloaded: {save_path}")
            
            if filename.endswith('.gz'):
                import gzip
                import shutil
                extracted_path = save_path[:-3]
                with gzip.open(save_path, 'rb') as f_in:
                    with open(extracted_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                self.processing_log.append(f"Extracted: {extracted_path}")
                save_path = extracted_path
            
            if "count" in filename.lower() or "matrix" in filename.lower() or "expression" in filename.lower():
                reply = QMessageBox.question(
                    self, "Download Complete",
                    f"Downloaded: {filename}\n\nThis looks like a count matrix. Load into DE Analysis?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                )
                if reply == QMessageBox.StandardButton.Yes:
                    self.load_count_matrix_from_path(save_path)
            else:
                QMessageBox.information(self, "Download Complete", f"Downloaded: {save_path}")
                
        except Exception as e:
            self.processing_log.append(f"❌ Download failed: {e}")
            QMessageBox.warning(self, "Download Failed", str(e))
    
    def create_de_analysis_tab(self):
        """Create the Differential Expression Analysis tab."""
        tab = QWidget()
        layout = QVBoxLayout(tab)
        
        if not HAS_DESEQ or not HAS_MATPLOTLIB:
            missing = []
            if not HAS_DESEQ:
                missing.append("pydeseq2, pandas, numpy, seaborn")
            if not HAS_MATPLOTLIB:
                missing.append("matplotlib")
            label = QLabel(f"Install required packages for DE analysis:\npip install {' '.join(missing)}")
            label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            layout.addWidget(label)
            self.tabs.addTab(tab, "🧬 DE Analysis")
            return
        
        # Controls
        control_group = QGroupBox("Differential Expression Analysis")
        control_layout = QVBoxLayout(control_group)
        
        info_label = QLabel("""
        <b>Differential Expression Analysis with pyDESeq2</b><br><br>
        This demonstrates RNA-seq differential expression analysis:<br>
        • Load or generate count matrix data<br>
        • Run DESeq2 statistical analysis<br>
        • Visualize results with volcano plots and heatmaps<br>
        • Export significant differentially expressed genes
        """)
        control_layout.addWidget(info_label)
        
        btn_row1 = QHBoxLayout()
        
        self.load_counts_btn = QPushButton("📂 Load Count Matrix")
        self.load_counts_btn.clicked.connect(self.load_count_matrix)
        btn_row1.addWidget(self.load_counts_btn)
        
        self.generate_demo_btn = QPushButton("🎲 Generate Demo Data")
        self.generate_demo_btn.clicked.connect(self.generate_demo_data)
        btn_row1.addWidget(self.generate_demo_btn)
        
        self.run_deseq_btn = QPushButton("▶ Run DESeq2 Analysis")
        self.run_deseq_btn.clicked.connect(self.run_deseq_analysis)
        self.run_deseq_btn.setEnabled(False)
        btn_row1.addWidget(self.run_deseq_btn)
        
        btn_row1.addStretch()
        control_layout.addLayout(btn_row1)
        
        layout.addWidget(control_group)
        
        # Results area with splitter
        splitter = QSplitter(Qt.Orientation.Horizontal)
        
        # Left: Results table
        table_group = QGroupBox("DE Results")
        table_layout = QVBoxLayout(table_group)
        
        self.de_table = QTableWidget()
        self.de_table.setColumnCount(6)
        self.de_table.setHorizontalHeaderLabels(["Gene", "BaseMean", "Log2FC", "lfcSE", "P-value", "Adj. P-value"])
        header = self.de_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for i in range(1, 6):
            header.setSectionResizeMode(i, QHeaderView.ResizeMode.ResizeToContents)
        self.de_table.setAlternatingRowColors(True)
        table_layout.addWidget(self.de_table)
        
        # Export button
        export_layout = QHBoxLayout()
        self.export_de_btn = QPushButton("💾 Export Results")
        self.export_de_btn.clicked.connect(self.export_de_results)
        self.export_de_btn.setEnabled(False)
        export_layout.addWidget(self.export_de_btn)
        
        self.de_summary_label = QLabel("Load data to begin analysis")
        export_layout.addWidget(self.de_summary_label)
        export_layout.addStretch()
        table_layout.addLayout(export_layout)
        
        splitter.addWidget(table_group)
        
        # Right: Visualization
        viz_group = QGroupBox("Visualization")
        viz_layout = QVBoxLayout(viz_group)
        
        viz_controls = QHBoxLayout()
        viz_controls.addWidget(QLabel("Plot:"))
        self.de_plot_combo = QComboBox()
        self.de_plot_combo.addItems(["Volcano Plot", "MA Plot", "Heatmap (Top 50)", "PCA"])
        self.de_plot_combo.currentIndexChanged.connect(self.update_de_plot)
        viz_controls.addWidget(self.de_plot_combo)
        viz_controls.addStretch()
        viz_layout.addLayout(viz_controls)
        
        self.de_figure = Figure(figsize=(8, 6), dpi=100)
        self.de_canvas = FigureCanvas(self.de_figure)
        viz_layout.addWidget(self.de_canvas)
        
        splitter.addWidget(viz_group)
        splitter.setSizes([400, 600])
        
        layout.addWidget(splitter)
        
        # Initialize data storage
        self.count_matrix = None
        self.sample_metadata = None
        self.de_results = None
        
        self.tabs.addTab(tab, "🧬 DE Analysis")
    
    def generate_demo_data(self):
        """Generate synthetic RNA-seq count data for demonstration."""
        np.random.seed(42)
        
        n_genes = 1000
        n_samples = 12
        
        gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
        sample_names = [f"Sample_{i+1}" for i in range(n_samples)]
        
        conditions = ["Control"] * 6 + ["Treatment"] * 6
        
        base_expression = np.random.negative_binomial(n=5, p=0.1, size=(n_genes, n_samples))
        
        n_de_genes = 100
        de_indices = np.random.choice(n_genes, n_de_genes, replace=False)
        
        for idx in de_indices[:50]:
            base_expression[idx, 6:] = (base_expression[idx, 6:] * np.random.uniform(2, 5)).astype(int)
        for idx in de_indices[50:]:
            base_expression[idx, 6:] = (base_expression[idx, 6:] * np.random.uniform(0.2, 0.5)).astype(int)
        
        base_expression = np.maximum(base_expression, 0)
        
        self.count_matrix = pd.DataFrame(base_expression, index=gene_names, columns=sample_names)
        self.sample_metadata = pd.DataFrame({
            "sample": sample_names,
            "condition": conditions
        }).set_index("sample")
        
        self.de_summary_label.setText(f"Demo data: {n_genes} genes, {n_samples} samples (6 Control, 6 Treatment)")
        self.run_deseq_btn.setEnabled(True)
        
        QMessageBox.information(self, "Demo Data Generated", 
            f"Generated synthetic count matrix:\n"
            f"• {n_genes} genes\n"
            f"• {n_samples} samples (6 Control, 6 Treatment)\n"
            f"• {n_de_genes} differentially expressed genes embedded\n\n"
            f"Click 'Run DESeq2 Analysis' to analyze.")
    
    def load_count_matrix(self):
        """Load count matrix from CSV file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load Count Matrix", "", "CSV Files (*.csv);;TSV Files (*.tsv);;All Files (*)"
        )
        if not file_path:
            return
        
        try:
            sep = '\t' if file_path.endswith('.tsv') else ','
            self.count_matrix = pd.read_csv(file_path, index_col=0, sep=sep)
            
            samples = self.count_matrix.columns.tolist()
            n_half = len(samples) // 2
            conditions = ["Control"] * n_half + ["Treatment"] * (len(samples) - n_half)
            
            self.sample_metadata = pd.DataFrame({
                "sample": samples,
                "condition": conditions
            }).set_index("sample")
            
            self.de_summary_label.setText(f"Loaded: {len(self.count_matrix)} genes, {len(samples)} samples")
            self.run_deseq_btn.setEnabled(True)
            
            QMessageBox.information(self, "Data Loaded", 
                f"Loaded count matrix:\n"
                f"• {len(self.count_matrix)} genes\n"
                f"• {len(samples)} samples\n\n"
                f"Samples auto-assigned to Control/Treatment groups.\n"
                f"Click 'Run DESeq2 Analysis' to analyze.")
            
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to load file: {e}")
    
    def run_deseq_analysis(self):
        """Run DESeq2 differential expression analysis."""
        if self.count_matrix is None:
            QMessageBox.warning(self, "No Data", "Please load or generate count data first.")
            return
        
        self.de_summary_label.setText("Running DESeq2 analysis...")
        QApplication.processEvents()
        
        try:
            counts = self.count_matrix.T
            
            row_sums = counts.sum(axis=0)
            counts = counts.loc[:, row_sums > 10]
            
            dds = DeseqDataSet(
                counts=counts,
                metadata=self.sample_metadata,
                design="~condition",
                refit_cooks=False
            )
            
            dds.deseq2()
            
            stat_res = DeseqStats(dds, contrast=["condition", "Treatment", "Control"])
            stat_res.run_wald_test()
            stat_res.summary()
            
            self.de_results = stat_res.results_df.copy()
            self.de_results = self.de_results.sort_values("padj")
            
            self.dds = dds
            
            self.populate_de_table()
            self.update_de_plot()
            
            sig_genes = (self.de_results["padj"] < 0.05).sum()
            up_genes = ((self.de_results["padj"] < 0.05) & (self.de_results["log2FoldChange"] > 0)).sum()
            down_genes = ((self.de_results["padj"] < 0.05) & (self.de_results["log2FoldChange"] < 0)).sum()
            
            self.de_summary_label.setText(
                f"Analysis complete: {sig_genes} significant genes (↑{up_genes} ↓{down_genes}) at FDR < 0.05"
            )
            self.export_de_btn.setEnabled(True)
            
        except Exception as e:
            QMessageBox.warning(self, "Analysis Error", f"DESeq2 analysis failed: {e}")
            self.de_summary_label.setText("Analysis failed")
            import traceback
            traceback.print_exc()
    
    def populate_de_table(self):
        """Populate the DE results table."""
        if self.de_results is None:
            return
        
        df = self.de_results.head(500)
        self.de_table.setRowCount(len(df))
        
        for i, (gene, row) in enumerate(df.iterrows()):
            self.de_table.setItem(i, 0, QTableWidgetItem(str(gene)))
            self.de_table.setItem(i, 1, QTableWidgetItem(f"{row.get('baseMean', 0):.2f}"))
            self.de_table.setItem(i, 2, QTableWidgetItem(f"{row.get('log2FoldChange', 0):.3f}"))
            self.de_table.setItem(i, 3, QTableWidgetItem(f"{row.get('lfcSE', 0):.3f}"))
            self.de_table.setItem(i, 4, QTableWidgetItem(f"{row.get('pvalue', 1):.2e}"))
            self.de_table.setItem(i, 5, QTableWidgetItem(f"{row.get('padj', 1):.2e}"))
            
            if row.get('padj', 1) < 0.05:
                color = QColor("#2d5a3d") if row.get('log2FoldChange', 0) > 0 else QColor("#5a2d2d")
                for j in range(6):
                    item = self.de_table.item(i, j)
                    if item:
                        item.setBackground(color)
    
    def update_de_plot(self):
        """Update the DE visualization plot."""
        if self.de_results is None:
            return
        
        self.de_figure.clear()
        ax = self.de_figure.add_subplot(111)
        
        plot_type = self.de_plot_combo.currentText()
        
        if plot_type == "Volcano Plot":
            self.plot_volcano(ax)
        elif plot_type == "MA Plot":
            self.plot_ma(ax)
        elif plot_type == "Heatmap (Top 50)":
            self.plot_heatmap(ax)
        elif plot_type == "PCA":
            self.plot_pca(ax)
        
        self.de_figure.tight_layout()
        self.de_canvas.draw()
    
    def plot_volcano(self, ax):
        """Create volcano plot."""
        df = self.de_results.dropna(subset=["log2FoldChange", "padj"])
        
        df["-log10(padj)"] = -np.log10(df["padj"].clip(lower=1e-300))
        
        not_sig = df["padj"] >= 0.05
        sig_up = (df["padj"] < 0.05) & (df["log2FoldChange"] > 1)
        sig_down = (df["padj"] < 0.05) & (df["log2FoldChange"] < -1)
        sig_mid = (df["padj"] < 0.05) & (df["log2FoldChange"].abs() <= 1)
        
        ax.scatter(df.loc[not_sig, "log2FoldChange"], df.loc[not_sig, "-log10(padj)"], 
                   c="#888888", alpha=0.5, s=10, label="Not significant")
        ax.scatter(df.loc[sig_mid, "log2FoldChange"], df.loc[sig_mid, "-log10(padj)"], 
                   c="#f4a460", alpha=0.7, s=15, label="Sig. (|LFC|≤1)")
        ax.scatter(df.loc[sig_up, "log2FoldChange"], df.loc[sig_up, "-log10(padj)"], 
                   c="#4ade80", alpha=0.8, s=20, label="Up-regulated")
        ax.scatter(df.loc[sig_down, "log2FoldChange"], df.loc[sig_down, "-log10(padj)"], 
                   c="#f87171", alpha=0.8, s=20, label="Down-regulated")
        
        ax.axhline(y=-np.log10(0.05), color="gray", linestyle="--", alpha=0.5)
        ax.axvline(x=1, color="gray", linestyle="--", alpha=0.5)
        ax.axvline(x=-1, color="gray", linestyle="--", alpha=0.5)
        
        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("-Log10(Adjusted P-value)")
        ax.set_title("Volcano Plot")
        ax.legend(loc="upper right", fontsize=8)
        ax.set_facecolor("#1a1a2e")
        self.de_figure.set_facecolor("#1a1a2e")
        ax.tick_params(colors="#eaeaea")
        ax.xaxis.label.set_color("#eaeaea")
        ax.yaxis.label.set_color("#eaeaea")
        ax.title.set_color("#eaeaea")
    
    def plot_ma(self, ax):
        """Create MA plot."""
        df = self.de_results.dropna(subset=["baseMean", "log2FoldChange", "padj"])
        
        df["log10_baseMean"] = np.log10(df["baseMean"].clip(lower=1))
        
        not_sig = df["padj"] >= 0.05
        sig = df["padj"] < 0.05
        
        ax.scatter(df.loc[not_sig, "log10_baseMean"], df.loc[not_sig, "log2FoldChange"],
                   c="#888888", alpha=0.5, s=10, label="Not significant")
        ax.scatter(df.loc[sig, "log10_baseMean"], df.loc[sig, "log2FoldChange"],
                   c="#4ade80", alpha=0.7, s=15, label="Significant")
        
        ax.axhline(y=0, color="gray", linestyle="-", alpha=0.5)
        
        ax.set_xlabel("Log10(Mean Expression)")
        ax.set_ylabel("Log2 Fold Change")
        ax.set_title("MA Plot")
        ax.legend(loc="upper right", fontsize=8)
        ax.set_facecolor("#1a1a2e")
        self.de_figure.set_facecolor("#1a1a2e")
        ax.tick_params(colors="#eaeaea")
        ax.xaxis.label.set_color("#eaeaea")
        ax.yaxis.label.set_color("#eaeaea")
        ax.title.set_color("#eaeaea")
    
    def plot_heatmap(self, ax):
        """Create heatmap of top DE genes."""
        top_genes = self.de_results.dropna(subset=["padj"]).head(50).index.tolist()
        
        if not top_genes or self.count_matrix is None:
            ax.text(0.5, 0.5, "No significant genes", ha="center", va="center", color="#eaeaea")
            return
        
        expr_data = self.count_matrix.loc[self.count_matrix.index.isin(top_genes)]
        expr_data = np.log2(expr_data + 1)
        expr_data = (expr_data.T - expr_data.mean(axis=1)).T
        
        sns.heatmap(expr_data, ax=ax, cmap="RdBu_r", center=0, 
                    xticklabels=True, yticklabels=False, cbar_kws={"shrink": 0.5})
        ax.set_title("Top 50 DE Genes (Z-scored)")
        ax.set_xlabel("Samples")
        ax.set_ylabel("Genes")
        ax.set_facecolor("#1a1a2e")
        self.de_figure.set_facecolor("#1a1a2e")
        ax.tick_params(colors="#eaeaea")
        ax.xaxis.label.set_color("#eaeaea")
        ax.yaxis.label.set_color("#eaeaea")
        ax.title.set_color("#eaeaea")
    
    def plot_pca(self, ax):
        """Create PCA plot."""
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        
        if self.count_matrix is None:
            return
        
        expr_data = np.log2(self.count_matrix + 1).T
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(expr_data)
        
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(scaled_data)
        
        colors = ["#4ade80" if c == "Control" else "#f87171" for c in self.sample_metadata["condition"]]
        
        ax.scatter(pca_result[:, 0], pca_result[:, 1], c=colors, s=100, alpha=0.8)
        
        for i, sample in enumerate(self.count_matrix.columns):
            ax.annotate(sample, (pca_result[i, 0], pca_result[i, 1]), fontsize=8, color="#eaeaea")
        
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
        ax.set_title("PCA Plot")
        ax.set_facecolor("#1a1a2e")
        self.de_figure.set_facecolor("#1a1a2e")
        ax.tick_params(colors="#eaeaea")
        ax.xaxis.label.set_color("#eaeaea")
        ax.yaxis.label.set_color("#eaeaea")
        ax.title.set_color("#eaeaea")
        
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor="#4ade80", label="Control"),
                          Patch(facecolor="#f87171", label="Treatment")]
        ax.legend(handles=legend_elements, loc="upper right")
    
    def export_de_results(self):
        """Export DE results to CSV."""
        if self.de_results is None:
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export DE Results", "de_results.csv", "CSV Files (*.csv)"
        )
        if file_path:
            self.de_results.to_csv(file_path)
            self.status_bar.showMessage(f"Exported to {file_path}", 5000)
    
    def apply_dark_theme(self):
        """Apply dark theme to the application."""
        self.setStyleSheet("""
            QMainWindow, QWidget {
                background-color: #1a1a2e;
                color: #eaeaea;
            }
            QGroupBox {
                border: 1px solid #3d3d5c;
                border-radius: 5px;
                margin-top: 10px;
                padding-top: 10px;
                font-weight: bold;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
            }
            QTableWidget {
                background-color: #16213e;
                alternate-background-color: #1a1a2e;
                gridline-color: #3d3d5c;
                border: 1px solid #3d3d5c;
                border-radius: 5px;
            }
            QTableWidget::item {
                padding: 5px;
            }
            QHeaderView::section {
                background-color: #0f3460;
                padding: 5px;
                border: 1px solid #3d3d5c;
                font-weight: bold;
            }
            QPushButton {
                background-color: #0f3460;
                border: 1px solid #3d3d5c;
                border-radius: 5px;
                padding: 8px 15px;
                min-width: 80px;
            }
            QPushButton:hover {
                background-color: #1a4b8c;
            }
            QPushButton:pressed {
                background-color: #0a2540;
            }
            QPushButton:disabled {
                background-color: #2d2d4a;
                color: #666;
            }
            QLineEdit, QTextEdit, QComboBox, QSpinBox {
                background-color: #16213e;
                border: 1px solid #3d3d5c;
                border-radius: 5px;
                padding: 5px;
            }
            QProgressBar {
                border: 1px solid #3d3d5c;
                border-radius: 5px;
                text-align: center;
                background-color: #16213e;
            }
            QProgressBar::chunk {
                background-color: #4ade80;
                border-radius: 4px;
            }
            QTabWidget::pane {
                border: 1px solid #3d3d5c;
                border-radius: 5px;
            }
            QTabBar::tab {
                background-color: #16213e;
                border: 1px solid #3d3d5c;
                padding: 10px 20px;
                margin-right: 2px;
            }
            QTabBar::tab:selected {
                background-color: #0f3460;
                border-bottom: 2px solid #4ade80;
            }
            QStatusBar {
                background-color: #0f3460;
                border-top: 1px solid #3d3d5c;
            }
            QMenuBar {
                background-color: #0f3460;
            }
            QMenuBar::item:selected {
                background-color: #1a4b8c;
            }
            QMenu {
                background-color: #16213e;
                border: 1px solid #3d3d5c;
            }
            QMenu::item:selected {
                background-color: #0f3460;
            }
        """)
    
    def load_default_experiments(self):
        """Load the default GSE experiments."""
        default_gses = ["GSE192721", "GSE244263", "GSE227835", "GSE218936", "GSE200743"]
        for gse in default_gses:
            exp = SRAExperiment(gse_id=gse)
            self.experiments.append(exp)
        self.update_table()
    
    def add_experiment(self):
        """Add a single experiment."""
        gse_id = self.gse_input.text().strip().upper()
        if not gse_id:
            return
        
        if not gse_id.startswith("GSE"):
            QMessageBox.warning(self, "Invalid ID", "Please enter a valid GSE ID (e.g., GSE192721)")
            return
        
        if any(exp.gse_id == gse_id for exp in self.experiments):
            QMessageBox.information(self, "Duplicate", f"{gse_id} already exists")
            return
        
        exp = SRAExperiment(gse_id=gse_id)
        self.experiments.append(exp)
        self.update_table()
        self.gse_input.clear()
        
        # Auto-fetch info
        self.fetch_experiment_info(exp)
    
    def add_batch_experiments(self):
        """Add multiple experiments at once."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Add Multiple GSE IDs")
        dialog.setMinimumWidth(400)
        
        layout = QVBoxLayout(dialog)
        
        label = QLabel("Enter GSE IDs (one per line):")
        layout.addWidget(label)
        
        text_edit = QTextEdit()
        text_edit.setPlaceholderText("GSE192721\nGSE244263\nGSE227835")
        layout.addWidget(text_edit)
        
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        
        if dialog.exec() == QDialog.DialogCode.Accepted:
            lines = text_edit.toPlainText().strip().split('\n')
            for line in lines:
                gse_id = line.strip().upper()
                if gse_id.startswith("GSE") and not any(exp.gse_id == gse_id for exp in self.experiments):
                    exp = SRAExperiment(gse_id=gse_id)
                    self.experiments.append(exp)
            self.update_table()
    
    def fetch_experiment_info(self, experiment: SRAExperiment):
        """Fetch info for a single experiment."""
        worker = NCBIFetchWorker(experiment.gse_id)
        worker.progress.connect(lambda msg: self.status_bar.showMessage(msg))
        worker.finished_fetch.connect(lambda exp, w=worker: self.on_experiment_fetched(experiment, exp, w))
        worker.error.connect(lambda err, w=worker: self.on_fetch_error(experiment, err, w))
        worker.start()
        self.workers.append(worker)
    
    def on_experiment_fetched(self, old_exp: SRAExperiment, new_exp: SRAExperiment, worker=None):
        """Handle fetched experiment data."""
        try:
            idx = self.experiments.index(old_exp)
            self.experiments[idx] = new_exp
            self.update_table()
            self.viz_widget.set_experiments(self.experiments)
            self.status_bar.showMessage(f"Fetched {new_exp.gse_id}: {new_exp.total_runs} runs", 5000)
        except ValueError:
            pass
        finally:
            if worker and worker in self.workers:
                self.workers.remove(worker)
    
    def on_fetch_error(self, experiment: SRAExperiment, error: str, worker=None):
        """Handle fetch error."""
        experiment.status = "error"
        self.update_table()
        self.status_bar.showMessage(f"Error fetching {experiment.gse_id}: {error}", 5000)
        if worker and worker in self.workers:
            self.workers.remove(worker)
    
    def fetch_all_experiments(self):
        """Fetch info for all experiments with staggered starts."""
        pending = [exp for exp in self.experiments if exp.status == "pending" or exp.total_runs == 0]
        
        for i, exp in enumerate(pending):
            # Stagger requests by 500ms to avoid NCBI rate limiting
            QTimer.singleShot(i * 500, lambda e=exp: self.fetch_experiment_info(e))
        
        if pending:
            self.status_bar.showMessage(f"Fetching {len(pending)} experiments...")
    
    def clear_experiments(self):
        """Clear all experiments."""
        self.experiments.clear()
        self.update_table()
        self.viz_widget.set_experiments([])
    
    def update_table(self):
        """Update the experiments table."""
        self.exp_table.setRowCount(len(self.experiments))
        
        for row, exp in enumerate(self.experiments):
            self.exp_table.setItem(row, 0, QTableWidgetItem(exp.gse_id))
            self.exp_table.setItem(row, 1, QTableWidgetItem(exp.title[:50] + "..." if len(exp.title) > 50 else exp.title if exp.title else "Pending..."))
            self.exp_table.setItem(row, 2, QTableWidgetItem(exp.bioproject if exp.bioproject else "—"))
            self.exp_table.setItem(row, 3, QTableWidgetItem(exp.platform[:30] if exp.platform else "—"))
            self.exp_table.setItem(row, 4, QTableWidgetItem(str(exp.num_samples) if exp.num_samples else "—"))
            self.exp_table.setItem(row, 5, QTableWidgetItem(str(exp.total_runs) if exp.total_runs else "—"))
            
            status_item = QTableWidgetItem(exp.status)
            if exp.status == "ready":
                status_item.setBackground(QColor("#22c55e"))
            elif exp.status == "error" or exp.status == "no_sra_data":
                status_item.setBackground(QColor("#ef4444"))
            elif exp.status == "pending":
                status_item.setBackground(QColor("#f59e0b"))
            self.exp_table.setItem(row, 6, status_item)
            
            # Actions
            action_btn = QPushButton("Fetch")
            action_btn.setStyleSheet("""
                QPushButton {
                    background-color: #2563eb;
                    color: #ffffff;
                    font-weight: bold;
                    border: 1px solid #3b82f6;
                    border-radius: 4px;
                    padding: 5px 10px;
                }
                QPushButton:hover {
                    background-color: #3b82f6;
                }
            """)
            action_btn.clicked.connect(lambda checked, e=exp: self.fetch_experiment_info(e))
            self.exp_table.setCellWidget(row, 7, action_btn)
        
        total_runs = sum(exp.total_runs for exp in self.experiments)
        self.exp_count_label.setText(f"Experiments: {len(self.experiments)} | Runs: {total_runs}")
    
    def browse_output_dir(self):
        """Browse for output directory."""
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory", self.output_dir)
        if dir_path:
            self.output_dir = dir_path
            self.output_dir_edit.setText(dir_path)
    
    def start_download(self):
        """Start downloading FASTQ files."""
        ready_experiments = [exp for exp in self.experiments if exp.status == "ready" and exp.run_accessions]
        
        if not ready_experiments:
            QMessageBox.warning(self, "No Data", "No experiments ready for download. Fetch experiment info first.")
            return
        
        all_runs = []
        for exp in ready_experiments:
            all_runs.extend(exp.run_accessions)
        
        # Apply max runs limit
        max_runs = self.max_runs_spin.value()
        if max_runs > 0 and len(all_runs) > max_runs:
            all_runs = all_runs[:max_runs]
        
        self.output_dir = self.output_dir_edit.text()
        fastq_dir = os.path.join(self.output_dir, "fastq")
        os.makedirs(fastq_dir, exist_ok=True)
        
        # Skip already downloaded if checkbox is checked
        if self.skip_existing_check.isChecked():
            runs_to_download = []
            for srr in all_runs:
                # Check if FASTQ files exist
                fastq_pattern = os.path.join(fastq_dir, f"{srr}*.fastq*")
                import glob
                existing = glob.glob(fastq_pattern)
                if existing:
                    self.download_log.append(f"⏭️ Skipping {srr} (already downloaded)")
                else:
                    runs_to_download.append(srr)
            all_runs = runs_to_download
        
        if not all_runs:
            QMessageBox.information(self, "Complete", "All requested runs have already been downloaded!")
            return
        
        self.download_log.clear()
        self.download_log.append(f"Starting download of {len(all_runs)} runs to {self.output_dir}\n")
        
        self.download_worker = DownloadWorker(all_runs, self.output_dir, self.sra_bin_path)
        self.download_worker.progress.connect(self.on_download_progress)
        self.download_worker.finished_download.connect(self.on_run_downloaded)
        self.download_worker.log.connect(lambda msg: self.download_log.append(msg))
        self.download_worker.finished.connect(self.on_download_finished)
        self.download_worker.start()
        
        self.start_download_btn.setEnabled(False)
        self.stop_download_btn.setEnabled(True)
    
    def on_download_progress(self, message: str, percentage: int):
        """Handle download progress."""
        self.download_progress.setValue(percentage)
        self.current_download_label.setText(message)
    
    def on_run_downloaded(self, srr: str, success: bool):
        """Handle individual run download completion."""
        status = "✅" if success else "❌"
        self.download_log.append(f"{status} {srr}")
    
    def on_download_finished(self):
        """Handle download completion."""
        self.start_download_btn.setEnabled(True)
        self.stop_download_btn.setEnabled(False)
        self.status_bar.showMessage("Download complete", 5000)
    
    def stop_download(self):
        """Stop the current download."""
        if hasattr(self, 'download_worker'):
            self.download_worker.cancel()
            self.download_log.append("\n⚠️ Download cancelled by user")
    
    def show_downloaded_files(self):
        """Show dialog with list of downloaded files."""
        import glob
        
        fastq_dir = os.path.join(self.output_dir_edit.text(), "fastq")
        sra_cache = os.path.join(os.path.expanduser("~"), "ncbi", "sra")
        
        dialog = QDialog(self)
        dialog.setWindowTitle("Downloaded Files")
        dialog.setMinimumSize(600, 400)
        
        layout = QVBoxLayout(dialog)
        
        # FASTQ files
        fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq*"))
        sra_files = glob.glob(os.path.join(sra_cache, "*.sra"))
        
        info = QTextEdit()
        info.setReadOnly(True)
        info.setFont(QFont("Consolas", 10))
        
        text = f"<h3>📁 FASTQ Files ({len(fastq_files)})</h3>"
        text += f"<p>Location: {fastq_dir}</p>"
        
        total_size = 0
        for f in sorted(fastq_files):
            size = os.path.getsize(f)
            total_size += size
            size_gb = size / (1024**3)
            text += f"<br>• {os.path.basename(f)} ({size_gb:.2f} GB)"
        
        text += f"<br><br><b>Total FASTQ: {total_size / (1024**3):.2f} GB</b>"
        
        text += f"<h3>📦 SRA Cache ({len(sra_files)})</h3>"
        text += f"<p>Location: {sra_cache}</p>"
        
        sra_total = 0
        for f in sorted(sra_files):
            size = os.path.getsize(f)
            sra_total += size
            size_gb = size / (1024**3)
            text += f"<br>• {os.path.basename(f)} ({size_gb:.2f} GB)"
        
        text += f"<br><br><b>Total SRA: {sra_total / (1024**3):.2f} GB</b>"
        
        info.setHtml(text)
        layout.addWidget(info)
        
        # Buttons
        btn_layout = QHBoxLayout()
        
        open_fastq_btn = QPushButton("Open FASTQ Folder")
        open_fastq_btn.clicked.connect(lambda: os.startfile(fastq_dir) if os.path.exists(fastq_dir) else None)
        btn_layout.addWidget(open_fastq_btn)
        
        clear_sra_btn = QPushButton("Clear SRA Cache")
        clear_sra_btn.clicked.connect(lambda: self.clear_sra_cache(sra_cache, dialog))
        btn_layout.addWidget(clear_sra_btn)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        btn_layout.addWidget(close_btn)
        
        layout.addLayout(btn_layout)
        dialog.exec()
    
    def clear_sra_cache(self, sra_cache: str, dialog: QDialog):
        """Clear SRA cache to free disk space."""
        import shutil
        reply = QMessageBox.question(
            self, "Clear Cache",
            "This will delete all SRA files to free disk space.\nFASTQ files will be kept.\n\nContinue?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
        )
        if reply == QMessageBox.StandardButton.Yes:
            try:
                for f in os.listdir(sra_cache):
                    if f.endswith('.sra'):
                        os.remove(os.path.join(sra_cache, f))
                QMessageBox.information(self, "Done", "SRA cache cleared!")
                dialog.close()
                self.show_downloaded_files()
            except Exception as e:
                QMessageBox.warning(self, "Error", f"Failed to clear cache: {e}")
    
    def run_dagster_job(self):
        """Run a Dagster job."""
        job_name = self.job_combo.currentText()
        pipeline_path = self.pipeline_path_edit.text()
        
        self.dagster_output.clear()
        self.dagster_output.append(f"Starting Dagster job: {job_name}\n")
        
        self.dagster_worker = DagsterWorker(pipeline_path, job_name)
        self.dagster_worker.output.connect(lambda msg: self.dagster_output.append(msg))
        self.dagster_worker.finished_dagster.connect(self.on_dagster_finished)
        self.dagster_worker.start()
        
        self.run_dagster_btn.setEnabled(False)
    
    def on_dagster_finished(self, success: bool):
        """Handle Dagster job completion."""
        status = "✅ Job completed successfully" if success else "❌ Job failed"
        self.dagster_output.append(f"\n{status}")
        self.run_dagster_btn.setEnabled(True)
    
    def launch_dagster_ui(self):
        """Launch Dagster web UI."""
        pipeline_path = self.pipeline_path_edit.text()
        subprocess.Popen([sys.executable, "-m", "dagster", "dev", "-f", pipeline_path])
        QMessageBox.information(self, "Dagster UI", "Dagster UI launching at http://localhost:3000")
    
    def set_sra_path(self):
        """Set SRA Toolkit path."""
        dir_path = QFileDialog.getExistingDirectory(self, "Select SRA Toolkit bin Directory", self.sra_bin_path)
        if dir_path:
            self.sra_bin_path = dir_path
    
    def set_output_dir(self):
        """Set output directory."""
        self.browse_output_dir()
    
    def export_experiments(self):
        """Export experiments to JSON."""
        file_path, _ = QFileDialog.getSaveFileName(self, "Export Experiments", "", "JSON Files (*.json)")
        if file_path:
            data = [asdict(exp) for exp in self.experiments]
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=2)
            self.status_bar.showMessage(f"Exported to {file_path}", 5000)
    
    def import_experiments(self):
        """Import experiments from JSON."""
        file_path, _ = QFileDialog.getOpenFileName(self, "Import Experiments", "", "JSON Files (*.json)")
        if file_path:
            with open(file_path, 'r') as f:
                data = json.load(f)
            self.experiments = [SRAExperiment(**exp) for exp in data]
            self.update_table()
            self.viz_widget.set_experiments(self.experiments)
            self.status_bar.showMessage(f"Imported {len(self.experiments)} experiments", 5000)
    
    def show_about(self):
        """Show about dialog."""
        QMessageBox.about(self, "About", """
        <h2>NCBI RNA-seq FASTQ Downloader</h2>
        <p>Version 1.0</p>
        <p>A Qt6 GUI application for downloading RNA-seq FASTQ files from NCBI GEO/SRA.</p>
        <p>Features:</p>
        <ul>
            <li>Fetch experiment metadata from NCBI</li>
            <li>Download FASTQ files using SRA Toolkit</li>
            <li>Visualize experiment statistics</li>
            <li>Dagster pipeline integration</li>
        </ul>
        """)


# ============================================================================
# Entry Point
# ============================================================================

def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    
    window = NCBIRNASeqGUI()
    window.show()
    
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
