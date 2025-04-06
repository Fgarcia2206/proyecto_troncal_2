# 02_blast.py
import subprocess
from Bio.Blast import NCBIXML
import pandas as pd
from pathlib import Path
import logging

# Configuración
RESULTS_PATH = Path(r"C:\Users\fgarc\OneDrive\Escritorio\Doctorado\Ramos\1° Semestre\Troncal\proyecto_troncal2\results\2025-03-30_FG\identificacion_patogeno")
BLAST_DB = r"C:\blast_db\nr"  # Ajusta esta ruta

def setup_logging():
    RESULTS_PATH.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(RESULTS_PATH / "blast_analysis.log"),
            logging.StreamHandler()
        ]
    )

def run_blast(query: Path, db: str) -> Path:
    """Ejecuta BLASTn y devuelve ruta de resultados"""
    output_file = RESULTS_PATH / "blast_results.xml"
    
    cmd = [
        "blastn",
        "-query", str(query),
        "-db", db,
        "-out", str(output_file),
        "-outfmt", "5",
        "-evalue", "1e-5",
        "-max_target_seqs", "10"
    ]
    
    subprocess.run(cmd, check=True)
    return output_file

def parse_blast_results(xml_path: Path) -> pd.DataFrame:
    """Convierte resultados BLAST a DataFrame"""
    hits = []
    with open(xml_path) as f:
        for record in NCBIXML.parse(f):
            for align in record.alignments:
                hits.append({
                    "query": record.query[:50],
                    "hit": align.hit_def.split("|")[0],
                    "evalue": min(hsp.expect for hsp in align.hsps)
                })
    
    return pd.DataFrame(hits)

if __name__ == "__main__":
    setup_logging()
    try:
        query_seq = RESULTS_PATH / "cleaned_sequences.fasta"
        blast_results = run_blast(query_seq, BLAST_DB)
        
        df = parse_blast_results(blast_results)
        csv_path = RESULTS_PATH / "blast_summary.csv"
        df.to_csv(csv_path, index=False)
        logging.info(f"Resultados BLAST guardados en: {csv_path}")
        
    except Exception as e:
        logging.error(f"Error en análisis BLAST: {e}", exc_info=True)