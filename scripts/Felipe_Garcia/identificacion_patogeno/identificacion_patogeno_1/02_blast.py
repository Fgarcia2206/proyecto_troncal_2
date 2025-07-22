from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import pandas as pd
from pathlib import Path
import logging
import re
from typing import Optional

# Configuración básica
RESULTS_PATH = Path("results/blast_analysis")
DEFAULT_PARAMS = {
    "program": "blastn",
    "database": "refseq_viral",
    "e_value": 1e-10,
    "hitlist_size": 50,
    "word_size": 11,
    "gapcosts": "5 2"  # Solo parámetros soportados por qblast
}

def setup_logging():
    """Configura el sistema de logging"""
    RESULTS_PATH.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(RESULTS_PATH / "blast_analysis.log"),
            logging.StreamHandler()
        ]
    )

def run_remote_blast(
    query_file: Path,
    program: str = DEFAULT_PARAMS["program"],
    database: str = DEFAULT_PARAMS["database"],
    e_value: float = DEFAULT_PARAMS["e_value"],
    hitlist_size: int = DEFAULT_PARAMS["hitlist_size"],
    word_size: int = DEFAULT_PARAMS["word_size"],
    gapcosts: str = DEFAULT_PARAMS["gapcosts"],
) -> Path:
    """
    Ejecuta BLAST remoto con parámetros personalizables
    
    Args:
        query_file: Archivo FASTA con secuencias
        program: blastn, blastp, etc.
        database: Base de datos a usar
        e_value: Umbral de significancia
        hitlist_size: Número máximo de resultados
        word_size: Tamaño de palabra (7-15 para blastn)
        gapcosts: Costos de gaps ("abrir extender")
        
    Returns:
        Ruta al archivo XML con resultados
    """
    output_file = RESULTS_PATH / "blast_results.xml"
    
    records = list(SeqIO.parse(query_file, "fasta"))
    if not records:
        raise ValueError("No se encontraron secuencias en el archivo de consulta")
    
    logging.info(f"Iniciando BLAST remoto con {len(records)} secuencias")
    logging.info(f"Parámetros: word_size={word_size}, gapcosts={gapcosts}, evalue={e_value}")

    with open(query_file) as f:
        fasta_data = f.read()
    
    try:
        # Solo parámetros soportados por NCBIWWW.qblast()
        blast_params = {
            "program": program,
            "database": database,
            "sequence": fasta_data,
            "expect": e_value,
            "hitlist_size": hitlist_size,
            "format_type": "XML",
            "word_size": word_size,
            "gapcosts": gapcosts,
        }
        
        handle = NCBIWWW.qblast(**blast_params)
        
        with open(output_file, "w") as f:
            blast_results = handle.read()
            f.write(blast_results)
        
        logging.info(f"Resultados guardados en: {output_file}")
        return output_file
        
    except Exception as e:
        logging.error(f"Error en BLAST remoto: {e}")
        raise

def parse_blast_results(xml_path: Path) -> pd.DataFrame:
    """Procesa resultados XML de BLAST"""
    hits = []
    
    try:
        with open(xml_path, 'r') as f:
            blast_records = NCBIXML.parse(f)
            
            for record in blast_records:
                if not record.alignments:
                    continue
                    
                for align in record.alignments:
                    for hsp in align.hsps:
                        hit_def = align.hit_def
                        organism = re.search(r'\[(.*?)\]', hit_def)
                        organism = organism.group(1) if organism else "Desconocido"
                        
                        hits.append({
                            "query_id": record.query,
                            "query_length": record.query_length,
                            "hit_accession": align.accession,
                            "hit_definition": hit_def,
                            "organism": organism,
                            "evalue": hsp.expect,
                            "bit_score": hsp.bits,
                            "alignment_length": hsp.align_length,
                            "percent_identity": hsp.identities / hsp.align_length * 100,
                            "query_start": hsp.query_start,
                            "query_end": hsp.query_end,
                            "hit_start": hsp.sbjct_start,
                            "hit_end": hsp.sbjct_end,
                            "alignment": f"{hsp.query[0:50]}..." if len(hsp.query) > 50 else hsp.query
                        })
    except Exception as e:
        logging.error(f"Error al parsear resultados: {str(e)}")
        raise
    
    df = pd.DataFrame(hits)
    
    if not df.empty:
        df = df.sort_values(["evalue", "bit_score"], ascending=[True, False])
        df = df.drop_duplicates(subset=["query_id", "hit_accession"], keep="first")
    
    return df

def save_results(df: pd.DataFrame, params: dict):
    """Guarda los resultados en archivos CSV"""
    csv_path = RESULTS_PATH / "blast_results_detailed.csv"
    df.to_csv(csv_path, index=False)
    logging.info(f"Resultados detallados guardados en: {csv_path}")
    
    if not df.empty:
        summary = df.groupby("organism").agg({
            "evalue": "min",
            "bit_score": "max",
            "percent_identity": "mean",
            "query_id": "count",
            "hit_accession": lambda x: ", ".join(sorted(set(x)))
        }).sort_values(["evalue", "bit_score"], ascending=[True, False])
        
        summary_path = RESULTS_PATH / "blast_results_summary.csv"
        summary.to_csv(summary_path)
        
        params_path = RESULTS_PATH / "blast_parameters.txt"
        with open(params_path, 'w') as f:
            for key, value in params.items():
                f.write(f"{key}: {value}\n")
        
        logging.info(f"Resumen de resultados guardado en: {summary_path}")
        logging.info("\nMejores hits encontrados:\n" + str(summary.head(10)))
    else:
        logging.warning("No se encontraron hits significativos")

if __name__ == "__main__":
    setup_logging()
    
    try:
        query_seq = Path(r"C:\Users\fgarc\OneDrive\Escritorio\Doctorado\Ramos\1° Semestre\Troncal\proyecto_troncal2\results\2025-03-30_FG\identificacion_patogeno\cleaned_sequences.fasta")
        
        if not query_seq.exists():
            raise FileNotFoundError(f"No se encontró el archivo de consulta: {query_seq}")
        
        # Parámetros personalizados (solo los soportados por qblast)
        custom_params = {
            "program": "blastn",
            "database": "refseq_viral",
            "e_value": 1e-10,
            "hitlist_size": 50,
            "word_size": 11,  # Ajustable entre 7-15
            "gapcosts": "5 2"  # Costos de gaps
        }
        
        blast_results = run_remote_blast(query_seq, **custom_params)
        df = parse_blast_results(blast_results)
        save_results(df, custom_params)
        
    except Exception as e:
        logging.error(f"Error inesperado: {str(e)}", exc_info=True)