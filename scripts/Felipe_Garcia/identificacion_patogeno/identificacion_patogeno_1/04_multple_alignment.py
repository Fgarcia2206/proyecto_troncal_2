from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from pathlib import Path
import logging
import subprocess

# Configuración
CLUSTAL_BIN = r"C:\Program Files\Clustal\Omega\clustalo.exe"
RESULTS_PATH = Path(r"C:\Users\fgarc\OneDrive\Escritorio\Doctorado\Ramos\1° Semestre\Troncal\proyecto_troncal2\results\2025-03-30_FG\identificacion_patogeno")

def setup_logging():
    """Configura el sistema de logging"""
    RESULTS_PATH.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(RESULTS_PATH / "alignment.log"),
            logging.StreamHandler()
        ]
    )

def run_clustal(input_fasta: Path) -> Path:
    """
    Ejecuta Clustal Omega para alineamiento múltiple
    
    Args:
        input_fasta: Archivo FASTA con secuencias a alinear
        
    Returns:
        Ruta al archivo de alineamiento en formato CLUSTAL
    """
    output_aln = RESULTS_PATH / "alignment.clustal"
    
    # Configurar línea de comando
    clustal_cmd = ClustalOmegaCommandline(
        cmd=CLUSTAL_BIN,
        infile=str(input_fasta),
        outfile=str(output_aln),
        outfmt="clustal",
        verbose=True,
        auto=True,
        threads=4
    )
    
    logging.info(f"Ejecutando Clustal Omega con: {str(clustal_cmd)}")
    
    try:
        stdout, stderr = clustal_cmd()
        logging.info(f"Clustal Omega completado. Salida:\n{stdout}")
        
        if stderr:
            logging.warning(f"Advertencias de Clustal:\n{stderr}")
            
        return output_aln
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Error en Clustal Omega:\n{e.stderr}")
        raise RuntimeError(f"Falló el alineamiento: {e.stderr}")

def analyze_alignment(aln_file: Path) -> dict:
    """
    Analiza un alineamiento múltiple y calcula métricas básicas
    
    Args:
        aln_file: Archivo de alineamiento
        
    Returns:
        Diccionario con métricas del alineamiento
    """
    try:
        alignment = AlignIO.read(aln_file, "clustal")
        
        metrics = {
            "num_sequences": len(alignment),
            "alignment_length": alignment.get_alignment_length(),
            "average_identity": None,
            "conserved_positions": None
        }
        
        # Calcular identidad promedio (simplificado)
        identities = []
        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            most_common = max(set(column), key=column.count)
            identities.append(sum(1 for base in column if base == most_common) / len(column))
        
        metrics["average_identity"] = sum(identities) / len(identities) * 100
        metrics["conserved_positions"] = sum(1 for x in identities if x == 1.0)
        
        return metrics
        
    except Exception as e:
        logging.error(f"Error al analizar alineamiento: {e}")
        raise

if __name__ == "__main__":
    setup_logging()
    try:
        # Ejecutar Clustal con las secuencias limpias
        input_seqs = RESULTS_PATH / "cleaned_sequences.fasta"
        aln_file = run_clustal(input_seqs)
        
        # Analizar resultados
        metrics = analyze_alignment(aln_file)
        logging.info(f"Métricas del alineamiento:\n{metrics}")
        
    except Exception as e:
        logging.error(f"Error en alineamiento: {e}", exc_info=True)