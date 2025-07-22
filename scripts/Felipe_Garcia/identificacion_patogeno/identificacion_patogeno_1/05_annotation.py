from Bio import Entrez, SeqIO
from pathlib import Path
import logging
import pandas as pd
import time

# Configuración
Entrez.email = "tu.email@institucion.edu"  # Obligatorio para usar NCBI
RESULTS_PATH = Path(r"C:\Users\fgarc\OneDrive\Escritorio\Doctorado\Ramos\1° Semestre\Troncal\proyecto_troncal2\results\2025-03-30_FG\identificacion_patogeno")
NCBI_DB = "nucleotide"
MAX_RETRIES = 3
RETRY_DELAY = 5  # segundos

def setup_logging():
    """Configura el sistema de logging"""
    RESULTS_PATH.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(RESULTS_PATH / "annotation.log"),
            logging.StreamHandler()
        ]
    )

def fetch_genome(accession: str) -> Path:
    """
    Descarga un genoma completo desde NCBI
    
    Args:
        accession: Número de acceso del genoma
        
    Returns:
        Ruta al archivo FASTA descargado
    """
    output_file = RESULTS_PATH / f"{accession}.fasta"
    
    if output_file.exists():
        logging.info(f"Archivo ya existe: {output_file}")
        return output_file
    
    for attempt in range(MAX_RETRIES):
        try:
            logging.info(f"Descargando genoma {accession} (intento {attempt+1})...")
            
            # Obtener el genoma
            handle = Entrez.efetch(
                db=NCBI_DB,
                id=accession,
                rettype="fasta",
                retmode="text"
            )
            
            # Guardar en archivo
            with open(output_file, "w") as f:
                f.write(handle.read())
            
            handle.close()
            logging.info(f"Genoma descargado en: {output_file}")
            return output_file
            
        except Exception as e:
            logging.warning(f"Error en intento {attempt+1}: {e}")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
    
    raise RuntimeError(f"No se pudo descargar el genoma {accession} después de {MAX_RETRIES} intentos")

def annotate_genome(fasta_file: Path) -> pd.DataFrame:
    """
    Realiza una anotación básica del genoma usando NCBI
    
    Args:
        fasta_file: Archivo FASTA con el genoma
        
    Returns:
        DataFrame con las características anotadas
    """
    try:
        # Leer el genoma
        record = next(SeqIO.parse(fasta_file, "fasta"))
        
        # Obtener anotación desde NCBI
        handle = Entrez.efetch(
            db=NCBI_DB,
            id=record.id,
            rettype="gb",
            retmode="text"
        )
        
        gb_record = SeqIO.read(handle, "genbank")
        handle.close()
        
        # Extraer características
        features = []
        for feature in gb_record.features:
            if feature.type not in ["source", "gene"]:
                continue
                
            qualifiers = feature.qualifiers
            start = feature.location.start
            end = feature.location.end
            strand = feature.location.strand
            
            features.append({
                "type": feature.type,
                "start": start,
                "end": end,
                "strand": strand,
                "product": qualifiers.get("product", [""])[0],
                "gene": qualifiers.get("gene", [""])[0],
                "protein_id": qualifiers.get("protein_id", [""])[0],
                "note": qualifiers.get("note", [""])[0]
            })
        
        return pd.DataFrame(features)
        
    except Exception as e:
        logging.error(f"Error en anotación: {e}")
        raise

if __name__ == "__main__":
    setup_logging()
    try:
        # Obtener el mejor hit de BLAST
        blast_summary = pd.read_csv(RESULTS_PATH / "blast_results_summary.csv")
        best_hit = blast_summary.iloc[0]
        
        # Extraer número de acceso (simplificado)
        accession = best_hit["accessions"].split(",")[0].strip()
        logging.info(f"Mejor hit encontrado: {accession}")
        
        # Descargar genoma
        genome_file = fetch_genome(accession)
        
        # Anotar genoma
        annotations = annotate_genome(genome_file)
        
        # Guardar anotaciones
        anno_path = RESULTS_PATH / f"{accession}_annotations.csv"
        annotations.to_csv(anno_path, index=False)
        logging.info(f"Anotaciones guardadas en: {anno_path}")
        
    except Exception as e:
        logging.error(f"Error en anotación genómica: {e}", exc_info=True)