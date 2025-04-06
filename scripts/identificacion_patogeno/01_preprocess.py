from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
from pathlib import Path
import logging

# Configuración de rutas
INPUT_PATH = Path(r"C:\Users\fgarc\OneDrive\Escritorio\Doctorado\Ramos\1° Semestre\Troncal\proyecto_troncal2\data\secuencias.fasta.txt")
RESULTS_PATH = Path(r"C:\Users\fgarc\OneDrive\Escritorio\Doctorado\Ramos\1° Semestre\Troncal\proyecto_troncal2\results\2025-03-30_FG\identificacion_patogeno")

def setup_logging():
    """Configura el sistema de logging"""
    RESULTS_PATH.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(RESULTS_PATH / "preprocess.log"),
            logging.StreamHandler()
        ]
    )

def clean_sequences(input_path: Path, output_path: Path) -> None:
    """
    Limpia y filtra secuencias FASTA:
    - Elimina caracteres no estándar
    - Convierte a mayúsculas
    - Filtra secuencias demasiado cortas
    - Separa secuencias concatenadas con 'xxxxx'
    """
    clean_seqs = []
    
    for record in SeqIO.parse(input_path, "fasta"):
        # Procesar cada secuencia
        seq = str(record.seq).upper()
        
        # Separar secuencias concatenadas con 'xxxxx'
        sub_seqs = re.split(r'[XxNn]{4,}', seq)
        
        for i, sub_seq in enumerate(sub_seqs):
            # Limpiar caracteres no estándar
            clean_seq = re.sub(r'[^ATCG]', '', sub_seq)
            
            if len(clean_seq) >= 100:  # Filtro por longitud mínima
                # Crear nuevo registro con ID único
                new_id = f"{record.id}_part{i+1}" if len(sub_seqs) > 1 else record.id
                clean_record = SeqRecord(
                    Seq(clean_seq),
                    id=new_id,
                    description=f"Cleaned version of {record.id}, length={len(clean_seq)}"
                )
                clean_seqs.append(clean_record)
    
    # Guardar resultados en archivo FASTA normal (sin comprimir)
    with open(output_path, 'w') as f:
        SeqIO.write(clean_seqs, f, "fasta")
    
    logging.info(f"Procesadas {len(clean_seqs)} secuencias. Guardadas en: {output_path}")
    logging.info(f"Puedes abrir el archivo en un editor de texto: {output_path}")

if __name__ == "__main__":
    setup_logging()
    try:
        output_file = RESULTS_PATH / "cleaned_sequences.fasta"
        clean_sequences(INPUT_PATH, output_file)
    except Exception as e:
        logging.error(f"Error en preprocesamiento: {e}", exc_info=True)