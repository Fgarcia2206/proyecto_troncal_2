# 03_visualize.py
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import logging

# Configuración
RESULTS_PATH = Path(r"C:\Users\fgarc\OneDrive\Escritorio\Doctorado\Ramos\1° Semestre\Troncal\proyecto_troncal2\results\2025-03-30_FG\identificacion_patogeno")

def setup_logging():
    RESULTS_PATH.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(RESULTS_PATH / "visualization.log"),
            logging.StreamHandler()
        ]
    )

def generate_blast_plot():
    """Genera gráfico de resultados BLAST"""
    csv_path = RESULTS_PATH / "blast_summary.csv"
    df = pd.read_csv(csv_path)
    
    top_hits = df.nsmallest(10, "evalue")
    
    plt.figure(figsize=(10, 6))
    plt.barh(
        top_hits["hit"].str[:40] + "...",
        -top_hits["evalue"].apply(lambda x: 0 if x == 0 else -x),
        color='#1f77b4'
    )
    
    plt.title("Top 10 Hits por E-value")
    plt.xlabel("-log10(E-value)")
    plt.tight_layout()
    
    plot_path = RESULTS_PATH / "blast_results_plot.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    return plot_path

if __name__ == "__main__":
    setup_logging()
    try:
        plot_file = generate_blast_plot()
        logging.info(f"Gráfico guardado en: {plot_file}")
    except Exception as e:
        logging.error(f"Error en visualización: {e}", exc_info=True)