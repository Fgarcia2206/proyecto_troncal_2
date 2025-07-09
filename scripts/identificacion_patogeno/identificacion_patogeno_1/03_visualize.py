import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo import draw
from Bio import AlignIO, Phylo
import pylab

# Configuración
RESULTS_PATH = Path(r"C:\Users\fgarc\OneDrive\Escritorio\Doctorado\Ramos\1° Semestre\Troncal\proyecto_troncal2\results\2025-03-30_FG\identificacion_patogeno")
PLOT_STYLE = "seaborn-v0_8-poster"
COLOR_PALETTE = "viridis"

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

def plot_blast_results():
    """Genera múltiples visualizaciones de los resultados BLAST"""
    # Configurar estilo
    plt.style.use(PLOT_STYLE)
    sns.set_palette(COLOR_PALETTE)
    
    # Cargar datos
    csv_path = RESULTS_PATH / "blast_results_detailed.csv"
    df = pd.read_csv(csv_path)
    
    if df.empty:
        logging.warning("No hay datos para visualizar")
        return []
    
    # Filtrar los mejores hits por query
    top_hits = df.loc[df.groupby("query_id")["evalue"].idxmin()]
    
    # 1. Gráfico de los mejores hits por E-value
    plt.figure(figsize=(12, 8))
    top_organisms = top_hits.nsmallest(20, "evalue")
    # Calcular -log10(evalue) para mejor visualización
    top_organisms["-log_evalue"] = -np.log10(top_organisms["evalue"].astype(float) + 1e-300)  # Evitar log(0)
    sns.barplot(
        data=top_organisms,
        y="organism",
        x="-log_evalue",
        estimator=np.median,
        errorbar=None
    )
    plt.title("Top 20 Organismos por E-value (mejores hits)")
    plt.xlabel("-log10(E-value)")
    plt.ylabel("Organismo")
    plt.tight_layout()
    plot1 = RESULTS_PATH / "blast_top_organisms.png"
    plt.savefig(plot1, dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Heatmap de identidad de secuencia
    plt.figure(figsize=(12, 8))
    pivot_data = top_hits.pivot_table(
        index="organism",
        columns="query_id",
        values="percent_identity",
        aggfunc="mean"
    ).fillna(0)
    sns.heatmap(
        pivot_data,
        cmap="YlOrRd",
        annot=True,
        fmt=".1f",
        linewidths=.5
    )
    plt.title("Porcentaje de Identidad por Organismo y Query")
    plt.tight_layout()
    plot2 = RESULTS_PATH / "blast_identity_heatmap.png"
    plt.savefig(plot2, dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Gráfico de distribución de E-values
    plt.figure(figsize=(10, 6))
    sns.histplot(
        data=top_hits,
        x="evalue",
        bins=50,
        log_scale=True
    )
    plt.title("Distribución de E-values de los mejores hits")
    plt.xlabel("E-value (log scale)")
    plt.ylabel("Frecuencia")
    plt.tight_layout()
    plot3 = RESULTS_PATH / "blast_evalue_distribution.png"
    plt.savefig(plot3, dpi=300, bbox_inches='tight')
    plt.close()
    
    logging.info(f"Gráficos guardados en:\n- {plot1}\n- {plot2}\n- {plot3}")
    return [plot1, plot2, plot3]

def generate_phylogenetic_tree():
    """Genera un árbol filogenético basado en los alineamientos"""
    try:
        # Esto requeriría un archivo de alineamiento múltiple (ej. de CLUSTAL)
        aln_file = RESULTS_PATH / "alignment.clustal"
        
        if not aln_file.exists():
            logging.warning("Archivo de alineamiento no encontrado. Ejecuta CLUSTAL primero.")
            return None
            
        alignment = AlignIO.read(aln_file, "clustal")
        
        # Calcular matriz de distancia
        calculator = DistanceTreeConstructor()
        dm = calculator.get_distance(alignment)
        
        # Construir árbol (método UPGMA)
        tree = calculator.upgma(dm)
        
        # Dibujar árbol
        plt.figure(figsize=(15, 10))
        Phylo.draw(tree, do_show=False)
        plt.title("Árbol Filogenético basado en alineamiento múltiple")
        
        plot_path = RESULTS_PATH / "phylogenetic_tree.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"Árbol filogenético guardado en: {plot_path}")
        return plot_path
        
    except Exception as e:
        logging.error(f"Error al generar árbol filogenético: {e}")
        return None

if __name__ == "__main__":
    setup_logging()
    try:
        plots = plot_blast_results()
        tree_plot = generate_phylogenetic_tree()
        
        if tree_plot:
            plots.append(tree_plot)
            
    except Exception as e:
        logging.error(f"Error en visualización: {e}", exc_info=True)