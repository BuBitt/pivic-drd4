import os
import logging
import shutil
from datetime import datetime

# Configuração inicial
EMAIL = "seu_email@example.com"  # Substitua pelo seu e-mail para usar a API do NCBI

# Expandir o ~ no caminho do ClustalW
CLUSTALW_PATH = os.path.expanduser("~/Downloads/clustalw-2.1/clustalw2")

# Diretório base para o projeto
BASE_DIR = "drd4_analysis"

# Diretório compartilhado para sequências (independente da ferramenta)
SHARED_SEQUENCES_DIR = os.path.join(BASE_DIR, "sequences")
SHARED_REFERENCE_DIR = os.path.join(BASE_DIR, "reference")
SHARED_CACHE_DIR = os.path.join(BASE_DIR, "cache")

# Variáveis globais para armazenar os diretórios específicos da ferramenta
OUTPUT_DIR = None
REFERENCE_DIR = None
SEQUENCES_DIR = None
ALIGNMENTS_DIR = None
REPORTS_DIR = None
LOGS_DIR = None
CACHE_DIR = None
VISUALIZATION_DIR = None

# Parâmetros configuráveis
MAX_SEQUENCES = 200
USE_CACHE = True
CACHE_TTL = 86400 * 7
MIN_SEQ_LENGTH = 100

# Termos de busca simplificados
SEARCH_TERMS = {
    "ADHD": ["ADHD"],
    "Autismo": ["autism"],
}

# Parâmetros específicos para ClustalW por condição
CLUSTALW_PARAMS = {
    "ADHD": {
        "GAPOPEN": 8,
        "GAPEXT": 0.15,
        "PWDNAMATRIX": "IUB",
        "DNAMATRIX": "IUB",
    },
    "Autismo": {
        "GAPOPEN": 6,
        "GAPEXT": 0.3,
        "PWDNAMATRIX": "IUB",
        "DNAMATRIX": "IUB",
    },
    "ADHD_Variantes": {
        "GAPOPEN": 5,
        "GAPEXT": 0.2,
        "PWDNAMATRIX": "IUB",
        "DNAMATRIX": "IUB",
        "KIMURA": "ON",
    },
}

# Configurações de detecção de polimorfismos
POLYMORPHISM_CONFIG = {
    "ADHD": {
        "min_vntr_length": 24,
        "max_gap_ratio": 0.7,
        "gap_threshold": 0.25,
        "snp_min_coverage": 0.3,
    },
    "Autismo": {
        "min_vntr_length": 24,
        "max_gap_ratio": 0.65,
        "gap_threshold": 0.2,
        "snp_min_coverage": 0.4,
    },
    "ADHD_Variantes": {
        "min_vntr_length": 20,
        "max_gap_ratio": 0.8,
        "gap_threshold": 0.3,
        "snp_min_coverage": 0.25,
    },
}

# Parâmetros específicos para MAFFT por condição
MAFFT_PARAMS = {
    "ADHD": {
        "algorithm": "--auto",  # Seleção automática de algoritmo
        "maxiterate": 2,  # Número de iterações de refinamento
        "reorder": True,  # Reordenar sequências
        "output_format": "clustal",  # Formato de saída
    },
    "Autismo": {
        "algorithm": "--genafpair",  # Algoritmo de precisão para sequências divergentes
        "maxiterate": 2,  # Número de iterações de refinamento
        "reorder": True,  # Reordenar sequências
        "output_format": "clustal",  # Formato de saída
    },
    "ADHD_Variantes": {
        "algorithm": "--localpair",  # Algoritmo de alta precisão para variantes
        "maxiterate": 3,  # Mais iterações para variantes
        "reorder": True,  # Reordenar sequências
        "output_format": "clustal",  # Formato de saída
        "adjustdirection": True,  # Ajustar direção da sequência (para variantes)
    },
}

# Descrições dos algoritmos para reports
ALIGNMENT_DESCRIPTIONS = {
    "clustalw": "ClustalW: Algoritmo de alinhamento progressivo que constrói o alinhamento "
    "múltiplo incrementalmente, começando com o alinhamento das sequências mais "
    "similares. Usa uma matriz de substituição para calcular as pontuações.",
    "mafft": "MAFFT: Algoritmo de alinhamento múltiplo rápido que utiliza transformadas de "
    "Fourier para identificar homologia entre sequências. Oferece vários modos de operação "
    "que balanceiam velocidade e precisão.",
    "mafft_auto": "MAFFT (--auto): Seleção automática de algoritmo baseada nas características das sequências.",
    "mafft_genafpair": "MAFFT (G-INS-i): Método iterativo global para sequências com regiões homólogas longas.",
    "mafft_localpair": "MAFFT (L-INS-i): Método de alinhamento iterativo local, ideal para sequências "
    "com múltiplos domínios conservados e regiões altamente variáveis.",
}


def setup_directories(output_dir=None):
    """Cria os diretórios necessários para o projeto."""
    global \
        OUTPUT_DIR, \
        REFERENCE_DIR, \
        SEQUENCES_DIR, \
        ALIGNMENTS_DIR, \
        REPORTS_DIR, \
        LOGS_DIR, \
        CACHE_DIR, \
        VISUALIZATION_DIR

    # Criar diretórios compartilhados (independente da ferramenta)
    os.makedirs(SHARED_SEQUENCES_DIR, exist_ok=True)
    os.makedirs(SHARED_REFERENCE_DIR, exist_ok=True)
    os.makedirs(SHARED_CACHE_DIR, exist_ok=True)

    # Atualizar o diretório de saída específico da ferramenta
    if output_dir:
        OUTPUT_DIR = output_dir
    else:
        OUTPUT_DIR = os.path.join(BASE_DIR, "results")

    # Definir caminhos dos diretórios específicos da ferramenta
    REFERENCE_DIR = SHARED_REFERENCE_DIR  # Usar diretório compartilhado
    SEQUENCES_DIR = SHARED_SEQUENCES_DIR  # Usar diretório compartilhado
    ALIGNMENTS_DIR = os.path.join(OUTPUT_DIR, "alignments")
    REPORTS_DIR = os.path.join(OUTPUT_DIR, "reports")
    LOGS_DIR = os.path.join(OUTPUT_DIR, "logs")
    CACHE_DIR = SHARED_CACHE_DIR  # Usar diretório compartilhado
    VISUALIZATION_DIR = os.path.join(OUTPUT_DIR, "visualization")

    # Criar os diretórios específicos da ferramenta
    for directory in [
        OUTPUT_DIR,
        ALIGNMENTS_DIR,
        REPORTS_DIR,
        LOGS_DIR,
        VISUALIZATION_DIR,
    ]:
        os.makedirs(directory, exist_ok=True)


def setup_logging(verbose=False, output_dir=None):
    """Configura o sistema de logging."""
    global LOGS_DIR

    # Atualizar diretório de logs se output_dir foi especificado
    if output_dir:
        LOGS_DIR = os.path.join(output_dir, "logs")
        os.makedirs(LOGS_DIR, exist_ok=True)

    log_file = os.path.join(
        LOGS_DIR, f"drd4_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    )

    # Configuração do logger principal
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),  # Para saída no console também
        ],
    )

    # Reduzir verbosidade de bibliotecas externas
    logging.getLogger("Bio").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)

    logging.info("Sistema de logging configurado")
