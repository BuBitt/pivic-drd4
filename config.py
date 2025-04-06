import os
import logging
from datetime import datetime

# Configuração inicial
EMAIL = "seu_email@example.com"  # Substitua pelo seu e-mail para usar a API do NCBI
OUTPUT_DIR = "drd4_analysis"

# Caminho para o executável do ClustalW
CLUSTALW_PATH = "/home/bruno/Downloads/clustalw-2.1/clustalw2"

# Definir os diretórios para a estrutura do projeto
REFERENCE_DIR = os.path.join(OUTPUT_DIR, "reference")
SEQUENCES_DIR = os.path.join(OUTPUT_DIR, "sequences")
ALIGNMENTS_DIR = os.path.join(OUTPUT_DIR, "alignments")
REPORTS_DIR = os.path.join(OUTPUT_DIR, "reports")
LOGS_DIR = os.path.join(OUTPUT_DIR, "logs")
CACHE_DIR = os.path.join(OUTPUT_DIR, "cache")
VISUALIZATION_DIR = os.path.join(OUTPUT_DIR, "visualization")

# Parâmetros configuráveis
MAX_SEQUENCES = 200  # Número máximo de sequências a baixar
USE_CACHE = True  # Se True, usa cache para evitar refazer downloads
CACHE_TTL = 86400 * 7  # Tempo de vida do cache em segundos (1 semana)
MIN_SEQ_LENGTH = 100  # Comprimento mínimo de sequência para consideração

# Termos de busca expandidos para melhorar a cobertura
SEARCH_TERMS = {
    "ADHD": ["ADHD"],
    "Autismo": ["autism"],
}

# Parâmetros específicos para ClustalW por condição
CLUSTALW_PARAMS = {
    "ADHD": {
        "GAPOPEN": 8,  # Ajustado para ser menos penalizado em abrir gaps
        "GAPEXT": 0.15,  # Valor médio para extensão de gaps
        "PWDNAMATRIX": "IUB",  # Matriz de DNA para sequências mais próximas
        "DNAMATRIX": "IUB",
    },
    "Autismo": {
        "GAPOPEN": 6,  # Menos penalidade para abertura de gaps
        "GAPEXT": 0.3,  # Mais penalidade para extensão de gaps
        "PWDNAMATRIX": "IUB",  # Matriz de substituição de DNA
        "DNAMATRIX": "IUB",  # Matriz para DNA, não proteína
    },
    "Comorbidade": {
        "GAPOPEN": 7,  # Valor intermediário entre ADHD e Autismo
        "GAPEXT": 0.2,  # Valor intermediário entre ADHD e Autismo
        "PWDNAMATRIX": "IUB",
        "DNAMATRIX": "IUB",
    },
    "ADHD_Variantes": {
        "GAPOPEN": 5,  # Mais flexível para variantes divergentes
        "GAPEXT": 0.2,  # Valor intermediário
        "PWDNAMATRIX": "IUB",
        "DNAMATRIX": "IUB",
        "KIMURA": "ON",  # Usar correção de Kimura para distâncias
    },
}

# Configurações adicionais para buscas
SEARCH_CONFIG = {
    "Autismo": {
        "min_identity": 70,  # Identidade mínima para considerar sequências relevantes
        "expand_results": True,  # Buscar sequências relacionadas às encontradas
        "include_variants": True,  # Incluir variantes do gene
        "use_protein_search": True,  # Buscar também sequências de proteínas
    },
    "ADHD": {
        "min_identity": 80,
        "expand_results": False,
        "include_variants": False,
        "use_protein_search": False,
    },
}

# Configurações de detecção de polimorfismos
POLYMORPHISM_CONFIG = {
    "ADHD": {
        "min_vntr_length": 24,  # Menor tamanho mínimo para detectar VNTRs menores
        "max_gap_ratio": 0.7,  # Maior tolerância a gaps
        "gap_threshold": 0.25,  # Limiar de gaps para considerar um VNTR
        "snp_min_coverage": 0.3,  # Cobertura mínima para considerar um SNP
    },
    "Autismo": {
        "min_vntr_length": 24,  # Menor para detectar VNTRs menores
        "max_gap_ratio": 0.65,  # Ajustado para autismo
        "gap_threshold": 0.2,  # Limiar de gaps mais restritivo
        "snp_min_coverage": 0.4,  # Cobertura mínima mais alta para autismo
    },
    "Comorbidade": {
        "min_vntr_length": 20,  # Menor ainda para melhor detecção em casos comórbidos
        "max_gap_ratio": 0.68,  # Intermediário entre ADHD e Autismo
        "gap_threshold": 0.22,  # Intermediário
        "snp_min_coverage": 0.35,  # Intermediário
    },
    "ADHD_Variantes": {
        "min_vntr_length": 20,  # Mais permissivo para variantes divergentes
        "max_gap_ratio": 0.8,  # Maior tolerância a gaps
        "gap_threshold": 0.3,  # Maior limiar para considerar gaps
        "snp_min_coverage": 0.25,  # Menor cobertura devido à maior divergência
    },
}


def setup_directories():
    """Cria os diretórios necessários para o projeto."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(REFERENCE_DIR, exist_ok=True)
    os.makedirs(SEQUENCES_DIR, exist_ok=True)
    os.makedirs(ALIGNMENTS_DIR, exist_ok=True)
    os.makedirs(REPORTS_DIR, exist_ok=True)
    os.makedirs(LOGS_DIR, exist_ok=True)
    os.makedirs(CACHE_DIR, exist_ok=True)
    os.makedirs(VISUALIZATION_DIR, exist_ok=True)


def setup_logging(verbose=False):
    """Configura o sistema de logging."""
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
