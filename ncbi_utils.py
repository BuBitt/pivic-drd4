from Bio import Entrez, SeqIO
import os
from config import EMAIL

# Configuração do email para o NCBI
Entrez.email = EMAIL


def fetch_reference_sequence(reference_dir=None):
    """
    Baixa a sequência de referência do gene DRD4.

    Args:
        reference_dir: Diretório onde salvar a referência. Se None, importa de config.
    """
    # Importar REFERENCE_DIR dinâmicamente para garantir que foi inicializado
    if reference_dir is None:
        from config import REFERENCE_DIR

        reference_dir = REFERENCE_DIR

    if reference_dir is None:
        raise ValueError(
            "Diretório de referência não inicializado. Execute setup_directories() primeiro."
        )

    print("Baixando a sequência de referência do gene DRD4...")
    search_term = "DRD4[Gene] AND Homo sapiens[Organism] AND RefSeq"
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        raise ValueError(
            "Não foi possível encontrar a sequência de referência do DRD4."
        )

    ref_id = record["IdList"][0]
    handle = Entrez.efetch(db="nucleotide", id=ref_id, rettype="fasta", retmode="text")
    ref_sequence = SeqIO.read(handle, "fasta")
    handle.close()

    ref_file = os.path.join(reference_dir, "drd4_reference.fasta")
    os.makedirs(reference_dir, exist_ok=True)  # Garantir que o diretório existe
    SeqIO.write([ref_sequence], ref_file, "fasta")
    print(f"Sequência de referência salva no arquivo: {ref_file}")
    return ref_file
