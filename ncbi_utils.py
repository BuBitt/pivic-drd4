from Bio import Entrez, SeqIO
import os
from config import REFERENCE_DIR, EMAIL

# Configuração do email para o NCBI
Entrez.email = EMAIL


def fetch_reference_sequence():
    """Baixa a sequência de referência do gene DRD4."""
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

    ref_file = os.path.join(REFERENCE_DIR, "drd4_reference.fasta")
    SeqIO.write([ref_sequence], ref_file, "fasta")
    print(f"Sequência de referência salva no arquivo: {ref_file}")
    return ref_file
