from Bio import SeqIO
import os
from config import SEQUENCES_DIR, OUTPUT_DIR


def extract_drd4_region(fasta_file):
    """Extrai a região de interesse do gene DRD4 (exon 3) das sequências."""
    print("Extraindo a região de interesse do gene DRD4 (exon 3)...")
    extracted_sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Ajustar o intervalo para corresponder ao tamanho das sequências baixadas
        start, end = 0, 250  # Intervalo ajustado
        if len(record.seq) > end:
            record.seq = record.seq[start:end]
            extracted_sequences.append(record)
            print(f"Extraída região de {start} a {end} para sequência: {record.id}")
        else:
            extracted_sequences.append(record)  # Usar a sequência completa se for curta
            print(f"Usando sequência completa para: {record.id} ({len(record.seq)} bp)")

    extracted_file = os.path.join(OUTPUT_DIR, "drd4_extracted.fasta")
    SeqIO.write(extracted_sequences, extracted_file, "fasta")
    print(f"Regiões extraídas salvas no arquivo: {extracted_file}")
    return extracted_file


def extract_coding_region(reference_file, fasta_file):
    """Extrai a região codificante do gene DRD4 com base na sequência de referência."""
    print(
        "Extraindo a região codificante do gene DRD4 com base na sequência de referência..."
    )

    # Definir os limites da região codificante (ajuste conforme necessário)
    coding_start, coding_end = (
        0,
        250,
    )  # Exemplo: ajustar conforme a anotação da referência

    # Carregar a sequência de referência
    reference_record = SeqIO.read(reference_file, "fasta")
    reference_coding_region = reference_record.seq[coding_start:coding_end]
    print(f"Região codificante esperada (posição {coding_start} a {coding_end}):")
    print(reference_coding_region)

    extracted_sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Verificar se a sequência contém a região codificante (parcial ou completa)
        if reference_coding_region in record.seq:
            start = record.seq.find(reference_coding_region)
            end = start + len(reference_coding_region)
            record.seq = record.seq[start:end]
            extracted_sequences.append(record)
            print(f"Extraída região codificante completa para sequência: {record.id}")
        else:
            # Permitir correspondências parciais
            overlap_start = max(0, record.seq.find(reference_coding_region[:50]))
            overlap_end = min(
                len(record.seq), overlap_start + len(reference_coding_region)
            )
            if (
                overlap_end - overlap_start > 50
            ):  # Exigir pelo menos 50 bases de sobreposição
                record.seq = record.seq[overlap_start:overlap_end]
                extracted_sequences.append(record)
                print(
                    f"Extraída região codificante parcial para sequência: {record.id}"
                )
            else:
                print(
                    f"Sequência {record.id} ignorada (não contém a região codificante)."
                )

    extracted_file = os.path.join(OUTPUT_DIR, "drd4_coding_extracted.fasta")
    SeqIO.write(extracted_sequences, extracted_file, "fasta")
    print(f"Regiões codificantes extraídas salvas no arquivo: {extracted_file}")
    return extracted_file
