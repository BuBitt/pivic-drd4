from Bio import SeqIO
import os
import subprocess
import logging
from concurrent.futures import ThreadPoolExecutor
from config import ALIGNMENTS_DIR, CLUSTALW_PATH


def run_alignment_in_thread(clustalw_cmd):
    """Executa o alinhamento em uma thread separada."""
    try:
        subprocess.run(clustalw_cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Erro no processo ClustalW: {e}")
        return False


def align_sequences_for_condition(
    reference_file,
    fasta_file,
    condition,
    clustalw_params=None,
    separate_by_condition=True,
):
    """
    Realiza alinhamento múltiplo entre a sequência de referência e sequências de uma condição específica.

    Args:
        reference_file: Arquivo com sequência(s) de referência
        fasta_file: Arquivo com sequências da condição específica
        condition: Nome da condição (ADHD, Autismo, etc)
        clustalw_params: Parâmetros específicos para o ClustalW
        separate_by_condition: Se True, garante que sequências sejam agrupadas por condição
    """
    print(f"Iniciando alinhamento múltiplo para a condição: {condition}...")
    combined_file = os.path.join(ALIGNMENTS_DIR, f"drd4_combined_{condition}.fasta")

    # Verificar arquivos de entrada
    if not os.path.exists(reference_file):
        logging.error(f"Arquivo de referência não encontrado: {reference_file}")
        return None

    if not os.path.exists(fasta_file):
        logging.error(
            f"Arquivo de sequências não encontrado para {condition}: {fasta_file}"
        )
        return None

    if os.path.getsize(fasta_file) == 0:
        logging.error(f"Arquivo de sequências vazio para {condition}: {fasta_file}")
        return None

    # Ler a sequência de referência
    try:
        reference_sequences = list(SeqIO.parse(reference_file, "fasta"))
        if not reference_sequences:
            logging.error(
                f"Nenhuma sequência de referência encontrada em {reference_file}"
            )
            return None

        # Prefixar sequências de referência
        for ref_seq in reference_sequences:
            original_id = ref_seq.id
            ref_seq.id = f"REF_{original_id}"
            ref_seq.description = f"Reference {ref_seq.description}"
    except Exception as e:
        logging.error(f"Erro ao processar sequência de referência: {str(e)}")
        return None

    # Ler as sequências da condição específica
    try:
        condition_sequences = list(SeqIO.parse(fasta_file, "fasta"))
        if not condition_sequences:
            logging.error(
                f"Nenhuma sequência encontrada para {condition} em {fasta_file}"
            )
            return None
    except Exception as e:
        logging.error(f"Erro ao processar sequências de {condition}: {str(e)}")
        return None

    # Combinar sequências e salvar no arquivo
    try:
        combined_sequences = reference_sequences + condition_sequences
        SeqIO.write(combined_sequences, combined_file, "fasta")

        print(
            f"Combinadas {len(combined_sequences)} sequências para alinhamento "
            f"({len(reference_sequences)} referências + {len(condition_sequences)} de {condition})"
        )
    except Exception as e:
        logging.error(f"Erro ao combinar sequências para alinhamento: {str(e)}")
        return None

    # Definir arquivo de saída do alinhamento
    aligned_file = os.path.join(ALIGNMENTS_DIR, f"drd4_aligned_{condition}.aln")

    # Configurar parâmetros do ClustalW
    if clustalw_params is None:
        clustalw_params = {}

    # Parâmetros padrão
    gap_open = clustalw_params.get("GAPOPEN", 10)
    gap_ext = clustalw_params.get("GAPEXT", 0.1)

    # Verificar o executável ClustalW
    if not os.path.exists(CLUSTALW_PATH):
        logging.error(f"Executável ClustalW não encontrado: {CLUSTALW_PATH}")
        return None

    # Construir o comando ClustalW
    clustalw_cmd = [
        CLUSTALW_PATH,
        "-INFILE=" + combined_file,
        "-OUTFILE=" + aligned_file,
        "-OUTPUT=CLUSTAL",
        f"-GAPOPEN={gap_open}",
        f"-GAPEXT={gap_ext}",
    ]

    # Adicionar parâmetros extras
    for param, value in clustalw_params.items():
        if param not in ["GAPOPEN", "GAPEXT"]:
            clustalw_cmd.append(f"-{param}={value}")

    # Executar o alinhamento
    success = False
    try:
        with ThreadPoolExecutor() as executor:
            future = executor.submit(run_alignment_in_thread, clustalw_cmd)
            success = future.result()

        if (
            success
            and os.path.exists(aligned_file)
            and os.path.getsize(aligned_file) > 0
        ):
            print(
                f"Alinhamento concluído para {condition}. Arquivo gerado: {aligned_file}"
            )
            return aligned_file
        else:
            logging.error(
                f"Falha no alinhamento para {condition}. Verificar: {aligned_file}"
            )
            return None
    except Exception as e:
        logging.error(f"Erro durante o alinhamento para {condition}: {str(e)}")
        return None


def align_all_sequences(all_sequences_file, reference_file):
    """
    Realiza alinhamento de todas as sequências DRD4 disponíveis com a referência,
    para comparação geral (sem foco em patologias específicas).
    """
    print(f"Iniciando alinhamento geral de sequências DRD4...")
    combined_file = os.path.join(ALIGNMENTS_DIR, "drd4_combined_all.fasta")

    # Ler a sequência de referência
    reference_sequences = list(SeqIO.parse(reference_file, "fasta"))
    for i, ref_seq in enumerate(reference_sequences):
        original_id = ref_seq.id
        ref_seq.id = f"REF_{original_id}"

    # Ler todas as sequências disponíveis
    all_sequences = list(SeqIO.parse(all_sequences_file, "fasta"))

    # Marcar como sequências genéricas sem condição específica
    for i, seq in enumerate(all_sequences):
        if not (seq.id.startswith("ADHD_") or seq.id.startswith("Autismo_")):
            seq.id = f"Generic_{seq.id}"

    # Combinar para alinhamento geral
    combined_sequences = reference_sequences + all_sequences
    SeqIO.write(combined_sequences, combined_file, "fasta")

    aligned_file = os.path.join(ALIGNMENTS_DIR, "drd4_aligned_all.aln")

    clustalw_cmd = [
        CLUSTALW_PATH,
        "-INFILE=" + combined_file,
        "-OUTFILE=" + aligned_file,
        "-OUTPUT=CLUSTAL",
        "-GAPOPEN=10",
        "-GAPEXT=0.1",
    ]

    with ThreadPoolExecutor() as executor:
        future = executor.submit(run_alignment_in_thread, clustalw_cmd)
        future.result()

    print(f"Alinhamento geral concluído. Arquivo gerado: {aligned_file}")
    return aligned_file


def align_adhd_sequences_only(fasta_file):
    """Realiza alinhamento múltiplo apenas entre as sequências relacionadas ao TDAH."""
    print(
        "Iniciando alinhamento múltiplo apenas entre as sequências relacionadas ao TDAH..."
    )
    aligned_file = os.path.join(ALIGNMENTS_DIR, "drd4_adhd_aligned.aln")

    clustalw_cmd = [
        CLUSTALW_PATH,
        "-INFILE=" + fasta_file,
        "-OUTFILE=" + aligned_file,
        "-OUTPUT=CLUSTAL",
        "-GAPOPEN=10",  # Ajuste para penalizar abertura de gaps
        "-GAPEXT=0.1",  # Penalidade menor para extensão de gaps
    ]

    try:
        subprocess.run(clustalw_cmd, check=True)
        print(f"Alinhamento concluído. Arquivo gerado: {aligned_file}")
    except subprocess.CalledProcessError as e:
        print(f"Erro ao executar ClustalW: {e}")
        raise

    return aligned_file
