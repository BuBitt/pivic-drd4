from Bio import SeqIO
import os
import subprocess
import logging
import shutil
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


def check_alignment_tools():
    """Verifica se as ferramentas de alinhamento estão disponíveis."""
    tools_status = {}

    # Verificar MAFFT
    try:
        result = subprocess.run(
            ["mafft", "--help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        tools_status["mafft"] = True
    except (FileNotFoundError, subprocess.SubprocessError):
        tools_status["mafft"] = False
        logging.warning(
            "Ferramenta MAFFT não encontrada no PATH. Por favor, instale-a ou verifique sua instalação."
        )

    # Verificar ClustalW
    clustalw_path = os.path.expanduser(CLUSTALW_PATH)  # Expandir o ~ no caminho

    if os.path.exists(clustalw_path) and os.access(clustalw_path, os.X_OK):
        tools_status["clustalw"] = True
    else:
        # Tente encontrar o comando no PATH
        clustalw_in_path = shutil.which("clustalw") or shutil.which("clustalw2")
        if clustalw_in_path:
            tools_status["clustalw"] = True
            # Não tentamos modificar a variável global aqui
        else:
            tools_status["clustalw"] = False
            logging.warning(
                f"ClustalW não encontrado no caminho configurado: {CLUSTALW_PATH}. "
                f"Por favor, verifique o caminho em config.py ou instale o ClustalW."
            )

    return tools_status


def align_sequences_for_condition(
    reference_file,
    fasta_file,
    condition,
    clustalw_params=None,
    alignment_tool="clustalw",
    separate_by_condition=True,
):
    """
    Realiza o alinhamento de sequências usando ClustalW ou MAFFT.

    Args:
        reference_file: Caminho para o arquivo de sequência de referência.
        fasta_file: Caminho para o arquivo FASTA com sequências a serem alinhadas.
        condition: Nome da condição (ADHD, Autismo, etc).
        clustalw_params: Parâmetros específicos para ClustalW.
        alignment_tool: Ferramenta de alinhamento local ('clustalw' ou 'mafft').
        separate_by_condition: Se True, garante separação por condição.

    Returns:
        Caminho para o arquivo de alinhamento gerado.
    """
    # Usar o diretório de alinhamentos definido em config
    from config import ALIGNMENTS_DIR

    combined_file = os.path.join(ALIGNMENTS_DIR, f"drd4_combined_{condition}.fasta")
    aligned_file = os.path.join(ALIGNMENTS_DIR, f"drd4_aligned_{condition}.aln")

    # Combinar a sequência de referência com as sequências da condição
    with open(combined_file, "w") as outfile:
        with open(reference_file, "r") as ref:
            outfile.write(ref.read())
        with open(fasta_file, "r") as seqs:
            outfile.write(seqs.read())

    # Escolher a ferramenta de alinhamento
    if alignment_tool == "clustalw":
        # Expandir o ~ no caminho
        clustalw_path = os.path.expanduser(CLUSTALW_PATH)

        # Verificar ClustalW antes de executar
        if not os.path.exists(clustalw_path):
            clustalw_cmd = shutil.which("clustalw") or shutil.which("clustalw2")
            if not clustalw_cmd:
                raise FileNotFoundError(
                    f"ClustalW não encontrado em: {CLUSTALW_PATH} ou no PATH do sistema."
                )
            clustalw_path = clustalw_cmd

        cmd = [
            clustalw_path,  # Usar o caminho expandido ou encontrado no PATH
            "-INFILE=" + combined_file,
            "-OUTFILE=" + aligned_file,
            "-OUTPUT=CLUSTAL",
        ]
        if clustalw_params:
            for key, value in clustalw_params.items():
                cmd.append(f"-{key.upper()}={value}")

        # Executar o comando
        subprocess.run(cmd, check=True)

    elif alignment_tool == "mafft":
        # Obter parâmetros configurados para MAFFT
        from config import MAFFT_PARAMS

        mafft_params = MAFFT_PARAMS.get(condition, MAFFT_PARAMS.get("ADHD", {}))

        # Montar o comando com base nos parâmetros
        cmd = ["mafft", "--clustalout"]  # Saída em formato CLUSTAL

        # Adicionar algoritmo específico
        algorithm = mafft_params.get("algorithm")
        if algorithm:
            cmd.append(algorithm)

        # Configurar iterações
        max_iterate = mafft_params.get("maxiterate")
        if max_iterate is not None:
            cmd.extend(["--maxiterate", str(max_iterate)])

        # Opções adicionais
        if mafft_params.get("reorder", True):
            cmd.append("--reorder")

        if mafft_params.get("adjustdirection", False):
            cmd.append("--adjustdirection")

        cmd.append("--quiet")  # Reduzir a saída verbosa
        cmd.append(combined_file)  # Arquivo de entrada

        # Executar o comando
        with open(aligned_file, "w") as outfile:
            subprocess.run(cmd, stdout=outfile, check=True)

    else:
        raise ValueError(f"Ferramenta de alinhamento desconhecida: {alignment_tool}")

    return aligned_file


def align_all_sequences(all_sequences_file, reference_file, alignment_tool="clustalw"):
    """
    Realiza alinhamento de todas as sequências DRD4 disponíveis com a referência,
    para comparação geral (sem foco em patologias específicas).

    Args:
        all_sequences_file: Arquivo com todas as sequências a serem alinhadas.
        reference_file: Arquivo contendo a sequência de referência.
        alignment_tool: Ferramenta de alinhamento a ser usada ('clustalw' ou 'mafft').

    Returns:
        Caminho para o arquivo de alinhamento gerado.
    """
    # Importar o diretório de alinhamentos atualizado
    from config import ALIGNMENTS_DIR

    if alignment_tool == "mafft":
        logging.info(f"Iniciando alinhamento geral de sequências DRD4 com MAFFT...")
    else:
        logging.info(f"Iniciando alinhamento geral de sequências DRD4 com ClustalW...")

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

    if alignment_tool == "mafft":
        # Usar MAFFT para o alinhamento
        cmd = [
            "mafft",
            "--clustalout",  # Saída em formato CLUSTAL
            "--retree",
            "1",  # Opção rápida para conjuntos grandes
            "--quiet",  # Reduzir saída verbosa
            combined_file,
        ]

        with open(aligned_file, "w") as outfile:
            subprocess.run(cmd, stdout=outfile, check=True)
    else:
        # Usar ClustalW para o alinhamento (comportamento padrão anterior)
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

    logging.info(f"Alinhamento geral concluído. Arquivo gerado: {aligned_file}")
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
