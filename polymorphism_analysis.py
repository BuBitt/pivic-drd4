import logging
import time
import numpy as np
from Bio import AlignIO
from Bio.AlignIO import ClustalIO
from config import POLYMORPHISM_CONFIG


def detect_polymorphisms(aligned_file, min_vntr_length=48, max_gap_ratio=0.4):
    """Detecta VNTRs e SNPs no alinhamento com parâmetros ajustáveis."""
    logging.info(f"Analisando alinhamento para detectar polimorfismos: {aligned_file}")
    start_time = time.time()

    try:
        # Primeiro tentamos ler como formato Clustal padrão
        try:
            alignment = AlignIO.read(aligned_file, "clustal")
        except Exception as first_error:
            # Se falhar, podemos estar lidando com um formato MAFFT
            logging.warning(f"Erro na leitura padrão do alinhamento: {first_error}")
            logging.info(
                "Tentando ler o arquivo usando alternativas para suporte ao MAFFT..."
            )

            try:
                # Tentativa alternativa: ler linha por linha e parsear
                with open(aligned_file, "r") as f:
                    lines = f.readlines()

                # Remover linhas de cabeçalho MAFFT que poderiam estar causando problemas
                cleaned_lines = []
                started = False

                for line in lines:
                    # Pular linhas de cabeçalho do MAFFT até encontrar a primeira linha de sequência
                    if (
                        not started
                        and line.strip()
                        and not line.startswith("CLUSTAL")
                        and not line.startswith("#")
                    ):
                        started = True

                    if started:
                        cleaned_lines.append(line)
                    elif line.startswith("CLUSTAL"):
                        cleaned_lines.append(line)  # Manter o cabeçalho CLUSTAL

                # Escrever em um arquivo temporário
                import tempfile

                with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp:
                    temp.writelines(cleaned_lines)
                    temp_name = temp.name

                # Tentar ler o arquivo temporário
                alignment = AlignIO.read(temp_name, "clustal")

                # Remover o arquivo temporário
                import os

                os.unlink(temp_name)

                logging.info(
                    "Sucesso ao ler o arquivo de alinhamento com tratamento especial."
                )
            except Exception as e:
                raise RuntimeError(f"Falha completa ao ler arquivo de alinhamento: {e}")

    except Exception as e:
        logging.error(f"Erro ao ler arquivo de alinhamento: {e}")
        return {"VNTRs": [], "SNPs": {}}

    alignment_length = alignment.get_alignment_length()
    num_sequences = len(alignment)
    logging.info(f"Tamanho do alinhamento: {alignment_length} posições.")

    # Determinar condição com base no nome do arquivo
    condition = "ADHD" if "ADHD" in aligned_file else "Autismo"

    # Usar configurações específicas para cada condição
    config = POLYMORPHISM_CONFIG.get(condition, {})
    min_vntr_length = config.get("min_vntr_length", min_vntr_length)
    max_gap_ratio = config.get("max_gap_ratio", max_gap_ratio)
    gap_threshold = config.get("gap_threshold", 0.2)
    snp_min_coverage = config.get("snp_min_coverage", 0.3)

    logging.info(
        f"Usando parâmetros para {condition}: min_vntr_length={min_vntr_length}, "
        f"max_gap_ratio={max_gap_ratio}, gap_threshold={gap_threshold}"
    )

    # Inicializar estruturas de dados
    vntr_regions = []
    snps = {}
    current_vntr = None

    # Pré-processar para detectar regiões com alta densidade de gaps - possíveis VNTRs
    gap_counts = []
    for i in range(alignment_length):
        column = alignment[:, i]
        gap_count = column.count("-")
        gap_counts.append(gap_count / num_sequences)

    # Suavizar as contagens de gap para melhor detecção
    window_size = 5
    smoothed_gaps = np.convolve(
        gap_counts, np.ones(window_size) / window_size, mode="same"
    )

    # Analisar de forma mais precisa
    for i in range(alignment_length):
        column = alignment[:, i]
        gap_ratio = smoothed_gaps[i]

        # Descartar colunas com alinhamento de baixa qualidade (quase tudo gap)
        if gap_ratio > 0.9:
            continue

        # Extrair bases não-gap
        bases = [base for base in column if base != "-"]
        if not bases:  # Skip positions where all sequences have gaps
            continue

        unique_bases = set(bases)

        # Detectar SNPs com melhor critério
        base_counts = {}
        for base in bases:
            if base in base_counts:
                base_counts[base] += 1
            else:
                base_counts[base] = 1

        # SNP quando há mais de uma base com cobertura significativa
        significant_variants = []
        for base, count in base_counts.items():
            if count / len(bases) >= snp_min_coverage:
                significant_variants.append(base)

        if len(significant_variants) > 1:
            snps[i] = significant_variants

        # Detectar VNTR com sistema mais sensível
        if gap_threshold <= gap_ratio <= max_gap_ratio:
            # É um possível VNTR
            if current_vntr is None:
                current_vntr = [i]
            elif i - current_vntr[-1] <= 3:  # Permitir pequenas lacunas
                current_vntr.append(i)
            else:
                # Finalizar o VNTR anterior se for grande o suficiente
                if len(current_vntr) >= min_vntr_length:
                    vntr_regions.append(current_vntr)
                current_vntr = [i]
        elif current_vntr is not None:
            # Finalizar VNTR atual se ele for longo o suficiente
            if len(current_vntr) >= min_vntr_length:
                vntr_regions.append(current_vntr)
            current_vntr = None

    # Verificar último possível VNTR
    if current_vntr and len(current_vntr) >= min_vntr_length:
        vntr_regions.append(current_vntr)

    # Combinar pequenos VNTRs próximos
    combined_vntrs = []
    if vntr_regions:
        current_region = vntr_regions[0]
        for i in range(1, len(vntr_regions)):
            region = vntr_regions[i]
            # Se início da região atual está próximo do fim da anterior
            if min(region) - max(current_region) < 20:
                # Combinar regiões
                current_region.extend(region)
            else:
                combined_vntrs.append(current_region)
                current_region = region
        combined_vntrs.append(current_region)

    # Registrar resultados
    for region in combined_vntrs:
        if region:  # Verificar se a região não está vazia
            # Para exibição, converter para posição baseada em 1
            start_pos = min(region) + 1  # +1 para converter de base 0 para base 1
            end_pos = max(region) + 1  # +1 para converter de base 0 para base 1
            length = end_pos - start_pos + 1
            logging.info(
                f"VNTR detectado nas posições {start_pos} a {end_pos} (comprimento: {length} pb)"
            )

    # Eliminar SNPs dentro de VNTRs para evitar sobreposição
    filtered_snps = {}
    for pos, bases in snps.items():
        is_in_vntr = False
        for region in combined_vntrs:
            if min(region) <= pos <= max(region):
                is_in_vntr = True
                break
        if not is_in_vntr:
            filtered_snps[pos] = bases

    execution_time = time.time() - start_time
    logging.info(
        f"Detecção de polimorfismos concluída em {execution_time:.2f} segundos"
    )
    logging.info(
        f"Total de SNPs: {len(filtered_snps)}, Total de VNTRs: {len(combined_vntrs)}"
    )

    return {"VNTRs": combined_vntrs, "SNPs": filtered_snps}
