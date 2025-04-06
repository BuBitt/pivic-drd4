from datetime import datetime
import os
import logging
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np
from config import REPORTS_DIR, VISUALIZATION_DIR


def generate_report(polymorphisms, report_name, category):
    """Gera um relatório com os polimorfismos encontrados, categorizado por TDAH ou autismo."""
    logging.info(f"Gerando relatório dos polimorfismos detectados para {category}...")
    report_file = os.path.join(REPORTS_DIR, report_name)
    with open(report_file, "w") as f:
        f.write(f"Relatório de Polimorfismos no Gene DRD4 ({category})\n")
        f.write("=" * 60 + "\n\n")

        # Adicionar seção de glossário com explicações das siglas
        f.write("GLOSSÁRIO DE TERMOS E SIGLAS:\n")
        f.write("-" * 50 + "\n")
        f.write(
            "DRD4  - Receptor de Dopamina D4: Gene que codifica um dos cinco receptores de dopamina\n"
        )
        f.write(
            "         encontrados em humanos. Está associado à regulação de funções cognitivas e\n"
        )
        f.write("         comportamentais.\n\n")
        f.write(
            "SNP   - Polimorfismo de Nucleotídeo Único (Single Nucleotide Polymorphism): \n"
        )
        f.write(
            "         Variação em um único nucleotídeo que ocorre em uma posição específica do genoma.\n\n"
        )
        f.write(
            "VNTR  - Repetições em Tandem de Número Variável (Variable Number Tandem Repeat):\n"
        )
        f.write(
            "         Regiões onde sequências curtas de DNA se repetem um número variável de vezes.\n"
        )
        f.write(
            "         No DRD4, o VNTR no éxon 3 (variando de 2 a 11 repetições) tem sido associado\n"
        )
        f.write("         a diferenças em fenótipos comportamentais.\n\n")
        f.write(
            "ADHD  - Transtorno de Déficit de Atenção e Hiperatividade (TDAH): Condição\n"
        )
        f.write(
            "         neurocomportamental caracterizada por desatenção, hiperatividade e impulsividade.\n\n"
        )
        f.write(
            "ASD   - Transtorno do Espectro Autista (Autism Spectrum Disorder): Grupo de condições\n"
        )
        f.write(
            "         caracterizadas por desafios com interação social, comunicação e\n"
        )
        f.write("         comportamentos repetitivos.\n\n")
        f.write(
            "pb    - Pares de bases: Unidade de medida para comprimento de segmentos de DNA.\n\n"
        )
        f.write(
            "EXON  - Segmento do gene que codifica diretamente para a proteína.\n\n"
        )
        f.write(
            "ALELO - Variante de um gene em uma posição específica do cromossomo.\n"
        )
        f.write(
            "         Ex.: Alelo de 7 repetições (7R) do VNTR do éxon 3 do DRD4.\n\n"
        )

        if category == "ADHD":
            f.write("ESPECÍFICO PARA ADHD/TDAH:\n")
            f.write(
                "7R    - Alelo de 7 repetições do VNTR no éxon 3 do DRD4, frequentemente\n"
            )
            f.write("         associado ao TDAH.\n\n")
            f.write(
                "4R    - Alelo de 4 repetições do VNTR, considerado o mais comum na população.\n\n"
            )
        elif category == "Autismo":
            f.write("ESPECÍFICO PARA AUTISMO:\n")
            f.write(
                "TEA   - Transtorno do Espectro Autista, termo atual para designar o conjunto de\n"
            )
            f.write(
                "         condições anteriormente separadas como autismo, síndrome de Asperger, etc.\n\n"
            )
            f.write(
                "PDD-NOS - Transtorno Invasivo do Desenvolvimento Sem Outra Especificação, uma\n"
            )
            f.write("         classificação anterior dentro do espectro autista.\n\n")

        f.write("-" * 50 + "\n\n")

        # Adicionar informações sobre as variantes LC812 caso seja a categoria de variantes
        if "Variantes" in category:
            f.write("NOTA SOBRE SEQUÊNCIAS VARIANTES:\n")
            f.write("-" * 50 + "\n")
            f.write(
                "As sequências com prefixo 'LC812' representam variantes divergentes do gene DRD4 em\n"
            )
            f.write(
                "Homo sapiens. Estas sequências apresentam diferenças significativas em comparação\n"
            )
            f.write(
                "com as sequências típicas de DRD4, mas foram mantidas na análise por sua relevância\n"
            )
            f.write(
                "potencial para o estudo de ADHD. Estas variantes podem representar:\n\n"
            )
            f.write(
                "1. Isoformas alternativas do gene DRD4 resultantes de splicing alternativo\n"
            )
            f.write(
                "2. Variantes populacionais específicas com alterações estruturais significativas\n"
            )
            f.write(
                "3. Regiões diferentes do gene DRD4 com relevância funcional distinta\n"
            )
            f.write("4. Polimorfismos raros com possíveis implicações fenotípicas\n\n")
            f.write(
                "A análise destas variantes separadamente permite uma melhor compreensão da\n"
            )
            f.write("diversidade genética do DRD4 em casos de ADHD.\n\n")

        f.write(f"Data da análise: {datetime.now().strftime('%d/%m/%Y %H:%M')}\n\n")
        f.write(f"Organismo: Homo sapiens (humano)\n\n")

        # Resumo
        total_snps = len(polymorphisms["SNPs"])
        total_vntrs = len(polymorphisms["VNTRs"])
        f.write(f"RESUMO:\n")
        f.write(f"- Total de SNPs detectados: {total_snps}\n")
        f.write(f"- Total de VNTRs detectados: {total_vntrs}\n\n")

        # Informações sobre VNTRs
        if polymorphisms["VNTRs"]:
            f.write("\nVNTRs Detectados:\n")
            f.write("-" * 40 + "\n")

            vntr_lengths = []
            for i, region in enumerate(polymorphisms["VNTRs"]):
                start, end = min(region), max(region)
                length = end - start + 1
                vntr_lengths.append(length)

                f.write(
                    f"VNTR #{i + 1}: Posições {start} a {end} (comprimento: {length} pb)\n"
                )

                # Adicionar interpretação biológica
                if length >= 2000:
                    f.write(
                        f"  * Este é um VNTR extenso que pode afetar a estrutura da proteína\n"
                    )
                elif length >= 400:
                    f.write(f"  * Pode corresponder à região VNTR de 48-bp no éxon 3\n")
                elif length >= 150:
                    f.write(
                        f"  * Tamanho sugestivo do VNTR associado a fenótipos comportamentais\n"
                    )
                elif length >= 48:
                    f.write(
                        f"  * Potencial variação no número de repetições (ex.: 2R, 4R, 7R)\n"
                    )

            # Estatísticas sobre VNTRs
            if vntr_lengths:
                avg_length = sum(vntr_lengths) / len(vntr_lengths)
                f.write(f"\nComprimento médio dos VNTRs: {avg_length:.2f} pb\n")
                f.write(f"Maior VNTR: {max(vntr_lengths)} pb\n")
                f.write(f"Menor VNTR: {min(vntr_lengths)} pb\n")
        else:
            f.write("\nVNTRs: Não detectados ou inconclusivos.\n")

        # Informações sobre SNPs
        f.write("\nSNPs Detectados:\n")
        f.write("-" * 40 + "\n")

        if polymorphisms["SNPs"]:
            # Agrupar SNPs próximos para melhor visualização
            positions = sorted(polymorphisms["SNPs"].keys())
            snp_clusters = []
            current_cluster = [positions[0]]

            for i in range(1, len(positions)):
                if (
                    positions[i] - positions[i - 1] <= 10
                ):  # SNPs dentro de 10 bases são agrupados
                    current_cluster.append(positions[i])
                else:
                    snp_clusters.append(current_cluster)
                    current_cluster = [positions[i]]

            if current_cluster:
                snp_clusters.append(current_cluster)

            # Mostrar SNPs por clusters
            for i, cluster in enumerate(snp_clusters):
                f.write(f"\nCluster de SNPs #{i + 1} ({len(cluster)} variantes):\n")
                for pos in cluster:
                    bases = polymorphisms["SNPs"][pos]
                    f.write(f"  - Posição {pos}: {', '.join(bases)}\n")

            # Adicionar informações sobre regiões mais variáveis
            hot_regions = [cluster[0] for cluster in snp_clusters if len(cluster) >= 3]
            if hot_regions:
                f.write("\nRegiões de alta variabilidade (hotspots):\n")
                for pos in hot_regions:
                    f.write(f"  - Região a partir da posição {pos}\n")
        else:
            f.write("Nenhum SNP identificado.\n")

        # Adicionar notas interpretativas específicas para cada condição
        f.write("\nNOTAS INTERPRETATIVAS:\n")
        f.write("-" * 40 + "\n")

        if category == "ADHD":
            f.write(
                "- O alelo de 7 repetições (7R) do VNTR de 48-pb no éxon 3 do DRD4 "
            )
            f.write("tem sido associado a risco aumentado para TDAH.\n")
            f.write(
                "- SNPs nas regiões promotoras podem afetar a expressão do receptor D4 "
            )
            f.write("e modular a resposta a medicamentos como o metilfenidato.\n")
        elif category == "Autismo":
            f.write(
                "- Variações no DRD4 podem contribuir para aspectos comportamentais "
            )
            f.write(
                "do espectro autista, especialmente relacionados à atenção e controle motor.\n"
            )
            f.write(
                "- A relação entre polimorfismos do DRD4 e autismo é menos direta que no TDAH, "
            )
            f.write("possivelmente como parte de um perfil genético mais complexo.\n")

    logging.info(f"Relatório gerado no arquivo: {report_file}")
    return report_file


def analyze_polymorphism_report(report_file):
    """Realiza uma análise estatística sobre o relatório de polimorfismos."""
    logging.info("Analisando estatisticamente o relatório de polimorfismos...")
    total_snps = 0
    snp_positions = []
    total_vntrs = 0
    vntr_lengths = []
    hotspots = 0

    with open(report_file, "r") as f:
        lines = f.readlines()

    reading_vntr = False
    reading_snp = False

    for line in lines:
        # Identificar seção
        if "VNTRs Detectados:" in line:
            reading_vntr = True
            reading_snp = False
            continue
        elif "SNPs Detectados:" in line:
            reading_vntr = False
            reading_snp = True
            continue

        # Processar VNTRs
        if reading_vntr and line.startswith("VNTR #"):
            total_vntrs += 1
            try:
                # Extrair as posições de início e fim e o comprimento do VNTR
                parts = line.split("(comprimento: ")
                if len(parts) > 1:
                    length_part = parts[1].split(" pb")[0]
                    length = int(length_part)
                    vntr_lengths.append(length)
            except (ValueError, IndexError) as e:
                logging.error(f"Erro ao processar linha de VNTR: {line.strip()} ({e})")

        # Processar SNPs
        elif reading_snp and "Posição" in line:
            total_snps += 1
            try:
                pos_part = line.split("Posição ")[1].split(":")[0]
                position = int(pos_part)
                snp_positions.append(position)
            except (ValueError, IndexError) as e:
                logging.error(f"Erro ao processar linha de SNP: {line.strip()} ({e})")

        # Identificar hotspots
        elif reading_snp and "Regiões de alta variabilidade" in line:
            hotspots += 1

    # Calcular estatísticas adicionais
    avg_vntr_length = sum(vntr_lengths) / len(vntr_lengths) if vntr_lengths else 0

    # Analisar distribuição de SNPs
    snp_density = (
        len(snp_positions) / (max(snp_positions) - min(snp_positions) + 1)
        if snp_positions
        else 0
    )

    # Estatísticas de SNPs
    logging.info(f"Total de SNPs detectados: {total_snps}")
    if snp_positions:
        logging.info(f"Posições dos SNPs: {snp_positions}")
        logging.info(f"Densidade de SNPs: {snp_density:.6f} SNPs/pb")

    # Estatísticas de VNTRs
    logging.info(f"Total de VNTRs detectados: {total_vntrs}")
    if vntr_lengths:
        logging.info(f"Comprimento médio dos VNTRs: {avg_vntr_length:.2f} pb")
        logging.info(f"Comprimentos individuais: {vntr_lengths}")

    # Construir o resultado
    result = {
        "total_snps": total_snps,
        "snp_positions": snp_positions,
        "snp_density": snp_density if snp_positions else 0,
        "total_vntrs": total_vntrs,
        "vntr_lengths": vntr_lengths,
        "avg_vntr_length": avg_vntr_length,
        "hotspots": hotspots,
    }

    return result


def create_summary_visualization(polymorphisms, condition):
    """Cria visualizações para os polimorfismos encontrados."""
    try:
        # Criar diretório para visualizações se não existir
        os.makedirs(VISUALIZATION_DIR, exist_ok=True)

        # Preparar dados
        snp_positions = sorted(list(polymorphisms["SNPs"].keys()))
        vntr_regions = []

        for region in polymorphisms["VNTRs"]:
            start, end = min(region), max(region)
            vntr_regions.append((start, end))

        # Calcular intervalo total
        all_positions = snp_positions + [
            pos for region in vntr_regions for pos in region
        ]
        if not all_positions:
            logging.warning(
                f"Sem dados suficientes para criar visualização para {condition}"
            )
            return

        min_pos = max(0, min(all_positions) - 50)
        max_pos = max(all_positions) + 50

        # Criar figura
        plt.figure(figsize=(12, 6))

        # Plotar SNPs como pontos
        if snp_positions:
            plt.scatter(
                snp_positions,
                [1] * len(snp_positions),
                color="red",
                marker="o",
                s=50,
                label="SNPs",
            )

        # Plotar VNTRs como regiões
        for i, (start, end) in enumerate(vntr_regions):
            plt.axvspan(
                start, end, alpha=0.3, color="blue", label="VNTR" if i == 0 else ""
            )

        # Adicionar anotações
        plt.title(f"Polimorfismos do DRD4 - {condition}")
        plt.xlabel("Posição na sequência")
        plt.yticks([])  # Remover ticks do eixo Y
        plt.xlim(min_pos, max_pos)

        # Adicionar legenda e informações
        if snp_positions or vntr_regions:
            plt.legend()

        # Salvar figura
        output_file = os.path.join(
            VISUALIZATION_DIR, f"drd4_polymorphisms_{condition.lower()}.png"
        )
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()

        logging.info(f"Visualização salva em: {output_file}")

        # Criar visualização da distribuição de comprimentos de VNTR
        if len(polymorphisms["VNTRs"]) > 0:
            plt.figure(figsize=(8, 5))

            vntr_lengths = []
            for region in polymorphisms["VNTRs"]:
                start, end = min(region), max(region)
                vntr_lengths.append(end - start + 1)

            plt.hist(
                vntr_lengths,
                bins=min(10, len(vntr_lengths)),
                alpha=0.7,
                color="blue",
                edgecolor="black",
            )
            plt.title(f"Distribuição de Comprimentos de VNTR - {condition}")
            plt.xlabel("Comprimento (pb)")
            plt.ylabel("Frequência")

            output_file = os.path.join(
                VISUALIZATION_DIR, f"vntr_distribution_{condition.lower()}.png"
            )
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            plt.close()

            logging.info(
                f"Visualização de distribuição de VNTRs salva em: {output_file}"
            )

    except Exception as e:
        logging.error(f"Erro ao criar visualização para {condition}: {str(e)}")
