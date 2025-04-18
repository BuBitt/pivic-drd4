from datetime import datetime
import os
import logging
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np


def generate_report(
    polymorphisms,
    report_name,
    category,
    reports_dir=None,
    alignment_tool="mafft",
    reference_seq_id="NM_000797.4",
):
    """Gera um relatório em formato Markdown com os polimorfismos encontrados."""
    # Obter o diretório de relatórios atualizado
    if reports_dir is None:
        from config import REPORTS_DIR

        reports_dir = REPORTS_DIR

    if not reports_dir:
        raise ValueError(
            "O diretório de relatórios (REPORTS_DIR) não foi configurado corretamente"
        )

    # Garantir que o diretório existe
    os.makedirs(reports_dir, exist_ok=True)

    # Mudar a extensão para .md
    if report_name.endswith(".txt"):
        report_name = report_name.replace(".txt", ".md")
    if not report_name.endswith(".md"):
        report_name = report_name + ".md"

    logging.info(
        f"Gerando relatório Markdown dos polimorfismos detectados para {category}..."
    )
    report_file = os.path.join(reports_dir, report_name)

    with open(report_file, "w") as f:
        f.write(f"# Relatório de Polimorfismos no Gene DRD4 ({category})\n\n")

        # Adicionar informação sobre a ferramenta de alinhamento usada
        from config import ALIGNMENT_DESCRIPTIONS

        f.write("## FERRAMENTA DE ALINHAMENTO\n\n")

        if alignment_tool == "mafft":
            from config import MAFFT_PARAMS

            params = MAFFT_PARAMS.get(category, MAFFT_PARAMS.get("ADHD", {}))
            algorithm = params.get("algorithm", "--auto").strip("-")

            if algorithm == "auto":
                f.write(f"{ALIGNMENT_DESCRIPTIONS['mafft_auto']}\n\n")
            elif algorithm == "genafpair":
                f.write(f"{ALIGNMENT_DESCRIPTIONS['mafft_genafpair']}\n\n")
            elif algorithm == "localpair":
                f.write(f"{ALIGNMENT_DESCRIPTIONS['mafft_localpair']}\n\n")
            else:
                f.write(f"{ALIGNMENT_DESCRIPTIONS['mafft']}\n\n")

            f.write(f"**Parâmetros utilizados:**\n\n")
            f.write(f"- **Algoritmo**: {algorithm}\n")
            f.write(f"- **Iterações de refinamento**: {params.get('maxiterate', 2)}\n")
            if params.get("reorder", False):
                f.write(f"- **Reordenação de sequências**: Sim\n")
            if params.get("adjustdirection", False):
                f.write(f"- **Ajuste automático de direção**: Sim\n")
            f.write("\n")
        else:
            f.write(f"{ALIGNMENT_DESCRIPTIONS['clustalw']}\n\n")

            # Adicionar parâmetros do ClustalW se necessário
            from config import CLUSTALW_PARAMS

            params = CLUSTALW_PARAMS.get(category, {})
            if params:
                f.write(f"**Parâmetros utilizados:**\n\n")
                f.write(f"- **Gap Open**: {params.get('GAPOPEN', 'Padrão')}\n")
                f.write(f"- **Gap Extension**: {params.get('GAPEXT', 'Padrão')}\n")
                f.write(f"- **Matriz DNA**: {params.get('DNAMATRIX', 'Padrão')}\n\n")

        # Adicionar seção de glossário com explicações das siglas
        f.write("## GLOSSÁRIO DE TERMOS E SIGLAS\n\n")
        f.write(
            "**DRD4**  - Receptor de Dopamina D4: Gene que codifica um dos cinco receptores de dopamina "
        )
        f.write(
            "encontrados em humanos. Está associado à regulação de funções cognitivas e "
        )
        f.write("comportamentais.\n\n")

        f.write(
            "**SNP**   - Polimorfismo de Nucleotídeo Único (Single Nucleotide Polymorphism): "
        )
        f.write(
            "Variação em um único nucleotídeo que ocorre em uma posição específica do genoma.\n\n"
        )

        f.write(
            "**VNTR**  - Repetições em Tandem de Número Variável (Variable Number Tandem Repeat): "
        )
        f.write(
            "Regiões onde sequências curtas de DNA se repetem um número variável de vezes.\n"
        )
        f.write(
            "No DRD4, o VNTR no éxon 3 (variando de 2 a 11 repetições) tem sido associado "
        )
        f.write("a diferenças em fenótipos comportamentais.\n\n")

        f.write(
            "**ADHD**  - Transtorno de Déficit de Atenção e Hiperatividade (TDAH): Condição "
        )
        f.write(
            "neurocomportamental caracterizada por desatenção, hiperatividade e impulsividade.\n\n"
        )

        f.write(
            "**ASD**   - Transtorno do Espectro Autista (Autism Spectrum Disorder): Grupo de condições "
        )
        f.write("caracterizadas por desafios com interação social, comunicação e ")
        f.write("comportamentos repetitivos.\n\n")

        f.write(
            "**pb**    - Pares de bases: Unidade de medida para comprimento de segmentos de DNA.\n\n"
        )

        f.write(
            "**EXON**  - Segmento do gene que codifica diretamente para a proteína.\n\n"
        )

        f.write(
            "**ALELO** - Variante de um gene em uma posição específica do cromossomo.\n"
        )
        f.write("Ex.: Alelo de 7 repetições (7R) do VNTR do éxon 3 do DRD4.\n\n")

        if category == "ADHD":
            f.write("### ESPECÍFICO PARA ADHD/TDAH\n\n")
            f.write(
                "**7R**    - Alelo de 7 repetições do VNTR no éxon 3 do DRD4, frequentemente "
            )
            f.write("associado ao TDAH.\n\n")
            f.write(
                "**4R**    - Alelo de 4 repetições do VNTR, considerado o mais comum na população.\n\n"
            )
        elif category == "Autismo":
            f.write("### ESPECÍFICO PARA AUTISMO\n\n")
            f.write(
                "**TEA**   - Transtorno do Espectro Autista, termo atual para designar o conjunto de "
            )
            f.write(
                "condições anteriormente separadas como autismo, síndrome de Asperger, etc.\n\n"
            )
            f.write(
                "**PDD-NOS** - Transtorno Invasivo do Desenvolvimento Sem Outra Especificação, uma "
            )
            f.write("classificação anterior dentro do espectro autista.\n\n")

        # Adicionar informações sobre as variantes LC812 caso seja a categoria de variantes
        if "Variantes" in category:
            f.write("## NOTA SOBRE SEQUÊNCIAS VARIANTES\n\n")
            f.write(
                "As sequências classificadas como variantes do gene DRD4 (especialmente as com prefixo LC812) "
            )
            f.write(
                "representam variantes estruturalmente divergentes em Homo sapiens. Estas sequências foram "
            )
            f.write(
                "selecionadas através de um processo de identificação específico:\n\n"
            )

            f.write(
                "- **Identificação por prefixos específicos**: Seleção exclusiva das sequências LC812 e outras LC* "
            )
            f.write(
                "que representam variantes moleculares oficialmente catalogadas na Library of Congress\n\n"
            )

            f.write(
                "Esta abordagem precisa garante a identificação das variantes funcionalmente mais "
            )
            f.write(
                "relevantes para a compreensão da diversidade do gene DRD4 em casos de ADHD.\n\n"
            )

        f.write(f"**Data da análise**: {datetime.now().strftime('%d/%m/%Y %H:%M')}\n\n")
        f.write(f"**Organismo**: Homo sapiens (humano)\n\n")
        f.write(
            f"**Sequência de referência**: Gene DRD4 humano, {reference_seq_id} (NCBI Reference Sequence)\n\n"
        )

        # Resumo
        total_snps = len(polymorphisms["SNPs"])
        total_vntrs = len(polymorphisms["VNTRs"])
        f.write("## RESUMO\n\n")
        f.write(f"- Total de SNPs detectados: **{total_snps}**\n")
        f.write(f"- Total de VNTRs detectados: **{total_vntrs}**\n\n")

        # Informações sobre VNTRs
        if polymorphisms["VNTRs"]:
            f.write("## VNTRs Detectados\n\n")

            vntr_lengths = []
            for i, region in enumerate(polymorphisms["VNTRs"]):
                start, end = min(region), max(region)
                length = end - start + 1
                vntr_lengths.append(length)

                f.write(
                    f"### VNTR #{i + 1}: Posições {start} a {end} (comprimento: {length} pb)\n\n"
                )

                # Adicionar interpretação biológica
                if length >= 2000:
                    f.write(
                        f"  * Este é um VNTR extenso que pode afetar a estrutura da proteína\n\n"
                    )
                elif length >= 400:
                    f.write(
                        f"  * Pode corresponder à região VNTR de 48-bp no éxon 3\n\n"
                    )
                elif length >= 150:
                    f.write(
                        f"  * Tamanho sugestivo do VNTR associado a fenótipos comportamentais\n\n"
                    )
                elif length >= 48:
                    f.write(
                        f"  * Potencial variação no número de repetições (ex.: 2R, 4R, 7R)\n\n"
                    )

            # Estatísticas sobre VNTRs
            if vntr_lengths:
                avg_length = sum(vntr_lengths) / len(vntr_lengths)
                f.write(f"**Comprimento médio dos VNTRs**: {avg_length:.2f} pb\n\n")
                f.write(f"**Maior VNTR**: {max(vntr_lengths)} pb\n\n")
                f.write(f"**Menor VNTR**: {min(vntr_lengths)} pb\n\n")
        else:
            f.write("## VNTRs\n\nNão detectados ou inconclusivos.\n\n")

        # Informações sobre SNPs
        f.write("## SNPs Detectados\n\n")

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

            # Adicionar nota sobre convenção de posição
            f.write(
                "> **NOTA:** Todas as posições seguem a convenção biológica e começam em 1 (não em 0).\n\n"
            )

            # Adicionar nota esclarecendo a posição relativa à referência
            f.write(
                "> **NOTA IMPORTANTE:** Todas as posições seguem a convenção biológica (começam em 1, não em 0) "
            )
            f.write(
                f"e são relativas à sequência de referência do gene DRD4 ({reference_seq_id}). "
            )
            f.write(
                "Os SNPs indicam locais onde as sequências analisadas divergem da referência.\n\n"
            )

            # Mostrar SNPs por clusters
            for i, cluster in enumerate(snp_clusters):
                f.write(f"### Cluster de SNPs #{i + 1} ({len(cluster)} variantes)\n\n")
                f.write("| Posição (ref) | Ref | Variantes |\n")
                f.write("|---------------|-----|----------|\n")
                for pos in cluster:
                    bases = polymorphisms["SNPs"][pos]

                    # Tentativa de identificar a base de referência (geralmente a primeira na lista)
                    ref_base = "N"  # Default se não puder determinar

                    # Tentar obter a base de referência (isso depende da estrutura dos seus dados)
                    if (
                        hasattr(polymorphisms, "reference_seq")
                        and len(polymorphisms["reference_seq"]) > pos
                    ):
                        ref_base = polymorphisms["reference_seq"][pos].lower()
                    elif len(bases) > 0:
                        # Se não temos a referência, assumimos a primeira base como referência
                        ref_base = bases[0]

                    # Exibir a posição com base 1 e base de referência
                    f.write(f"| {pos + 1} | {ref_base} | {', '.join(bases)} |\n")
                f.write("\n")

            # Adicionar informações sobre regiões mais variáveis
            hot_regions = [cluster[0] for cluster in snp_clusters if len(cluster) >= 3]
            if hot_regions:
                f.write("### Regiões de alta variabilidade (hotspots)\n\n")
                for pos in hot_regions:
                    f.write(f"- Região a partir da posição {pos}\n")
                f.write("\n")
        else:
            f.write("Nenhum SNP identificado.\n\n")

        # Adicionar notas interpretativas específicas para cada condição
        f.write("## NOTAS INTERPRETATIVAS\n\n")

        if category == "ADHD":
            f.write(
                "- O alelo de 7 repetições (7R) do VNTR de 48-pb no éxon 3 do DRD4 "
            )
            f.write("tem sido associado a risco aumentado para TDAH.\n\n")
            f.write(
                "- SNPs nas regiões promotoras podem afetar a expressão do receptor D4 "
            )
            f.write("e modular a resposta a medicamentos como o metilfenidato.\n")
        elif category == "Autismo":
            f.write(
                "- Variações no DRD4 podem contribuir para aspectos comportamentais "
            )
            f.write(
                "do espectro autista, especialmente relacionados à atenção e controle motor.\n\n"
            )
            f.write(
                "- A relação entre polimorfismos do DRD4 e autismo é menos direta que no TDAH, "
            )
            f.write("possivelmente como parte de um perfil genético mais complexo.\n")

    logging.info(f"Relatório Markdown gerado no arquivo: {report_file}")
    return report_file


def analyze_polymorphism_report(report_file):
    """
    Analisa um relatório de polimorfismos para extrair estatísticas.

    Args:
        report_file: Caminho para o arquivo de relatório

    Returns:
        Dicionário com estatísticas sobre os polimorfismos
    """
    logging.info(f"Analisando estatisticamente o relatório de polimorfismos...")

    # Inicializar estatísticas
    stats = {
        "total_snps": 0,
        "snp_positions": [],
        "snp_density": 0,
        "total_vntrs": 0,
        "vntr_lengths": [],
        "avg_vntr_length": 0,
        "hotspots": 0,
    }

    try:
        # Verificar se o relatório existe
        if not os.path.exists(report_file):
            logging.warning(f"Arquivo de relatório não encontrado: {report_file}")
            return stats

        # CORREÇÃO: Extrair diretamente dos polimorfismos brutos em vez de tentar analisar o relatório
        # Usar o mesmo caminho do relatório para determinar o caminho do arquivo de alinhamento
        alignment_dir = os.path.join(
            os.path.dirname(os.path.dirname(report_file)), "alignments"
        )
        condition = os.path.basename(report_file).split("_")[-1].split(".")[0]
        alignment_file = os.path.join(alignment_dir, f"drd4_aligned_{condition}.aln")

        # Se o arquivo de alinhamento existir, analisar para obter os polimorfismos brutos
        if os.path.exists(alignment_file):
            from polymorphism_analysis import detect_polymorphisms
            from config import POLYMORPHISM_CONFIG

            # Obter configuração específica
            config = POLYMORPHISM_CONFIG.get(condition.capitalize(), {})

            # Detectar polimorfismos diretamente
            logging.info(
                f"Re-analisando polimorfismos a partir do alinhamento: {alignment_file}"
            )
            polymorphisms = detect_polymorphisms(
                alignment_file,
                min_vntr_length=config.get("min_vntr_length", 24),
                max_gap_ratio=config.get("max_gap_ratio", 0.7),
                gap_threshold=config.get("gap_threshold", 0.25),
            )

            # Extrair estatísticas dos polimorfismos detectados
            if polymorphisms and "SNPs" in polymorphisms:
                stats["total_snps"] = len(polymorphisms["SNPs"])
                stats["snp_positions"] = list(polymorphisms["SNPs"].keys())

                # Calcular densidade de SNPs se houver alinhamento
                alignment_length = 0
                try:
                    from Bio import AlignIO

                    alignment = AlignIO.read(alignment_file, "clustal")
                    alignment_length = alignment.get_alignment_length()
                    if alignment_length > 0 and stats["total_snps"] > 0:
                        stats["snp_density"] = stats["total_snps"] / alignment_length
                except Exception as e:
                    logging.warning(f"Erro ao calcular densidade de SNPs: {e}")

                # Identificar hotspots (regiões com 3+ SNPs próximos)
                positions = sorted(stats["snp_positions"])
                hotspot_count = 0
                i = 0
                while i < len(positions) - 2:
                    if (
                        positions[i + 2] - positions[i] <= 10
                    ):  # 3 SNPs dentro de 10 bases
                        hotspot_count += 1
                        # Pular para depois deste hotspot
                        while (
                            i < len(positions) - 1
                            and positions[i + 1] - positions[i] <= 10
                        ):
                            i += 1
                    i += 1
                stats["hotspots"] = hotspot_count

            if polymorphisms and "VNTRs" in polymorphisms:
                stats["total_vntrs"] = len(polymorphisms["VNTRs"])
                vntr_lengths = []

                for region in polymorphisms["VNTRs"]:
                    if region:  # Verificar se a região não está vazia
                        length = max(region) - min(region) + 1
                        vntr_lengths.append(length)

                stats["vntr_lengths"] = vntr_lengths

                # Calcular comprimento médio de VNTRs se houver algum
                if vntr_lengths:
                    stats["avg_vntr_length"] = sum(vntr_lengths) / len(vntr_lengths)

            # Log dos resultados obtidos
            logging.info(f"Total de SNPs detectados: {stats['total_snps']}")
            logging.info(f"Total de VNTRs detectados: {stats['total_vntrs']}")
        else:
            logging.warning(f"Arquivo de alinhamento não encontrado: {alignment_file}")

    except Exception as e:
        logging.error(
            f"Erro ao analisar relatório de polimorfismos: {str(e)}", exc_info=True
        )

    return stats


def create_summary_visualization(polymorphisms, condition, visualization_dir=None):
    """Cria visualizações para os polimorfismos encontrados."""
    # Obter o diretório de visualização atualizado
    if visualization_dir is None:
        from config import VISUALIZATION_DIR

        visualization_dir = VISUALIZATION_DIR

    if not visualization_dir:
        raise ValueError(
            "O diretório de visualização (VISUALIZATION_DIR) não foi configurado corretamente"
        )

    try:
        # Criar diretório para visualizações se não existir
        os.makedirs(visualization_dir, exist_ok=True)

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
            visualization_dir, f"drd4_polymorphisms_{condition.lower()}.png"
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
                visualization_dir, f"vntr_distribution_{condition.lower()}.png"
            )
            plt.savefig(output_file, dpi=300, bbox_inches="tight")
            plt.close()

            logging.info(
                f"Visualização de distribuição de VNTRs salva em: {output_file}"
            )

    except Exception as e:
        logging.error(f"Erro ao criar visualização para {condition}: {str(e)}")
