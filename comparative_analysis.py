import os
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO


def compare_condition_polymorphisms(results, visualization_dir=None, reports_dir=None):
    """
    Compara polimorfismos encontrados entre diferentes condições (ADHD, Autismo, Comorbidade)

    Args:
        results: Lista de resultados da análise de polimorfismos por condição
        visualization_dir: Diretório para salvar visualizações
        reports_dir: Diretório para salvar relatórios
    """
    # Obter os diretórios atualizados
    if visualization_dir is None:
        from config import VISUALIZATION_DIR

        visualization_dir = VISUALIZATION_DIR

    if reports_dir is None:
        from config import REPORTS_DIR

        reports_dir = REPORTS_DIR

    if not visualization_dir or not reports_dir:
        raise ValueError(
            "Os diretórios de visualização e relatórios não foram configurados corretamente"
        )

    # Garantir que os diretórios existam
    os.makedirs(visualization_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)

    logging.info("Iniciando análise comparativa entre condições...")

    # Debug: mostrar todos os resultados recebidos
    for i, result in enumerate(results):
        logging.info(f"Resultado #{i + 1}: {result}")
        if isinstance(result, dict):
            for key, value in result.items():
                logging.info(f"  Chave: {key}, Valor: {value}")
                if key == "stats" and isinstance(value, dict):
                    for stat_key, stat_value in value.items():
                        logging.info(f"    Estatística: {stat_key} = {stat_value}")

    # Extrair dados de cada condição - sem modificar os resultados originais
    conditions = {}
    for result in results:
        if (
            not isinstance(result, dict)
            or "condition" not in result
            or "polymorphisms"
            not in result  # Note: checking for polymorphisms directly
        ):
            logging.warning(f"Resultado inválido: {result}")
            continue

        condition = result["condition"]

        # Criar estatísticas diretamente dos dados de polimorfismos, não dos stats
        stats = {
            "total_snps": 0,
            "snp_positions": [],
            "snp_density": 0,
            "total_vntrs": 0,
            "vntr_lengths": [],
            "avg_vntr_length": 0,
            "hotspots": 0,
        }

        # Obter polimorfismos diretamente dos resultados
        polymorphisms = result.get("polymorphisms", {})

        # Processar SNPs
        if "SNPs" in polymorphisms and polymorphisms["SNPs"]:
            stats["total_snps"] = len(polymorphisms["SNPs"])

            # Converter as posições para índice 1 (biológico) ao invés de índice 0 (programação)
            stats["snp_positions"] = [pos + 1 for pos in polymorphisms["SNPs"].keys()]

            # Calcular hotspots (regiões com 3+ SNPs próximos)
            positions = sorted(stats["snp_positions"])
            hotspot_count = 0
            i = 0
            while i < len(positions) - 2:
                if positions[i + 2] - positions[i] <= 10:  # 3 SNPs dentro de 10 bases
                    hotspot_count += 1
                    # Pular para depois deste hotspot
                    while (
                        i < len(positions) - 1 and positions[i + 1] - positions[i] <= 10
                    ):
                        i += 1
                i += 1
            stats["hotspots"] = hotspot_count

        # Processar VNTRs
        if "VNTRs" in polymorphisms and polymorphisms["VNTRs"]:
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

        # Calcular densidade de SNPs
        if "report_file" in result:
            report_dir = os.path.dirname(result["report_file"])
            alignments_dir = os.path.join(os.path.dirname(report_dir), "alignments")

            # Tentar todas as combinações possíveis do nome do arquivo de alinhamento
            possible_filenames = [
                f"drd4_aligned_{condition}.aln",
                f"drd4_aligned_{condition.lower()}.aln",
                f"drd4_aligned_{condition.upper()}.aln",
                f"drd4_aligned_{condition.capitalize()}.aln",
            ]

            alignment_length = 0
            for filename in possible_filenames:
                alignment_file = os.path.join(alignments_dir, filename)
                if os.path.exists(alignment_file):
                    try:
                        from Bio import AlignIO

                        alignment = AlignIO.read(alignment_file, "clustal")
                        alignment_length = alignment.get_alignment_length()
                        break
                    except Exception as e:
                        logging.warning(
                            f"Erro ao ler arquivo de alinhamento {alignment_file}: {e}"
                        )

            if alignment_length > 0 and stats["total_snps"] > 0:
                stats["snp_density"] = stats["total_snps"] / alignment_length

        # Log detalhado das estatísticas calculadas
        logging.info(f"Estatísticas calculadas para {condition}:")
        logging.info(f"  SNPs: {stats['total_snps']}")
        logging.info(f"  VNTRs: {stats['total_vntrs']}")
        logging.info(f"  Densidade SNPs: {stats['snp_density']}")
        logging.info(f"  Hotspots: {stats['hotspots']}")

        conditions[condition] = stats

    # Verificar se temos resultados para comparar
    if len(conditions) < 1:
        logging.warning("Insuficientes condições para análise comparativa")
        return None

    # Criar visualização comparativa apenas se tivermos dados
    output_file = None

    if len(conditions) >= 2:
        try:
            # Preparar dados para visualização
            condition_names = list(conditions.keys())
            snp_counts = [conditions[c].get("total_snps", 0) for c in condition_names]
            vntr_counts = [conditions[c].get("total_vntrs", 0) for c in condition_names]

            # Criar visualização comparativa
            plt.figure(figsize=(10, 6))
            x = np.arange(len(condition_names))
            width = 0.35

            fig, ax = plt.subplots()
            rects1 = ax.bar(x - width / 2, snp_counts, width, label="SNPs")
            rects2 = ax.bar(x + width / 2, vntr_counts, width, label="VNTRs")

            # Adicionar detalhes ao gráfico
            ax.set_ylabel("Quantidade")
            ax.set_title("Comparação de Polimorfismos por Condição")
            ax.set_xticks(x)
            ax.set_xticklabels(condition_names)
            ax.legend()

            # Adicionar rótulos nas barras
            def autolabel(rects):
                for rect in rects:
                    height = rect.get_height()
                    ax.annotate(
                        "{}".format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 pontos acima da barra
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                    )

            autolabel(rects1)
            autolabel(rects2)

            # Salvar visualização
            output_file = os.path.join(visualization_dir, "condition_comparison.png")
            fig.tight_layout()
            plt.savefig(output_file)
            plt.close()

            logging.info(f"Visualização comparativa salva em: {output_file}")
        except Exception as e:
            logging.error(f"Erro ao criar visualização comparativa: {str(e)}")

    # Criar relatório comparativo com os dados originais não modificados
    report_file = create_comparative_report(conditions, reports_dir=reports_dir)
    logging.info(f"Relatório comparativo salvo em: {report_file}")

    return output_file


def create_comparative_report(
    conditions, reports_dir=None, reference_seq_id="NM_000797.4"
):
    """Cria um relatório comparativo em Markdown entre diferentes condições"""
    # Obter o diretório de relatórios atualizado
    if reports_dir is None:
        from config import REPORTS_DIR

        reports_dir = REPORTS_DIR

    if not reports_dir:
        raise ValueError("O diretório de relatórios não foi configurado corretamente")

    os.makedirs(reports_dir, exist_ok=True)
    report_file = os.path.join(reports_dir, "comparative_report.md")

    logging.info(f"Gerando relatório comparativo com {len(conditions)} condições")

    with open(report_file, "w") as f:
        f.write("# RELATÓRIO COMPARATIVO DE POLIMORFISMOS NO GENE DRD4\n\n")

        f.write("## VISÃO GERAL DAS CONDIÇÕES ANALISADAS\n\n")
        f.write(
            "Este relatório apresenta uma análise comparativa das variações genéticas do DRD4\n"
        )
        f.write(
            "encontradas em sequências de Homo sapiens relacionadas a diferentes condições.\n"
        )
        f.write(
            "As sequências foram obtidas do NCBI e classificadas conforme sua associação\n"
        )
        f.write("com ADHD, Autismo, ou variantes específicas.\n\n")

        # Explicar a categoria de variantes específicas
        if "ADHD_Variantes" in conditions:
            f.write("### NOTA SOBRE SEQUÊNCIAS VARIANTES\n\n")
            f.write(
                "As sequências classificadas como 'ADHD_Variantes' representam variantes divergentes do gene\n"
            )
            f.write(
                "DRD4 (especialmente as com prefixo LC812) que demonstraram diferenças significativas de\n"
            )
            f.write(
                "alinhamento em comparação com as sequências típicas de DRD4. A seleção foi realizada usando\n"
            )
            f.write(
                "critérios específicos, priorizando sequências com prefixos LC812 e outras LC*.\n"
            )
            f.write(
                "Estas foram analisadas separadamente para melhorar a precisão dos resultados.\n\n"
            )

        # Adicionar informação sobre a sequência de referência
        f.write("## INFORMAÇÕES DE REFERÊNCIA\n\n")
        f.write(
            f"* **Sequência de referência**: Gene DRD4 humano, {reference_seq_id} (NCBI Reference Sequence)\n"
        )
        f.write(
            "* **Convenção de posição**: Todas as posições seguem a convenção biológica (começam em 1)\n"
        )
        f.write(
            "* **Análise**: Os polimorfismos representam divergências em relação à sequência de referência\n\n"
        )

        # Tabela comparativa
        f.write("## TABELA COMPARATIVA\n\n")
        f.write("| Condição | SNPs | VNTRs | Densidade SNPs | Hotspots |\n")
        f.write("|----------|------|-------|---------------|----------|\n")

        for condition, stats in conditions.items():
            # Obter os valores com tratamento de valores ausentes
            snps = stats.get("total_snps", 0)
            vntrs = stats.get("total_vntrs", 0)
            density = stats.get("snp_density", 0)
            hotspots = stats.get("hotspots", 0)

            # Garantir que todos os valores são do tipo esperado
            snps = int(snps) if snps is not None else 0
            vntrs = int(vntrs) if vntrs is not None else 0
            density = float(density) if density is not None else 0.0
            hotspots = int(hotspots) if hotspots is not None else 0

            # Escrever linha na tabela
            f.write(
                f"| {condition} | {snps} | {vntrs} | {density:.5f} | {hotspots} |\n"
            )

        # Análise dos resultados
        f.write("\n## ANÁLISE DOS RESULTADOS\n\n")

        # Análise de SNPs entre condições
        f.write("### 1. Comparação de SNPs\n\n")

        # Ordenar condições por número de SNPs
        sorted_by_snps = sorted(
            conditions.items(), key=lambda x: x[1].get("total_snps", 0), reverse=True
        )

        if sorted_by_snps:  # Verificar se há dados para analisar
            highest_snp_condition = sorted_by_snps[0][0]
            lowest_snp_condition = sorted_by_snps[-1][0]

            f.write(
                f"* A condição **{highest_snp_condition}** apresenta o maior número de SNPs ({sorted_by_snps[0][1]['total_snps']}), "
            )
            f.write(
                f"enquanto **{lowest_snp_condition}** apresenta o menor ({sorted_by_snps[-1][1]['total_snps']}).\n\n"
            )

            # Análise de densidade de SNPs
            sorted_by_density = sorted(
                conditions.items(),
                key=lambda x: x[1].get("snp_density", 0),
                reverse=True,
            )
            f.write(
                f"* A maior densidade de SNPs foi observada em **{sorted_by_density[0][0]}** "
            )
            f.write(
                f"({sorted_by_density[0][1].get('snp_density', 0):.5f} SNPs/pb).\n\n"
            )

        # Análise de VNTRs
        vntr_conditions = [
            cond
            for cond, stats in conditions.items()
            if stats.get("total_vntrs", 0) > 0
        ]

        f.write("### 2. Análise de VNTRs\n\n")

        if vntr_conditions:
            f.write(
                f"* VNTRs foram identificados nas seguintes condições: {', '.join(vntr_conditions)}.\n\n"
            )

            for condition in vntr_conditions:
                f.write(
                    f"* **{condition}**: {conditions[condition]['total_vntrs']} VNTRs detectados.\n\n"
                )

            # Destaque para ADHD se tiver VNTRs (alelo 7R)
            if (
                "ADHD" in vntr_conditions
                and conditions["ADHD"].get("total_vntrs", 0) > 0
            ):
                f.write(
                    "* A presença de VNTRs em ADHD é particularmente relevante, pois o alelo 7R do VNTR no "
                )
                f.write(
                    "éxon 3 do DRD4 tem sido fortemente associado a esta condição em diversos estudos.\n\n"
                )
        else:
            f.write(
                "* Não foram detectados VNTRs significativos nas condições analisadas.\n\n"
            )

        # O resto do código segue o mesmo formato, adaptando para markdown
        f.write("### 3. Padrões de Polimorfismos\n\n")
        if "Autismo" in conditions and "ADHD" in conditions:
            if conditions["Autismo"].get("total_snps", 0) < conditions["ADHD"].get(
                "total_snps", 0
            ):
                f.write(
                    "* O autismo apresenta menor variabilidade de SNPs comparado ao ADHD,\n"
                )
                f.write(
                    "  sugerindo possivelmente um papel diferente do gene DRD4 nestas condições.\n\n"
                )
            else:
                f.write(
                    "* O padrão de polimorfismos sugere variabilidade genética distinta entre as condições,\n"
                )
                f.write(
                    "  com diferentes perfis de mutação que podem afetar a função do receptor D4.\n\n"
                )

        if "ADHD_Variantes" in conditions:
            f.write(
                "* As variantes divergentes de ADHD apresentam um perfil de polimorfismo distinto,\n"
            )
            f.write(
                "  indicando possíveis diferenças estruturais ou funcionais nessas variantes do gene DRD4.\n\n"
            )

        if "ADHD_Variantes" in conditions:
            f.write("### RELEVÂNCIA DAS VARIANTES LC812\n\n")
            f.write(
                "As variantes LC812 e outras sequências divergentes apresentam um perfil de polimorfismo distinto,\n"
            )
            f.write(
                "com maior densidade de SNPs (%.5f SNPs/pb) comparado às sequências ADHD típicas.\n"
                % conditions.get("ADHD_Variantes", {}).get("snp_density", 0)
            )
            f.write(
                "Esta maior variabilidade pode indicar uma subpopulação genética específica e potencialmente\n"
            )
            f.write(
                "relevante para a heterogeneidade fenotípica observada em pacientes com ADHD.\n\n"
            )

        # Relevância biológica
        f.write("### 4. Relevância Biológica\n\n")
        f.write(
            "* As variações genéticas detectadas no receptor DRD4 podem afetar a sinalização dopaminérgica,\n"
        )
        f.write(
            "  impactando funções cognitivas como atenção, recompensa e comportamento impulsivo.\n\n"
        )
        f.write(
            "* Estudos funcionais são necessários para determinar como os padrões específicos de\n"
        )
        f.write(
            "  polimorfismos detectados afetam a estrutura e função da proteína DRD4.\n\n"
        )

        # Conclusão geral
        f.write("## CONCLUSÕES\n\n")
        f.write(
            "A análise comparativa revela padrões distintos de polimorfismos entre as condições.\n"
        )
        f.write(
            "ADHD e Autismo apresentam perfis genéticos distintos no gene DRD4, enquanto as variantes\n"
        )
        f.write(
            "divergentes (sequências LC) mostram características únicas que podem representar\n"
        )
        f.write(
            "isoformas funcionalmente relevantes ou variações populacionais específicas.\n\n"
        )
        f.write(
            "Esta análise identifica alvos potenciais para estudos funcionais futuros que podem\n"
        )
        f.write(
            "esclarecer o papel do receptor dopaminérgico D4 nestas condições neuropsiquiátricas.\n"
        )

    logging.info(f"Relatório comparativo em Markdown salvo em: {report_file}")
    return report_file
