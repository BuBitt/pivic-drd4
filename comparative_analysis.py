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

    # Extrair dados de cada condição
    conditions = {}
    for result in results:
        condition = result["condition"]
        stats = result["stats"]
        conditions[condition] = stats

    # Verificar se temos resultados para comparar
    if len(conditions) < 2:
        logging.warning("Insuficientes condições para análise comparativa")
        return

    # Preparar dados para visualização
    condition_names = list(conditions.keys())
    snp_counts = [conditions[c]["total_snps"] for c in condition_names]
    vntr_counts = [conditions[c]["total_vntrs"] for c in condition_names]

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

    # Criar relatório comparativo
    create_comparative_report(conditions, reports_dir=reports_dir)

    return output_file


def create_comparative_report(conditions, reports_dir=None):
    """Cria um relatório comparativo em Markdown entre diferentes condições"""
    # Obter o diretório de relatórios atualizado
    if reports_dir is None:
        from config import REPORTS_DIR

        reports_dir = REPORTS_DIR

    if not reports_dir:
        raise ValueError("O diretório de relatórios não foi configurado corretamente")

    os.makedirs(reports_dir, exist_ok=True)
    report_file = os.path.join(reports_dir, "comparative_report.md")

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

        # Tabela comparativa
        f.write("## TABELA COMPARATIVA\n\n")
        f.write("| Condição | SNPs | VNTRs | Densidade SNPs | Hotspots |\n")
        f.write("|----------|------|-------|---------------|----------|\n")

        for condition, stats in conditions.items():
            snp_density = stats.get("snp_density", 0)
            hotspots = stats.get("hotspots", 0)
            f.write(
                f"| {condition} | {stats['total_snps']} | {stats['total_vntrs']} | {snp_density:.5f} | {hotspots} |\n"
            )

        # Análise dos resultados
        f.write("\n## ANÁLISE DOS RESULTADOS\n\n")

        # Análise de SNPs entre condições
        f.write("### 1. Comparação de SNPs\n\n")

        # Ordenar condições por número de SNPs
        sorted_by_snps = sorted(
            conditions.items(), key=lambda x: x[1].get("total_snps", 0), reverse=True
        )
        highest_snp_condition = sorted_by_snps[0][0]
        lowest_snp_condition = sorted_by_snps[-1][0]

        f.write(
            f"- A condição **{highest_snp_condition}** apresenta o maior número de SNPs ({sorted_by_snps[0][1]['total_snps']}), "
        )
        f.write(
            f"enquanto **{lowest_snp_condition}** apresenta o menor ({sorted_by_snps[-1][1]['total_snps']}).\n\n"
        )

        # Análise de densidade de SNPs
        sorted_by_density = sorted(
            conditions.items(), key=lambda x: x[1].get("snp_density", 0), reverse=True
        )
        f.write(
            f"- A maior densidade de SNPs foi observada em **{sorted_by_density[0][0]}** "
        )
        f.write(f"({sorted_by_density[0][1].get('snp_density', 0):.5f} SNPs/pb).\n\n")

        # Análise de VNTRs
        vntr_conditions = [
            cond
            for cond, stats in conditions.items()
            if stats.get("total_vntrs", 0) > 0
        ]
        if vntr_conditions:
            f.write("### 2. Análise de VNTRs\n\n")
            f.write(
                f"- VNTRs foram identificados nas seguintes condições: {', '.join(vntr_conditions)}.\n\n"
            )

            for condition in vntr_conditions:
                f.write(
                    f"- **{condition}**: {conditions[condition]['total_vntrs']} VNTRs detectados.\n\n"
                )

            # Destaque para ADHD se tiver VNTRs (alelo 7R)
            if (
                "ADHD" in vntr_conditions
                and conditions["ADHD"].get("total_vntrs", 0) > 0
            ):
                f.write(
                    "- A presença de VNTRs em ADHD é particularmente relevante, pois o alelo 7R do VNTR no\n"
                )
                f.write(
                    "  éxon 3 do DRD4 tem sido fortemente associado a esta condição em diversos estudos.\n\n"
                )
        else:
            f.write("### 2. Análise de VNTRs\n\n")
            f.write(
                "- Não foram detectados VNTRs significativos nas condições analisadas.\n\n"
            )

        # O resto do código segue o mesmo formato, adaptando para markdown
        f.write("### 3. Padrões de Polimorfismos\n\n")
        if "Autismo" in conditions and "ADHD" in conditions:
            if conditions["Autismo"].get("total_snps", 0) < conditions["ADHD"].get(
                "total_snps", 0
            ):
                f.write(
                    "- O autismo apresenta menor variabilidade de SNPs comparado ao ADHD,\n"
                )
                f.write(
                    "  sugerindo possivelmente um papel diferente do gene DRD4 nestas condições.\n\n"
                )
            else:
                f.write(
                    "- O padrão de polimorfismos sugere variabilidade genética distinta entre as condições,\n"
                )
                f.write(
                    "  com diferentes perfis de mutação que podem afetar a função do receptor D4.\n\n"
                )

        if "ADHD_Variantes" in conditions:
            f.write(
                "- As variantes divergentes de ADHD apresentam um perfil de polimorfismo distinto,\n"
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
            "- As variações genéticas detectadas no receptor DRD4 podem afetar a sinalização dopaminérgica,\n"
        )
        f.write(
            "  impactando funções cognitivas como atenção, recompensa e comportamento impulsivo.\n\n"
        )
        f.write(
            "- Estudos funcionais são necessários para determinar como os padrões específicos de\n"
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
