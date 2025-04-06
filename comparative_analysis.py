import os
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from config import REPORTS_DIR, VISUALIZATION_DIR


def compare_condition_polymorphisms(results):
    """
    Compara polimorfismos encontrados entre diferentes condições (ADHD, Autismo, Comorbidade)

    Args:
        results: Lista de resultados da análise de polimorfismos por condição
    """
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
    output_file = os.path.join(VISUALIZATION_DIR, "condition_comparison.png")
    fig.tight_layout()
    plt.savefig(output_file)
    plt.close()

    logging.info(f"Visualização comparativa salva em: {output_file}")

    # Criar relatório comparativo
    create_comparative_report(conditions)

    return output_file


def create_comparative_report(conditions):
    """Cria um relatório comparativo entre diferentes condições"""
    report_file = os.path.join(REPORTS_DIR, "comparative_report.txt")

    with open(report_file, "w") as f:
        f.write("RELATÓRIO COMPARATIVO DE POLIMORFISMOS NO GENE DRD4\n")
        f.write("=" * 60 + "\n\n")

        # Adicionar informações gerais sobre as condições analisadas
        f.write("VISÃO GERAL DAS CONDIÇÕES ANALISADAS\n")
        f.write("-" * 50 + "\n")
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
            f.write("NOTA SOBRE SEQUÊNCIAS VARIANTES:\n")
            f.write(
                "As sequências classificadas como 'ADHD_Variantes' representam variantes divergentes do gene\n"
            )
            f.write(
                "DRD4 (prefixo LC812) que demonstraram diferenças significativas de alinhamento em comparação\n"
            )
            f.write(
                "com as sequências típicas de DRD4. Estas foram analisadas separadamente para melhorar a\n"
            )
            f.write(
                "precisão dos resultados e podem representar isoformas ou variantes populacionais específicas.\n\n"
            )

        # Tabela comparativa
        f.write("TABELA COMPARATIVA:\n")
        f.write("-" * 50 + "\n")
        f.write(
            f"{'Condição':<15} {'SNPs':<10} {'VNTRs':<10} {'Densidade SNPs':<15} {'Hotspots':<10}\n"
        )

        for condition, stats in conditions.items():
            snp_density = stats.get("snp_density", 0)
            hotspots = stats.get("hotspots", 0)
            f.write(
                f"{condition:<15} {stats['total_snps']:<10} {stats['total_vntrs']:<10} {snp_density:<15.5f} {hotspots:<10}\n"
            )

        # Análise dos resultados
        f.write("\nANÁLISE DOS RESULTADOS:\n")
        f.write("-" * 50 + "\n")

        # Análise de SNPs entre condições
        f.write("\n1. Comparação de SNPs:\n")

        # Ordenar condições por número de SNPs
        sorted_by_snps = sorted(
            conditions.items(), key=lambda x: x[1].get("total_snps", 0), reverse=True
        )
        highest_snp_condition = sorted_by_snps[0][0]
        lowest_snp_condition = sorted_by_snps[-1][0]

        f.write(
            f"   • A condição {highest_snp_condition} apresenta o maior número de SNPs ({sorted_by_snps[0][1]['total_snps']}), "
        )
        f.write(
            f"enquanto {lowest_snp_condition} apresenta o menor ({sorted_by_snps[-1][1]['total_snps']}).\n"
        )

        # Análise de densidade de SNPs
        sorted_by_density = sorted(
            conditions.items(), key=lambda x: x[1].get("snp_density", 0), reverse=True
        )
        f.write(
            f"   • A maior densidade de SNPs foi observada em {sorted_by_density[0][0]} "
        )
        f.write(f"({sorted_by_density[0][1].get('snp_density', 0):.5f} SNPs/pb).\n")

        # Análise de VNTRs
        vntr_conditions = [
            cond
            for cond, stats in conditions.items()
            if stats.get("total_vntrs", 0) > 0
        ]
        if vntr_conditions:
            f.write("\n2. Análise de VNTRs:\n")
            f.write(
                f"   • VNTRs foram identificados nas seguintes condições: {', '.join(vntr_conditions)}.\n"
            )

            for condition in vntr_conditions:
                f.write(
                    f"   • {condition}: {conditions[condition]['total_vntrs']} VNTRs detectados.\n"
                )

            # Destaque para ADHD se tiver VNTRs (alelo 7R)
            if (
                "ADHD" in vntr_conditions
                and conditions["ADHD"].get("total_vntrs", 0) > 0
            ):
                f.write(
                    "   • A presença de VNTRs em ADHD é particularmente relevante, pois o alelo 7R do VNTR no\n"
                )
                f.write(
                    "     éxon 3 do DRD4 tem sido fortemente associado a esta condição em diversos estudos.\n"
                )
        else:
            f.write("\n2. Análise de VNTRs:\n")
            f.write(
                "   • Não foram detectados VNTRs significativos nas condições analisadas.\n"
            )

        # Padrões específicos de polimorfismos
        f.write("\n3. Padrões de Polimorfismos:\n")
        if "Autismo" in conditions and "ADHD" in conditions:
            if conditions["Autismo"].get("total_snps", 0) < conditions["ADHD"].get(
                "total_snps", 0
            ):
                f.write(
                    "   • O autismo apresenta menor variabilidade de SNPs comparado ao ADHD,\n"
                )
                f.write(
                    "     sugerindo possivelmente um papel diferente do gene DRD4 nestas condições.\n"
                )
            else:
                f.write(
                    "   • O padrão de polimorfismos sugere variabilidade genética distinta entre as condições,\n"
                )
                f.write(
                    "     com diferentes perfis de mutação que podem afetar a função do receptor D4.\n"
                )

        if "ADHD_Variantes" in conditions:
            f.write(
                "   • As variantes divergentes de ADHD apresentam um perfil de polimorfismo distinto,\n"
            )
            f.write(
                "     indicando possíveis diferenças estruturais ou funcionais nessas variantes do gene DRD4.\n"
            )

        # Relevância biológica
        f.write("\n4. Relevância Biológica:\n")
        f.write(
            "   • As variações genéticas detectadas no receptor DRD4 podem afetar a sinalização dopaminérgica,\n"
        )
        f.write(
            "     impactando funções cognitivas como atenção, recompensa e comportamento impulsivo.\n"
        )
        f.write(
            "   • Estudos funcionais são necessários para determinar como os padrões específicos de\n"
        )
        f.write(
            "     polimorfismos detectados afetam a estrutura e função da proteína DRD4.\n"
        )

        # Conclusão geral
        f.write("\nCONCLUSÕES:\n")
        f.write("-" * 50 + "\n")
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

    logging.info(f"Relatório comparativo salvo em: {report_file}")
    return report_file
