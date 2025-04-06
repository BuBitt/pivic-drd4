import os
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def compare_alignment_results(
    primary_results, secondary_results, primary_tool, secondary_tool, reports_dir=None
):
    """
    Compara os resultados de duas ferramentas de alinhamento diferentes e gera relatório em Markdown.

    Args:
        primary_results: Resultados da ferramenta preferencial
        secondary_results: Resultados da ferramenta secundária
        primary_tool: Nome da ferramenta principal (ex: 'mafft')
        secondary_tool: Nome da ferramenta secundária (ex: 'clustalw')
        reports_dir: Diretório para salvar relatórios

    Returns:
        Caminho do relatório de comparação
    """
    if reports_dir is None:
        from config import REPORTS_DIR

        reports_dir = REPORTS_DIR

    # Garantir que o diretório existe
    os.makedirs(reports_dir, exist_ok=True)

    # Criar mapeamento de resultados por condição
    primary_by_condition = {r["condition"]: r for r in primary_results}
    secondary_by_condition = {r["condition"]: r for r in secondary_results}

    # Condições disponíveis em ambos
    common_conditions = set(primary_by_condition.keys()) & set(
        secondary_by_condition.keys()
    )

    if not common_conditions:
        logging.warning(
            "Não há condições comuns para comparar entre as ferramentas de alinhamento"
        )
        return None

    # Relatório de saída
    report_path = os.path.join(reports_dir, f"alignment_tools_comparison.md")

    with open(report_path, "w") as f:
        f.write("# RELATÓRIO COMPARATIVO DE FERRAMENTAS DE ALINHAMENTO\n\n")

        f.write(f"**Ferramenta principal**: {primary_tool.upper()}\n\n")
        f.write(f"**Ferramenta secundária**: {secondary_tool.upper()}\n\n")

        f.write("## OBSERVAÇÕES IMPORTANTES\n\n")
        f.write(
            "Diferenças entre as ferramentas são esperadas devido aos diferentes algoritmos usados.\n"
        )
        f.write(
            "O MAFFT geralmente proporciona alinhamentos mais precisos para sequências diversas.\n"
        )
        f.write(
            "O ClustalW é mais tradicional e pode ser preferível para sequências altamente similares.\n\n"
        )

        f.write("## SUMÁRIO COMPARATIVO\n\n")

        # Usar tabelas HTML para maior compatibilidade quando a tabela precisa ser complexa
        f.write("<table>\n")
        f.write("<tr>\n")
        f.write("  <th>Condição</th>\n")
        f.write("  <th>Ferramenta</th>\n")
        f.write("  <th>SNPs</th>\n")
        f.write("  <th>VNTRs</th>\n")
        f.write("  <th>Diferença SNPs</th>\n")
        f.write("  <th>Diferença VNTRs</th>\n")
        f.write("</tr>\n")

        total_snp_diff = 0
        total_vntr_diff = 0

        for condition in sorted(common_conditions):
            p_stats = primary_by_condition[condition]["stats"]
            s_stats = secondary_by_condition[condition]["stats"]

            snp_diff = p_stats["total_snps"] - s_stats["total_snps"]
            vntr_diff = p_stats["total_vntrs"] - s_stats["total_vntrs"]

            total_snp_diff += abs(snp_diff)
            total_vntr_diff += abs(vntr_diff)

            # Linha para ferramenta primária
            f.write("<tr>\n")
            f.write(f"  <td rowspan='2'>{condition}</td>\n")
            f.write(f"  <td><b>{primary_tool}</b></td>\n")
            f.write(f"  <td>{p_stats['total_snps']}</td>\n")
            f.write(f"  <td>{p_stats['total_vntrs']}</td>\n")
            f.write(f"  <td></td>\n")
            f.write(f"  <td></td>\n")
            f.write("</tr>\n")

            # Linha para ferramenta secundária
            f.write("<tr>\n")
            f.write(f"  <td><b>{secondary_tool}</b></td>\n")
            f.write(f"  <td>{s_stats['total_snps']}</td>\n")
            f.write(f"  <td>{s_stats['total_vntrs']}</td>\n")
            f.write(f"  <td>{'+' if snp_diff > 0 else ''}{snp_diff}</td>\n")
            f.write(f"  <td>{'+' if vntr_diff > 0 else ''}{vntr_diff}</td>\n")
            f.write("</tr>\n")

        f.write("</table>\n\n")

        # O resto do relatório segue o mesmo padrão
        f.write("\n## ANÁLISE DETALHADA\n\n")

        avg_snp_diff = total_snp_diff / len(common_conditions)
        avg_vntr_diff = total_vntr_diff / len(common_conditions)

        f.write(f"Diferença média absoluta em SNPs: {avg_snp_diff:.2f}\n")
        f.write(f"Diferença média absoluta em VNTRs: {avg_vntr_diff:.2f}\n\n")

        for condition in sorted(common_conditions):
            f.write(f"\n### Análise da condição: {condition}\n\n")

            p_stats = primary_by_condition[condition]["stats"]
            s_stats = secondary_by_condition[condition]["stats"]

            # Comparar posições de SNPs
            p_snps = set(p_stats["snp_positions"])
            s_snps = set(s_stats["snp_positions"])
            common_snps = p_snps & s_snps
            p_unique_snps = p_snps - s_snps
            s_unique_snps = s_snps - p_snps

            f.write(f"SNPs encontrados por ambas ferramentas: {len(common_snps)}\n")
            if common_snps:
                f.write(f"  Posições: {sorted(common_snps)}\n")

            f.write(f"SNPs exclusivos do {primary_tool}: {len(p_unique_snps)}\n")
            if p_unique_snps and len(p_unique_snps) < 20:
                f.write(f"  Posições: {sorted(p_unique_snps)}\n")

            f.write(f"SNPs exclusivos do {secondary_tool}: {len(s_unique_snps)}\n")
            if s_unique_snps and len(s_unique_snps) < 20:
                f.write(f"  Posições: {sorted(s_unique_snps)}\n")

            # Comparar VNTRs (mais complexo devido à natureza dos intervalos)
            f.write(
                f"\nVNTRs detectados por {primary_tool}: {p_stats['total_vntrs']}\n"
            )
            if p_stats["vntr_lengths"]:
                f.write(f"  Comprimentos: {p_stats['vntr_lengths']}\n")

            f.write(
                f"VNTRs detectados por {secondary_tool}: {s_stats['total_vntrs']}\n"
            )
            if s_stats["vntr_lengths"]:
                f.write(f"  Comprimentos: {s_stats['vntr_lengths']}\n")

        # Recomendação
        f.write("\n## RECOMENDAÇÕES\n\n")
        f.write(
            "Com base na literatura científica atual e nas diferenças observadas:\n\n"
        )
        f.write(
            "1. O MAFFT é geralmente recomendado para análises de polimorfismos em genes como o DRD4 devido à:\n"
        )
        f.write("   - Maior precisão no alinhamento de sequências divergentes\n")
        f.write("   - Melhor detecção de regiões de alta variabilidade como VNTRs\n")
        f.write("   - Algoritmos mais modernos e eficientes\n\n")

        f.write("2. Para maior confiabilidade:\n")
        f.write(
            "   - Considere os SNPs detectados por ambas as ferramentas como altamente confiáveis\n"
        )
        f.write(
            "   - Examine VNTRs detectados por pelo menos uma ferramenta para confirmação visual\n"
        )
        f.write("   - Documente a ferramenta usada ao relatar resultados\n\n")

        f.write("3. Para regiões com alta divergência (VNTRs):\n")
        f.write("   - O MAFFT tende a produzir alinhamentos mais precisos\n")
        f.write(
            "   - Validação adicional pode ser necessária para polimorfismos críticos\n"
        )

    return report_path
