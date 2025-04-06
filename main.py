import asyncio
import os
import logging
import time
import concurrent.futures
import argparse
from Bio import Entrez, SeqIO

# Importar funções de módulos específicos
from config import (
    setup_directories,
    setup_logging,
    CLUSTALW_PARAMS,
    POLYMORPHISM_CONFIG,
)
from ncbi_utils import fetch_reference_sequence
from alignment_utils import align_sequences_for_condition
from polymorphism_analysis import detect_polymorphisms
from report_utils import (
    generate_report,
    analyze_polymorphism_report,
)
from async_utils import run_all_async_tasks


def process_condition(
    reference_file,
    filtered_file,
    condition,
    alignment_tool="clustalw",
    reports_dir=None,
    visualization_dir=None,
):
    """
    Processa uma condição específica (ADHD ou autism) - função isolada para permitir
    processamento paralelo.
    """
    try:
        if os.path.getsize(filtered_file) == 0:
            logging.warning(
                f"Nenhuma sequência relacionada a {condition} foi encontrada para alinhamento."
            )
            return None

        logging.info(
            f"Processando alinhamento e análise para {condition} usando {alignment_tool}..."
        )
        start_time = time.time()

        # Realizar alinhamento com base na ferramenta selecionada
        aligned_file = align_sequences_for_condition(
            reference_file,
            filtered_file,
            condition,
            CLUSTALW_PARAMS.get(condition, {}),
            alignment_tool=alignment_tool,  # Passar a ferramenta selecionada
            separate_by_condition=True,  # Garantir separação por condição
        )

        # Obter configuração específica para detecção de polimorfismos
        config = POLYMORPHISM_CONFIG.get(condition, {})

        # Detectar polimorfismos com parâmetros ajustados
        polymorphisms = detect_polymorphisms(
            aligned_file,
            min_vntr_length=config.get("min_vntr_length", 48),
            max_gap_ratio=config.get("max_gap_ratio", 0.4),
        )

        # Gerar relatório
        report_name = f"report_with_reference_{condition.lower()}.txt"
        report_file = generate_report(
            polymorphisms, report_name, condition, reports_dir=reports_dir
        )

        # Analisar estatisticamente
        stats = analyze_polymorphism_report(report_file)

        elapsed_time = time.time() - start_time
        logging.info(
            f"Análise de {condition} concluída em {elapsed_time:.2f} segundos. Estatísticas: {stats}"
        )

        return {"condition": condition, "stats": stats, "report_file": report_file}
    except Exception as e:
        logging.error(
            f"Erro ao processar condição {condition}: {str(e)}", exc_info=True
        )
        return None


def main(args=None):
    """Função principal do programa."""
    # Iniciar cronômetro
    start_time = time.time()

    # Processar argumentos de linha de comando
    parser = argparse.ArgumentParser(
        description="Análise de polimorfismos do gene DRD4"
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Pular o download de sequências se já existirem",
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Mostrar logs detalhados"
    )
    parser.add_argument(
        "--force-recalc",
        action="store_true",
        help="Forçar recálculo mesmo com cache disponível",
    )
    parser.add_argument(
        "--max-seq", type=int, default=250, help="Número máximo de sequências a baixar"
    )
    parser.add_argument(
        "--alignment-tool",
        choices=["clustalw", "mafft", "both"],
        default="mafft",  # Mudando o padrão para MAFFT
        help="Escolha a ferramenta de alinhamento: 'clustalw', 'mafft' ou 'both' (ambos para comparação)",
    )
    parser.add_argument(
        "--preferred-tool",
        choices=["clustalw", "mafft"],
        default="mafft",
        help="Em caso de 'both', qual ferramenta será usada para o relatório final",
    )

    if args is None:
        args = parser.parse_args()

    # Definir o diretório de saída
    if args.alignment_tool == "both":
        output_dir = f"drd4_analysis_{args.preferred_tool}_validated"
        secondary_output_dir = f"drd4_analysis_{('clustalw' if args.preferred_tool == 'mafft' else 'mafft')}"
    else:
        output_dir = f"drd4_analysis_{args.alignment_tool}"
        secondary_output_dir = None

    # Configurar diretórios e logging ANTES de qualquer operação
    setup_directories(output_dir=output_dir)
    setup_logging(verbose=args.verbose, output_dir=output_dir)

    # Verificar se as ferramentas necessárias estão disponíveis
    from alignment_utils import check_alignment_tools

    tools_status = check_alignment_tools()

    # Verificar se a ferramenta selecionada está disponível
    if not tools_status.get(args.alignment_tool, False):
        if args.alignment_tool == "clustalw":
            logging.warning(
                "ClustalW não está disponível. Por favor, verifique o caminho em config.py."
            )
            logging.info("Vai tentar usar MAFFT como alternativa...")
            args.alignment_tool = "mafft"
            if not tools_status.get("mafft", False):
                logging.error(
                    "Nem MAFFT está disponível. Por favor, instale pelo menos uma ferramenta de alinhamento."
                )
                return 1
        elif args.alignment_tool == "mafft":
            logging.error(
                "MAFFT não está disponível. Por favor, instale-o usando 'sudo apt-get install mafft' ou similar."
            )
            return 1

    # Obter as variáveis globais atualizadas após setup_directories
    from config import (
        REFERENCE_DIR,
        SEQUENCES_DIR,
        ALIGNMENTS_DIR,
        REPORTS_DIR,
        CACHE_DIR,
        VISUALIZATION_DIR,
    )

    logging.info(
        f"Iniciando análise de polimorfismos do DRD4 com {args.alignment_tool.upper()}..."
    )
    logging.info(f"Arquivos serão salvos no diretório: {output_dir}")

    try:
        # Baixar a sequência de referência passando o diretório configurado
        reference_file = fetch_reference_sequence(reference_dir=REFERENCE_DIR)

        # Executar todas as tarefas assíncronas em um único loop
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        try:
            logging.info(
                f"Buscando sequências com limite máximo de {args.max_seq} por condição..."
            )
            adhd_filtered_file, autism_filtered_file, all_sequences_file = (
                loop.run_until_complete(
                    run_all_async_tasks(
                        skip_download=args.skip_download,
                        force_recalc=args.force_recalc,
                        max_sequences=args.max_seq,
                        sequences_dir=SEQUENCES_DIR,
                        cache_dir=CACHE_DIR,
                    )
                )
            )
        finally:
            loop.close()

        logging.info(
            f"Todas as sequências DRD4 processadas. Arquivos criados em: {os.path.dirname(all_sequences_file)}"
        )

        # Processar as condições em paralelo para melhor eficiência
        conditions = [("ADHD", adhd_filtered_file), ("Autismo", autism_filtered_file)]

        # Adicionar variantes como uma condição separada, se existirem
        adhd_variants_file = os.path.join(
            os.path.dirname(adhd_filtered_file), "drd4_adhd_variants.fasta"
        )
        if (
            os.path.exists(adhd_variants_file)
            and os.path.getsize(adhd_variants_file) > 0
        ):
            conditions.append(("ADHD_Variantes", adhd_variants_file))
            logging.info(
                f"Incluindo análise de ADHD_Variantes do arquivo: {adhd_variants_file}"
            )
        else:
            logging.info(
                "Nenhuma variante ADHD específica encontrada para análise separada"
            )
            # Vamos tentar gerar um arquivo de variantes caso ainda não exista
            if not os.path.exists(adhd_variants_file):
                from async_utils import extract_adhd_variants

                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
                try:
                    variants_file = loop.run_until_complete(
                        extract_adhd_variants(adhd_filtered_file, SEQUENCES_DIR)
                    )
                    if variants_file and os.path.getsize(variants_file) > 0:
                        conditions.append(("ADHD_Variantes", variants_file))
                        logging.info(
                            f"Geradas variantes ADHD para análise: {variants_file}"
                        )
                finally:
                    loop.close()

        results = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
            # Submeter as tarefas de alinhamento e análise em paralelo
            future_to_condition = {
                executor.submit(
                    process_condition,
                    reference_file,
                    filtered_file,
                    condition,
                    args.alignment_tool,  # Passar a ferramenta selecionada
                    REPORTS_DIR,  # Passar o diretório de relatórios
                    VISUALIZATION_DIR,  # Passar o diretório de visualização
                ): condition
                for condition, filtered_file in conditions
            }

            for future in concurrent.futures.as_completed(future_to_condition):
                condition = future_to_condition[future]
                try:
                    result = future.result()
                    if result:
                        results.append(result)
                except Exception as e:
                    logging.error(f"Erro ao processar {condition}: {str(e)}")

        # Se a análise for com ambas as ferramentas
        if args.alignment_tool == "both":
            # Rodar primeiro com a ferramenta secundária
            secondary_tool = "clustalw" if args.preferred_tool == "mafft" else "mafft"
            logging.info(
                f"Realizando análise inicial com {secondary_tool.upper()} para validação cruzada..."
            )

            # Salvar o alignment_tool atual
            primary_tool = args.preferred_tool

            # Substituir temporariamente pelo secundário
            args.alignment_tool = secondary_tool

            # Configurar diretório para a ferramenta secundária
            setup_directories(output_dir=secondary_output_dir)

            # Executar análise com ferramenta secundária (simplificada)
            secondary_results = run_simplified_analysis(
                args, reference_file, adhd_filtered_file, autism_filtered_file
            )

            # Restaurar a configuração para continuar com a ferramenta preferida
            args.alignment_tool = primary_tool
            setup_directories(output_dir=output_dir)

            # Continuar com a análise principal
            # ...continuar o fluxo existente...

            # Depois de obter os resultados principais, comparar:
            if results and secondary_results:
                try:
                    from validation_utils import compare_alignment_results

                    comparison_report = compare_alignment_results(
                        results,
                        secondary_results,
                        args.preferred_tool,
                        secondary_tool,
                        REPORTS_DIR,
                    )

                    logging.info(
                        f"Relatório de comparação entre {args.preferred_tool} e {secondary_tool} "
                        f"gerado em: {comparison_report}"
                    )
                except Exception as e:
                    logging.error(
                        f"Erro ao comparar resultados das ferramentas: {str(e)}"
                    )

        # Adicionar alinhamento de sequências genéricas (não associadas a patologias)
        logging.info(
            f"Iniciando alinhamento de sequências genéricas usando {args.alignment_tool.upper()}..."
        )
        try:
            from alignment_utils import align_all_sequences

            align_all_sequences(
                all_sequences_file, reference_file, alignment_tool=args.alignment_tool
            )
            logging.info("Alinhamento de sequências genéricas concluído.")
        except Exception as e:
            logging.error(f"Erro ao alinhar sequências genéricas: {str(e)}")

        # Resumo dos resultados
        logging.info(f"Análise completa para {len(results)} condições")
        for result in results:
            logging.info(
                f"Condição: {result['condition']} - SNPs: {result['stats']['total_snps']}, "
                f"VNTRs: {result['stats']['total_vntrs']}"
            )

        # Adicionar análise comparativa se tivermos múltiplas condições
        if len(results) > 1:
            try:
                from comparative_analysis import compare_condition_polymorphisms

                compare_result = compare_condition_polymorphisms(
                    results,
                    visualization_dir=VISUALIZATION_DIR,
                    reports_dir=REPORTS_DIR,
                )
                if compare_result:
                    logging.info(
                        f"Análise comparativa concluída e salva em: {compare_result}"
                    )
            except Exception as e:
                logging.error(f"Erro na análise comparativa: {str(e)}")

        total_time = time.time() - start_time
        logging.info(f"Processo completo executado em {total_time:.2f} segundos")

    except Exception as e:
        logging.error(f"Erro durante a execução: {str(e)}", exc_info=True)
        return 1

    return 0


def run_simplified_analysis(args, reference_file, adhd_file, autism_file):
    """Executa uma versão simplificada da análise para fins de validação cruzada."""
    from config import REPORTS_DIR, VISUALIZATION_DIR

    results = []
    conditions = [("ADHD", adhd_file), ("Autismo", autism_file)]

    for condition, filtered_file in conditions:
        result = process_condition(
            reference_file,
            filtered_file,
            condition,
            args.alignment_tool,
            REPORTS_DIR,
            VISUALIZATION_DIR,
        )
        if result:
            results.append(result)

    return results


if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)
