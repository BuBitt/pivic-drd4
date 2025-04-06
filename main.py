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


def process_condition(reference_file, filtered_file, condition):
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

        logging.info(f"Processando alinhamento e análise para {condition}...")
        start_time = time.time()

        # Realizar alinhamento - apenas entre a sequência de referência e a condição específica
        aligned_file = align_sequences_for_condition(
            reference_file,
            filtered_file,
            condition,
            CLUSTALW_PARAMS.get(condition, {}),
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
        report_file = generate_report(polymorphisms, report_name, condition)

        # Analisar estatisticamente
        stats = analyze_polymorphism_report(report_file)

        # Removida a chamada para visualização
        # create_summary_visualization(polymorphisms, condition)

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

    if args is None:
        args = parser.parse_args()

    # Configurar diretórios e logging
    setup_directories()
    setup_logging(verbose=args.verbose)

    logging.info(
        f"Iniciando análise de polimorfismos do DRD4 com ClustalW (modo {'verboso' if args.verbose else 'normal'})..."
    )

    try:
        # Baixar a sequência de referência
        reference_file = fetch_reference_sequence()

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

        # Separar as sequências LC em um grupo distinto para análise de forma mais confiável
        adhd_variants_file = os.path.join(
            os.path.dirname(adhd_filtered_file), "drd4_adhd_variants.fasta"
        )
        try:
            # Verificar que os arquivos existem antes de processá-los
            if (
                os.path.exists(adhd_filtered_file)
                and os.path.getsize(adhd_filtered_file) > 0
            ):
                # Separar as sequências LC (variantes divergentes) das outras sequências ADHD
                regular_adhd_seqs = []
                variant_adhd_seqs = []

                for seq in SeqIO.parse(adhd_filtered_file, "fasta"):
                    if "LC812" in seq.id:  # Identificar as sequências LC divergentes
                        # Criar uma cópia da sequência para não modificar a original
                        variant_seq = seq[:]
                        # Adicionar 'Variante' ao ID para marcação, mas manter estrutura original
                        variant_seq.id = f"ADHD_Variante_{seq.id.replace('ADHD_', '')}"
                        variant_adhd_seqs.append(variant_seq)
                    else:
                        regular_adhd_seqs.append(seq)

                # Salvar as sequências variantes em arquivo separado se houver alguma
                if variant_adhd_seqs:
                    SeqIO.write(variant_adhd_seqs, adhd_variants_file, "fasta")
                    logging.info(
                        f"Separadas {len(variant_adhd_seqs)} sequências variantes LC em {adhd_variants_file}"
                    )

                    # Adicionar variantes como condição apenas se tiver sequências
                    conditions.append(("ADHD_Variantes", adhd_variants_file))

                # Reescrever o arquivo ADHD original apenas com sequências regulares se necessário
                if len(regular_adhd_seqs) != len(regular_adhd_seqs) + len(
                    variant_adhd_seqs
                ):
                    try:
                        temp_regular_file = adhd_filtered_file + ".tmp"
                        SeqIO.write(regular_adhd_seqs, temp_regular_file, "fasta")
                        os.replace(temp_regular_file, adhd_filtered_file)
                        logging.info(
                            f"Arquivo ADHD original atualizado com {len(regular_adhd_seqs)} sequências regulares"
                        )
                    except Exception as file_error:
                        logging.error(
                            f"Erro ao atualizar arquivo ADHD: {str(file_error)}"
                        )
        except Exception as e:
            logging.error(f"Erro ao processar sequências variantes: {str(e)}")
            # Continuar com as condições originais em caso de erro

        results = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
            # Submeter as tarefas de alinhamento e análise em paralelo
            future_to_condition = {
                executor.submit(
                    process_condition, reference_file, filtered_file, condition
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

        # Adicionar alinhamento de sequências genéricas (não associadas a patologias)
        logging.info(
            "Iniciando alinhamento de sequências genéricas (não patológicas)..."
        )
        try:
            from alignment_utils import align_all_sequences

            align_all_sequences(all_sequences_file, reference_file)
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

                compare_result = compare_condition_polymorphisms(results)
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


if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)
