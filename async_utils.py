import asyncio
import logging
import os
import time
import pickle
import hashlib
from Bio import Entrez, SeqIO
from tqdm import tqdm
from config import (
    SEQUENCES_DIR,
    EMAIL,
    CACHE_DIR,
    SEARCH_TERMS,
    USE_CACHE,
    CACHE_TTL,
    MIN_SEQ_LENGTH,
)

# Configurar email para NCBI
Entrez.email = EMAIL


def get_cache_path(search_term):
    """Gera um nome de arquivo de cache baseado no termo de busca"""
    # Criar um hash do termo para evitar nomes de arquivos inválidos
    hash_obj = hashlib.md5(search_term.encode())
    safe_term = hash_obj.hexdigest()
    return os.path.join(CACHE_DIR, f"{safe_term}.pickle")


async def fetch_sequence_async(seq_ids, webenv, query_key, start, end):
    """Função assíncrona para buscar sequências do NCBI."""
    try:
        logging.info(
            f"Buscando {len(seq_ids)} sequências a partir da posição {start}..."
        )
        handle = Entrez.efetch(
            db="nucleotide",
            rettype="fasta",
            retmode="text",
            webenv=webenv,
            query_key=query_key,
            retstart=start,
            retmax=end - start,
        )
        sequences = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return sequences
    except Exception as e:
        logging.error(
            f"Erro ao buscar sequências a partir da posição {start}: {str(e)}"
        )
        return []


async def fetch_condition_sequences_async(
    condition, max_sequences=200, force_recalc=False
):
    """Busca sequências diretamente relacionadas a uma condição específica."""
    if condition not in SEARCH_TERMS:
        logging.warning(f"Condição não reconhecida: {condition}")
        return []

    all_sequences = []

    # Buscar usando múltiplos termos de busca para cada condição
    for term in SEARCH_TERMS[condition]:
        cache_path = get_cache_path(term)

        # Verificar se há cache disponível
        if USE_CACHE and os.path.exists(cache_path) and not force_recalc:
            try:
                logging.info(f"Carregando sequências para '{term}' do cache...")
                with open(cache_path, "rb") as f:
                    cache_data = pickle.load(f)
                    cache_time = cache_data.get("time", 0)

                    # Usar cache se estiver dentro do TTL
                    if time.time() - cache_time < CACHE_TTL:
                        logging.info(f"Cache válido encontrado para '{term}'")
                        cached_sequences = cache_data.get("sequences", [])
                        all_sequences.extend(cached_sequences)
                        logging.info(
                            f"Carregadas {len(cached_sequences)} sequências do cache para '{term}'"
                        )
                        continue
                    else:
                        logging.info(
                            f"Cache expirado para '{term}', baixando novamente..."
                        )
            except Exception as e:
                logging.warning(f"Erro ao ler cache para '{term}': {e}")

        # Se não houver cache ou estiver expirado, buscar do NCBI
        logging.info(
            f"Buscando sequências relacionadas a '{term}' diretamente do NCBI..."
        )
        search_term = f"DRD4[Gene] AND Homo sapiens[Organism] AND {term}"

        try:
            handle = Entrez.esearch(
                db="nucleotide", term=search_term, retmax=max_sequences, usehistory="y"
            )
            record = Entrez.read(handle)
            handle.close()

            webenv = record["WebEnv"]
            query_key = record["QueryKey"]
            ids = record["IdList"]

            if not ids:
                logging.info(f"Nenhuma sequência para '{term}' encontrada.")
                continue

            logging.info(
                f"Encontradas {len(ids)} sequências para '{term}'. Iniciando download..."
            )

            # Usado await para buscar sequências de forma assíncrona
            sequences = await fetch_all_sequences_async(ids, webenv, query_key)

            # Filtrar sequências para garantir que sejam do gene DRD4 e tenham tamanho suficiente
            filtered_sequences = [
                seq
                for seq in sequences
                if "DRD4" in seq.description and len(seq.seq) > MIN_SEQ_LENGTH
            ]

            logging.info(
                f"Obtidas {len(filtered_sequences)} sequências válidas para '{term}'"
            )

            # Save results to cache
            if USE_CACHE:
                try:
                    with open(cache_path, "wb") as f:
                        pickle.dump(
                            {"sequences": filtered_sequences, "time": time.time()}, f
                        )
                        logging.info(f"Cache salvo para '{term}'")
                except Exception as e:
                    logging.warning(f"Erro ao salvar cache para '{term}': {e}")

            all_sequences.extend(filtered_sequences)

        except Exception as e:
            logging.error(f"Erro ao buscar sequências para '{term}': {str(e)}")

    # Remover duplicatas
    unique_sequences = []
    unique_ids = set()

    for seq in all_sequences:
        if seq.id not in unique_ids:
            unique_ids.add(seq.id)
            unique_sequences.append(seq)

    logging.info(f"Total de {len(unique_sequences)} sequências únicas para {condition}")

    return unique_sequences


async def fetch_all_sequences_async(ids, webenv, query_key, batch_size=5):
    """Busca todas as sequências do NCBI de forma assíncrona."""
    tasks = []
    for start in range(0, len(ids), batch_size):
        end = min(start + batch_size, len(ids))
        tasks.append(
            fetch_sequence_async(ids[start:end], webenv, query_key, start, end)
        )

    results = []
    for task_batch in asyncio.as_completed(tasks):
        batch_results = await task_batch
        results.extend(batch_results)

    return results


async def get_all_conditions_sequences_async(max_sequences=200, force_recalc=False):
    """Função assíncrona para buscar sequências de ADHD e autismo."""
    logging.info("Buscando todas as sequências relacionadas a ADHD e autismo...")

    # Buscar sequências de ADHD
    adhd_unique_sequences = await fetch_condition_sequences_async(
        "ADHD", max_sequences, force_recalc
    )

    # Buscar sequências de autismo
    autism_unique_sequences = await fetch_condition_sequences_async(
        "Autismo", max_sequences, force_recalc
    )

    logging.info(
        f"Total de sequências únicas relacionadas ao ADHD: {len(adhd_unique_sequences)}"
    )
    logging.info(
        f"Total de sequências únicas relacionadas ao autismo: {len(autism_unique_sequences)}"
    )

    # Salvar as sequências
    adhd_filtered_file = os.path.join(SEQUENCES_DIR, "drd4_adhd_sequences.fasta")
    autism_filtered_file = os.path.join(SEQUENCES_DIR, "drd4_autism_sequences.fasta")

    SeqIO.write(adhd_unique_sequences, adhd_filtered_file, "fasta")
    SeqIO.write(autism_unique_sequences, autism_filtered_file, "fasta")

    logging.info(
        f"Sequências relacionadas ao ADHD salvas no arquivo: {adhd_filtered_file}"
    )
    logging.info(
        f"Sequências relacionadas ao autismo salvas no arquivo: {autism_filtered_file}"
    )

    return adhd_filtered_file, autism_filtered_file


async def get_all_drd4_sequences_async(max_sequences=500, force_recalc=False):
    """Função assíncrona para buscar todas as sequências do gene DRD4."""
    cache_path = os.path.join(CACHE_DIR, "all_drd4_sequences.pickle")

    # Verificar cache
    if USE_CACHE and os.path.exists(cache_path) and not force_recalc:
        try:
            logging.info("Verificando cache para todas as sequências DRD4...")
            with open(cache_path, "rb") as f:
                cache_data = pickle.load(f)
                cache_time = cache_data.get("time", 0)

                if time.time() - cache_time < CACHE_TTL:
                    logging.info("Usando sequências DRD4 em cache")
                    sequences = cache_data.get("sequences", [])
                    fasta_file = os.path.join(SEQUENCES_DIR, "all_drd4_sequences.fasta")
                    SeqIO.write(sequences, fasta_file, "fasta")
                    return fasta_file
                else:
                    logging.info(
                        "Cache de sequências DRD4 expirado, baixando novamente..."
                    )
        except Exception as e:
            logging.warning(f"Erro ao ler cache: {e}")

    logging.info("Buscando todas as sequências relacionadas ao gene DRD4 no NCBI...")
    search_term = "DRD4[Gene] AND Homo sapiens[Organism]"

    try:
        handle = Entrez.esearch(
            db="nucleotide", term=search_term, retmax=max_sequences, usehistory="y"
        )
        record = Entrez.read(handle)
        handle.close()

        webenv = record["WebEnv"]
        query_key = record["QueryKey"]
        ids = record["IdList"]
        logging.info(
            f"Encontradas {len(ids)} sequências. Iniciando download em lotes..."
        )

        sequences = await fetch_all_sequences_async(
            ids, webenv, query_key, batch_size=10
        )

        # Filtrar sequências para garantir que sejam do gene DRD4
        filtered_sequences = [
            seq
            for seq in sequences
            if "DRD4" in seq.description and len(seq.seq) > MIN_SEQ_LENGTH
        ]

        # Mostrar resumo
        logging.info(
            f"Processadas {len(filtered_sequences)} sequências válidas do DRD4"
        )

        # Salvar no cache
        if USE_CACHE:
            try:
                with open(cache_path, "wb") as f:
                    pickle.dump(
                        {"sequences": filtered_sequences, "time": time.time()}, f
                    )
                    logging.info("Cache de todas as sequências DRD4 salvo")
            except Exception as e:
                logging.warning(f"Erro ao salvar cache: {e}")

        fasta_file = os.path.join(SEQUENCES_DIR, "all_drd4_sequences.fasta")
        SeqIO.write(filtered_sequences, fasta_file, "fasta")
        logging.info(
            f"Todas as sequências relacionadas ao DRD4 salvas no arquivo: {fasta_file}"
        )
        return fasta_file

    except Exception as e:
        logging.error(f"Erro ao buscar sequências DRD4: {str(e)}")
        # Retornar um arquivo vazio em caso de erro
        empty_file = os.path.join(SEQUENCES_DIR, "all_drd4_sequences.fasta")
        with open(empty_file, "w") as f:
            pass
        return empty_file


async def run_all_async_tasks(
    skip_download=False, force_recalc=False, max_sequences=200
):
    """Executa todas as tarefas assíncronas em uma única função."""
    # Verificar se arquivos já existem quando skip_download=True
    adhd_filtered_file = os.path.join(SEQUENCES_DIR, "drd4_adhd_sequences.fasta")
    autism_filtered_file = os.path.join(SEQUENCES_DIR, "drd4_autism_sequences.fasta")
    all_sequences_file = os.path.join(SEQUENCES_DIR, "all_drd4_sequences.fasta")

    files_exist = (
        os.path.exists(adhd_filtered_file)
        and os.path.exists(autism_filtered_file)
        and os.path.exists(all_sequences_file)
    )

    if skip_download and files_exist:
        logging.info("Pulando download de sequências, usando arquivos existentes")
        return adhd_filtered_file, autism_filtered_file, all_sequences_file

    # Obter sequências diretamente relacionadas ao ADHD e autismo
    adhd_filtered_file, autism_filtered_file = await get_all_conditions_sequences_async(
        max_sequences=max_sequences, force_recalc=force_recalc
    )

    # Buscar todas as sequências do gene DRD4 para referência
    all_sequences_file = await get_all_drd4_sequences_async(
        max_sequences=max_sequences * 2,  # Mais sequências para referência geral
        force_recalc=force_recalc,
    )

    return adhd_filtered_file, autism_filtered_file, all_sequences_file
