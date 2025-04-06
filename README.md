# Análise de Polimorfismos do Gene DRD4

## Visão Geral
Este projeto analisa polimorfismos no gene do receptor de dopamina D4 (DRD4), com foco em sequências associadas ao TDAH (ADHD) e transtorno do espectro autista. O software realiza alinhamento de sequências, detecta polimorfismos (SNPs e VNTRs) e gera relatórios comparativos entre as condições.

## Funcionalidades
- Recuperação automática de sequências do gene DRD4 do NCBI
- Alinhamento de sequências usando ClustalW
- Detecção de Polimorfismos de Nucleotídeo Único (SNPs)
- Identificação de Repetições em Tandem de Número Variável (VNTRs)
- Análise comparativa entre sequências associadas ao TDAH e autismo
- Tratamento especial de variantes de sequência divergentes

## Instalação

### Pré-requisitos
- Python 3.8+
- BioPython
- ClustalW (deve estar instalado e acessível no PATH)
- Matplotlib (opcional, para visualizações)

### Configuração
1. Clone este repositório
```bash
git clone https://github.com/seuusuario/analise-polimorfismos-drd4.git
cd analise-polimorfismos-drd4
```

2. Crie e ative um ambiente virtual (opcional, mas recomendado)
```bash
python -m venv venv
source venv/bin/activate  # No Windows: venv\Scripts\activate
```

3. Instale as dependências
```bash
pip install -r requirements.txt
```

## Uso

### Uso Básico
Execute o script principal para realizar uma análise completa:
```bash
python main.py
```

### Argumentos de Linha de Comando
- `--skip-download`: Pula o download de sequências se elas já existirem localmente
- `--verbose`: Habilita logs detalhados
- `--force-recalc`: Força o recálculo mesmo que existam resultados em cache
- `--max-seq N`: Limita o número máximo de sequências a serem baixadas por condição (padrão: 250)

Exemplo:
```bash
python main.py --verbose --max-seq 100
```

## Estrutura do Projeto

- `main.py`: Ponto de entrada da aplicação
- `config.py`: Configurações e parâmetros
- `ncbi_utils.py`: Funções para recuperar sequências do NCBI
- `alignment_utils.py`: Funções de alinhamento de sequências usando ClustalW
- `polymorphism_analysis.py`: Detecção de SNPs e VNTRs
- `report_utils.py`: Geração de relatórios de análise
- `comparative_analysis.py`: Análise comparativa entre diferentes condições
- `async_utils.py`: Utilitários assíncronos para processamento paralelo

### Diretórios de Dados
- `drd4_analysis/`: Diretório raiz para arquivos de análise
  - `sequences/`: Sequências baixadas e filtradas
  - `alignments/`: Saída de alinhamentos de sequências
  - `reports/`: Relatórios gerados
  - `visualizations/`: Visualizações estatísticas (se habilitadas)

## Como Funciona

1. **Recuperação de Sequências**: Recupera sequências DRD4 do NCBI, filtrando aquelas associadas ao TDAH e autismo
2. **Tratamento de Variantes**: Identifica e separa variantes de sequência divergentes (série LC812)
3. **Alinhamento**: Alinha sequências usando ClustalW com parâmetros personalizados por condição
4. **Detecção de Polimorfismos**:
   - SNPs: Identifica posições com variações de nucleotídeos
   - VNTRs: Detecta repetições em tandem de número variável
5. **Geração de Relatórios**: Cria relatórios detalhados dos achados para cada condição
6. **Análise Comparativa**: Se várias condições forem analisadas, compara padrões entre elas

## Contexto Científico

O gene DRD4 codifica o receptor de dopamina D4, que tem sido implicado em várias condições neuropsiquiátricas. De particular interesse é o VNTR de 48 pb no éxon 3, com o alelo de 7 repetições associado ao TDAH em múltiplos estudos. Esta ferramenta ajuda pesquisadores a analisar esses polimorfismos em dados de sequência.

## Licença

Este projeto é licenciado sob a GNU General Public License v3.0 (GPL-3.0)

Permissões:
- Uso comercial
- Modificação
- Distribuição
- Uso em patentes
- Uso privado

Condições:
- Divulgação da fonte (source code)
- Manutenção da mesma licença
- Declaração de mudanças
- Disponibilização do código-fonte

Para mais informações, consulte o arquivo LICENSE ou visite:
[https://www.gnu.org/licenses/gpl-3.0.html](https://www.gnu.org/licenses/gpl-3.0.html)

## Agradecimentos

Este projeto foi desenvolvido como parte de um programa PIVIC (Programa Institucional Voluntário de Iniciação Científica).
