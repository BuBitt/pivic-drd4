# Análise de Polimorfismos do Gene DRD4

## Visão Geral
Este projeto analisa polimorfismos no gene do receptor de dopamina D4 (DRD4), com foco em sequências associadas ao TDAH (ADHD) e transtorno do espectro autista. O software realiza alinhamento de sequências, detecta polimorfismos (SNPs e VNTRs) e gera relatórios comparativos entre as condições.

## Funcionalidades
- Recuperação automática de sequências do gene DRD4 do NCBI
- Alinhamento de sequências usando ClustalW ou MAFFT (configurável)
- Detecção de Polimorfismos de Nucleotídeo Único (SNPs)
- Identificação de Repetições em Tandem de Número Variável (VNTRs)
- Análise comparativa entre sequências associadas ao TDAH e autismo
- Tratamento especial de variantes LC (Library of Congress) do DRD4
- Comparação de resultados entre diferentes ferramentas de alinhamento (opcional)

## Instalação

### Pré-requisitos
- Python 3.8+
- BioPython
- ClustalW e/ou MAFFT (devem estar instalados e acessíveis no PATH)
- Matplotlib (opcional, para visualizações)
- NumPy

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

4. Instale as ferramentas de alinhamento
```bash
# Para MAFFT (recomendado)
sudo apt-get install mafft  # Em sistemas baseados em Debian/Ubuntu

# Para ClustalW (opcional)
# Baixe do site oficial e configure o caminho em config.py
```

## Uso

### Uso Básico
Execute o script principal para realizar uma análise completa:
```bash
# Usando MAFFT (recomendado)
python main.py

# Usando ClustalW
python main.py --alignment-tool clustalw

# Usando ambas as ferramentas para comparação
python main.py --alignment-tool both
```

### Argumentos de Linha de Comando
- `--skip-download`: Pula o download de sequências se elas já existirem localmente
- `--verbose`: Habilita logs detalhados
- `--force-recalc`: Força o recálculo mesmo que existam resultados em cache
- `--max-seq N`: Limita o número máximo de sequências a serem baixadas por condição (padrão: 250)
- `--alignment-tool {mafft,clustalw,both}`: Escolhe a ferramenta de alinhamento (padrão: mafft)
- `--preferred-tool {mafft,clustalw}`: Em caso de usar 'both', define qual ferramenta será usada para o relatório final

Exemplo:
```bash
python main.py --verbose --max-seq 100 --alignment-tool mafft
```

## Algoritmos de Alinhamento

### MAFFT (Recomendado)
O MAFFT (Multiple Alignment using Fast Fourier Transform) é uma ferramenta moderna para alinhamento de sequências com várias estratégias:
- **--auto**: Seleciona automaticamente o algoritmo baseado nas características das sequências (padrão para ADHD)
- **--genafpair**: G-INS-i, método iterativo global para sequências com regiões homólogas longas (usado para Autismo)
- **--localpair**: L-INS-i, método iterativo local para sequências com regiões variáveis (usado para variantes LC)

### ClustalW
Algoritmo tradicional de alinhamento progressivo, configurado com diferentes parâmetros para cada condição:
- **ADHD**: GAPOPEN=8, GAPEXT=0.15
- **Autismo**: GAPOPEN=6, GAPEXT=0.3
- **ADHD_Variantes**: GAPOPEN=5, GAPEXT=0.2, KIMURA=ON

## Estrutura do Projeto

- `main.py`: Ponto de entrada da aplicação
- `config.py`: Configurações e parâmetros
- `ncbi_utils.py`: Funções para recuperar sequências do NCBI
- `alignment_utils.py`: Funções de alinhamento usando ClustalW ou MAFFT
- `polymorphism_analysis.py`: Detecção de SNPs e VNTRs
- `report_utils.py`: Geração de relatórios de análise
- `comparative_analysis.py`: Análise comparativa entre diferentes condições
- `async_utils.py`: Utilitários assíncronos para processamento paralelo
- `validation_utils.py`: Comparação entre resultados de diferentes ferramentas

### Diretórios de Dados
- `drd4_analysis/`: Diretório base compartilhado
  - `sequences/`: Sequências baixadas (compartilhado entre ferramentas)
  - `reference/`: Sequências de referência (compartilhado entre ferramentas)
  - `cache/`: Cache de dados (compartilhado entre ferramentas)
- `drd4_analysis_mafft/`: Resultados da análise com MAFFT
  - `alignments/`: Saída de alinhamentos MAFFT
  - `reports/`: Relatórios gerados pelo MAFFT
  - `visualization/`: Visualizações estatísticas MAFFT
- `drd4_analysis_clustalw/`: Resultados da análise com ClustalW
  - `alignments/`: Saída de alinhamentos ClustalW
  - `reports/`: Relatórios gerados pelo ClustalW
  - `visualization/`: Visualizações estatísticas ClustalW

## Variantes LC (Library of Congress)

As sequências LC (especialmente LC812) representam variantes catalogadas do gene DRD4 que são oficialmente reconhecidas. O sistema identifica exclusivamente estas sequências como "ADHD_Variantes", priorizando a precisão científica sobre outras variações. Este tratamento específico permite:

1. Alinhamento de alta precisão com parâmetros especializados (--localpair no MAFFT)
2. Detecção de polimorfismos específicos às variantes oficialmente reconhecidas
3. Comparação estruturada entre os polimorfismos das sequências padrão e variantes
4. Documentação em relatórios Markdown que especifica claramente a fonte da variante

Este enfoque é particularmente valioso para estudos que buscam caracterizar as diferenças estruturais entre as variantes fundamentais do gene DRD4 e sequências típicas.

## Relatórios em Markdown

Todos os relatórios gerados pelo sistema estão em formato Markdown, proporcionando:

1. **Melhor formatação e leitura** - hierarquia clara com cabeçalhos, listas e tabelas
2. **Facilidade de conversão** - possibilidade de converter para HTML, PDF ou outros formatos
3. **Suporte a tabelas** - visualização tabular de SNPs, VNTRs e comparações
4. **Compatibilidade com GitHub/GitLab** - visualização direta em repositórios de código
5. **Portabilidade** - fácil incorporação em documentos científicos e publicações

Os relatórios incluem:
- Detecção de polimorfismos por condição
- Análises comparativas entre condições
- Comparação entre ferramentas de alinhamento (quando aplicável)

## Comparação de Ferramentas de Alinhamento

Quando o modo de comparação é utilizado (`--alignment-tool both`), o sistema:

1. Executa as análises com ambas as ferramentas (ClustalW e MAFFT)
2. Gera um relatório comparativo detalhado que inclui:
   - Contagem de SNPs e VNTRs detectados por cada ferramenta
   - Análise de sobreposição de polimorfismos 
   - Avaliação de divergências entre as ferramentas
   - Recomendações baseadas nos resultados

### Interpretando as Diferenças

As diferenças entre os resultados de ClustalW e MAFFT ocorrem principalmente devido a:

- **Estratégias de alinhamento**: MAFFT utiliza transformadas de Fourier e é geralmente mais preciso com sequências divergentes
- **Tratamento de gaps**: As diferentes penalidades para abertura e extensão de gaps podem resultar em diferentes detecções de regiões VNTR
- **Abordagem estatística**: Os algoritmos empregam diferentes métricas para calcular a similaridade entre sequências

**Recomendação**: Para sequências altamente variáveis como as do gene DRD4, os polimorfismos detectados consistentemente por ambas as ferramentas têm maior probabilidade de representar variações biológicas reais, enquanto discrepâncias podem exigir validação adicional.

## Como Funciona

1. **Recuperação de Sequências**: Recupera sequências DRD4 do NCBI, filtrando aquelas associadas ao TDAH e autismo
2. **Tratamento de Variantes**: Identifica sequências LC e as separa para análise específica
3. **Alinhamento**: Alinha sequências usando MAFFT ou ClustalW com parâmetros personalizados por condição
4. **Detecção de Polimorfismos**:
   - SNPs: Identifica posições com variações de nucleotídeos
   - VNTRs: Detecta repetições em tandem de número variável
5. **Geração de Relatórios**: Cria relatórios detalhados dos achados para cada condição
6. **Análise Comparativa**: Se várias condições forem analisadas, compara padrões entre elas

## Contexto Científico

O gene DRD4 codifica o receptor de dopamina D4, que tem sido implicado em várias condições neuropsiquiátricas. De particular interesse é o VNTR de 48 pb no éxon 3, com o alelo de 7 repetições associado ao TDAH em múltiplos estudos.

As variantes LC (especialmente LC812) representam sequências oficialmente catalogadas na Library of Congress e são particularmente valiosas para pesquisa por fornecerem referências estáveis para o estudo de variações nas populações. A detecção precisa destas variantes e seus polimorfismos pode contribuir para a compreensão da heterogeneidade genética subjacente ao TDAH.

Esta ferramenta ajuda pesquisadores a analisar sistematicamente esses polimorfismos em dados de sequência, permitindo análises comparativas robustas entre condições neuropsiquiátricas e variantes genéticas.

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
