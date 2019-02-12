# deg_to_path.r
Рим Губаев, 2017

Скрипт deg_to_path.r предназначен для того, чтобы распределить дифференциально экспрессирующиеся гены (ДЭГ) *Lentinus edodes* по метаболическим путям близкородственного вида (*Agaricus bisphorus*), а также всем известным метаболическим путям KEGG, определить пути, достоверно представленые активируемыми/репрессируемыми генами. Кроме того, с помощью данного скрипта ДЭГ *L. edodes* можно присвоить идентификационные номера генов *Penicillium rubens*, что необходимо для проведения функционального анализа с помощью ресурса BioCyc.
Для того, чтобы произвести вышеописанные манипуляции с помощью данного скрипта для других видов нужно подготовить входные данные в соответствии с примерами, представленными в директории ./deg_to_path и описанными ниже.

Входные данные:

1) Таблицы Led_Abi_pathways.csv, Led_KEGG_pathways.csv, в которых собранным транскриптам *L. edodes* присвоены номера путей *A. bisporus*, а также всех известных путей базы данных KEGG. Таблицы были получены с помощью скрипта path_annotation.r.

2) Таблицы abi_pathway_name,  KEGG_pathway_name, содержащие названия и номера путей A. bisporus, а также всех путей базы данных KEGG. Получены c помощью KEGG API (http://rest.kegg.jp/list/pathway/abv, http://rest.kegg.jp/list/pathway).

3) Таблица gene_exp.diff, содержащая информацию о дифференциально экспрессирующихся генах *L. edodes*. Таблица была получена с помощью программы Cuffdiff.

4) Таблица Pru_transcript_to_loci.csv, в которой собранным транскриптам *L. edodes* присвоены идентификационные номера генов *P. rubens*. Таблица была получена с помощью скрипта path_annotation.r.

Выполняемые операции:

На первом этапе ДЭГ *L. edodes* присвиваются номера метаболических путей *A.bisporus*. Затем ДЭГ с присвоенными номерами метаболических путей *A.bisporus* делятся на две группы: 1) достоверно активируемые, 2) достоверено репрессируемые. Затем, в каждой группе подсчитывается количество ДЭГ находящихся в каждом пути. Подсчитывается достоверность представленности пути в каждой группе по сравнению с представленностью данного пути в общей совокопности ДЭГ с помощью теста Фишера. Для каждого пути подсчитывается соотношение шансов.
К ДЭГ *L. edodes* дополнительно присваюваются идентификационные номера *A. bisporus* для проведения функционального анализа с помощью ресурса BioCyc.

Все манипуляции, производимые с входными данными, описаны непосредственно в R скрипте, в коментариях к командам (строки, начинающиеся с символа "#").

Выходные данные:

1) Таблица LED_UP_pathways.csv, содержащая информацию о метаболических путях *A.bisporus*, насыщенных активируемыми генами (количество генов в пути, количество активируемых генов  в пути, уровень достоверности представленности данного пути, соотношение шансов).

2) Таблица LED_DOWN_pathways.csv, содержащая информацию о метаболических путях *A.bisporus*, насыщенных репрессируемымми генами (количество генов в пути, количество активируемых генов  в пути, уровень достоверности представленности данного пути, соотношение шансов).

3) Таблица KEGG_UP_pathways.csv, содержащая информацию о всех известных метаболических путях KEGG, насыщенных активируемыми генами (количество генов в пути, количество активируемых генов  в пути, уровень достоверности представленности данного пути, соотношение шансов).

4) Таблица KEGG_DOWN_pathways.csv, содержащая информацию о всех известных метаболических путях KEGG, насыщенных репрессируемымми генами (количество генов в пути, количество активируемых генов  в пути, уровень достоверности представленности данного пути, соотношение шансов).

5) Таблицы Led_to_biocyc_down.csv и Led_to_biocyc_down.csv, содержащие идентификационные номера активируемых и репрессируемых генов *P. rubens* с соответствующими уровнями экспрессии ДЭГ *L. edodes*.

Email: rimgubaev@gmail.com

Rim Gubaev, 2017

The deg_to_path.r script is designed to classify differentially expressed genes (DEG) of *Lentinus edodes* into the closely related species (*Agaricus bisphorus*) metabolic pathways, as well as all known KEGG metabolic pathways, to identify pathways that are significantly represented by activated/repressed genes. In addition, using this DEG script *L. edodes*, you can assign identification numbers of *Penicillium rubens* genes, which is necessary for carrying out a functional analysis using the BioCyc resource. In order to make the above manipulations using this script for other types, you need to prepare the input data in accordance with the examples presented in the ./deg_to_path directory and described below.

Input data:

1) Tables Led_Abi_pathways.csv, Led_KEGG_pathways.csv, in which the collected L. edodes transcripts are assigned the path numbers of A. bisporus, as well as all known KEGG database paths. Tables were generated using the path_annotation.r script.

2) The abi_pathway_name, KEGG_pathway_name tables contain the names and numbers of the A. bisporus paths, as well as all the paths of the KEGG database. Obtained using the KEGG API (http://rest.kegg.jp/list/pathway/abv, http://rest.kegg.jp/list/pathway).

3) The gene_exp.diff table, containing information about differentially expressed L. edodes genes. The table was obtained using the Cuffdiff program.

4) The Pru_transcript_to_loci.csv table, in which the identification numbers of P. rubens genes are assigned to the assembled L. edodes transcripts. The table was obtained using the path_annotation.r script.

Operations performed:

In the first stage of DEG L. edodes, the numbers of the metabolic pathways of A.bisporus are assigned. Then DEG with assigned numbers of the metabolic pathways of A.bisporus are divided into two groups: 1) reliably activated, 2) reliably repressed. Then, in each group, the number of DEG in each path is calculated. The accuracy of the representation of the path in each group is calculated compared to the representation of this path in the overall sosocopism of DEG using the Fisher test. For each path is calculated odds ratio. DEG L. edodes are additionally assigned identification numbers of A. bisporus for functional analysis using the BioCyc resource.

All manipulations made with the input data are described directly in the R script, in the comments to the commands (lines starting with the "#" character).

Output:

1) Table LED_UP_pathways.csv, containing information on the A.bisporus metabolic pathways saturated with activated genes (the number of genes in the path, the number of activated genes in the path, the level of confidence in the representation of this path, the odds ratio).

2) Table LED_DOWN_pathways.csv, containing information on the metabolic pathways of A.bisporus saturated with repressed genes (the number of genes in the path, the number of activated genes in the path, the level of confidence in the representation of this path, the odds ratio).

3) Table KEGG_UP_pathways.csv, containing information on all known KEGG metabolic pathways saturated with activated genes (the number of genes in the path, the number of activated genes in the path, the confidence level of the representation of this path, the odds ratio).

4) Table KEGG_DOWN_pathways.csv, containing information on all known KEGG metabolic pathways saturated with repressed genes (the number of genes in the path, the number of activated genes in the path, the confidence level of the representation of this path, the odds ratio).

5) Tables Led_to_biocyc_down.csv and Led_to_biocyc_down.csv, containing the identification numbers of the activated and repressed P. rubens genes with the corresponding expression levels of DEG L. edodes.

Email: rimgubaev@gmail.com

