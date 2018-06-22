# deg_to_path
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
