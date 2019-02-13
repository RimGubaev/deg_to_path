# Устанавливаем рабочую директорию
# Set the working directory
setwd("~/your/working/directory/")

# Загружаем библиотеку "plyr"
# Load the "plyr" library
library("plyr")

############################################
# 1 Выявляем пути A.bispodus, достоверно насыщенные активируемыми/репрессируемыми генами L. edodes
# 1 Find A.bispodus pathways significantly enriched by activated/repressedgenes of L. edodes
############################################

# Загружаем таблицу c транскриптами L.edodes (Led), распределенными по функциональным категориям KEGG для A.bishorus (Abi)
# Upload table with L.edodes (Led) transcripts, combined into KEGG functional categories of A.bishorus (Abi)
Led_path <- read.csv(file = "Led_Abi_pathways.csv", header = T)
path_names <- read.csv(file = "abi_pathway_name", sep = "\t", header = F)
path_names <- data.frame(path_num=path_names$V1, path_name=path_names$V2)

# Распределяем дифференциально экспрессирующиеся гены по метаболическим путям
# Classify DEGs according to metabolic pathways 
diff_exp <- read.csv(file = "gene_exp.diff", sep = "\t")
diff_exp <- diff_exp[ ,c(4,11,13,14,15)]
diff_exp_path <- merge(x = diff_exp, y = Led_path, by.x = "locus", by.y = "assembled_transcript")
diff_exp_path <- subset(diff_exp_path, diff_exp_path$p_value < 0.05)

# Из таблицы с аннотированными ДЭГ извлекаем достоверно активируемые
# Extracting significantly activated DEGs from the table
diff_exp_up <- subset(diff_exp_path, diff_exp_path$log2.fold_change. > 0)

# Из таблицы с аннотированными ДЭГ извлекаем достоверно репрессируемые
# Extracting significantly repressed DEGs from the table
diff_exp_down <- subset(diff_exp_path, diff_exp_path$log2.fold_change. < 0)

# Подсчитываем общее количество генов для каждого пути
# Count the total ammount of genes for each pathway
pathway_table <- as.data.frame(count(Led_path$path_num))

# Подсчитываем количество активируемых генов для каждого пути
# Count the total ammount of activated genes for each pathway
pathway_table_up <- as.data.frame(count(diff_exp_up$path_num))
colnames(pathway_table_up)[2] <- "up_in_pathway"


# Подсчитываем количество репрессируемых генов для каждого пути
# Count the total ammount of repressed genes for each pathway
pathway_table_down <- as.data.frame(count(diff_exp_down$path_num))
colnames(pathway_table_down)[2] <- "down_in_pathway"

# Объединяем посчитанные количество активируемых/репрессируемых генов для всех путей
# Merge counted amount of activated/repressed genes for all pathways
UP_DOWN_pathways <- Reduce(function(x, y) merge(x, y, all=TRUE), list(pathway_table, pathway_table_up, pathway_table_down))
UP_DOWN_pathways[is.na(UP_DOWN_pathways )] <- 0 
colnames(UP_DOWN_pathways)[1] <- "path_num"
colnames(UP_DOWN_pathways)[2] <- "genes_in_path"


# Подсчет количества активируемых и репрессируемых генов, всех генов, вошедших в пути
# Count the aounts of activated and repressed genes, all genes within pathways  
up_DEG <- nrow(diff_exp_up)
down_DEG <- nrow(diff_exp_down)
all_Genes_in_path <- nrow(Led_path)

#Тест Фишера для активируемых генов
# Fisher test for activated genes
UP_pathways <- UP_DOWN_pathways[ ,c(-4)]
UP_pathways$up_not_in_pathway <- up_DEG - UP_pathways$up_in_pathway
UP_pathways$no_up_in_pathway <- UP_pathways$genes_in_path - UP_pathways$up_in_pathway
UP_pathways$no_up_not_in_pathway <- ((all_Genes_in_path - UP_pathways$up_not_in_pathway) - UP_pathways$genes_in_path)
n <- nrow(UP_pathways)
for(i in 1:n) {
  print(i)
  dat <- c(UP_pathways$up_in_pathway[i], UP_pathways$up_not_in_pathway[i], UP_pathways$no_up_in_pathway[i], UP_pathways$no_up_not_in_pathway[i])
  dfxm <- matrix(dat, nrow = 2)
  UP_pathways$p.value[i] <- fisher.test(dfxm, alternative = "greater")$p.value
  # Подсчет соотношения шансов для каждого пути
  # Count the odds ration for each pathway
  UP_pathways$odds.ratio[i] <- fisher.test(dfxm, alternative = "greater")$estimate
  print (UP_pathways$p.value[i])
}

# Делаем поправку на множественное тестирование (False Discovery Rate)
UP_pathways$q.value <- (p.adjust(UP_pathways$p.value, method = "fdr"))
# Присваеваем путям названия
# Adding a names to the pathway                  
UP_pathways <- merge(x = UP_pathways[,-c(4:6)], y = path_names, by = "path_num")

# Тест Фишера для репрессируемых генов
# Fisher test for activated genes
DOWN_pathways <- UP_DOWN_pathways[ ,-3]
DOWN_pathways$down_not_in_pathway <- down_DEG - DOWN_pathways$down_in_pathway
DOWN_pathways$no_down_in_pathway <- DOWN_pathways$genes_in_path - DOWN_pathways$down_in_pathway
DOWN_pathways$no_down_not_in_pathway <- ((all_Genes_in_path - DOWN_pathways$down_not_in_pathway) - DOWN_pathways$genes_in_path)

n <- nrow(DOWN_pathways)
for(i in 1:n) {
  print(i)
  dat <- c(DOWN_pathways$down_in_pathway[i], DOWN_pathways$down_not_in_pathway[i], DOWN_pathways$no_down_in_pathway[i], DOWN_pathways$no_down_not_in_pathway[i])
  dfxm <- matrix(dat, nrow = 2)
  DOWN_pathways$p.value[i] <- fisher.test(dfxm, alternative = "greater")$p.value
  # Подсчет соотношения шансов для каждого пути
  # Count the odds ration for each pathway
  DOWN_pathways$odds.ratio[i] <- fisher.test(dfxm, alternative = "greater")$estimate
  print (DOWN_pathways$p.value[i])
}

# Делаем поправку на множественное тестирование (False Discovery Rate)
# Make an FDR correction
DOWN_pathways$q.value <- (p.adjust(DOWN_pathways$p.value, method = "fdr"))
# Присваеваем путям названия
# Adding a names to the pathways
DOWN_pathways <- merge(x = DOWN_pathways[,-c(4:6)], y = path_names, by = "path_num")

# Записываем полученные результаты в таблицу
# Save the results into a table
write.table(x = UP_pathways, sep = "\t", file = "LED_UP_pathways.csv", row.names = F, quote=FALSE)
write.table(x = DOWN_pathways, sep = "\t", file = "LED_DOWN_pathways.csv", row.names = F, quote=FALSE)

# очищаем среду
# clean the environment 
rm(list = ls())

############################################
# 2 Выявляем пути KEGG, достоверно насыщенные активируемыми/репрессируемыми генами L. edodes
# 2 Find KEGG pathwys significantly enriched by activated/repressed genes of L. edodes
############################################

# Загружаем таблицу c транскриптами L.edodes (Led), распределенными по функциональным категориям KEGG
# Upload table with L.edodes (Led) transcripts, classified into KEGG functional categories                   
Led_KEGG_path <- read.csv(file = "Led_KEGG_pathways.csv", header = T)
path_names <- read.csv(file = "KEGG_pathway_name", sep = "\t", header = F)
path_names <- data.frame(path_num=path_names$V1, path_name=path_names$V2)

# Распределяем дифференциально экспрессирующиеся гены по метаболическим путям
# Classify DEGs into matabolic pathways
diff_exp <- read.csv(file = "gene_exp.diff", sep = "\t")
diff_exp <- diff_exp[ ,c(4,11,13,14,15)]
diff_exp_path <- merge(x = diff_exp, y = Led_KEGG_path, by.x = "locus")
diff_exp_path <- subset(diff_exp_path, diff_exp_path$p_value < 0.05)

# Из таблицы с аннотированными ДЭГ извлекаем достоверно активируемые
# Extract the significantly activated genes from table with annotated genes
diff_exp_up <- subset(diff_exp_path, diff_exp_path$log2.fold_change. > 0)

# Из таблицы с аннотированными ДЭГ извлекаем достоверно репрессируемые
# Extract the significantly repressed genes from table with annotated genes
diff_exp_down <- subset(diff_exp_path, diff_exp_path$log2.fold_change. < 0)

# Подсчитываем общее количество генов для каждого пути
# Count the total ammount of genes for each pathway                   
pathway_table <- as.data.frame(count(Led_KEGG_path$map))

# Подсчитываем количество активируемых генов для каждого пути
# Count the total ammount of activated genes for each pathway
pathway_table_up <- as.data.frame(count(diff_exp_up$map))
colnames(pathway_table_up)[2] <- "up_in_pathway"

# Подсчитываем количество репрессируемых генов для каждого пути
# Count the total ammount of repressed genes for each pathway
pathway_table_down <- as.data.frame(count(diff_exp_down$map))
colnames(pathway_table_down)[2] <- "down_in_pathway"

# Объединяем посчитанные количество активируемых/репрессируемых генов для всех путей
# Merge counted amount of activated/repressed genes for all pathways
UP_DOWN_pathways <- Reduce(function(x, y) merge(x, y, all=TRUE), list(pathway_table, pathway_table_up, pathway_table_down))
UP_DOWN_pathways[is.na(UP_DOWN_pathways )] <- 0 
colnames(UP_DOWN_pathways)[1] <- "path_num"
colnames(UP_DOWN_pathways)[2] <- "genes_in_path"

# подсчет количества активируемых и репрессируемых генов, всех генов, вошедших в пути
# Count the aounts of activated and repressed genes, all genes within pathways  
up_DEG <- nrow(diff_exp_up)
down_DEG <- nrow(diff_exp_down)
all_Genes_in_path <- nrow(Led_KEGG_path)

#Тест Фишера для активируемых генов
# Fisher test for activated genes
UP_pathways <- UP_DOWN_pathways[ ,c(-4)]
UP_pathways$up_not_in_pathway <- up_DEG - UP_pathways$up_in_pathway
UP_pathways$no_up_in_pathway <- UP_pathways$genes_in_path - UP_pathways$up_in_pathway
UP_pathways$no_up_not_in_pathway <- ((all_Genes_in_path - UP_pathways$up_not_in_pathway) - UP_pathways$genes_in_path)
n <- nrow(UP_pathways)
for(i in 1:n) {
  print(i)
  dat <- c(UP_pathways$up_in_pathway[i], UP_pathways$up_not_in_pathway[i], UP_pathways$no_up_in_pathway[i], UP_pathways$no_up_not_in_pathway[i])
  dfxm <- matrix(dat, nrow = 2)
  UP_pathways$p.value[i] <- fisher.test(dfxm, alternative = "greater")$p.value  
  # Подсчет соотношения шансов для каждого пути
  # Count the odds ration for each pathway
  UP_pathways$odds.ratio[i] <- fisher.test(dfxm, alternative = "greater")$estimate
  print (UP_pathways$p.value[i])
}
                                               
# Делаем поправку на множественное тестирование (False Discovery Rate)
# Make an FDR correction
UP_pathways$q.value <- (p.adjust(UP_pathways$p.value, method = "fdr"))
# Присваеваем путям названия
# Adding a names to the pathways
UP_pathways <- merge(x = UP_pathways[,-c(4:6)], y = path_names, by = "path_num")

# Тест Фишера для репрессируемых генов 
# Fisher test for repressed genes
DOWN_pathways <- UP_DOWN_pathways[ ,-3]
DOWN_pathways$down_not_in_pathway <- down_DEG - DOWN_pathways$down_in_pathway
DOWN_pathways$no_down_in_pathway <- DOWN_pathways$genes_in_path - DOWN_pathways$down_in_pathway
DOWN_pathways$no_down_not_in_pathway <- ((all_Genes_in_path - DOWN_pathways$down_not_in_pathway) - DOWN_pathways$genes_in_path)

n <- nrow(DOWN_pathways)
for(i in 1:n) {
  print(i)
  dat <- c(DOWN_pathways$down_in_pathway[i], DOWN_pathways$down_not_in_pathway[i], DOWN_pathways$no_down_in_pathway[i], DOWN_pathways$no_down_not_in_pathway[i])
  dfxm <- matrix(dat, nrow = 2)
  DOWN_pathways$p.value[i] <- fisher.test(dfxm, alternative = "greater")$p.value
  # Подсчет соотношения шансов для каждого пути
  # Count the odds ration for each pathway
  DOWN_pathways$odds.ratio[i] <- fisher.test(dfxm, alternative = "greater")$estimate
  print (DOWN_pathways$p.value[i])
}

# Делаем поправку на множественное тестирование (False Discovery Rate)
# Make an FDR correction                           
DOWN_pathways$q.value <- (p.adjust(DOWN_pathways$p.value, method = "fdr"))

# Присваеваем путям названия
# Adding a names to the pathways
DOWN_pathways <- merge(x = DOWN_pathways[,-c(4:6)], y = path_names, by = "path_num")

# Записываем полученные результаты в таблицу
# Save the results into a table  
write.table(x = UP_pathways, sep = "\t", file = "KEGG_UP_pathways.csv", row.names = F, quote=FALSE)
write.table(x = DOWN_pathways, sep = "\t", file = "KEGG_DOWN_pathways.csv", row.names = F, quote=FALSE)

############################################
# 3 Формируем таблицу входных данных для BioCyc
# 3 Create the table of input data for BioCyc
############################################

# Создание таблиц, необходимых для визуализации ДЭГ в BioCyc
# Create the tables for DEG visualization in BioCyc
Pru_transcript_to_loci <- read.csv(file = "Pru_transcript_to_loci.csv")
Pru_transcript_to_loci <- merge(x = Pru_transcript_to_loci, y = diff_exp, by.x = "assembled_transcript", by.y = "locus")
Pru_transcript_to_loci_up <- subset(Pru_transcript_to_loci, Pru_transcript_to_loci$significant == "yes" & Pru_transcript_to_loci$log2.fold_change. > 0)
Pru_transcript_to_loci_down <- subset(Pru_transcript_to_loci, Pru_transcript_to_loci$significant == "yes" & Pru_transcript_to_loci$log2.fold_change. < 0)
write.table(x = Pru_transcript_to_loci_up[,c(2,3)], sep = "\t", file = "Led_to_biocyc_up.csv", row.names = F, quote=FALSE)
write.table(x = Pru_transcript_to_loci_down[,c(2,3)], sep = "\t", file = "Led_to_biocyc_down.csv", row.names = F, quote=FALSE)
