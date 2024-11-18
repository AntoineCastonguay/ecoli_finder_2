ecoli <- read.table("analyse/output.txt", header = T)

ecoli_BW25113_2 <- read.table("data/All-genes-of-E.-coli-K-12-substr.-BW25113.txt", header = T)
ecoli_BW25113_2 <- ecoli_BW25113_2[order(ecoli_BW25113_2$Left.End.Position),1:4]
rownames(ecoli_BW25113_2) <- 1:nrow(ecoli_BW25113_2)

ecoli_BW25113 <- read.table("data/BW25113_features_supplemented.bed")
ecoli_BW25113 <- ecoli_BW25113[,c(4,2,3,6)]
colnames(ecoli_BW25113) <- c("Gene_Name","Left.End.Position", "Right.End.Position", "sens")

# ajoute longeur à la base de donné de référence
gene_length <- c()
for (i in 1:nrow(ecoli_BW25113)) {
  res = ecoli_BW25113$Right.End.Position[i] - ecoli_BW25113$Left.End.Position[i]
  gene_length <- append(gene_length, res)
}
ecoli_BW25113 <- cbind(ecoli_BW25113, length = gene_length)

# Trouve la position du debut de gene selon la position du primer donc position initialer de primer left +47 car 3 nucléotide de codon d'intiation dans le primer.
# Pour le primer right, position initial du primer + 20 car 20 nucléotide du gene dans le primer
#ecoli$first_pos <- ifelse(ecoli$flag == 99 | ecoli$flag == 97 | ecoli$flag == 163 | ecoli$flag == 161,
#                          ecoli$first_pos + 47,
#                          ifelse(ecoli$flag == 147 | ecoli$flag == 145 | ecoli$flag == 83 | ecoli$flag == 81,
#                                 ecoli$first_pos + 20, ecoli$first_pos))
#ecoli$second_pos <- ifelse(ecoli$flag == 99 | ecoli$flag == 97 | ecoli$flag == 163 | ecoli$flag == 161,
#                           ecoli$second_pos + 20,
#                           ifelse(ecoli$flag == 147 | ecoli$flag == 145 | ecoli$flag == 83 | ecoli$flag == 81,
#                                  ecoli$second_pos + 47, ecoli$second_pos))

for (i in 1:length(ecoli$first_pos)) {
  if (ecoli$init[i] == 'f') {
    if (ecoli$flag[i] == 99 | ecoli$flag[i] == 97 | ecoli$flag[i] == 81 | ecoli$flag[i] == 83) {
      ecoli$first_pos[i] <- ecoli$first_pos[i] + 47
      ecoli$second_pos[i] <- ecoli$second_pos[i] + 20
    }else if(ecoli$flag[i] == 147 | ecoli$flag[i] == 145 | ecoli$flag[i] == 163 | ecoli$flag[i] == 161){
      ecoli$first_pos[i] <- ecoli$first_pos[i] + 29
      ecoli$second_pos[i] <- ecoli$second_pos[i] + 2
    }
  }else{
    if (ecoli$flag[i] == 99 | ecoli$flag[i] == 97 | ecoli$flag[i] == 83 | ecoli$flag[i] == 81) {
      ecoli$first_pos[i] <- ecoli$first_pos[i] + 29
      ecoli$second_pos[i] <- ecoli$second_pos[i] +2
    }else if(ecoli$flag[i] == 147 | ecoli$flag[i] == 145 | ecoli$flag[i] == 163 | ecoli$flag[i] == 161){
      ecoli$first_pos[i] <- ecoli$first_pos[i] + 47
      ecoli$second_pos[i] <- ecoli$second_pos[i] + 20
    }
  }
}

ecoli2 <- read.table("analyse/output.txt", header = T)
for (i in 1:length(ecoli2$first_pos)) {
    if (ecoli2$flag[i] == 99 | ecoli2$flag[i] == 97 | ecoli2$flag[i] == 83 | ecoli2$flag[i] == 81) {
      ecoli2$first_pos[i] <- ecoli2$first_pos[i] + 50
      ecoli2$second_pos[i] <- ecoli2$second_pos[i]
    }else if(ecoli2$flag[i] == 147 | ecoli2$flag[i] == 145 | ecoli2$flag[i] == 163 | ecoli2$flag[i] == 161){
      ecoli2$first_pos[i] <- ecoli2$first_pos[i] 
      ecoli2$second_pos[i] <- ecoli2$second_pos[i] + 50
    }
    ecoli2$length[i] <- ecoli2$second_pos[i] - ecoli2$first_pos[i]
}
ecoli_positif2 <- subset(ecoli2,ecoli2$length > 0 & ecoli2$length < 10000)
write.csv(ecoli_positif2, file = "../keio_long_read/data/ecoli_bw25113_deletion.csv")


# 77 car 47 pour le primer left qui est a l'exterieur du gene et 30 pour le right. 
for (i in 1:nrow(ecoli)) {
  if (ecoli$length[i] > 0) {
    ecoli$length[i] = ecoli$length[i] - 77
  }else if(ecoli$length[i] < 0){
    ecoli$length[i] = ecoli$length[i] + 77
  }
}

huge_gene <- subset(ecoli,ecoli$length > 100000)
# Il y a 10 gene dont la longueur est aberente dont 3 presume essentiel.

other_flag <- subset(ecoli, ecoli$flag %in% c(65, 129, 81, 161, 97, 145, 113, 177))
good_flag <- subset(ecoli, ecoli$flag %in% c(83,99,147,163))
table(ecoli$flag)

# Les flag
# Mapped within the insert size and in correct orientation
# 99 = read paired, read mapped in proper pair, *mate* reverse strand et *first* in pair.
# 147 = read paired, read mapped in proper pair, *read* reverse strand et *second* in pair.
# 83 = read paired, read mapped in proper pair, *read* reverse strand et *first* in pair.
# 163 = read paired, read mapped in proper pair, *mate* reverse strand et *second* in pair.
# Mapped uniquely, but with wrong insert size
# 65 = read paired et first in pair
# 129 = read paired et second in pair
# 81 = read paired, read reverse strand et first in pair
# 161 = read paired, mate reverse strand et second in pair.
# 97 = read paired, mate reverse strand et first in pair.
# 145 = read paired, read reverse strand et second in pair.
# 113 = read paired, read reverse strand, mate reverse strand et first in pair
# 177 = read paired, read reverse strand, mate reverse strand et second in pair

# 113 et 177 problematique puisque aligner sur le meme brin.
# 65, 129, 81, 161, 97 et 145 insert > 3550 inclus les pairs aberentes.

quality <- subset(ecoli,ecoli$quality != '50M')
table(ecoli$quality)
# 26 donc la qualite n'est pas 50M sur 8576.
# Il y a 25 non-essentiel et 1 presume essentiel.

ecoli_positif <- subset(ecoli,ecoli$length > 0 & ecoli$length < 10000)
rownames(ecoli_positif) <- 1:nrow(ecoli_positif)
ecoli_new_gene <- ecoli_positif[,c(1,2,4,5)]
# garde une copie par gene et enleve les longeur aberente
new_gene_f <- c()
new_gene_s <- c()
new_gene_f_pos <- data.frame()
new_gene_s_pos <- data.frame()
diff_f <- c()
diff_s <- c()
for (i in 1:nrow(ecoli_new_gene)) {
  valeur_first <- ecoli_new_gene$first_pos[i]
  valeur_second <- ecoli_new_gene$second_pos[i]
  # Calculer la correspondance la plus proche
  index_proche_first <- which.min(abs(ecoli_BW25113$Left.End.Position - valeur_first))
  index_proche_second <- which.min(abs(ecoli_BW25113$Right.End.Position - valeur_second))
  
  new_gene_f[i] <- ecoli_BW25113$Gene_Name[index_proche_first]
  new_gene_s[i] <- ecoli_BW25113$Gene_Name[index_proche_second]
  new_gene_f_pos[i,1] <- ecoli_BW25113$Left.End.Position[index_proche_first]
  new_gene_f_pos[i,2] <- ecoli_BW25113$Right.End.Position[index_proche_first]
  new_gene_s_pos[i,1] <- ecoli_BW25113$Right.End.Position[index_proche_second]
  new_gene_s_pos[i,2] <- ecoli_BW25113$Left.End.Position[index_proche_second]
  diff_f[i] <- ecoli_new_gene$first_pos[i] - new_gene_f_pos[i,1]
  diff_s[i] <- ecoli_new_gene$second_pos[i] - new_gene_s_pos[i,1]
}
ecoli_new_gene <- cbind(ecoli_new_gene, 
                       new_gene_f = new_gene_f,
                       new_gene_s = new_gene_s,
                       new_gene_f_pos = new_gene_f_pos[,1],
                       new_gene_f_s_pos = new_gene_f_pos[,2],
                       new_gene_s_pos = new_gene_s_pos[,1],
                       new_gene_s_f_pos = new_gene_s_pos[,2],
                       diff_f = diff_f,
                       diff_s = diff_s,
                       fiable = (new_gene_f == new_gene_s))

ecoli_positif <- cbind(ecoli_positif,
                       new_gene = new_gene_f,
                       fiable = ecoli_new_gene$new_gene_f == ecoli_new_gene$new_gene_s,
                       modif = ecoli_new_gene$new_gene_f != ecoli_new_gene$gene | ecoli_new_gene$new_gene_s != ecoli_new_gene$gene)

ecoli_positif <- ecoli_positif[,c(1,2,9,10,11,12,3,4,5,6,7,8)]
ecoli_new_gene <- ecoli_new_gene[,c('id','gene','fiable','new_gene_f','new_gene_s','first_pos','new_gene_f_pos','new_gene_f_s_pos','diff_f','second_pos','new_gene_s_pos','new_gene_s_f_pos','diff_s')]

table(ecoli_new_gene$new_gene_f == ecoli_new_gene$new_gene_s)
sub_new_gene_f_s <- subset(ecoli_new_gene,ecoli_new_gene$new_gene_f != ecoli_new_gene$new_gene_s)
sub_new_gene_old_new <- subset(ecoli_new_gene, ecoli_new_gene$new_gene_f != ecoli_new_gene$gene | ecoli_new_gene$new_gene_s != ecoli_new_gene$gene)

write.csv(sub_new_gene_old_new, file = "../keio_long_read/data/ecoli_bw25113_gene_name_modif.csv")

correspondance_l <- c()
correspondance_r <- c()
c_l <- c()
c_r <- c()
c_diff_l <- c()
c_diff_r <- c()
for (i in 1:nrow(ecoli_positif)) {
  res <- ecoli_positif$first_pos[i] %in% ecoli_BW25113$Left.End.Position
  res2 <- ecoli_positif$second_pos[i] %in% ecoli_BW25113$Right.End.Position
  correspondance_l <- append(correspondance_l, res)
  correspondance_r <- append(correspondance_r, res2)
  
  if (ecoli_positif$new_gene[i] != 'none') {
    result <- subset(ecoli_BW25113, Gene_Name == ecoli_positif$new_gene[i])
    if (nrow(result) == 1){
      if(ecoli_positif$first_pos[i] == result$Left.End.Position){
        c_l <- append(c_l,'T')
        c_diff_l <- append(c_diff_l,NA)
      }else{
        c_l <- append(c_l,result$Left.End.Position)
        c_diff_l <- append(c_diff_l, ecoli_positif$first_pos[i] - result$Left.End.Position)
      }
      
      if(ecoli_positif$second_pos[i] == result$Right.End.Position){
        c_r <- append(c_r,'T')
        c_diff_r <- append(c_diff_r,NA)
      }else{
        c_r <- append(c_r,result$Right.End.Position)
        c_diff_r <- append(c_diff_r, ecoli_positif$second_pos[i] - result$Right.End.Position)
      }
    }else if(nrow(result) == 0){
      c_l <- append(c_l,NA)
      c_diff_l <- append(c_diff_l,NA)
      c_r <- append(c_r,NA)
      c_diff_r <- append(c_diff_r,NA)
    }else{
      c_l <- append(c_l,'>2')
      c_diff_l <- append(c_diff_l,NA)
      c_r <- append(c_r,'>2')
      c_diff_r <- append(c_diff_r,NA)
    }
  }else{
    c_l <- append(c_l,NA)
    c_diff_l <- append(c_diff_l,NA)
    c_r <- append(c_r,NA)
    c_diff_r <- append(c_diff_r,NA)
  }
  
}
ecoli_positif <- cbind(ecoli_positif, correspondance_l = correspondance_l, c_gene_l = c_l, c_diff_l = c_diff_l, correspondance_r = correspondance_r, c_gene_r = c_r,c_diff_r=c_diff_r)

table(ecoli_positif$correspondance_l)
table(ecoli_positif$correspondance_r)

ecoli_positif_first <- ecoli_positif[,c('id','gene','init','first_pos','correspondance_l','c_gene_l','c_diff_l','second_pos','correspondance_r','c_gene_r','c_diff_r')]

write.csv(ecoli_positif, file = "../keio_long_read/data/ecoli_bw25113_gene.csv")
