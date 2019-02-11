#!/usr/bin/env/ Rscript

library(ROCR)
library(ggplot2)
library(plotrix)
library(gridExtra)

# set working directory
setwd("../")

## function added to use the CAGI answer key, it changes predictions row to template order according to ID column...
row.corrector <- function(template = temp, sub.lines = lines){
  sub.table <- strsplit(sub.lines[-1], "\t")
  pred.ID <- sapply(1:length(sub.table), function(x){sub.table[[x]][1]})
  correct.order <- sapply(1:length(sub.table), function(x){grep(paste("^", template[x], "$", sep=""), pred.ID)})
  sub.lines <- sub.lines[-1] # remove header
  lines.ordered <- sub.lines[correct.order]
  return(lines.ordered) 
}


couple.color = c("blue", "yellow")

diseases <- vector(mode="list", length=7)
names(diseases) <- c("ID", "ASD.Autistic.traits", "Epilepsy", "Microcephaly", "Macrocephaly", "Hypotonia", "Ataxia")
diseases[[1]] <- 2 #7
diseases[[2]] <- 3 #8
diseases[[3]] <- 4 #9
diseases[[4]] <- 5 #10
diseases[[5]] <- 6 #11
diseases[[6]] <- 7 #12
diseases[[7]] <- 8 #13

#Import real values
Answers.table  <- read.table('./data/experimental_value/CAGI_answer_key_V02.txt', sep = '\t', quote = "", header = TRUE, stringsAsFactors = FALSE) #file with experimental data

exp.val <- cbind(Answers.table$Name, Answers.table$ID, Answers.table$ASD.Autistic.traits, Answers.table$Epilepsy,
                 Answers.table$Microcephaly, Answers.table$Macrocephaly, Answers.table$Hypotonia, Answers.table$Ataxia)
colnames(exp.val) <- c( "Name", "ID", "ASD.Autistic.traits", "Epilepsy", "Microcephaly", "Macrocephaly", "Hypotonia", "Ataxia")
exp.val <- data.frame(exp.val, stringsAsFactors = F)

#Import submission
path.prediction <- './data/submission/'
submission.files <- sort(list.files(path = path.prediction, pattern = "Group", full.names = TRUE))
submission.files.name <- gsub("-prediction_file-", ".", gsub(".txt$", "", gsub("^Group_", "", basename(submission.files))))
if(any(grepl("-late", submission.files.name))){
  submission.files.name <- gsub("-late", "", submission.files.name)
}

# extract patient ID from template file 
sub.template <- './data/template/template.txt'  # submission template
temp <- strsplit(scan(file = sub.template, what = character(), sep = '\n', quiet = T), '\t')
temp <- sapply(2:length(temp), function(x) {temp[[x]][1]})

# extract experimental classes
exp.val.classes <- exp.val[, 2:ncol(exp.val)]

# summary of patient classes
tab.class <- sapply(1:nrow(exp.val.classes), function(x) {paste(colnames(exp.val.classes)[which(exp.val.classes[x, ]==1)], collapse = ", ")})
clas.s.t <-unique(tab.class) # type of elements
counts.s.t <- sapply(1:length(clas.s.t), function(x){length(grep(paste("^", clas.s.t[x], "$", sep = ""), tab.class))})
tabel.1 <- cbind(clas.s.t, counts.s.t)
colnames(tabel.1) <- c("class detected", "number of patients")
write.table(tabel.1[order(-counts.s.t), ], file = './results/class_phenotype_patients_detected.txt', sep = "\t", quote = FALSE, row.names = F)

default.group.colors <- c("black", "blue", "magenta", "red", "orange", "gold", "pink", "yellow", "green", "grey", "brown")
group.colors <- c(default.group.colors[1], rep(default.group.colors[2], 6), rep(default.group.colors[3], 3), rep(default.group.colors[5], 3), default.group.colors[7])



# pie chart of patients with and without variants
pv <- sapply(1:nrow(Answers.table), function(x) {all(is.na(Answers.table[x, c(10:12)]))})
p.w.v <- which(pv==T)
p.v <- which(pv==F)


##################
# figure A #######
##################

if(TRUE){
  png("./results/figure_A.png", width = 8, height = 6, units = 'in', res = 300)
  pat.no.var <- length(p.w.v)/150*100
  pat.var <- length(p.v)/150*100
  slices <- c(pat.no.var, pat.var)
  lbls <- c(paste("No associated variants ", round(pat.no.var, digits = 0), "%, ", length(p.w.v), sep=""), paste("At least one noted variant ", round(pat.var, digits = 0),"%, ", length(p.v), sep = "" ))
  pie(slices, labels = lbls, col=c('red', 'blue'), main="Presence of annotated variants in patients of the ID challange")
  dev.off()
}

###########################################
# extract information about variants#######
###########################################


# create causing.variants.table from CAGI_ID_document
Pat.caus.var <- which(is.na(Answers.table$Disease.causing) == F) # check ID of patients with causative variants
caus.var.table <- cbind(Answers.table$Name[Pat.caus.var], Answers.table$Disease.causing[Pat.caus.var])
# check for patients with multiple causative variants
multiplevariants <- sapply(1:length(caus.var.table[, 2]), function(x){length(strsplit(caus.var.table[x, 2], ',')[[1]])})
# create new table with one row for each variant
causing.variants.table <- c()
# patient with one variant

for(i in (1:nrow(caus.var.table))){
  if(multiplevariants[i] == 1){
    causing.variants.table <- rbind(causing.variants.table, caus.var.table[i,])
  } else {
    patient.multiple.variants <- strsplit(caus.var.table[i, 2], ',')[[1]]
    additional.lines <- c()
    for (j in (1:multiplevariants[i])){
      additional.lines <- rbind(additional.lines, c(caus.var.table[i,1], patient.multiple.variants[j] ))
    }
    causing.variants.table <- rbind(causing.variants.table, additional.lines)
  }
}

colnames(causing.variants.table) <- c("Patient.Code", "Chr")
causing.variants.table <- data.frame(causing.variants.table, stringsAsFactors = F)
causing.variants <- as.vector(causing.variants.table$Chr)
names(causing.variants) <- causing.variants.table$Patient.Code

# get putative causative variants
Pat.pc.var <- which(is.na(Answers.table$Putative) == F) # check ID of patients with causative variants
pc.var.table <- cbind(Answers.table$Name[Pat.pc.var], Answers.table$Putative[Pat.pc.var])
# check for patients with multiple causative variants
multiplevariants <- sapply(1:length(pc.var.table[, 2]), function(x){length(strsplit(pc.var.table[x, 2], ',')[[1]])})
# create new table with one row for each variant
pc.variants.table <- c()
# patient with one variant

for(i in (1:nrow(pc.var.table))){
  if(multiplevariants[i] == 1){
    pc.variants.table <- rbind(pc.variants.table, pc.var.table[i,])
  } else {
    patient.multiple.variants <- strsplit(pc.var.table[i, 2], ',')[[1]]
    additional.lines <- c()
    for (j in (1:multiplevariants[i])){
      additional.lines <- rbind(additional.lines, c(pc.var.table[i,1], patient.multiple.variants[j] ))
    }
    pc.variants.table <- rbind(pc.variants.table, additional.lines)
  }
}
colnames(pc.variants.table) <- c("Patient.Code", "Chr")
putative.causing.variants.table <- data.frame(pc.variants.table, stringsAsFactors = F)

putative.causing.variants <- as.vector(putative.causing.variants.table$Chr)
names(putative.causing.variants) <- putative.causing.variants.table$Patient.Code

# get contributing factor variants
Cont.fact.var <- which(is.na(Answers.table$Contributing.factor) == F) # check ID of patients with causative variants
Cont.fact.var.table <- cbind(Answers.table$Name[Cont.fact.var], Answers.table$Contributing.factor[Cont.fact.var])
# check for patients with multiple causative variants
multiplevariants <- sapply(1:length(Cont.fact.var.table[, 2]), function(x){length(strsplit(Cont.fact.var.table[x, 2], ',')[[1]])})
# create new table with one row for each variant
cf.variants.table <- c()
# patient with one variant

for(i in (1:nrow(Cont.fact.var.table))){
  if(multiplevariants[i] == 1){
    cf.variants.table <- rbind(cf.variants.table, Cont.fact.var.table[i,])
  } else {
    patient.multiple.variants <- strsplit(Cont.fact.var.table[i, 2], ',')[[1]]
    additional.lines <- c()
    for (j in (1:multiplevariants[i])){
      additional.lines <- rbind(additional.lines, c(Cont.fact.var.table[i,1], patient.multiple.variants[j] ))
    }
    cf.variants.table <- rbind(cf.variants.table, additional.lines)
  }
}

colnames(cf.variants.table) <- c("Patient.Code", "Chr")
contributing.factor.variants.table <- data.frame(cf.variants.table, stringsAsFactors = F)

contributing.factor.variants <- as.vector(contributing.factor.variants.table$Chr)
names(contributing.factor.variants) <- contributing.factor.variants.table$Patient.Code

##################
# figure B #######
##################

if(TRUE)
{
  png('./results/figure_B.png', width = 3, height = 4, units = 'in', res = 300)
  data <- c(length(unique(causing.variants.table$Chr)), length(unique(putative.causing.variants.table$Chr)), 
            length(unique(contributing.factor.variants.table$Chr)))
  plt <- barplot(data, col=c(couple.color, "purple"), xaxt="n", ylim=c(0, 35))
  text(mean(plt), 37, labels = "Known variants statistic", xpd = TRUE, cex=0.8, font=2) 
  text(plt, data+1, labels = data, xpd = TRUE, cex=0.8, col = "black")
  text(plt, -1, labels = c("Causative", "Putative", "Contributing"), xpd = TRUE, cex=0.6) 
  dev.off()
}

if(TRUE){
  #############
  # Table 4 #
  ############
  
  # make variant statistics:
  matr.cv <- matrix(data = "empty", nrow = length(submission.files), ncol = 2)
  matr.pcv <- matrix(data = "empty", nrow = length(submission.files), ncol = 2)
  matr.cf <- matrix(data = "empty", nrow = length(submission.files), ncol = 2)
  
  for(j in 1:length(submission.files)) {
    CurrentSubmission <- submission.files[j]
    #print("***************************************************************")
    #print(paste(basename(CurrentSubmission), "->", j))
    lines <- scan(file = CurrentSubmission ,what = character(), sep = "\n", quiet = T)
    header <- strsplit(lines[1], "\t")[[1]]
    sub.list <- strsplit(row.corrector(template = temp, sub.lines = lines), '\t')
    ID <- as.vector(sapply(1:length(sub.list), function(x) suppressWarnings(as.character(sub.list[[x]][1]))))
    col_of_V <- grepl('-V$', header)
    V <- t(as.data.frame(sapply(1:length(sub.list), function(x) suppressWarnings(as.character(sub.list[[x]][col_of_V])))))
    rownames(V) <- ID
    # causing variants match
    match.var <- sapply(1:length(causing.variants), function(x){length(grep(causing.variants[x], V[which(rownames(V)==names(causing.variants)[x]), ]))})
    match.var <-  length(which( match.var > 0))
    matr.cv[j, ] <- c(gsub("-prediction_file-", ".", gsub(".txt$", "", gsub("^Group_", "", basename(submission.files[j])))), match.var)
    # putative causing variants match
    match.var <- sapply(1:length(putative.causing.variants), function(x){length(grep(putative.causing.variants[x], V[which(rownames(V)==names(putative.causing.variants)[x]), ]))})
    match.var <-  length(which( match.var > 0))
    matr.pcv[j, ] <- c(gsub("-prediction_file-", ".", gsub(".txt$", "", gsub("^Group_", "", basename(submission.files[j])))), match.var)
    # causing factors match
    match.var <- sapply(1:length(contributing.factor.variants), function(x){length(grep(contributing.factor.variants[x], V[which(rownames(V)==names(contributing.factor.variants)[x]), ]))})
    match.var <-  length(which( match.var > 0))
    matr.cf[j, ] <- c(gsub("-prediction_file-", ".", gsub(".txt$", "", gsub("^Group_", "", basename(submission.files[j])))), match.var)
  }
  
  colnames(matr.cv) <- c("Submission ID", "# causative variant found")
  colnames(matr.pcv) <- c("Submission ID", "# putative causative variant found")
  colnames(matr.cf) <- c("Submission ID", "# contributing factor variant found")
  # create causative and putative table
  c.pc.matrix <- cbind(matr.cv, matr.pcv[, 2])
  # create causative, putative table and contributing factor table
  c.pc.cf.matrix <- cbind(c.pc.matrix, matr.cf[, 2])
  row.sum <- as.numeric(matr.cv[, 2])+as.numeric(matr.pcv[, 2])+as.numeric(matr.cf[, 2])
  total.vars.results <- cbind(c.pc.cf.matrix, row.sum, rank(-row.sum, ties.method = "min"))
  colnames(total.vars.results) <- c("Submission ID", "# causative variants", "# putative causative variants", "# contributing factor variants", "total predicted variants", "rank")

  # find total number of variants for each prediction
  V.list <-  vector(length = length(submission.files))
  ID <- rep("empty", length(submission.files))
  
  for (i in 1:length(submission.files)) {
    prediction <- strsplit(scan(file = submission.files[i], what = character(), sep = "\n", quiet = T) , "\t")
    #Sub.name <- strsplit(submission.files.name[i],"(?=[.])", perl = TRUE)[[1]][1]
    Sub.name <- strsplit(submission.files.name[i],"(?=[-._])", perl = TRUE)[[1]]
    ID[i] <- paste(Sub.name[grep('[0-9]', Sub.name)], collapse  = ".")
    header <- prediction[[1]]
    col_of_V <- grepl('-V', header)
    V <- sapply(2:length(prediction), function(x) prediction[[x]][col_of_V])
    V <- t(V)
    for(k in 1:nrow(V)){
      for(j in 1:ncol(V)){
        V[k,j] <- paste(unique(strsplit(V[k,j], ",")[[1]]), collapse = ",") # someone put the same variant multiple times in the same cell
      }
    }
    V <- as.character(V) #transpose P
    V <- unique(unlist(sapply(1:length(V), function(x){unlist(strsplit(V[x], ","))})))
    V.list[i] <- length(which(V!="*"))
  }
  
  ast <- t(rbind(ID, V.list))
  table.four <- cbind(ast[, 1], total.vars.results[, 5], ast[, 2], round(as.numeric(total.vars.results[, 5])/56, digits = 2), round(as.numeric(total.vars.results[, 5])/V.list, digits = 2))
  colnames(table.four) <- c("Submission ID", "Correctly predicted variants", "Total predicted variants", "Correctly predicted variants/Experimental variants", "Correctly predicted variants/Predicted Variants")
  write.table(table.four, file = './results/table4.txt', sep = "\t", quote = FALSE, row.names = F, col.names = T)

  #############
  # Figure 5  #
  #############
  
  # create histogram of all variants
  table.var.pred <- total.vars.results[, c(1:(ncol(total.vars.results)-2))]
  table.var.histog <- matrix(NA, nrow = 1, ncol = 3)
  for (i in 1:nrow(table.var.pred)) {
    table.var.histog <- rbind(table.var.histog, cbind(table.var.pred[i, 1], table.var.pred[i, 2], "causative"))
    table.var.histog <- rbind(table.var.histog, cbind(table.var.pred[i, 1], table.var.pred[i, 3], "putative"))
    table.var.histog <- rbind(table.var.histog, cbind(table.var.pred[i, 1], table.var.pred[i, 4], "contributing factor"))
  }
  table.var.histog <- table.var.histog[c(2:nrow(table.var.histog)), ] # remove empty row
  table.var.histog[, 1] <- gsub("-late", "", table.var.histog[, 1])
  colnames(table.var.histog) <- c("Submission ID", "Number", "Variant")
  rownames(table.var.histog) <- c()
  table.var.histog <- rbind(table.var.histog, rbind(c("EV",length(unique(causing.variants.table$Chr)), "causative"), c("EV", length(unique(putative.causing.variants.table$Chr)), "putative"), 
                                                    c("EV", length(unique(contributing.factor.variants.table$Chr)), "contributing factor") ))
  table.var.histog <- data.frame(table.var.histog, stringsAsFactors = F)
  table.var.histog$Submission.ID <- as.factor(table.var.histog$Submission.ID)
  table.var.histog$Number <- as.numeric(table.var.histog$Number)
  table.var.histog$Variant <- as.factor(table.var.histog$Variant)
  table.var.histog$Variant <- ordered(table.var.histog$Variant, levels = c("causative", "putative", "contributing factor"))
  png("./results/figure_5.png", width = 8, height = 4, units = 'in', res = 300)
  p <- ggplot(table.var.histog, aes(x = Submission.ID, y = Number, fill = Variant)) +
    ggtitle("Predicted variants barplot") +
    geom_bar(stat = "identity", position = position_dodge(0.9))+theme_bw() + scale_fill_manual("Variant", values = c("causative" = "blue", "putative" = "yellow", "contributing factor" = "purple")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.grid.major = element_blank(),
                                                                                                                                                                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                                                                                                                             plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
                                                                                                                                                                                                             #axis.text.y = element_text(size = 14),
                                                                                                                                                                                                             axis.title.y = element_text(size = 14, face = "bold"),
                                                                                                                                                                                                             #axis.text.x = element_text(size = 14),
                                                                                                                                                                                                             axis.title.x = element_text(size = 14, face = "bold"),
                                                                                                                                                                                                             axis.text = element_text(size = 8))
  p <- p + scale_x_discrete(limits=c("EV", gsub("-late", "", table.var.pred[order(-as.numeric(table.var.pred[, 2])), 1]))) + theme(axis.text = element_text(size = 12)) + labs(x = "Submission ID", y = "Number of predicted variants") + guides(fill=guide_legend(title="Type of variant"))
  print(p)
  dev.off()
}
table.caus <- cbind(contributing.factor.variants.table, rep("Contributing factor", nrow(contributing.factor.variants.table)))
colnames(table.caus) <- c("ID", "Chr", "Type")
table.dis <- cbind(causing.variants.table, rep("Disease causing", nrow(causing.variants.table)))
colnames(table.dis) <- c("ID", "Chr", "Type")
table.put <- cbind(putative.causing.variants.table, rep("Putative", nrow(putative.causing.variants.table)))
colnames(table.put) <- c("ID", "Chr", "Type")
var.general.table <- rbind(table.caus, table.dis, table.put)

#################
# table S1 ######
#################


if(TRUE){
  
  V.list <-  vector(mode = 'list', length = length(submission.files))
  ID <- rep("empty", length(submission.files))
  for (i in 1:length(submission.files)) {
    prediction <- strsplit(scan(file = submission.files[i], what = character(), sep = "\n", quiet = T) , "\t")
    Sub.name <- strsplit(submission.files.name[i],"(?=[-._])", perl = TRUE)[[1]]
    ID[i] <- paste(Sub.name[grep('[0-9]', Sub.name)], collapse  = ".")
    header <- prediction[[1]]
    col_of_V <- grepl('-V', header)
    V <- sapply(2:length(prediction), function(x) prediction[[x]][col_of_V])
    V <- t(V)
    for(k in 1:nrow(V)){
      for(j in 1:ncol(V)){
        V[k,j] <- paste(unique(strsplit(V[k,j], ",")[[1]]), collapse = ",") # someone put the same variant multiple times in the same cell
      }
    }
    V <- as.character(V)
    V <- unlist(sapply(1:length(V), function(x){unlist(strsplit(V[x], ","))}))
    V.list[[i]] <- V
  }
  Variant.Matrix <- matrix(data = 0, ncol =length(submission.files), nrow = nrow(var.general.table) )
  colnames(Variant.Matrix) <- ID
  for (i in 1:nrow(var.general.table)){
    for(j in 1:length(submission.files)){
      Variant.Matrix[i, j] <- length(grep(paste("^", var.general.table[i, 2], "$",  sep = ""), V.list[[j]]))
    }
  }
  number.of.pred.predicting <- sapply(1:nrow(Variant.Matrix), function(x){length(which(Variant.Matrix[x, ]>0))})
  id.of.predictors <- sapply(1:nrow(Variant.Matrix), function(x){paste(colnames(Variant.Matrix)[which(Variant.Matrix[x, ]>0)], collapse =", ")})
  # find number of group correctly predicting
  group.var.Matrix <- cbind(Variant.Matrix[,1], rowSums(Variant.Matrix[, c(2:7)]), rowSums(Variant.Matrix[, c(8:10)]), rowSums(Variant.Matrix[, c(11:13)]), Variant.Matrix[, 14]) 
  colnames(group.var.Matrix) <- seq(1:5)
  number.of.group.predicting <- sapply(1:nrow(group.var.Matrix), function(x){length(which(group.var.Matrix[x, ]>0))})
  id.of.groups.predicting <- sapply(1:nrow(group.var.Matrix), function(x){paste(colnames(group.var.Matrix)[which(group.var.Matrix[x, ]>0)], collapse = ", ")})
  total.by.variant <- rowSums(Variant.Matrix)
  total.by.group <- colSums(Variant.Matrix)
  Results.matrix <- cbind(var.general.table[, c(1,3,2)], number.of.group.predicting, id.of.groups.predicting, number.of.pred.predicting, id.of.predictors)
  colnames(Results.matrix) <- c("Name", "Class of variant", "Chr", "Number of Groups", "ID of Groups", "Number of predictions", "Predictions ID")
  Results.matrix <- Results.matrix[, c(1, 2, 3, 4, 6, 5, 7)]
  t.cf <- Results.matrix[grepl("Contributing factor", Results.matrix$`Class of variant` ),] [order(-Results.matrix$`Number of predictions`[grepl("Contributing factor", Results.matrix$`Class of variant` )]), ]
  t.dc <- Results.matrix[grepl("Disease causing", Results.matrix$`Class of variant` ),] [order(-Results.matrix$`Number of predictions`[grepl("Disease causing", Results.matrix$`Class of variant` )]), ]
  t.p <- Results.matrix[grepl("Putative", Results.matrix$`Class of variant` ),] [order(-Results.matrix$`Number of predictions`[grepl("Putative", Results.matrix$`Class of variant` )]), ]
  Results.matrix <- rbind(t.cf, t.dc, t.p)
  write.table(file = './results/tableS1.txt', Results.matrix, sep = '\t', row.names = F, col.names = T, quote = F)
  
  #################
  # figure 6 ######
  #################
  
  Results.matrix <- rbind(cbind(rep("Contributing", length(summary(as.factor(t.cf$`Number of Groups`)))), summary(as.factor(t.cf$`Number of Groups`)), names(summary(as.factor(t.cf$`Number of Groups`)))),
        cbind(rep("Causative", length(summary(as.factor(t.dc$`Number of Groups`)))), summary(as.factor(t.dc$`Number of Groups`)), names(summary(as.factor(t.dc$`Number of Groups`)))),
        cbind(rep("Putative", length(summary(as.factor(t.p$`Number of Groups`)))), summary(as.factor(t.p$`Number of Groups`)), names(summary(as.factor(t.p$`Number of Groups`)))))
  Results.matrix <- suppressWarnings(data.frame(Results.matrix, stringsAsFactors = F))
  colnames(Results.matrix) <- c("variant", "count", "NG")
  Results.matrix$variant <- factor(Results.matrix$variant, levels = c("Causative", "Putative", "Contributing"))
  Results.matrix$NG <- as.factor(Results.matrix$NG)
  Results.matrix$count <- as.numeric(Results.matrix$count)
  
  png("./results/figure_6.png", width = 8, height = 6, units = 'in', res = 300)
  p <- ggplot(Results.matrix, aes(x = variant, y = count, fill = NG)) + geom_bar(stat = "identity") + xlab("") + ggtitle("Number of groups correctly predicting variants") + 
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 8)) + 
    guides(fill=guide_legend(title="Number of Groups")) +
    ylab("Number of variants")
    print(p+scale_fill_brewer(palette = "Spectral"))
  dev.off()    

}

##################
# figure C #######
##################

# Statistical Validation
if(TRUE)
{ 
    statistical.data.complete = matrix(NA, 4, length(diseases))
    rownames(statistical.data.complete) <- c("Yes", "No", "NA", "Total")
    colnames(statistical.data.complete) <- c("ID", "ASD Autistic traits", "Epilepsy", "Microcephaly", "Macrocephaly", "Hypotonia", "Ataxia")
    for(column.disease.actual in (1:length(diseases)))
    {
      id <- column.disease.actual
      statistical.data.complete[1, id] = length(which(exp.val[, (column.disease.actual+1)]=="1"))
      statistical.data.complete[2, id] = length(which(exp.val[, (column.disease.actual+1)]=="0"))
      statistical.data.complete[3, id] = length(which(is.na(exp.val[, (column.disease.actual+1)])))
      statistical.data.complete[4, id] = sum(statistical.data.complete[1:2, id])
    }
    statistical.data <- statistical.data.complete[4, ]
    png('./results/figure_C.png', width = 5, height = 4, units = 'in', res = 300)
    plt <- barplot(statistical.data, col='steelblue', xaxt="n", ylim=c(0, 170))
    text(mean(plt), (max(statistical.data) + max(statistical.data)/4), labels = "Patients with known features", xpd = TRUE, cex=1.2, font=2) 
    xleft <- plt-0.5
    ybottom <- 0
    xright <- plt+0.5
    ytop1 <- statistical.data.complete[1,]
    ytop2 <- statistical.data.complete[2,]
    rect(xleft,ybottom,xright,ytop1,col=couple.color[1]) 
    rect(xleft,ytop1,xright,ytop1+ytop2,col=couple.color[2])
    text(plt, statistical.data+5, labels = statistical.data, xpd = TRUE, cex=0.8, col = "black" )
    text(plt, par("usr")[3], labels = colnames(statistical.data.complete), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=0.6) 
    legend("topright", legend=c("Yes", "No"), fill=couple.color, cex=0.6)
    dev.off()
    write.table(statistical.data.complete, file = './results/figure_C_statistics.txt', sep = "\t", quote = FALSE)
}

################################################
## Supplementary table 2, ROC plots ############ 
################################################

############### ROC BY DISEASE FR CORRECTED ##############
if(TRUE){
  matrix.threshold <- matrix(NA, nrow = length(submission.files), ncol = length(diseases) ) # create matirx of threshold (one for each submission and class), used later for single patients analysis
  occurrence.by.disease <- matrix(NA, length(submission.files), length(diseases))
  heat.map.performances <- matrix(NA, length(submission.files), length(diseases))
  rownames(heat.map.performances) <- paste("Submission", submission.files.name)
  colnames(heat.map.performances) <- names(diseases)
  scores <- c("AUC", "MCC", "ACC", "F1", "TPR", "PPV", "TNR", "NPV", "FNR", "TP", "TN", "FN", "FP")
  performances <- matrix(NA, length(submission.files), length(scores))
  rownames(performances) <- paste("Submission", submission.files.name)
  colnames(performances) <- scores
  color.group = c(rep("red", 1), rep("blue", 6), rep("purple", 3), rep("black", 3), "pink")
  #png("results/AUC_corrected_best_del.png", width = 5, height = 10, units = 'in', res = 300)
  #par(mfrow = c(4, 2)) #, pty='s') # decide number of plots in each row of the picture and total number of lines
  
  # prepare supplementary table 2
  sup.tab2 <- submission.files.name
  
  for(d in 1:length(exp.val.classes))
  {
    print(paste("DISEASE:", names(diseases)[d]))
     exp.val.classes.labels <- matrix(as.numeric(exp.val.classes[, d]))
    png(paste("results/AUC_corrected_", names(diseases)[d],".png", sep = ""), width = 5, height = 5, units = 'in', res = 300)
    pred <- prediction(seq(0.0,1.0,by=0.0001), c(1, rep(c(0,1), 10001/2)))
    roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
    plot(roc.perf, cex = 0.6)
    text(0.5, 1.1, labels = paste("ROC among submissions in", gsub("\\.", " ", names(diseases)[d])), xpd = TRUE, cex=0.8, font=2) 
    rocs <- c(roc.perf)
    abline(a=0, b=1, col="white", lty=2)
    color.group = c(rep("red", 1), rep("blue", 6), rep("purple", 3), rep("black", 3), "pink")
    legend("bottomright", legend = c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5"), fill = c('red', 'blue', 'purple', 'black', "pink"), cex = 0.6)
    for(j in 1:length(submission.files)) {
      CurrentSubmission <- submission.files[j]
      print(paste(basename(CurrentSubmission), "->", j))
      lines <- scan(file = CurrentSubmission ,what = character(), sep = "\n", quiet = T)
      header <- strsplit(lines[1], "\t")[[1]]
      sub.list <- strsplit(row.corrector(template = temp, sub.lines = lines), '\t')
      col_of_P <- grepl('-P$', header)
      P <- as.vector(sapply(1:length(sub.list), function(x) suppressWarnings(as.numeric(sub.list[[x]][col_of_P][d]))))
      labels <- exp.val.classes.labels[which(!is.na(exp.val.classes.labels))]
      values <- P[which(!is.na(exp.val.classes.labels))]
      occurrence.by.disease[j, d] <- length(which(!is.na(as.numeric(values))))
      values[which(is.na(values))] <- 0 # P values * are converted to 0
    
      if(length(unique(labels[which(!is.na(values))]))==2) 
      {
        pred <- prediction(values, labels)
        
        roc.perf <- performance(pred, measure = "tpr", x.measure = "fpr")
        auc.best.perf <- max(performance(pred, "auc")@y.values[[1]], na.rm = TRUE)
        
        mcc.perf <- performance(pred, "mat")
        mcc.best.perf.id <- which.max(mcc.perf@y.values[[1]])
        mcc.best.perf <- mcc.perf@y.values[[1]][mcc.best.perf.id]
        
        best.split = mcc.perf@x.values[[1]][mcc.best.perf.id]
        
        acc.perf <- performance(pred, "acc")
        acc.best.perf.by.mcc <- acc.perf@y.values[[1]][mcc.best.perf.id]
        
        f1.perf <- performance(pred, "f")
        f1.best.perf.by.mcc <- f1.perf@y.values[[1]][mcc.best.perf.id]
        
        tpr.perf <- performance(pred, "tpr")
        tpr.best.perf.by.mcc <- tpr.perf@y.values[[1]][mcc.best.perf.id]
        
        ppv.perf <- performance(pred, "ppv")
        ppv.best.perf.by.mcc <- ppv.perf@y.values[[1]][mcc.best.perf.id]
        
        tnr.perf <- performance(pred, "tnr")
        tnr.best.perf.by.mcc <- tnr.perf@y.values[[1]][mcc.best.perf.id]
        
        npv.perf <- performance(pred, "npv")
        npv.best.perf.by.mcc <- npv.perf@y.values[[1]][mcc.best.perf.id]
        
        fnr.perf <- performance(pred, "fnr")
        fnr.best.perf.by.mcc <- fnr.perf@y.values[[1]][mcc.best.perf.id]
        
        tp=pred@tp[[1]][mcc.best.perf.id]
        fn=pred@fn[[1]][mcc.best.perf.id]
        fp=pred@fp[[1]][mcc.best.perf.id]
        tn=pred@tn[[1]][mcc.best.perf.id]
        
        performances[j,1] <- auc.best.perf
        performances[j,2] <- mcc.best.perf
        performances[j,3] <- acc.best.perf.by.mcc
        performances[j,4] <- f1.best.perf.by.mcc
        performances[j,5] <- tpr.best.perf.by.mcc
        performances[j,6] <- ppv.best.perf.by.mcc
        performances[j,7] <- tnr.best.perf.by.mcc
        performances[j,8] <- npv.best.perf.by.mcc
        performances[j,9] <- fnr.best.perf.by.mcc
        performances[j,10] <- tp
        performances[j,11] <- tn
        performances[j,12] <- fn
        performances[j,13] <- fp
        
        heat.map.performances[j, d] <- auc.best.perf
        plot(roc.perf, col = color.group[j], add = TRUE)

        # store MCC threshold for other functions
        matrix.threshold[j, d] <- pred@cutoffs[[1]][mcc.best.perf.id] # extract threshold of best MCC
      }
      else
        print("Not enough classes")
    }
    dev.off()
    
    sup.tab2 <- cbind(sup.tab2, round(performances[, 1:4], 2))
    #write.table(performances, file = paste('./results/roc_performances_not_rounded_', names(diseases)[d],'.txt', sep =""), sep = "\t", quote = FALSE)
    write.table(round(performances, 2), file = paste('./results/roc_performances_', names(diseases)[d],'.txt', sep =""), sep = "\t", quote = FALSE)
    #write.table(round(performances[order(performances[, 1], decreasing=TRUE), ][1:3, ], 2), file = paste('./results/roc_performances_best_', names(diseases)[d],'.txt', sep =""), sep = "\t", quote = FALSE)
  }
  colnames(sup.tab2)[1] <- "Submission ID"
  sup.tab2 <- rbind(c("Disease", rep("ID", 4), rep("ASD Autistic traits", 4), rep("Epilepsy", 4), rep("Microcephaly", 4),
    rep("Macrocephaly", 4), rep("Hypotonia", 4), rep("Ataxia", 4)), colnames(sup.tab2), sup.tab2)
  write.table(sup.tab2, file = './results/supplementary_table_2.txt', sep='\t', quote = FALSE, col.names = F, row.names = F )
  #dev.off()
}

 # create a table showing group and predictors that have classified a patient as case (one table for each class)

if(TRUE){
  ID.col <- c("disease", submission.files.name)
  disease.matrix.list <- vector(mode = "list", length = 7)
  for(d in 1:length(exp.val.classes)) {
    exp.val.classes.labels <- matrix(as.numeric(exp.val.classes[, d]))
    matrix.disease  <- matrix(NA, nrow = length(exp.val.classes.labels), ncol = length(submission.files))
    for(j in 1:length(submission.files)) {
      CurrentSubmission <- submission.files[j]
      lines <- scan(file = CurrentSubmission ,what = character(), sep = "\n", quiet = T)
      header <- strsplit(lines[1], "\t")[[1]]
      sub.list <- strsplit(row.corrector(template = temp, sub.lines = lines), '\t')
      col_of_P <- grepl('-P$', header)
      P <- as.vector(sapply(1:length(sub.list), function(x) suppressWarnings(as.numeric(sub.list[[x]][col_of_P][d]))))
      values <- P
      values[which(is.na(values))] <- 0
      pred.class <- values
      pred.class[ which(values[] >= matrix.threshold[j, d]) ] <- 1
      pred.class[ which( values[] < matrix.threshold[j, d]) ] <- 0
      matrix.disease[, j] <- pred.class
    }
    matrix.disease <- cbind(exp.val.classes.labels, matrix.disease)
    colnames(matrix.disease) <- ID.col
    disease.matrix.list[[d]] <- matrix.disease 
    
  }
  
  for(i in 1:length(disease.matrix.list)){
    Variant.Matrix <- disease.matrix.list[[i]][, 2:ncol( disease.matrix.list[[i]])]
    number.of.pred.predicting <- sapply(1:nrow(Variant.Matrix), function(x){length(which(Variant.Matrix[x, ]>0))})
    id.of.predictors <- sapply(1:nrow(Variant.Matrix), function(x){paste(colnames(Variant.Matrix)[which(Variant.Matrix[x, ]>0)], collapse =", ")})
    # find number of group correctly predicting
    group.var.Matrix <- cbind(Variant.Matrix[,1], rowSums(Variant.Matrix[, c(2:7)]), rowSums(Variant.Matrix[, c(8:10)]), rowSums(Variant.Matrix[, c(11:13)]), Variant.Matrix[, 14]) 
    colnames(group.var.Matrix) <- seq(1:5)
    number.of.group.predicting <- sapply(1:nrow(group.var.Matrix), function(x){length(which(group.var.Matrix[x, ]>0))})
    id.of.groups.predicting <- sapply(1:nrow(group.var.Matrix), function(x){paste(colnames(group.var.Matrix)[which(group.var.Matrix[x, ]>0)], collapse = ", ")})
    
    total.by.variant <- rowSums(Variant.Matrix)
    total.by.group <- colSums(Variant.Matrix)
    Results.matrix <- cbind(disease.matrix.list[[i]][, 1 ], Variant.Matrix, total.by.variant, number.of.group.predicting, id.of.groups.predicting, number.of.pred.predicting, id.of.predictors)
    # add total by row to matrix
    colnames(Results.matrix)[c(1, (length(colnames(Results.matrix))-4):length(colnames(Results.matrix)))] <-c("Presence", "total predicted", "Number of Groups predicting", "ID of Groups predicting", "Total numbers of predictions", "Predictions ID")
    Results.matrix <- cbind(c(1:length(temp)), temp, Results.matrix[, c(1, 19, 17, 18, 20) ])
    colnames(Results.matrix)[c(1,2)] <- c("patient CAGI ID", "patient ID")
    write.table(file = paste('./results/match_', names(diseases)[i],'_FR.txt', sep = ""), Results.matrix, sep = '\t', row.names = F, col.names = T, quote = F)
  }
}

#############
# Figure 2 #
############

if(TRUE){
  #### histograms like Hopkins, with number of groups that have correctly predicted by at least one group etc....
  
  tables.disease <- list.files(path="./results/", pattern = "match", full.names = T)
  tab.hist <- matrix(c(NA,NA,NA,NA), 1, 4)
  for(i in (1:7)){
    tab.pheno <- read.table(tables.disease[grep(colnames(exp.val.classes)[i], tables.disease)], sep = "\t", header = T)
    tab.pheno <- tab.pheno[which(tab.pheno[, 3]==1), ] # consider only patients that were cases
    tab.temp <- summary(as.factor(tab.pheno$Number.of.Groups.predicting))
    tab.hist <- rbind(tab.hist, cbind(rep(colnames(exp.val.classes)[i], length(tab.temp)), names(tab.temp), as.numeric(tab.temp), rep(sum(as.numeric(tab.temp)), length(tab.temp))))
  }
  tab.hist <- tab.hist[which(is.na(tab.hist[, 1])!=T), ] # remove NA row
  tab.hist[, 1] <- gsub("ASD.Autistic.traits", "ASD Autistic traits", tab.hist[, 1])
  colnames(tab.hist)<- c("Submission.ID", "Variant", "Number", "total")
  tab.hist <- data.frame(tab.hist, stringsAsFactors = F)
  tab.hist$Submission.ID <- as.factor(tab.hist$Submission.ID)
  tab.hist$Number <- as.numeric(tab.hist$Number)
  tab.hist$Variant <- as.factor(tab.hist$Variant)
  tab.hist$total <- as.numeric(tab.hist$total)
  png("./results/figure_2.png", width = 8, height = 6, units = 'in', res = 300)
  p <- ggplot(tab.hist, aes(x = Submission.ID, y = Number, fill = Variant)) +
    ggtitle("Number of groups correctly predicting patient phenotype") + geom_bar(stat = "identity")+theme_bw()+geom_text(aes(Submission.ID, total, label = total, fill = NULL), vjust = - 0.8, data = tab.hist)+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=2), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          #axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 13, face = "bold"),
          #axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 8)) + labs(y = "Number of patients", x = "") + guides(fill=guide_legend(title="Number of groups"))
  print(p + scale_x_discrete(limits=c("ID", "ASD Autistic traits", "Epilepsy", "Microcephaly", "Macrocephaly", "Hypotonia", "Ataxia"))+scale_fill_brewer(palette = "Spectral"))
  dev.off()
  
}

###### P value and thresholds histograms TP, TN, FP, FN by submission 

####################################
# Figure 1 Supplementary Materials #
####################################

if(TRUE){
  MakeSupFig1 <- function(variants.only = F){
    Plot.list <- vector("list", length = 7)
    for(d in (1:7)){
      tab.hist <- matrix(c(NA,NA,NA), 1, 3)
      patient.pred <- c() # number of patients predicted
      if(variants.only){
        exp.val.classes.labels <- exp.val.classes[p.v, d]
      }else{
        exp.val.classes.labels <- exp.val.classes[, d]
      }
      for(j in 1:length(submission.files)) {
        CurrentSubmission <- submission.files[j]
        lines <- scan(file = CurrentSubmission ,what = character(), sep = "\n", quiet = T)
        header <- strsplit(lines[1], "\t")[[1]]
        sub.list <- strsplit(row.corrector(template = temp, sub.lines = lines), '\t')
        col_of_P <- grepl('-P$', header)
        P <- as.vector(sapply(1:length(sub.list), function(x) suppressWarnings(as.numeric(sub.list[[x]][col_of_P][d]))))
        if(variants.only){
          values <- P[p.v]
        }else{
          values <- P
        }
        # avid NA cells of experimental data in predicted and experimental P column
        values <- values[which(!is.na(exp.val.classes.labels))] 
        labels <- exp.val.classes.labels[which(!is.na(exp.val.classes.labels))]
        # avoid NA cells in predicted data in predicted and experimental P column
        labels <- labels[which(!is.na(values))]
        values <- values[which(!is.na(values))]
        # convert P vector in binary data
        values[ which(values[] >= matrix.threshold[j, d]) ] <- 1
        values[ which( values[] < matrix.threshold[j, d]) ] <- 0
        positive <- which(labels==1)
        negative <- which(labels==0)
        # compute contingency matrix 
        # positive
        TP <- length(which(values[positive] == labels[positive])) 
        FN <- length(which(values[positive] != labels[positive]))
        # negative
        TN <- length(which(values[negative] == labels[negative]))
        FP <- length(which(values[negative] != labels[negative]))
        
        tab.temp <- cbind(rep(submission.files.name[j], 2), c("TP", "TN", "FP", "FN"), c(TP, TN, FP, FN))
        tab.hist <- rbind(tab.hist, tab.temp)
        patient.pred <- rbind(patient.pred, rep(sum(sum(TP, TN), sum(FP, FN)), 3))
      }
      tab.hist <- tab.hist[which(is.na(tab.hist[, 1]) != T), ] # remove NA row
      # create row with total patients
      if(variants.only){
        indx.rem <- which(as.numeric(tab.hist[, 3]) < 1)
      } else {
        indx.rem <- which(as.numeric(tab.hist[, 3]) < 3)
      }
      second.ind <- tab.hist[, 3]
      second.ind[which(seq(1:nrow(tab.hist))%%2 == 0)] <-  as.numeric(second.ind[which(seq(1:nrow(tab.hist))%%2 == 0)]) + as.numeric(second.ind[which(seq(2:nrow(tab.hist))%%2 > 0)])
      TP.res <- second.ind[which(tab.hist[,2]=="TN")]
      second.ind[which(tab.hist[,2]=="FP")] <- as.numeric(tab.hist[, 3][which(tab.hist[,2]=="FP")]) + as.numeric(TP.res)
      second.ind[which(tab.hist[,2]=="FN")] <- as.numeric(tab.hist[, 3][which(tab.hist[,2]=="FN")]) + as.numeric(second.ind[which(tab.hist[,2]=="FP")])
      # put in correct order
      TPS <- as.numeric(tab.hist[, 3][which(tab.hist[, 2]=="TP")])+as.numeric(tab.hist[, 3][which(tab.hist[, 2]=="TN")])+as.numeric(tab.hist[, 3][which(tab.hist[, 2]=="FP")])+as.numeric(tab.hist[, 3][which(tab.hist[, 2]=="FN")])
      TNS <- as.numeric(tab.hist[, 3][which(tab.hist[, 2]=="TN")])+as.numeric(tab.hist[, 3][which(tab.hist[, 2]=="FP")])+as.numeric(tab.hist[, 3][which(tab.hist[, 2]=="FN")])
      FPS <- as.numeric(tab.hist[, 3][which(tab.hist[, 2]=="FP")])+as.numeric(tab.hist[, 3][which(tab.hist[, 2]=="FN")])
      FNS <- tab.hist[, 3][which(tab.hist[, 2]=="FN")]
      second.ind[which(tab.hist[,2]=="TP")] <- TPS
      second.ind[which(tab.hist[,2]=="TN")] <- TNS
      second.ind[which(tab.hist[,2]=="FP")] <- FPS
      second.ind[which(tab.hist[,2]=="FN")] <- FNS
      second.ind[indx.rem] <- 0
      tab.hist <- cbind(tab.hist, second.ind)
      # add experimental class
      Pos <- length(which(exp.val.classes.labels==1))
      Neg <- length(which(exp.val.classes.labels==0))
      exp.table.plot <- rbind(cbind("Exp","P", Pos, Pos),
                              cbind("Exp","N", Neg, sum(Pos,Neg)))
      colnames(exp.table.plot) <- NULL
      tab.hist <- rbind(tab.hist, exp.table.plot)
      
      colnames(tab.hist)<- c("Submissions", "Positive", "Patients", "label_pos")
      tab.hist <- data.frame(tab.hist, stringsAsFactors = F)
      tab.hist$Submissions <- as.factor(tab.hist$Submissions)
      tab.hist$Patients <- as.numeric(tab.hist$Patients)
      tab.hist$Positive <- factor(tab.hist$Positive, levels = c("TP", "TN", "FP", "FN", "N", "P"))
      tab.hist$label_pos <- as.numeric(tab.hist$label_pos)
      p <- ggplot(tab.hist, aes(x = Submissions, y = Patients, fill = Positive)) +
           ggtitle(paste("Coverage in ", gsub("\\.", " ", names(exp.val.classes)[d]), sep = "")) + geom_bar(stat = "identity")+theme_bw() + theme(legend.title = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=2), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                                                               plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
                                                                                                                                               #axis.text.y = element_text(size = 14),
                                                                                                                                               axis.title.y = element_text(size = 14, face = "bold"),
                                                                                                                                               #axis.text.x = element_text(size = 14),
                                                                                                                                               axis.title.x = element_text(size = 14, face = "bold"),
                                                                                                                                               axis.text.x = element_text(size = 8, color= group.colors))+
            geom_text(aes(y = label_pos, label = Patients), check_overlap = T, vjust = 1.1, 
                  color="white", size=3.5)+
            scale_fill_manual(values=c("firebrick2", "red4", "royalblue", "royalblue4", "yellow1", "blue"))
      Plot.list[[d]] <- p
    }    
    if(variants.only){
      png("./results/Supplementary_figure1_only_var.png", width = 20, height = 20, units = 'in', res = 300)
    } else {
      png("./results/Supplementary_figure1.png", width = 20, height = 20, units = 'in', res = 300)
    }
    numb.Na <- 7%%3 # fix remainder cells to NA (no plot)
    plot.matrix <- matrix(data = c(seq(1:7), rep(NA, 2)), nrow= 3, ncol = 3, byrow = T)    
    grid.arrange(
      grobs = Plot.list,
      layout_matrix = plot.matrix
    )
    dev.off()    
  }
  MakeSupFig1(variants.only = T)
  MakeSupFig1(variants.only = F)
}

if(TRUE){
  
  LoadSubmission <- function(path.prediction = './data/submission/'){
    # make a list to store predictions
    P.list <- vector(mode = 'list', length = length(submission.files))
    V.list <-  vector(mode = 'list', length = length(submission.files))
    for (i in 1:length(submission.files)) {
      CurrentSubmission <- submission.files[i]
      lines <- scan(file = CurrentSubmission ,what = character(), sep = "\n", quiet = T)
      header <- strsplit(lines[1], "\t")[[1]]
      sub.list <- strsplit(row.corrector(template = temp, sub.lines = lines), '\t')
      col_of_P <- grepl('-P$', header)
      P <- sapply(1:length(sub.list), function(x) {suppressWarnings(as.numeric(sub.list[[x]][col_of_P]))})
      values <- t(P)
      values[which(is.na(values))] <- 0
      pred.class <- values
      for(j in 1:7){
        pred.class[, j][ which(values[, j] >= matrix.threshold[i, j]) ] <- 1
        pred.class[, j][ which( values[, j] < matrix.threshold[i, j]) ] <- 0
      }
      P.list[[i]] <- pred.class
      col_of_V <- grepl('-V$', header)
      V <- sapply(1:length(sub.list), function(x) {sub.list[[x]][col_of_V]})
      V.list[[i]] <- t(V)
    }
    return(list("P" = P.list, "V" = V.list))
  }
  
  Sub.loaded <-  LoadSubmission()
  P.list <- Sub.loaded$P
  V.list <- Sub.loaded$V
  
  VariantComparison <- function(reference.var, predicted.var ){
    if(length(reference.var) < 0 ){
      return(0) 
    } # if no variant is reported in experimental data, return no match
    predicted.var <- unique(unlist(sapply(1:length(predicted.var), function(x){strsplit(predicted.var[x], ",")})))
    num.match <- sum(as.integer(reference.var%in%predicted.var))
    return(num.match)
  }
  
  # Create results matrix
  Index.score <- c('Patient', 'nClass', 'nC', 'nCV', 'Correct groups', 'Groups with correct variant', 'Correct predictions', 'Predictions with correct variant')
  result.matrix <- matrix('empty', nrow = nrow(exp.val.classes), ncol = length(Index.score))#length(submission.files.name))
  colnames(result.matrix) <- Index.score
  Group.names <- sapply(1:length(submission.files.name), function(x){strsplit(submission.files.name[x], "\\.")[[1]][1]})
  
  # compute statistics for each patient
  Hamming.matr <- matrix(data=NA, nrow = nrow(exp.val.classes), ncol=length(submission.files.name) )
  colnames(Hamming.matr) <- submission.files.name
  for(i in 1:nrow(exp.val.classes)){
    Ans_exp <- exp.val.classes[i, ] # correct answer patient i
    # set variables for statistics
    nClass <- 7 - length(which(is.na(Ans_exp))) # number of available experimental classes
    nC <- 0 #number of correct predictions
    nCV <- 0 # number of correct variants
    Corr.group <- c() # name of Group with correct prediction
    Corr.group.V <- c() # name of Group with correct variant
    Corr.pred <- c() # name of Group.prediction with correct prediction
    Corr.pred.V <- c() # name of Group.prediction with correct variant
    for( j in 1:length(P.list)){ # compute statistics for each submission
      pred.name <- submission.files.name[j]
      Group.id <- Group.names[j]
      P <- P.list[[j]][i, ]
      V <- V.list[[j]][i, ]
      if(all(P[which(is.na(Ans_exp)==F)] == Ans_exp[which(is.na(Ans_exp)==F)])){ # patients predicted with all equal disease P are discarded
        Corr.group <- c(Corr.group, Group.id )
        Corr.pred <- c(Corr.pred, pred.name)
        Hamming.matr[i, j] <- 1
      } else {
        Hamming.matr[i, j] <-  round(length(which(P[which(is.na(Ans_exp)==F)] == Ans_exp[which(is.na(Ans_exp)==F)])) / length(which(is.na(Ans_exp)==F)), digits = 1)
      }
      if(VariantComparison(reference.var = unlist(strsplit(paste(Answers.table[i, 10:12], collapse = ","), ",")), predicted.var = unique(V)) > 0){ # check if a submission correctly predicted at least one variant
        Corr.group.V <- c(Corr.group.V, Group.id )
        Corr.pred.V <- c(Corr.pred.V, pred.name)
      }
    }
    Patient.ID <- paste("P", i, sep = "")
    nC <- length(unique(Corr.group))
    nCV <- length(unique(Corr.group.V))
    result.matrix[i,] <- c(i, nClass, nC, nCV, paste(unique(Corr.group),collapse = ','), paste(unique(Corr.group.V), collapse = ','), paste(unique(Corr.pred), collapse = ','), paste(unique(Corr.pred.V), collapse = ","))
  }
  # create another result matrix to be printed
  result.matrix.corr <- result.matrix
  # change col (nCV) to NA when no experimental variant is present (index.pat.WV is a vector of patients ID withou experimental variant)
  result.matrix.corr[p.w.v, grepl("nCV", Index.score)] <- "NA"
  # change col 6 and 8 to NA when no experimental variant is present
  result.matrix.corr[p.w.v, 6] <- "NA"
  result.matrix.corr[p.w.v, 8] <- "NA"
  # convert "" in NA in result.matrix
  result.matrix.corr[which(result.matrix.corr == "")] <- "None"
  # add Hamming matrix
  result.matrix.corr <- cbind(result.matrix.corr, Hamming.matr)
  # print result.matrix
  write.table(result.matrix.corr, './results/general_statistics_by_patient_tab_4.txt', sep = "\t" , quote = F, row.names = F)
  
  # compute statistics on patients with variants
  tab.p.v <- result.matrix.corr[p.v, ]
  mat.p.v <- matrix(data = NA, ncol = 11, nrow= 14)
  colnames(mat.p.v) <- c("Prediction", "nC", "nCV", "nC_CV", seq(1:7))
  mat.p.v[, 1] <- submission.files.name
  for (i in (1:length(submission.files.name))) {
   nC <- length(grep(submission.files.name[i], tab.p.v[, 7], fixed = T))
   nCV <- length(grep(submission.files.name[i], tab.p.v[, 8], fixed = T))
   # compute intersection between nC and NCV: correct prediction and variant
   nCCV <- length(intersect(grep(submission.files.name[i], tab.p.v[, 7], fixed = T), grep(submission.files.name[i], tab.p.v[, 8], fixed = T)))
   # compute nC for class of patients
   nC.vector <- sapply(1:7, function(x){length(grep(submission.files.name[i], tab.p.v[which(tab.p.v[, 2] == as.character(x)), 7], fixed = T))})
   mat.p.v[i, 2:ncol(mat.p.v)] <- c(nC, nCV, nCCV, nC.vector)
  }
  write.table(mat.p.v, './results/general_table_pat_with_v.txt', sep = "\t" , quote = F, row.names = F)
  for (i in (1:length(submission.files.name))) {
    nC <- length(grep(submission.files.name[i], result.matrix.corr[, 7], fixed = T))
    nCV <- length(grep(submission.files.name[i], result.matrix.corr[, 8], fixed = T))
    # compute intersection between nC and NCV: correct prediction and variant
    nCCV <- length(intersect(grep(submission.files.name[i], result.matrix.corr[, 7], fixed = T), grep(submission.files.name[i], result.matrix.corr[, 8], fixed = T)))
    # compute nC for class of patients
    nC.vector <- sapply(1:7, function(x){length(grep(submission.files.name[i], result.matrix.corr[which(result.matrix.corr[, 2] == as.character(x)), 7], fixed = T))})
    mat.p.v[i, 2:ncol(mat.p.v)] <- c(nC, nCV, nCCV, nC.vector)
  }
  write.table(mat.p.v, './results/general_table_all_patients.txt', sep = "\t" , quote = F, row.names = F)
}

############CLASSIFY SUBMISSIONS#############

###################
# Figure, Table 3 #
###################

if(TRUE){
  heat.map.performances[which(is.na(heat.map.performances))] <- 0 # produced in ROC by disease script
  complete.score.matrix <- heat.map.performances
  complete.score.matrix.rank <- data.frame(sapply(1:ncol(complete.score.matrix), function(x) {rank(-complete.score.matrix[, x], ties.method = "average")}))
  rownames(complete.score.matrix.rank) <- c(rownames(complete.score.matrix))
  colnames(complete.score.matrix.rank) <- c(colnames(complete.score.matrix))
  complete.score.matrix.rank <- cbind(complete.score.matrix.rank, rep(NA, nrow(complete.score.matrix.rank)), rep(NA, nrow(complete.score.matrix.rank)))
  colnames(complete.score.matrix.rank) <- c(colnames(complete.score.matrix), "AVG", "Final")
  complete.score.matrix.rank$AVG <- rowMeans(complete.score.matrix.rank, na.rm = TRUE)
  complete.score.matrix.rank$Final <- rank(complete.score.matrix.rank$AVG, ties.method = "average")
  complete.score.matrix.rank.sorted <- complete.score.matrix.rank[order(complete.score.matrix.rank$Final), ]
  complete.score.matrix.rank.sorted$AVG <- round(complete.score.matrix.rank.sorted$AVG, digits = 2)
  complete.score.matrix.rank.sorted$Final <- round(complete.score.matrix.rank.sorted$Final, digits = 0)
  write.table(complete.score.matrix.rank.sorted, file='./results/table_3.txt', sep = "\t" , quote = F, col.names = NA)
  
  complete.score.matrix.sorted <- round(complete.score.matrix[order(complete.score.matrix.rank$Final), ], 2)
  cellcolors<-matrix(NA, ncol(complete.score.matrix.sorted), ncol(complete.score.matrix.sorted))
  threshold = 0.5
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted > threshold] <- color.scale(complete.score.matrix.sorted[complete.score.matrix.sorted > threshold & !is.na(complete.score.matrix.sorted)], cs1 = 0, cs2 = c(0.9, 0.4), cs3 = 0)
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted < threshold] <- color.scale(complete.score.matrix.sorted[complete.score.matrix.sorted < threshold & !is.na(complete.score.matrix.sorted)], cs1 = c(0.4, 0.9), cs2 = 0, cs3 = 0)
  cellcolors[!is.na(complete.score.matrix.sorted) & complete.score.matrix.sorted == threshold] <- "#FFFFFFFF"
  cellcolors[is.na(complete.score.matrix.sorted)] <- "#000000FF"
  
  png('./results/figure_3.png', width = 10, height = 10, units = 'in', res = 300)
  par(mar=c(0,nchar(max(rownames(complete.score.matrix.sorted)))/2,nchar(max(colnames(complete.score.matrix.sorted))),0)+.1)
  color2D.matplot(complete.score.matrix.sorted, show.values = 2, axes = FALSE, cellcolors=cellcolors, vcex = 1.6, xlab="", ylab="")
  axis(3,at=0.5:ncol(complete.score.matrix.sorted),las=2,labels=colnames(complete.score.matrix.sorted))
  axis(2,at=0.5:nrow(complete.score.matrix.sorted),las=2,labels=rev(gsub("_", ".", gsub("Submission_", "Submssion ", rownames(complete.score.matrix.sorted)))))
  ppi = 600
  dev.off()
}
