# Libraries ---------------------------------------------------------------
pkg <- c('dplyr','ggplot2','magrittr','stringr','reshape2','scales','ggthemes','ggrepel','data.table')
inst <- lapply(pkg,library,character.only =T)
library(microbenchmark)
# library(dplyr)
# library(ggplot2)
# library(magrittr)
# library(stringr)
# library(reshape2)
# library(scales)
# library(ggthemes)
# library(ggrepel)
# library(data.table)
# library(topGO)
rm(pkg,inst)
# Load and prepare database ---------------------------------------------------
load('./database/Human.info.Rdata')
load('./database/Mouse.info.Rdata')
load('./database/Contaminants.Rdata')
load('./database/Ortholog.pair.Rdata')
load('./database/Ortholog.hum.Rdata')
load('./database/Ortholog.mus.Rdata')
load('./database/ggtheme_xf.Rdata')
load('./database/hum.allppi.uni.Rdata') #ppi
load('./database/Hum_wikipath_list.Rdata')
load('./database/Mus_wikipath_list.Rdata')
load('./database/Hum_path2Ent_list.Rdata')
load('./database/Hum_Ent2GO_list.Rdata')
load('./database/Hum_Ent2path_list.Rdata')
load('./database/Hum_GO2Ent_list.Rdata')
load('./database/Mus_path2Ent_list.Rdata')
load('./database/Mus_Ent2path_list.Rdata')
load('./database/Mus_GO2Ent_list.Rdata')
load('./database/Mus_Ent2GO_list.Rdata') #GO and pathways.
load('./database/GOterm.Rdata')
load('./database/ReadMScsv.Rdata') #Function
load('./database/ReadCleanMSProRes.Rdata') #Function
load('./database/GenOrthpair.Rdata') #Function

ortholog.pairpro <- ortholog.pair[,c('Accession.hum','Accession.mus')] %>% unique()
orgholog.humpro <- ortholog.hum[,c('Accession.hum'),drop =F] %>% unique()
orgholog.muspro <- ortholog.mus[,c('Accession.mus'),drop =F] %>% unique() #Extract only Accession pairs.

setorder(hum.info,na.last = T)
setorder(mus.info,na.last = T)  
hum.info <- hum.info[!duplicated(Accession.hum) & !is.na(Accession.hum)]
mus.info <- mus.info[!duplicated(Accession.mus) & !is.na(Accession.mus)] # Sort and keep the first record of Accession.

# Define Clean Function (Do not run)---------------------------------------------------
ReadMScsv <- function(x,...){
  dt <- fread(x,...)
  col.name <- str_replace_all(names(dt),c(
    '# '='',
    'Abundances \\(Normalized\\):\\s*'='Abn_',
    'Abundance Ratio: '='Rat_',
    '\\s*/\\s*'='_',
    ' '= '_'))
  names(dt) <- col.name
  index <-  str_detect(col.name, paste(
    c('Accession',
      'Desc',
      'Rat_',
      'Abn_',
      'Unique',
      'Score',
      'PEP_Sequest'),
    collapse = '|'
  ))
  dt <- dt[,index, with = F]
  if ('Description' %in% col.name) {
    Genesym <- str_extract(dt$Description, "(?<=GN=)\\w*\\b") %>% coalesce('noname')
    dt[,Genesym:=Genesym]
  }
  if ('Master_Protein_Accessions' %in% col.name) {
    setnames(dt,'Master_Protein_Accessions','Accession')
  }
  dt
}
ReadCleanMSProRes <- function(x, bloodcontamine = F,commoncontamine=T,...){
  dt <- ReadMScsv(x,...)
  complete <- dt[complete.cases(dt),]
  uni2S10 <- complete[Unique_Peptides > 1 & Score_Sequest_HT >10]
  kera <- str_detect(uni2S10$Description,'^Keratin')
  contamin <- F
  Ig <- F
  HLA <- F
  Hemo <- F
  if(bloodcontamine == T){
    Ig <- str_detect(uni2S10$Description,'^Ig ')
    HLA <- str_detect(uni2S10$Description,'^HLA ')
    Hemo <- str_detect(uni2S10$Description,'^Hemoglobin')
  }
  if(commoncontamine ==T){
    contamin <- uni2S10$Accession %in% contaminants
  }
  final_t <- subset(uni2S10,!(contamin | Ig | HLA | Hemo | kera))
  # final_t <- subset(uni2S10,!(Ig | HLA | Hemo | kera))
  # final_t <- uni2S10[!(Ig | HLA | Hemo | kera)] #a little slower
  filter_out.pro <<- uni2S10$Accession[contamin | Ig | HLA | Hemo | kera]
  n_all <- nrow(table)
  n_uni2S10 <- nrow(uni2S10)
  n_com <- sum(contamin)
  n_kera <- sum(kera)
  n_Ig <- sum(Ig)
  n_HLA <- sum(HLA)
  n_final <- nrow(final_t)
  n_Hemo <- sum(Hemo)
  cat(paste('Total proteins detected:',n_all,'\n',
            '\bObservations with missing value were ommited.','\n',
            '\bProteins with Unique peptides > 1 and Score > 10:',n_uni2S10,'\n',
            '\bContaminant proteins filtered:',n_com,'\n',
            '\bKeratin proteins filtered:',n_kera,'\n',
            '\bIgG proteins filtered:',n_Ig,'\n',
            '\bHLA proteins filtered:',n_HLA,'\n',
            '\bHemoglobin proteins filtered:',n_Hemo,'\n',
            '\bFinal result was returned, protein:',n_final))
  final_t
}
GenOrthpair <- function(humcleanPD,muscleanPD,genclassified = T, report = T){
  humtitle <- deparse(substitute(humcleanPD))
  mustitle <- deparse(substitute(muscleanPD))
  completetable <- humcleanPD %>% left_join(ortholog.pairpro , by = c('Accession'= 'Accession.hum')) %>%
    inner_join(muscleanPD,by = c('Accession.mus' = 'Accession')) %>% filter(!duplicated(.)) %>% as.data.table()
  colname.humpart <- paste(names(humcleanPD), 'hum', sep = '.')
  colname.muspart <- paste(names(muscleanPD), 'mus', sep = '.')
  names(completetable) <-
    c(colname.humpart, colname.muspart) #generate matchtable.
  humcleanPD <-
    humcleanPD %>% mutate(pairstate = 'noorth') %>%  as.data.table()
  humcleanPD[Accession %in% ortholog.pairpro$Accession.hum, pairstate := 'nopair']
  humcleanPD[Accession %in% completetable$Accession.hum, pairstate := 'paired']
  muscleanPD <-
    muscleanPD %>% mutate(pairstate = 'noorth') %>%  as.data.table()
  muscleanPD[Accession %in% ortholog.pairpro$Accession.mus, pairstate := 'nopair']
  muscleanPD[Accession %in% completetable$Accession.mus, pairstate := 'paired']
  if(genclassified == T){
    assign(paste(humtitle, 'classified', sep = '.'),
           humcleanPD,
           envir = parent.frame())
    assign(paste(mustitle, 'classified', sep = '.'),
           muscleanPD,
           envir = parent.frame())
  }

  if (report == T) {
    hum.t0 <- nrow(humcleanPD)
    hum.t1 <- sum(humcleanPD$pairstate == 'paired')
    hum.t2 <- sum(humcleanPD$pairstate == 'nopair')
    hum.t3 <- sum(humcleanPD$pairstate == 'noorth')
    mus.t0 <- nrow(muscleanPD)
    mus.t1 <- sum(muscleanPD$pairstate == 'paired')
    mus.t2 <- sum(muscleanPD$pairstate == 'nopair')
    mus.t3 <- sum(muscleanPD$pairstate == 'noorth')
    n <- nrow(completetable)
    cat(paste(n,'ortholor pairs identified','\n',
              '--------human proteins--------','\n',
              hum.t0,'human input proteins.','\n',
              hum.t1, "mouse ortholog proteins paired.",'\n',
              hum.t2, 'mouse ortholog proteins no paired.','\n',
              hum.t3, 'have no mouse orthlog proteins.','\n',
              '--------mouse proteins--------','\n',
              mus.t0,'mouse input proteins.','\n',
              mus.t1, "human ortholog proteins paired.",'\n',
              mus.t2, 'human ortholog proteins not detected.','\n',
              mus.t3, 'have no human ortholog proteins.','\n'
    ))
  }

  completetable
} #Genereate .matched .nomatch tables.
save(ReadMScsv, file = 'ReadMScsv.Rdata')
save(ReadCleanMSProRes, file = 'ReadCleanMSProRes.Rdata')
save(GenOrthpair, file = 'GenOrthpair.Rdata')

# Define GGPLOT2 theme -(Commented, will not run)-----------------------------------------
# theme_xf <- theme(panel.grid.major = element_blank(),
#                   panel.background = element_blank(),
#                              panel.grid.minor = element_blank(),
#                              panel.border = element_blank(),
#                              plot.title = element_text(size = rel(1.5)),
#                              axis.line.x = element_line(colour='black',size=1),
#                              axis.line.y = element_line(colour='black',size=1),
#                              axis.ticks = element_line(colour = 'black',size =0.8),
#                              axis.ticks.length = unit(0.2,'cm'),
#                              axis.title = element_text(size = rel(1.2)),
#                              axis.text = element_text(size = rel(1.2), colour = 'black'),
#                              legend.title = element_blank(),
#                              legend.text = element_text(size = rel(1.2)),
#                              legend.key = element_blank(),
#                              legend.key.size = unit(0.8,'cm')
# )
# save(theme_xf,file = 'ggtheme_xf.Rdata') 
# Clean and Pair Plaque Data----------------------------------------------
plaque_hum <- ReadCleanMSProRes('AD_Plaque_Human_Lumos_pro.csv',bloodcontamine=F,commoncontamine = T)
plaque_mus <- ReadCleanMSProRes('AD_Plaque_Mouse_Lumos_pro.csv',bloodcontamine=F,commoncontamine = T)
plaque_mus <- plaque_mus[,`:=`(Abn_AD_NPR = (Abn_AD_NPR1 + Abn_AD_NPR2)/2,
                                    Abn_AD_SPR = (Abn_AD_SPR1 + Abn_AD_SPR2)/2,
                                    Rat_AD_SPR_AD_NPR = (Rat_AD_SPR1_AD_NPR1 + Rat_AD_SPR2_AD_NPR2)/2,
                                    Rat_AD_NPR_WT_Contol = (Rat_AD_NPR1_WT_Control + Rat_AD_NPR2_WT_Control)/2)]
plaque_hum.mus_orth <- GenOrthpair(plaque_hum,plaque_mus) 
# Only paired orthologs,and add suffixes '.hum' and '.mus' to column names.

temp1 <- plaque_hum %>% `names<-`(paste0(names(plaque_hum),'.hum'))
temp2 <- plaque_mus %>% `names<-`(paste0(names(plaque_mus),'.mus'))
plaque_hum.mus <- plaque_hum.mus_orth %>% full_join(temp1) %>% full_join(temp2) %>% as.data.table()
# Full ortholog pair information.

write.csv(hum.info[Accession.hum %in% plaque_hum$Accession],file = './figures/hum_info_plaque.csv',quote = F,row.names = F)
write.csv(mus.info[Accession.mus %in% plaque_mus$Accession],file = './figures/mus_info_plaque.csv',quote = F,row.names = F)
# ???

plaque_hum_all.pro <-  ReadMScsv('AD_Plaque_Human_Lumos_pro.csv')
plaque_mus_all.pro <-  ReadMScsv('AD_Plaque_Mouse_Lumos_pro.csv')
# Includes the proteins with low score, contaminations etc.

plaque_hum.mus_all <- GenOrthpair(plaque_hum_all.pro,plaque_mus_all.pro)
plaque_hum.classified[pairstate == 'nopair' & Accession %in% plaque_hum.mus_all$Accession.hum, pairstate := 'filtered']
plaque_mus.classified[pairstate == 'nopair' & Accession %in% plaque_hum.mus_all$Accession.mus, pairstate := 'filtered']
# Add filter information

plaque_hum_fea <- plaque_hum %>% merge(hum.info, by.x ='Accession' , by.y = 'Accession.hum', all.x = T)# Add protein feature.
plaque_hum_fea[,Symbol.hum:= coalesce(Symbol.hum,Genesym)]
plaque_mus_fea <- plaque_mus %>% merge(mus.info, by.x ='Accession' , by.y = 'Accession.mus', all.x = T)# Add protein feature.
plaque_mus_fea[,Symbol.mus:= coalesce(Symbol.mus,Genesym)]
# Simplified gene infomation.

rm(temp1, temp2,
   plaque_hum_all.pro,
   plaque_mus_all.pro,
   plaque_hum_all.pro.classified,
   plaque_mus_all.pro.classified,
   plaque_hum.mus_all)

# Figure 1B Venn Plot (20161230) ---------------------------------------
library(VennDiagram)
tiff(filename = 'Figure1B_venn.tif', width = 6, height = 6, units = 'in', res = 300)
Fig1B_venn <- draw.pairwise.venn(area1 = nrow(plaque_hum),
                                 area2 = nrow(plaque_mus),
                                 cross.area = nrow(plaque_hum.mus_orth), 
                                 category = c('5401\n                  human proteins', '     4550\n mouse proteins'),
                                 col = c('black','black'),
                                 alpha = c(0,0),
                                 cex = rep(1.5,3),
                                 fontfamily = rep('sans',3),
                                 cat.pos = c(-30,30),
                                 cat.cex = rep(1.5,2),
                                 cat.col = c('black','black'),
                                 cat.just = c(list(c(0.8,-0.2),list(0.3,-0.5))),
                                 cat.fontfamily = rep('sans',2),
                                 margin = 0.05)
grid.draw(Fig1B_venn)
dev.off()
rm(Fig1B_venn)
detach('package:VennDiagram', unload = T)
# Further polishment can be done in photoshop.

# Figure 1CD Density Plot ----------------------------------------------------
ProDensPlot <- function(classifiedPD,plottitle){
  abn_matrix <- select(classifiedPD, starts_with(('Abn')))
  ave_abn <-
    transmute(abn_matrix, ave = apply(abn_matrix, 1, mean)) %>% mutate(pairstate = classifiedPD$pairstate) %>% 
    filter(pairstate!='noorth')
  P <- ggplot(ave_abn, aes(x = log10(ave))) +
    geom_line(stat = 'density', aes(color = pairstate), size = 1.2) +
    scale_color_manual(values = c('#924c4c','#afe3b1','#bb8f8f')) +
    labs(
      title = plottitle,
      x = expression('Normalized Protein abundance (log'[10] * ')'),
      y = 'Density'
    ) +
    scale_y_continuous(breaks = seq(0,1,by=0.25))
  P  + theme_xf +
    theme(legend.position = c(.85, .92),
          plot.title = element_text(hjust = 0.5))
}

ProDensPlot(plaque_hum.classified,plottitle = 'Human proteins')
ggsave('Fig1C_Human_Density.tiff',width = 4, height = 3, units = 'in', scale = 2)

ProDensPlot(plaque_mus.classified,plottitle ='Mouse proteins')
ggsave('Fig1D_Mouse_Density.tiff', width = 3, height = 2.5, units = 'in', scale = 2)

rm(ProDensPlot)
# Figure 1E Correlation Plot  ----------------------------------------------
library(corrplot)
abnmatrix <- select(plaque_hum.mus_orth, starts_with('Abn')) %>% log10()
abnmatrix[,6:7] <- abnmatrix[,7:6]
abnmatrix <- abnmatrix[,1:9]
names(abnmatrix) <- c('AD_SPs', 'AD_Brain', 'NonAD_SPs', 'NonAD_Brain', 'AD_Mouse_SPs1', 
                      'AD_Mouse_Brain1', 'AD_Mouse_SPs2','AD_Mouse_Brain2', 'WT_Mouse_Brain')
M <- cor(abnmatrix)
abnmatrix.hum <- select(plaque_hum,starts_with('Abn')) %>% log10()
M.hum <- cor(abnmatrix.hum)
abnmatrix.mus <- select(plaque_mus,starts_with('Abn')) %>% log10()
abnmatrix.mus[,2:3] <- abnmatrix.mus[,3:2]
abnmatrix.mus <- abnmatrix.mus[,1:5]
M.mus <- cor(abnmatrix.mus)
M[1:4,1:4] <- M.hum
M[5:9,5:9] <- M.mus

tiff(filename = 'Fig1E_Plaque_cor.tiff', width = 6, height = 6, units = 'in',res = 300)
corrplot.mixed(M, lower = 'ellipse',upper = 'number',
               tl.pos = 'lt',
               number.digits = 2,
               tl.col ='black'
               # col = colorRampPalette(c("red", "white", "blue"))(100),
               # cl.lim = c(0,1)
)
dev.off()
rm(M.hum,M.mus,M,abnmatrix,abnmatrix.hum,abnmatrix.mus)
# Should be abandoned?
detach('package:corrplot',unload = T)
# Figure 2A Violin Plot. -------------------------------------------------------
plaque_hum.proabn.all <- plaque_hum.classified %>% 
  select(Accession,Description,Abundance = Abn_nonAD_NPR) %>% 
  mutate(Org = 'Human Proteome')
#prepare human table.
plaque_mus.proabn.all <- plaque_mus.classified %>% 
  select(Accession,Description,Abundance = Abn_WT_Control) %>% 
  mutate(Org = 'Mouse Proteome')
#prepare mouse table.
brainproteinabn.hummus <- bind_rows(plaque_hum.proabn.all,plaque_mus.proabn.all)
brainproteinabn.hummus$Org <- as.factor(brainproteinabn.hummus$Org)
vioplot <- ggplot(brainproteinabn.hummus, aes(x=Org, y = Abundance)) + 
  geom_violin(size = 1.2, aes(fill = Org, color = Org)) + 
  geom_boxplot(width = .15, fill = 'white', outlier.color = NULL, outlier.size = 0, size = 1) +
  scale_y_continuous(trans = 'log10',breaks = 10^seq(2:7),labels = comma) +
  scale_color_manual(values = c('#4682B4','#771F1F'), guide = F)+
  scale_fill_manual(values = c('#4682B4','#771F1F'), guide = F) + theme_xf +
  labs(x = NULL, y = 'Normalized Protein Abundance') +
  theme(axis.title.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(0.8)),
        axis.text.y = element_text(size = rel(0.8)))
vioplot
ggsave('Figure 2A (violin plot).tiff', height = 4, width = 5, units = 'in', dpi = 150)
# Should be abandoned.

# Figure 2A Violin Plot Top Part  -------------------------------------------
# plaque_hum.2part <- plaque_hum.classified
# plaque_hum.2part$pairstate[plaque_hum.classified$pairstate != 'matched'] <- 'nomatch'
top3perc.proteins.hum <- filter(plaque_hum.proabn.all,Abundance > quantile(plaque_hum.proabn.all[,3],0.97))
top3perc.proteins.mus <- filter(plaque_mus.proabn.all,Abundance > quantile(plaque_mus.proabn.all[,3],0.97))
top3perc.proteins <- bind_rows(top3perc.proteins.hum,
                               top3perc.proteins.mus)
top10 <- function(vec){
  top <-  vec > sort(vec,decreasing = T)[11]
  top
}

label.proteins <- top3perc.proteins %>% group_by(Org) %>% filter(top10(Abundance))
label.hum <- plaque_hum$Genesym[match(label.proteins$Accession,plaque_hum$Accession)]
label.mus <- plaque_mus$Genesym[match(label.proteins$Accession,plaque_mus$Accession)]
label.proteins$label <- coalesce(label.hum,label.mus)

dotplot <- ggplot(top3perc.proteins, aes(x = Org, y = Abundance)) +
  geom_dotplot(aes(color=Org,fill = Org), stackdir = 'center', binaxis = 'y', dotsize = 0.5, binwidth = 0.03) + 
  geom_text_repel(data = label.proteins,mapping = aes(label=label,color = Org,segment.color = Org),
                  fontface = 'bold',
                  size=3,
                  nudge_x = 0.3,
                  segment.size = 0.5)+
  scale_y_continuous(trans = 'log10', labels = c('100,000',rep('',3),'500,000',rep('',4),'1,000,000',rep('',3)),
                     breaks = seq(100000,1300000,by = 100000)) +
  scale_x_discrete(expand = c(0,0))+
  scale_color_manual(values =c('#4682B4','#771F1F'),guide =F)+
  scale_fill_manual(values =c('#4682B4','#771F1F'),guide=F)+
  labs(x=NULL, y = 'Normalized Protein Abundance') +
  theme_xf
dotplot
ggsave('Fig2B_Violin_Toppart.tiff',height = 4,width = 6,units = 'in',dpi = 300)
rm(plaque_hum.proabn.all,
   plaque_mus.proabn.all,
   brainproteinabn.hummus,
   vioplot,
   top3perc.proteins.hum,
   top3perc.proteins.mus,
   top3perc.proteins,
   label.proteins,
   label.hum,
   label.mus)
# Should be abandoned.

# Figure 2BC GO CC analysis (Part1 prepare)------------------------------------------------
library(topGO)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
cc.lv1 <- goterm2goid[c("cytoplasm", 
                        "plasma membrane", 
                        "nucleus", 
                        "extracellular region")]
cc.lv2 <- goterm2goid[c("cytosol",
                        "mitochondrion",
                        "cytoskeleton",
                        "endoplasmic reticulum",
                        "Golgi apparatus",
                        "endosome",
                        "ribosome",
                        "melanosome")]
cc.lv3 <- goterm2goid[c("synapse",
                        "axon",
                        "neuronal cell body",
                        "dendrite",
                        "growth cone")]
cc.table <- data.frame(id = unname(c(cc.lv1,cc.lv2,cc.lv3)),
                       term = c(names(cc.lv1),names(cc.lv2),names(cc.lv3)),
                       level = c(rep('lv1',length(cc.lv1)),
                                 rep('lv2',length(cc.lv2)),
                                 rep('lv3',length(cc.lv3))),
                       stringsAsFactors = F)
hum.gene2GO <- as.list(org.Hs.egGO) %>% lapply(names)
mus.gene2GO <- as.list(org.Mm.egGO) %>% lapply(names)#Generate go data.
# Figure 2BC GO CC analysis (Part2 Human.Mouse topGO cc generation.) -------
plaque_hum_fea <- merge(plaque_hum.classified,hum_feature.new,by.x = 'Accession',
                        by.y = 'Human.UniProt.SwissProt.Accession',all.x = T) # Add protein feature.
plaque_mus_fea <- merge(plaque_mus.classified,mus_feature.new,by.x = 'Accession',
                        by.y = 'Mouse.UniProt.SwissProt.Accession',all.x = T) # Add protein feature.

hum.mappedkeys <- mappedkeys(org.Hs.egGO) %>% 
  dplyr::intersect(hum_feature.new$Human.EntrezGene.ID) #Extract pros mapped to GO

humgen.plaque <- factor(as.integer(hum.mappedkeys %in% plaque_hum_fea$Human.EntrezGene.ID)) %>%
  `names<-`(hum.mappedkeys) # Generate data for topGO.

hum.GO.CC <- new('topGOdata',ontology = 'CC',
                 allGenes = humgen.plaque,
                 nodeSize = 10,
                 annot = annFUN.gene2GO,
                 gene2GO = hum.gene2GO)

mus.mappedkeys <- mappedkeys(org.Mm.egGO) %>% 
  dplyr::intersect(mus_feature.new$Mouse.EntrezGene.ID) #Extract pros mapped to GO

musgen.plaque <- factor(as.integer(mus.mappedkeys %in% plaque_mus_fea$Mouse.EntrezGene.ID)) %>%
  `names<-`(mus.mappedkeys) # Generate data for topGO.

mus.GO.CC <- new('topGOdata',ontology = 'CC',
                 allGenes = musgen.plaque,
                 nodeSize = 10,
                 annot = annFUN.gene2GO,
                 gene2GO = mus.gene2GO)
# Figure 2BC GO CC analysis (Part3 Table generation and plot) --------------
for(i in 1:nrow(cc.table)){
  hum.genes.termi <- unlist(genesInTerm(hum.GO.CC,whichGO = cc.table$id[i]))
  cc.table$hum.genes[i] <- length(hum.genes.termi)
  cc.table$hum.detected[i] <- length(intersect(plaque_hum_fea$Human.EntrezGene.ID,hum.genes.termi))
  mus.genes.termi <- unlist(genesInTerm(mus.GO.CC,whichGO = cc.table$id[i]))
  cc.table$mus.genes[i] <- length(mus.genes.termi)
  cc.table$mus.detected[i] <- length(intersect(plaque_mus_fea$Mouse.EntrezGene.ID,mus.genes.termi))
}
cc.table.full <- cc.table %>%
  mutate(fac_term = factor(term,levels = rev(term))) %>% mutate(hum.prot.perc = round(100*hum.detected/nrow(plaque_hum_fea),2),
                                                                mus.prot.perc = round(100*mus.detected/nrow(plaque_mus_fea))) %>%
  melt(measure.vars = c('hum.prot.perc','mus.prot.perc'))
cc.table.full$variable <- factor(cc.table.full$variable, levels = c('mus.prot.perc','hum.prot.perc'))
cc.plot <- ggplot(cc.table.full,aes(x=fac_term, y=value))
cc.plot + geom_bar(width = 0.8, aes(fill = variable),
                   position = 'dodge',
                   stat = 'identity') +
  # facet_wrap(~level,nrow = 1,scales  = 'free_x') +
  coord_flip() +
  labs(x=NULL,y=NULL) +
  theme_xf +
  scale_y_continuous(expand = c(0,0.2),breaks = seq(0,80,by =20),limits = c(0,80))+
  theme(axis.text.y = element_text(size = rel(1)),
        axis.text.x = element_text(size = rel(1.2)),
        legend.position = c(0.8,0.1))
ggsave('Fig2B_GOCC_sample.tiff',height = 9,width = 8,units = 'in',dpi = 150)


cc.table.full <- cc.table %>%
  mutate(fac_term = factor(term,levels = rev(term))) %>% mutate(hum.prot.perc = round(100*hum.detected/hum.genes,2),
                                                                mus.prot.perc = round(100*mus.detected/mus.genes)) %>%
  melt(measure.vars = c('hum.prot.perc','mus.prot.perc'))
cc.table.full$variable <- factor(cc.table.full$variable, levels = c('mus.prot.perc','hum.prot.perc'))
cc.plot <- ggplot(cc.table.full,aes(x=fac_term, y=value))
cc.plot + geom_bar(width = 0.8, aes(fill = variable),
                   position = 'dodge',
                   stat = 'identity') +
  # facet_wrap(~level,nrow = 1,scales  = 'free_x') +
  coord_flip() +
  labs(x=NULL,y=NULL) +
  theme_xf +
  scale_y_continuous(expand = c(0,0.2),breaks = seq(0,100,by =20),limits = c(0,90))+
  theme(axis.text.y = element_text(size = rel(1)),
        axis.text.x = element_text(size = rel(1.2)),
        legend.position = c(0.8,0.9))
ggsave('Fig2C_GOCC_term.tiff',height = 9,width = 8,units = 'in',dpi = 150)

# BP GO analysis functios (table and plot) (Updated 1208)-------------------------------------------------
gen.hum.table.GOBP.fis <- function(cleanPD,background=c('all','table'),backgroundtable=NULL){
  library(topGO)
  library(org.Hs.eg.db)
  hum.gene2GO <- as.list(org.Hs.egGO) %>% lapply(names)
  cleanPD <- as.data.table(cleanPD)
  cleanPD.fea <-cleanPD[,.(Accession.hum)] %>%
    merge(hum.info,by = 'Accession.hum',all.x =T) # Add protein feature.
  
  hum.mappedkeys <- mappedkeys(org.Hs.egGO) %>%
    dplyr::intersect(hum.info$Entrez.hum) #Extract pros mapped to GO
  
  if(background == 'all'){
    cleanPD.gens <-
      factor(as.integer(hum.mappedkeys %in% cleanPD.fea$Entrez.hum)) %>%
      `names<-`(hum.mappedkeys) # Generate data for topGO.
  } else {
    backgroundtable <- as.data.table(backgroundtable)
    background.fea <- backgroundtable[,.(Accession.hum)] %>% merge(hum.info, by= 'Accession.hum',all.x =T)
    backgroundkeys <- dplyr::intersect(hum.mappedkeys,background.fea$Entrez.hum)
    cleanPD.gens <- factor(as.integer(backgroundkeys %in% cleanPD.fea$Entrez.hum)) %>% `names<-`(backgroundkeys)
  }
  cleanPD.GO.BP <- new(
    'topGOdata',
    ontology = 'BP',
    allGenes = cleanPD.gens,
    nodeSize = 8,
    annot = annFUN.gene2GO,
    gene2GO = hum.gene2GO
  )
  cleanPD.res<- runTest(cleanPD.GO.BP,
                        algorithm = 'classic',
                        statistic = 'fisher')
  cleanPD.res.table <-
    GenTable(
      cleanPD.GO.BP,
      classicFisher = cleanPD.res,
      orderBy = 'classicFisher',
      ranksOf = 'classicFisher',
      topNodes = 200
    )
  cleanPD.res.table$classicFisher <-
    str_replace(cleanPD.res.table$classicFisher, '<', '') %>% as.numeric()
  cleanPD.res.table <- cleanPD.res.table %>% mutate(EnrichFold = Significant / Expected) %>%
    filter(EnrichFold >= 1.5,Significant >=5,classicFisher<0.01)
  for (i in 1:nrow(cleanPD.res.table)) {
    id.i <- cleanPD.res.table$GO.ID[i]
    allgene.i <- unlist(genesInTerm(cleanPD.GO.BP, whichGO = id.i))
    cleanPD.res.table$keep[i] <- 1
    if (i == 1)
      next()
    for (j in 1:(i - 1)) {
      if (cleanPD.res.table$keep[j] == 0) {
        next()
      } else {
        id.j <- cleanPD.res.table$GO.ID[j]
        allgene.j <- unlist(genesInTerm(cleanPD.GO.BP, whichGO = id.j))
        same.ij <- dplyr::intersect(allgene.j, allgene.i)
        len.i <- length(allgene.i)
        len.j <- length(allgene.j)
        len.ij <- length(same.ij)
        sim.i <- len.ij / len.i
        sim.j <- len.ij / len.j
        if (sim.j < 0.90 & sim.i < 0.90) {
          next()
        } else {
          mask <-
            cleanPD.res.table$EnrichFold[j] < cleanPD.res.table$EnrichFold[i]
          if (mask) {
            if (len.j > len.ij * 2) {
              next()
            }
            cleanPD.res.table$keep[j] <- 0
          } else {
            if (len.i > len.ij * 2) {
              next()
            }
            cleanPD.res.table$keep[i] <- 0
          }
        }
        
      }
    }
  }
  
  cleanPD.res.table <-
    cleanPD.res.table %>% filter(keep == 1, EnrichFold > 1.5) %>% arrange(desc(EnrichFold))
  result <- cleanPD.res.table[1:min(50,nrow(cleanPD.res.table)), ]
  allgene <- rep('0',nrow(result))
  siggene <- rep('0',nrow(result))
  for (i in 1:nrow(result)) {
    all <- unname(unlist(genesInTerm(cleanPD.GO.BP, whichGO = result[i,1])))
    all <- hum.info$Symbol.hum[hum.info$Entrez.hum %in% all]
    sig <- dplyr::intersect(cleanPD.fea$Symbol.hum,all)
    all <- paste(all,collapse = ',')
    sig <- paste(sig,collapse = ',')
    allgene[i] <- all
    siggene[i] <- sig
  }
  result <- result %>% cbind(allgene) %>% cbind(siggene)
  result
}
gen.mus.table.GOBP.fis <- function(cleanPD,background=c('all','table'),backgroundtable=NULL){
  library(topGO)
  library(org.Mm.eg.db)
  mus.gene2GO <- as.list(org.Mm.egGO) %>% lapply(names)
  cleanPD.fea <-
    merge(
      cleanPD,
      mus_feature.new,
      by.x = 'Accession',
      by.y = 'Mouse.UniProt.SwissProt.Accession',
      all.x = T
    ) # Add protein feature.
  
  mus.mappedkeys <- mappedkeys(org.Mm.egGO) %>%
    dplyr::intersect(mus_feature.new$Mouse.EntrezGene.ID) #Extract pros mapped to GO
  if(background == 'all'){
    cleanPD.gens <-
      factor(as.integer(mus.mappedkeys %in% cleanPD.fea$Mouse.EntrezGene.ID)) %>%
      `names<-`(mus.mappedkeys) # Generate data for topGO.
  } else {
    background.fea <-
      merge(
        backgroundtable,
        mus_feature.new,
        by.x = 'Accession',
        by.y = 'Mouse.UniProt.SwissProt.Accession',
        all.x = T
      )
    backgroundkeys <- dplyr::intersect(mus.mappedkeys,background.fea$Mouse.EntrezGene.ID)
    cleanPD.gens <- factor(as.integer(backgroundkeys %in% cleanPD.fea$Mouse.EntrezGene.ID)) %>% `names<-`(backgroundkeys)
  }
  
  cleanPD.GO.BP <- new(
    'topGOdata',
    ontology = 'BP',
    allGenes = cleanPD.gens,
    nodeSize = 8,
    annot = annFUN.gene2GO,
    gene2GO = mus.gene2GO
  )
  cleanPD.res<- runTest(cleanPD.GO.BP,
                        algorithm = 'classic',
                        statistic = 'fisher')
  cleanPD.res.table <-
    GenTable(
      cleanPD.GO.BP,
      classicFisher = cleanPD.res,
      orderBy = 'classicFisher',
      ranksOf = 'classicFisher',
      topNodes = 300
    )
  cleanPD.res.table$classicFisher <-
    str_replace(cleanPD.res.table$classicFisher, '<', '') %>% as.numeric()
  cleanPD.res.table <-  cleanPD.res.table %>% mutate(EnrichFold = Significant / Expected) %>%
    filter(EnrichFold >= 1.5,Significant >=5,classicFisher<0.01)
  for (i in 1:nrow(cleanPD.res.table)) {
    id.i <- cleanPD.res.table$GO.ID[i]
    allgene.i <- unlist(genesInTerm(cleanPD.GO.BP, whichGO = id.i))
    cleanPD.res.table$keep[i] <- 1
    if (i == 1)
      next()
    for (j in 1:(i - 1)) {
      if (cleanPD.res.table$keep[j] == 0) {
        next()
      } else {
        id.j <- cleanPD.res.table$GO.ID[j]
        allgene.j <- unlist(genesInTerm(cleanPD.GO.BP, whichGO = id.j))
        same.ij <- dplyr::intersect(allgene.j, allgene.i)
        len.i <- length(allgene.i)
        len.j <- length(allgene.j)
        len.ij <- length(same.ij)
        sim.i <- len.ij / len.i
        sim.j <- len.ij / len.j
        if (sim.j < 0.90 & sim.i < 0.90) {
          next()
        } else {
          mask <-
            cleanPD.res.table$EnrichFold[j] < cleanPD.res.table$EnrichFold[i]
          if (mask) {
            if (len.j > len.ij * 2) {
              next()
            }
            cleanPD.res.table$keep[j] <- 0
          } else {
            if (len.i > len.ij * 2) {
              next()
            }
            cleanPD.res.table$keep[i] <- 0
          }
        }
        
      }
    }
  }
  cleanPD.res.table <-
    cleanPD.res.table %>% filter(keep == 1, EnrichFold > 1.5) %>% arrange(desc(EnrichFold))
  result <- cleanPD.res.table[1:min(50,nrow(cleanPD.res.table)), ]
  allgene <- rep('0',nrow(result))
  siggene <- rep('0',nrow(result))
  for (i in 1:nrow(result)) {
    all <- unname(unlist(genesInTerm(cleanPD.GO.BP, whichGO = result[i,1])))
    all <- mus_feature.new$Mouse.Associated.Gene.Name[mus_feature.new$Mouse.EntrezGene.ID %in% all]
    sig <- dplyr::intersect(cleanPD.fea$Mouse.Associated.Gene.Name,all)
    all <- paste(all,collapse = ',')
    sig <- paste(sig,collapse = ',')
    allgene[i] <- all
    siggene[i] <- sig
  }
  result <- result %>% cbind(allgene) %>% cbind(siggene)
  result
}##Not Updated.

gen.bp.pointplot <- function(table.fis,label = 'Neural',keywords= c('neu','ves','mem','syn','exo','axon','vac','glutamate','dendr','nerv')){
  wordsdetect <- function(char){
    !all(!str_detect(char,keywords))
  }
  index <- sapply(table.fis$Term, wordsdetect)
  table.fis <- cbind(table.fis,index)
  
  plot <- ggplot(table.fis, aes(x = Significant, y = EnrichFold)) +
    geom_point(aes(fill = index),
               shape = 21,
               color = 'black',
               size = 3.5) +
    scale_fill_manual(values = c('white',muted('green'))) +
    scale_x_continuous(trans = 'log10',breaks = c(50,100,200,400))+
    labs(x = 'Number of genes in term',
         y = 'Enrichfold (Detected/Expected)') +
    expand_limits(y=1.5)+
    theme_xf+
    theme(legend.position = c(0.8,0.9))
  plot
}
gen.bp.pointplot.weak <- function(table.fis,label = 'Neural',keywords= c('neu','ves','mem','syn','exo','axon','vac','glutamate','dendr','nerv')){
  wordsdetect <- function(char){
    !all(!str_detect(char,keywords))
  }
  index <- sapply(table.fis$Term, wordsdetect)
  table.fis <- cbind(table.fis,index)
  
  plot <- ggplot(table.fis, aes(x = Significant, y = EnrichFold)) +
    geom_point(aes(fill = index),
               shape = 21,
               color = 'black',
               size = 3.5) +
    scale_fill_manual(values = c('white',muted('green'))) +
    labs(x = 'Number of genes in term',
         y = 'Enrichfold (Detected/Expected)') +
    expand_limits(y=1.5)+
    theme_xf+
    theme(legend.position = c(0.8,0.9))
  plot
}
# Figure 2D AD brain heatmap--------------------------------------------
library(gplots)
library(RColorBrewer)

GenOrderedTileplot <- function(table,breaks = breaks){
  table <- table %>% mutate(cluster = NULL) %>% 
    melt(id.vars ='order',variable.name='sample', value.name = 'ratio') %>% mutate(ratio = log2(ratio))
  plot <- ggplot(table,aes(sample,order))
  plot + geom_tile(aes(fill = ratio)) +
    scale_fill_gradientn(
      colours = c(
        'forestgreen',
        'forestgreen',
        'white',
        'firebrick4',
        'firebrick4'
      ),
      breaks = c(trunc(min(table$ratio, na.rm = T)*100)/100,
                 -breaks,
                 0,
                 breaks,
                 trunc(max(table$ratio, na.rm = T)*100)/100),
      values = rescale(c(min(table$ratio, na.rm = T),
                         -breaks,
                         0,
                         breaks,
                         max(table$ratio, na.rm = T))),
      labels = c(round(min(table$ratio, na.rm = T), 2),
                 as.character(-breaks),
                 '0',
                 as.character(breaks),
                 round(max(table$ratio, na.rm = T),2))) +
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    labs(x=NULL,y=NULL)+
    theme(legend.position = "top", 
          axis.ticks = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = rel(0.8), 
                                     angle = 330, 
                                     hjust = 0, 
                                     colour = "grey50"))
}
GenClusterMap <- function(matrix,clusternum,breaks=0.5,reorder=F, levels=NULL,name = T){
  library(cluster)
  matrix <- as.data.table(matrix)
  clarax <- clara(log2(matrix[,-1]),clusternum)
  if(name == F){
    matrix <- matrix[!duplicated(matrix[,1])]
    clarax <- clara(log2(matrix),clusternum)
  }
  matrix[,cluster := clarax$clustering]
  if(reorder ==T){
    matrix[,cluster:= factor(cluster,levels = levels)]
  }
  setorder(matrix,cluster)
  assign('clustered.table',matrix,envir = parent.frame())
  matrix <- matrix %>% mutate(order = seq(1,nrow(matrix)),cluster =NULL)
  if(name == T){
    matrix <- matrix[,-1] 
  } 
  GenOrderedTileplot(matrix,breaks = breaks)
}


brain.ratio <- plaque_hum.mus_orth[!duplicated(Accession.hum),.(Accession.hum,
                                 Rat_AD_NPR_nonAD_NPR.hum,
                                 Rat_AD_NPR_WT_Contol.mus)] 

GenClusterMap(brain.ratio,clusternum = 9,reorder = T,levels = c(6,5,4,3,1,2,7,8,9))
clustered.table[,sub:=4]

clustered.table[Rat_AD_NPR_nonAD_NPR.hum > 1.13 & Rat_AD_NPR_WT_Contol.mus > 1.13 & sub == 4, sub := 1]
# clustered.table[Rat_AD_NPR_nonAD_NPR.hum < 1.13 & Rat_AD_NPR_WT_Contol.mus > 1.5 , sub := 2]
clustered.table[Rat_AD_NPR_nonAD_NPR.hum > 1.13 & Rat_AD_NPR_WT_Contol.mus < 1.13 & sub == 4 , sub := 2]
clustered.table[Rat_AD_NPR_nonAD_NPR.hum < 1.13 & Rat_AD_NPR_WT_Contol.mus > 1.13 & sub == 4 , sub := 3]
clustered.table[Rat_AD_NPR_nonAD_NPR.hum < 0.88 & Rat_AD_NPR_WT_Contol.mus > 0.88 & sub == 4 , sub := 5]
clustered.table[Rat_AD_NPR_nonAD_NPR.hum > 0.88 & Rat_AD_NPR_WT_Contol.mus < 0.88 & sub == 4 , sub := 6]
clustered.table[Rat_AD_NPR_nonAD_NPR.hum < 0.88 & Rat_AD_NPR_WT_Contol.mus < 0.88 & sub == 4 , sub := 7]

temp <- list(NULL)
for(i in 1:7){
  sub <- clustered.table[sub == i] %>% as.data.frame()
  sub <- sub[sample(1:nrow(sub),nrow(sub)),]
  temp[[i]] <- sub
} #Randomlization
clustered.table <- rbindlist(temp) # Randomlization

p.table <- clustered.table %>% merge(plaque_hum.mus_orth[,.(Accession.hum,
                                                       Rat_AD_NPR1_WT_Control.mus,
                                                       Rat_AD_NPR2_WT_Control.mus)], by = 'Accession.hum')
setorder(p.table,-sub)
p.table <- p.table[,order := seq(.N)][,.(order,Rat_AD_NPR_nonAD_NPR.hum,Rat_AD_NPR1_WT_Control.mus,Rat_AD_NPR2_WT_Control.mus)]

GenOrderedTileplot(p.table,breaks = 0.5)
ggsave('./figures/Fig2D_ADheatmap_0108.tiff',width = 2,height = 8,units = 'in',dpi = 300)

R1 <- cor(brain.ratio$Rat_AD_NPR_nonAD_NPR.hum,brain.ratio$Rat_AD_NPR_WT_Contol.mus,method = 'spearman')
R2 <- cor(p.table$Rat_AD_NPR1_WT_Control.mus,p.table$Rat_AD_NPR2_WT_Control.mus,method = 'spearman')

rm(sub,clustered.table)
# Figure 2E Prepare AD up/down protein data -----------------------------------------------
# Genchangepro <- function(cleanPD, xcolnum, ycolnum, title, level = 0.05,trans = 'log10'){
#   library(car)
#   f <- paste(trans,'(',names(cleanPD)[ycolnum],')','~',trans,'(',names(cleanPD)[xcolnum],')')
#   linear <- lm(f, data = cleanPD)
#   bon.level <- 1- level/nrow(cleanPD)
#   pred <- predict(linear,cleanPD[,xcolnum,drop = F],interval= 'prediction', level = bon.level)
#   # outlier <- outlierTest(linear,n.max = nmax,cutoff = threshhold)
#   # up.index <- names(outlier$rstudent)[outlier$rstudent > 0] %>% as.numeric()
#   # down.index <- names(outlier$rstudent)[outlier$rstudent < 0] %>% as.numeric()
#   up.index <- log10(cleanPD[,ycolnum]) > pred[,3]
#   down.index <- log10(cleanPD[,ycolnum]) < pred[,2]
#   cleanPD$change <- 'nochange'
#   cleanPD$change[up.index] <- 'up'
#   cleanPD$change[down.index] <- 'down'
#   names(cleanPD)[ncol(cleanPD)] <- paste(title,'changed',sep = '.')
#   changed.table <- cleanPD[up.index|down.index,]
#   cleanPD[,ncol(cleanPD)] <- as.factor(cleanPD[,ncol(cleanPD)])
#   changed.table[,ncol(cleanPD)] <- as.factor(changed.table[,ncol(cleanPD)])
#   assign(paste(title,'changed',level,sep = '.'),changed.table,envir = parent.frame())
#   assign(paste(title,'all',sep = '.'),cleanPD,envir = parent.frame())
# }
Genchangelabelpt <- function(all.data,label.data, xcolname, ycolname,labelcolname,colorcolname, xlab, ylab,labelsize=3.5){
  # label.data <- as.data.table(label.data)
  # label.data[,colorcolname] <- factor(as.character(label.data[,colorcolname]),levels = c('down','nochange','up'))
  plot <-
    ggplot(all.data,
           aes_string(x = xcolname,
                      y = ycolname,
                      color = colorcolname)) +
    geom_point() + 
    geom_abline(slope = 1,
                intercept = 0,
                color = 'blue') +
    labs(x = xlab, 
         y = ylab) +
    geom_label_repel(
      data = label.data,
      aes_string(
        # color = colorcolname,
        fill = colorcolname,
        label = labelcolname
      ),
      fontface = 'bold',
      size = labelsize,
      label.size = 0,
      color = 'white',
      segment.color = 'black',
      box.padding = unit(0.05, 'lines'),
      nudge_x = ifelse(as.logical(label.data[,ncol(label.data),with = F]== 'up'), -0.4, 0.4)
    ) +
    theme(axis.title = element_text(size = rel(1.5))) +
    scale_x_continuous(
      trans = 'log10',
      breaks = c(10,10 ^ seq(2,6)),
      labels = comma,
      expand = c(0, 0),
      limits = c(10,1500000)
    ) +
    scale_y_continuous(
      trans = 'log10',
      breaks = c(10,10 ^ seq(2,6)),
      labels = comma,
      expand = c(0, 0),
      limits = c(10,1500000)
    ) +
    scale_color_manual(values = c('green4', 'black', 'red'), guide = F) +
    scale_fill_manual(values = c('forestgreen', 'firebrick4'),
                      guide = F) +
    theme_xf
  plot
}
Genstatplot <- function(value){
  table <- data.frame(value = value)
  plot <- ggplot(table, aes(x = value)) +
    geom_line(stat = 'density', size = 1.2)+
    scale_x_continuous(breaks = c(0.5,1.0,1.5),limits = c(0.5,1.5))+
    labs(x=NULL,y=NULL)+
    theme_xf+
    theme(axis.text=element_text(size=rel(1.8)))
  plot
}
f <- 'log10 (Abn_AD_NPR) ~ log10(Abn_nonAD_NPR)'
rep.linear <- lm(f, data = plaque_hum_fea)
RepLmTest <- function(cleanPD,lm = rep.linear,xcolname,ycolname,level=0.05,title){
  temp <- cleanPD[,xcolname,with = F] %>% setnames(xcolname,'Abn_nonAD_NPR')
  bon.level = 1- level/nrow(cleanPD)
  pred <- predict(rep.linear,temp,interval= 'prediction',level = bon.level)
  act <- log10(cleanPD[,ycolname,with = F])
  down <- act < pred[,2]
  up <- act > pred[,3]
  result <- rep('nochange',nrow(cleanPD))
  result[up] <- 'up'
  result[down] <- 'down'
  cleanPD$result <- result
  changed.table <- cleanPD[result!='nochange']
  setnames(cleanPD,'result',paste0(title,'.changed'))
  setnames(changed.table,'result',paste0(title,'.changed'))
  assign(paste0(title,'.all'),cleanPD,envir = parent.frame())
  assign(paste0(title,'.changed.',level),changed.table,envir = parent.frame())
}
# Figure 2E human ad change proteins. -------------------------------------
RepLmTest(cleanPD = plaque_hum_fea,xcolname = "Abn_nonAD_NPR",ycolname = "Abn_AD_NPR",level = 0.05,title = 'humad')
humad.plot <- Genchangelabelpt(all.data = humad.all,
                               label.data = humad.changed.0.05,
                               xcolname = 'Abn_nonAD_NPR',
                               ycolname = 'Abn_AD_NPR',
                               labelcolname = 'Symbol.hum',
                               colorcolname = 'humad.changed',
                               xlab = 'Abundance in nonAD NPR',
                               ylab = 'Abundance in AD NPR',
                               labelsize = 3.5)
humad.plot
ggsave('./figures/Figure2E_hum_brain_point_0108.tiff', width = 7, height = 6, units = 'in',dpi = 150)
Genstatplot(plaque_hum_fea$`Rat_AD_NPR/nonAD_NPR`)
ggsave('./figures/Fig2E_hum_brain_ratiodensity.tiff', width = 3,height = 3,units = 'in',dpi = 150)
rm(humad.plot)
# library(psych)
# describe(plaque_hum_fea$`Rat_AD_NPR/nonAD_NPR`)
# Figure 2F Mouse AD Proteins ---------------------------------------------
plaque_mus_fea[,Abn_AD_range := abs(Abn_AD_NPR1 - Abn_AD_NPR)/2]
RepLmTest(cleanPD = plaque_mus_fea,xcolname = 'Abn_WT_Control',ycolname = 'Abn_AD_NPR',level = 0.05,title = 'musad')
musad.plot <- Genchangelabelpt(all.data = musad.all,
                               label.data = musad.changed.0.05,
                               xcolname = 'Abn_WT_Control',
                               ycolname= 'Abn_AD_NPR',
                               labelcolname = "Symbol.mus",
                               colorcolname = 'musad.changed',
                               xlab = 'Abundance in WT mouse brain',
                               ylab = 'Abundance in APP/PS1 mouse NPR',
                               labelsize = 3.5)
musad.plot + geom_errorbar(data = musad.changed.0.05,aes(ymax = (Abn_AD_NPR+Abn_AD_range),
                                                         ymin = (Abn_AD_NPR- Abn_AD_range),
                                                         color = musad.changed,
                                                         width = 0.05))
ggsave('./figures/Figure2F_mus_brain_point_0108.tiff', width = 7, height = 6, units = 'in',dpi = 150)
rm(musad.plot)
Genstatplot(plaque_mus_fea$`Rat_AD_NPR/WT_Control`)
ggsave('Fig2F_mus_brain_ratiodensity.tiff', width = 3,height = 3,units = 'in',dpi = 150)
describe(plaque_mus_fea$`Rat_AD_NPR/WT_Control`)
# Figure 3A humadplaque Proteins----------------------------------------
RepLmTest(cleanPD = plaque_hum_fea,xcolname = 'Abn_AD_NPR',ycolname = 'Abn_AD_SPR',title = 'humadplaque')
humadplaque.plot <- Genchangelabelpt(all.data = humadplaque.all,
                                     label.data = humadplaque.changed.0.05,
                                     xcolname = 'Abn_AD_NPR',
                                     ycolname = 'Abn_AD_SPR',
                                     labelcolname = "Symbol.hum",
                                     colorcolname = "humadplaque.changed",
                                     xlab = 'Abundance in AD NPR',
                                     ylab = 'Abundance in AD SPR',
                                     labelsize = 3.5)
humadplaque.plot
ggsave('./figures/Figure3B_humadplaque_point_0108.tiff', width = 7, height = 6, units = 'in',dpi = 150)
rm(humadplaque.plot)
Genstatplot(plaque_hum_fea$`Rat_AD_SPR/AD_NPR`)
ggsave('Fig3B_humadplaque_ratiodensity.tiff', width = 3,height = 3,units = 'in',dpi = 150)
describe(plaque_hum_fea$`Rat_AD_SPR/AD_NPR`)
# Figure 3B human normal plaque -------------------------------------------
RepLmTest(cleanPD = plaque_hum_fea,xcolname = 'Abn_nonAD_NPR',ycolname = 'Abn_nonAD_SPR',title = 'humnormplaque')
humnormplaque.plot <- Genchangelabelpt(all.data = humnormplaque.all,
                                       label.data = humnormplaque.changed.0.05,
                                       xcolname = 'Abn_nonAD_NPR',
                                       ycolname = 'Abn_nonAD_SPR',
                                       labelcolname = "Symbol.hum",
                                       colorcolname = "humnormplaque.changed",
                                       xlab = 'Abundance in Non-AD NPR',
                                       ylab = 'Abundance in Non-AD SPR',
                                       labelsize = 3.5)
humnormplaque.plot
ggsave('./figures/Figure3C_humnormalplaque_point.tiff', width = 7, height = 6, units = 'in',dpi = 150)
rm(humnormplaque.plot)
Genstatplot(plaque_hum_fea$`Rat_nonAD_SPR/nonAD_NPR`)
ggsave('Fig3C_humnormalplaque_ratiodensity.tiff', width = 3,height = 3,units = 'in',dpi = 150)
describe(plaque_hum_fea$`Rat_nonAD_SPR/nonAD_NPR`)
# Figure 3C Mus AD plaque1-------------------------------------------------------
RepLmTest(cleanPD = plaque_mus_fea,xcolname = "Abn_AD_NPR1",ycolname = "Abn_AD_SPR1",title = 'musadplaque1')
musadplaque1.label.0.05 <- musadplaque1.changed.0.05[order(abs(log2(Rat_AD_SPR1_AD_NPR1)),decreasing = T)][1:40,]

musadplaque1.plot <- Genchangelabelpt(all.data = musadplaque1.all,
                                      label.data = musadplaque1.label.0.05,
                                      xcolname = "Abn_AD_NPR1",
                                      ycolname = "Abn_AD_SPR1",
                                      labelcolname = "Symbol.mus" ,
                                      colorcolname = "musadplaque1.changed" ,
                                      xlab = 'Abundance in APP/PS1 mouse NPR (1)',
                                      ylab = 'Abundance in APP/PS1 mosue SPR (1)')
musadplaque1.plot
ggsave('./figures/Figure3D_musadplaque1_point.tiff', width = 7, height = 6, units = 'in',dpi = 150)
rm(musadplaque1.plot)
Genstatplot(plaque_mus_fea$Rat_AD_SPR1_AD_NPR1)
ggsave('Fig3D_musadplaque1_ratiodensity.tiff', width = 3,height = 3,units = 'in',dpi = 150)
describe(plaque_mus_fea$Rat_AD_SPR1_AD_NPR1)
# Figure S3A Mus AD plaque2 -------------------------------------------------
RepLmTest(cleanPD = plaque_mus_fea,xcolname = "Abn_AD_NPR2",ycolname = "Abn_AD_SPR2",title = 'musadplaque2')
musadplaque2.changed.0.05$musadplaque2.changed <- factor(musadplaque2.changed.0.05$musadplaque2.changed,levels = c('down','up'))
musadplaque2.label.0.05 <- musadplaque2.changed.0.05[order(abs(log2(Rat_AD_SPR2_AD_NPR2)),decreasing = T)][1:40,]


musadplaque2.plot <- Genchangelabelpt(all.data = musadplaque2.all,
                                      label.data = musadplaque2.label.0.05,
                                      xcolname = "Abn_AD_NPR2",
                                      ycolname = "Abn_AD_SPR2",
                                      labelcolname = "Symbol.mus" ,
                                      colorcolname = "musadplaque2.changed" ,
                                      xlab = 'Abundance in APP/PS1 mouse NPR (2)',
                                      ylab = 'Abundance in APP/PS1 mosue SPR (2)')
musadplaque2.plot
ggsave('./figures/FigureS3A_musadplaque2_point.tiff', width = 7, height = 6, units = 'in',dpi = 150)
Genstatplot(plaque_mus_fea$Rat_AD_SPR2_AD_NPR2)
ggsave('FigS3A_musadplaque2_ratiodensity.tiff', width = 3,height = 3,units = 'in',dpi = 150)
describe(plaque_mus_fea$Rat_AD_SPR2_AD_NPR2)
# Figure 3D-E Generate AD changed protein table  ------------------------------------------
sprelated.hum <- union(humnormplaque.changed.0.05$Accession,humadplaque.changed.0.05$Accession)
adrelated.hum <- union(sprelated.hum,humad.changed.0.05$Accession)
sprelated.mus<- intersect(musadplaque2.changed.0.05$Accession,musadplaque1.changed.0.05$Accession) 
adrelated.mus <- union(sprelated.mus,musad.changed.0.05$Accession)
sprelated.mus.labeled <- union(musadplaque2.label.0.05$Accession,musadplaque1.label.0.05$Accession)

sp.core <- plaque_hum.mus_orth[Accession.hum %in% sprelated.hum & Accession.mus %in% sprelated.mus, Accession.hum]

sprelated.orth.hum <- ortholog.pairpro[Accession.mus %in% sprelated.mus]$Accession.hum %>% union(sprelated.hum)
adrelated.orth.hum <- ortholog.pairpro[Accession.mus %in% adrelated.mus]$Accession.hum %>% union(adrelated.hum)

# Figure 3D-E ADppi CPDB --------------------------------------------------------------
hum.sp.ppi <- hum.allppi.uni[Accession.hum.x %in% sprelated.hum & Accession.hum.y %in% sprelated.hum]
hum.sp.ppi.nodeinfo <- plaque_hum[Accession %in% sprelated.hum]
write.csv(hum.sp.ppi,file = 'hum.sp.ppi.csv',quote = F, row.names = F)
write.csv(hum.sp.ppi.nodeinfo,file = 'hum.sp.ppi.nodeinfo.csv',quote = F, row.names = F)



spppi <- hum.allppi.uni[Accession.hum.x %in% sprelated.orth.hum & Accession.hum.y %in% sprelated.orth.hum]

DegreeSum <- function(x,table){
  tb <- as.data.frame(table)
  res <- 0
  for(i in seq(ncol(tb))){
    res <- res + sum(tb[,i] == x)
  }
  res
}
spppi.node.degree <- sapply(sprelated.orth.hum,DegreeSum,table = spppi)
DegreeCat <- function(x){
  if(x > 30){
    res <- '30'
  } else if(x > 10){
    res <- '15'
  } else if(x > 5){
    res <- '5'
  } else {
    res <- '0'
  }
}
spppi.node.degree <-  data.table(Accession.hum = sprelated.orth.hum, Degreecat = sapply(spppi.node.degree,DegreeCat))


spppi.node.info <- data.table(Accession.hum = sprelated.orth.hum) %>% 
  merge(ortholog.pairpro, by = 'Accession.hum',all.x = T) %>% 
  # merge(humsp.all[,.(Accession,humsp.changed)],by.x = 'Accession.hum', by.y = 'Accession',all.x = T) %>% 
  merge(humadplaque.all[,.(Accession,humadplaque.changed)],by.x = 'Accession.hum', by.y = 'Accession',all.x = T) %>% 
  merge(humnormplaque.all[,.(Accession,humnormplaque.changed)],by.x = 'Accession.hum', by.y = 'Accession',all.x = T) %>% 
  # merge(musad.all[,.(Accession,musad.changed)],by.x = 'Accession.mus', by.y = 'Accession',all.x = T) %>% 
  merge(musadplaque1.all[,.(Accession,musadplaque1.changed)],by.x = 'Accession.mus', by.y = 'Accession',all.x = T) %>% 
  merge(musadplaque2.all[,.(Accession,musadplaque2.changed)],by.x = 'Accession.mus', by.y = 'Accession',all.x = T) %>% 
  merge(plaque_hum[,.(Accession,
                      Rat_AD_NPR_nonAD_NPR,
                      Rat_AD_SPR_AD_NPR,
                      Rat_nonAD_SPR_nonAD_NPR)],by.x = 'Accession.hum', by.y = 'Accession',all.x = T) %>% 
  merge(plaque_mus[,.(Accession,
                      Rat_AD_NPR_WT_Contol,
                      Rat_AD_SPR_AD_NPR)],by.x = 'Accession.mus', by.y = 'Accession',all.x = T,
        suffixes = c('.hum','.mus')) %>% 
  merge(hum.info[,.(Accession.hum,Symbol.hum)],by = 'Accession.hum',all.x = T) %>% 
  merge(mus.info[,.(Accession.mus,Symbol.mus)],by = 'Accession.mus',all.x = T) %>% 
  merge(spppi.node.degree, by = 'Accession.hum',all.x = T)
spppi.node.info <- spppi.node.info[!(is.na(Rat_AD_SPR_AD_NPR.hum) & is.na(Rat_AD_SPR_AD_NPR.mus))]

Na21Log2 <- function(vec){
  vec[is.na(vec)] <- 1
  vec <- log2(vec)
  vec
}

temp <- spppi.node.info[,lapply(.SD,Na21Log2),.SDcols =c('Rat_AD_NPR_nonAD_NPR',
                                                         'Rat_AD_SPR_AD_NPR.hum',
                                                         'Rat_nonAD_SPR_nonAD_NPR',
                                                         'Rat_AD_NPR_WT_Contol',
                                                         'Rat_AD_SPR_AD_NPR.mus')]
spppi.node.info[,c('Rat_AD_NPR_nonAD_NPR',
                   'Rat_AD_SPR_AD_NPR.hum',
                   'Rat_nonAD_SPR_nonAD_NPR',
                   'Rat_AD_NPR_WT_Contol',
                   'Rat_AD_SPR_AD_NPR.mus'):=as.list(temp)]

spppi.hum <- hum.allppi.uni[Accession.hum.x %in% sprelated.hum & Accession.hum.y %in% sprelated.hum] #Human only.
spppi.string <- fread('adrelated_from_string.tsv')
names(spppi.string)[1] <- 'node1'
spppi.hum.string <- spppi.string %>% merge(hum.info, by.x = 'node1',by.y = 'Symbol.hum') %>% 
  merge(hum.info,by.x = 'node2',by.y = 'Symbol.hum')
  
spppi.hum.string <-  spppi.hum.string[Accession.hum.x %in% sprelated.hum & Accession.hum.y %in% sprelated.hum,
                                      .(Accession.hum.x,Accession.hum.y)]
spppi.hum <- rbindlist(list(spppi.hum[,.(Accession.hum.x,Accession.hum.y)],spppi.hum.string)) %>% unique()


write.csv(spppi.hum,file = 'spppi.hum.csv',quote = F, row.names = F)
write.csv(spppi,file = 'spppi.csv',quote = F, row.names = F)
write.csv(spppi.node.info,file = 'spppi.node.csv',quote = F,row.names = F)

# Figure 3D-EMouse PPI complex function analysis -------------------------------------


goterm.vec <-  goterm$Name %>% `names<-`(goterm$ID)
MusSymComFuc <- function(vec){
    ent <- mus.info$Entrez.mus[mus.info$Symbol.mus %in% vec] %>% as.character()
    names(vec) <- ent
    gos <- mus.ent2go.list[ent] %>% lapply(unique) %>% unlist() %>% table() %>% as.data.table()
    names(gos) <- c('ID','N')
    co.go <- gos[N>1] %>% mutate(cat= 'GO') %>% as.data.table()
    go.syms <- sapply(co.go$ID,function(go){paste(vec[intersect(mus.go2ent[[go]], ent)], collapse = ',')})
    co.go <- merge(co.go,goterm, by= 'ID', all.x = T)[,.(Name,N,cat)]
    paths <- mus.ent2path.list[ent] %>% unlist() %>% table() %>% as.data.table()
    names(paths) <- c('Name','N')
    co.path <- paths[N>1]%>% mutate(cat= 'Path') %>% as.data.table()
    path.syms <- sapply(co.path$Name,function(path){paste(vec[intersect(mus.allpath.list[[path]],ent)], collapse = ',')})
    sym <- c(go.syms,path.syms)
    res <- rbindlist(list(co.go,co.path)) %>% mutate(Genes = sym) %>% as.data.table()
    res
}

vec <- c('Sfpq','Snrpg','Rpl21','U2af1','Psip1','Hnrnpa0','Rbmx','Aimp2','Ncstn','Ephx1')
vec <- c('Syt5','Syt13','Syt4','Nrxn','Sdcbp','Rundc3a')
vec <- c('Was','Wasf2','Hcls1','Abi3','Lcp1','Aif1','Cnn3','Myo1e','Actr1a','Rac2','Ncf1','Ptprc','Ptpn6','Vav1','Clcn7,Tmem63b','Ostm1')
vec <- c('Stx7','Stx8','Stx12','Vamp7','Ehd4','Bloc1s1','Ptn','Snapin','Vps41','Vps33a','Vps16','Vps39','Atp9a')
vec <- c('Eef1b','Igsf8','Eef1g','Flnc')
vec <- c('Apoe','Igtax','Gfap','Col1a2','Ctsb','Flt1','Col1a1',
         'Vtn','Anxa2','Vim','Icam1','Gpc1','Itgb5','Sdc4',
         'Itgam','Htra1','Gpc4','Vwf','Edil3','Mfge8','Slit2',
         'Slit1','Plek','Msn','Mdk','Tgm1','Ctgf','Nov','Apbb1ip','Gpc5','Clstn3','Aplp1','Aplp2','Spon1')
vec <- c('Rufy1','Rufy2','Anxa3','Anxa4','Capg','Mif','Spock2','Spock3','Csf1r','Inpp5d','Fcer1g','Hk2','Hexb','Mlf2')
vec <- c('Tpp1','Scarb2','Asah1','Arl8b','Rraga','Rragc','Lamtor3','Lamtor1','Lamtor2','Lamp1','Lamp2','Ncstn','Lamp2','Clcn7','Ostm1')
vec <- c('Itm2b','Necap2','Dus3l','Serpine2','Pi4k2a','Fermt3','Tbc1d9b','Tmem55b','Blmh','Clstn1','Clu','Asph')
vec <- c('Tom1','Tollip','Atg9a','Vps29','Abcb6','Psap','Itm2c','Bace1')
test <- MusSymComFuc(vec)


# Figure 3F Merge human and mouse string data.--------------------------------------
hum.interaction <- fread('./original_data//LCM_hum_interaction.csv') %>% `names<-`(c('node1','node2'))
mus.interaction <- fread('./original_data/LCM_mus_interaction.csv') %>% `names<-`(c('node1','node2'))
mus.interaction[,`:=`(node1 = toupper(node1),node2 = toupper(node2))]

hum.interaction <- merge(hum.interaction,hum.info,by.x = 'node1',by.y = 'Symbol.hum') %>% 
  merge(hum.info,by.x ='node2',by.y = 'Symbol.hum')
hum.interaction <- hum.interaction[,.(Accession.hum.x,Accession.hum.y)]

mus.interaction <- merge(mus.interaction,mus.info,by.x = 'node1',by.y = 'Symbol.mus') %>% 
  merge(mus.info,by.x ='node2',by.y = 'Symbol.mus')
mus.interaction <- merge(mus.interaction,ortholog.pairpro,by.x = 'Accession.mus.x',by.y = 'Accession.mus') %>% 
  merge(ortholog.pairpro,by.x='Accession.mus.y',by.y='Accession.mus')
mus.interaction <- mus.interaction[,.(Accession.hum.x,Accession.hum.y)]


UniDirec <- function(table){
  temp <- table %>% as.data.table() %>% unique() %>% t() %>% as.data.table()
  temp <- lapply(temp,function(vec)as.data.table(matrix(sort(vec),nrow = 1)))
  res <- rbindlist(temp) %>% unique()
  res
}

string.interaction <- rbind(hum.interaction,mus.interaction) %>% UniDirec() %>% 
  `names<-`(c('Accession.hum.x','Accession.hum.y'))
string.interaction <- string.interaction[Accession.hum.x %in% sp.core & Accession.hum.y %in% sp.core]
commonpath.interaction <- hum.allppi.uni[Accession.hum.x %in% sp.core & Accession.hum.y %in% sp.core,.(Accession.hum.x,Accession.hum.y)]
core.interaction <- rbind(commonpath.interaction,string.interaction) %>% UniDirec()
core.info <- plaque_hum.mus_orth[Accession.hum %in% sp.core][,Ratio:= log2((Rat_AD_SPR_AD_NPR.hum+Rat_nonAD_SPR_nonAD_NPR.hum+Rat_AD_SPR_AD_NPR.mus)/3)][,.(Accession.hum,Ratio)] %>% merge(hum.info,by = 'Accession.hum')

write.csv(core.info,file = './figures/Core_info.csv',quote = F,row.names = F)
write.csv(core.interaction,file = './figures/Core_interaction',quote = F,row.names = F)
rm(core.interaction,core.info,string.interaction,commonpath.interaction,hum.interaction,mus.interaction)
# Figure 3G ratio plot (follow figure 3F.s)----------------------------------------------
adrelated.table <- plaque_hum.mus[Accession.hum %in% adrelated.hum]
TilePlot <- function(table){
  
  plot.t <- melt(
    as.data.frame(table),
    id.vars = 1,
    value.name = 'Ratio',
    variable.name = 'Site') %>% as.data.table()
plot.t[,Ratio:=log2(Ratio)]
plot.t$Genesym.hum <- factor(plot.t$Genesym.hum, levels = rev(table$Genesym.hum))

plot <- ggplot(plot.t,aes(Site,Genesym.hum))
plot + geom_tile(aes(fill = Ratio)) +
  scale_fill_gradientn(
    colours = c(
      'forestgreen',
      'forestgreen',
      'white',
      'firebrick4',
      'firebrick4'),
    breaks = c(trunc(min(plot.t$Ratio, na.rm = T)*100)/100,
               -1,
               0,
               1,
               trunc(max(plot.t$Ratio, na.rm = T)*100)/100),
    values = scales::rescale(c(min(plot.t$Ratio, na.rm = T),
                       -1,
                       0,
                       1,
                       max(plot.t$Ratio, na.rm = T))),
    labels = c(trunc(min(plot.t$Ratio, na.rm = T)*100)/100,
               '-1',
               '0',
               '1',
               trunc(max(plot.t$Ratio, na.rm = T)*100)/100)) +
  # scale_color_manual(values = c(up = 'firebrick',down = 'forestgreen',nc=NA),na.value = NA,guide =F)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  labs(x=NULL,y=NULL)+
  theme(legend.position = "top", 
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = rel(0.8)),
        axis.text.x = element_text(size = rel(0.8), 
                                   angle = 330, 
                                   hjust = 0, 
                                   colour = "grey50"))
}

table.hum.index <- c("MDK","SMOC1","PTN","NTN1","TMEFF2","OLFML3","C4A","SLIT1","SDC4","SLIT2","FLT1","SPOCK2","SPOCK1","HTRA1","GPC1","SPON1","ITGB2","APP","COL1A1","COL1A2","GRN","GFAP","CLCN6","GPC5","ICAM1","APOE",'TMEFF1',"C4B","BLMH","RAN","CSRP1","TAGLN","ACSF3",'DSP',"GSTM5","ALDH2","DDX3Y","GPNMB","CTHRC1","CD163",  "CPS1", "CYBA", "SCIN", "HLA",  "AQP1", "CD44", "BPIFB1","SERPINB3","SPRR1B","AZGP1", "PRPH", "ECM1",'CDSN', "C4BPA","COL25A1","IGHV3","HRNR", "C1S",  "SNRPE","CHI3L1","DEFA1",'IGHG4','HBG1',"ALDH1A2","MT1G","TOP3B","FAM107A",  "CLDN10","PDLIM3","GSTM1","DCD","CD5L", "S100A7","MYH2", "CALML3")
adrelated.table <- adrelated.table[match(table.hum.index,adrelated.table$Genesym.hum),.(Genesym.hum,
                                                                                        Rat_AD_NPR_nonAD_NPR.hum,
                                                                                        Rat_AD_SPR_AD_NPR.hum,
                                                                                        Rat_nonAD_SPR_nonAD_NPR.hum,
                                                                                        Rat_AD_NPR_WT_Contol.mus,
                                                                                        Rat_AD_SPR1_AD_NPR1.mus,
                                                                                        Rat_AD_SPR2_AD_NPR2.mus)]
TilePlot(adrelated.table)
ggsave('./figures/Fig3E_HUMAN_adrelated_protein_ratio.tiff', width = 5, height = 15, units = 'in',dpi = 300)

temp <- plaque_hum.classified %>% merge(hum.info,by.x ='Accession',by.y='Accession.hum')
temp <- temp[Symbol.hum %in% table.hum.index,.(Symbol.hum,pairstate)]
temp[,Symbol.hum:=factor(Symbol.hum,levels = table.hum.index)]
setorder(temp,Symbol.hum)
rm(adrelated.table)
# Figure 3H analysis (Human)----------------------------------------------------------
plaque_hum.pep <- ReadMScsv('AD_Plaque_Human_Lumos_pep.csv') %>% 
  merge(hum.info,by.x = 'Accession',by.y = 'Accession.hum',all.x = T)
plaque_hum.pep <- plaque_hum.pep[Percolator_PEP_Sequest_HT < 0.01 & complete.cases(plaque_hum.pep)]

plaque_mus.pep <- ReadMScsv('AD_Plaque_Mouse_Lumos_pep.csv') %>% 
  merge(mus.info,by.x = 'Accession',by.y = 'Accession.mus',all.x = T)
plaque_mus.pep <- plaque_mus.pep[Percolator_PEP_Sequest_HT < 0.01 & complete.cases(plaque_mus.pep)]

temp <- plaque_hum.pep[Symbol.hum %in% c('APP','GAPDH','MAPT')]
plot <- ggplot(temp,aes(x=Symbol.hum, y = `Rat_AD_SPR_AD_NPR`))
plot + geom_jitter(width = 0.3,size=2) +
  geom_boxplot(outlier.size = 0,
               outlier.color = 'white',
               alpha = 0,
               # fill =NULL,
               size = 0.8,
               width = 0.45,
               color = 'forestgreen')+
  scale_y_continuous(trans = 'log2',breaks = c(0.5,1,2,4,8),limits = c(0.5,10))+
  theme_xf +
  scale_x_discrete(labels = c('APP','GAPDH','MAPT(tau)'))+
  labs(y = 'Enrichment Fold in AD SPR',
       x = NULL)
ggsave('./figures/Fig3H_peppoint.tiff',width = 5,height = 4,units = 'in',dpi = 150)
# Supplimentary table5, preparation prepare human peptides table.(Human)--------------------------------------
load('./database/ProteinSeq.Hum.Rdata')
load('./database/ProteinSeq.Mus.Rdata')

plaque_hum.pep.seqinfo <- plaque_hum.pep[Accession %in% plaque_hum$Accession,.(Abn_AD_SPR=sum(Abn_AD_SPR),
                                                                        Abn_AD_NPR=sum(Abn_AD_NPR),
                                                                        Abn_nonAD_SPR=sum(Abn_nonAD_SPR),
                                                                        Abn_nonAD_NPR=sum(Abn_nonAD_NPR)),
                                         by = .(Accession,Sequence,Symbol.hum)]
plaque_hum.pep.seqinfo[,`:=`(AD_SPR = Abn_AD_SPR/Abn_AD_NPR,
                             NonAD_SPR = Abn_nonAD_SPR/Abn_nonAD_NPR,
                             AD_NPR = Abn_AD_NPR/Abn_nonAD_NPR)] #Human

plaque_hum.pep.seqinfo <- merge(plaque_hum.pep.seqinfo,proteinseq.hum, by = 'Accession',suffix=c('.pep','.pro'))
temp <- t(plaque_hum.pep.seqinfo[,.(Sequence.pep,Sequence.pro)]) %>% as.data.table()
temp <- sapply(temp,function(vec)str_locate(vec[2],vec[1])) %>% t() %>% as.data.table() %>% `names<-`(c('start','end'))

plaque_hum.pep.seqinfo <- cbind(plaque_hum.pep.seqinfo, temp)
plaque_hum.pep.seqinfo[,`:=`(midpoint = (start + end)/2,
                             length = str_length(Sequence.pro))][,position := midpoint*100/length]

plaque_hum.peplist <- split(plaque_hum.pep.seqinfo,by='Accession')

PosCor <- function(table){
  if(nrow(table) < 3){
    p1 <- 1
    p2 <- 1
    p3 <- 1
  } else{
    fit1 <-lm(log2(AD_SPR) ~ position + I(position ^ 2),
              weights = Abn_AD_SPR,
              data = table)
    e1 <- summary(fit1)$coefficients
    if(nrow(e1) ==3){
      fit2 <-lm(log2(NonAD_SPR) ~ position + I(position ^ 2),
              weights = Abn_nonAD_SPR,
              data = table)
      fit3 <-lm(log2(AD_NPR) ~ position + I(position ^ 2),
              weights = Abn_AD_NPR,
              data = table)
      e2 <- summary(fit2)$coefficients
      e3 <- summary(fit3)$coefficients
      p1 <- min(e1[2:3,4])
      p2 <- min(e2[2:3,4])
      p3 <- min(e3[2:3,4])
    } else {
      p1 <- 1
      p2 <- 1
      p3 <- 1
    }
  }
  data.table(
      Accession = table$Accession[1],
      Genename = table$Symbol.hum[1],
      AD_SPR = p1,
      NonAD_SPR = p2,
      AD_NPR = p3
    )
}
part.enriched <- lapply(plaque_hum.peplist,PosCor) %>% rbindlist() # 30 seconds.

part.enriched.adsp <- plaque_hum.pep.seqinfo[AD_SPR >=2,Accession] %>% 
  intersect(part.enriched[AD_SPR <= 0.05,Accession])
part.enriched.nadsp <- plaque_hum.pep.seqinfo[NonAD_SPR >=2,Accession] %>%
  intersect(part.enriched[NonAD_SPR <= 0.05,Accession])
part.enriched.adbrain <- plaque_hum.pep.seqinfo[AD_NPR >=2,Accession] %>% 
  intersect(part.enriched[AD_NPR <= 0.05,Accession])

ExtremOutlier <- function(vec){
  q1 <- quantile(vec,0.25)
  q3 <- quantile(vec,0.75)
  llim <- q1-1.5*IQR(vec)
  hlim <- q3+1.5*IQR(vec)
  res <- vec > hlim | vec < llim
  res
}
AD_SPR_outlier <- plaque_hum.pep.seqinfo[ExtremOutlier(AD_SPR) & AD_SPR >=2,.SD, by = Accession]
nonAD_SPR_outlier <- plaque_hum.pep.seqinfo[ExtremOutlier(NonAD_SPR) & NonAD_SPR >= 2,.SD, by = Accession]
AD_NPR_outlier <- plaque_hum.pep.seqinfo[ExtremOutlier(AD_NPR) & AD_NPR >= 2,.SD, by = Accession]




AD_SPR_partial <- AD_SPR_outlier %>% 
  rbind(plaque_hum.pep.seqinfo[Accession %in% part.enriched.adsp & AD_SPR >= 2]) %>% unique()
nonAD_SPR_partial <- nonAD_SPR_outlier %>% 
  rbind(plaque_hum.pep.seqinfo[Accession %in% part.enriched.nadsp & NonAD_SPR >= 2]) %>% unique()
AD_NPR_partial <- AD_NPR_outlier %>% 
  rbind(plaque_hum.pep.seqinfo[Accession %in% part.enriched.adbrain & AD_NPR >= 2]) %>% unique()

write.csv(AD_SPR_partial,'./figures/AD_SPR_partial.csv')
write.csv(nonAD_SPR_outlier,'./figures/nonAD_SPR_partial.csv')
write.csv(AD_NPR_outlier,'./figures/AD_NPR_partial.csv')

plaque_mus.pep.seqinfo <- plaque_mus.pep[Accession %in% plaque_mus$Accession,.(Abn_AD_SPR1=sum(Abn_AD_SPR1),
                                                                               Abn_AD_SPR2=sum(Abn_AD_SPR2),
                                                                               Abn_AD_NPR1=sum(Abn_AD_NPR1),
                                                                               Abn_AD_NPR2=sum(Abn_AD_NPR2),
                                                                               Abn_WT_Control = sum(Abn_WT_Control)),
                                         by = .(Accession,Sequence,Symbol.mus)]
plaque_mus.pep.seqinfo[,`:=`(AD_SPR1 = Abn_AD_SPR1/Abn_AD_NPR1,
                             AD_SPR2 = Abn_AD_SPR2/Abn_AD_NPR2,
                             AD_NPR1 = Abn_AD_NPR1/Abn_WT_Control,
                             AD_NPR2 = Abn_AD_NPR2/Abn_WT_Control)] #Mouse

plaque_mus.pep.seqinfo <- merge(plaque_mus.pep.seqinfo,proteinseq.mus, by = 'Accession',suffix=c('.pep','.pro'))
temp <- t(plaque_mus.pep.seqinfo[,.(Sequence.pep,Sequence.pro)]) %>% as.data.table()
temp <- sapply(temp,function(vec)str_locate(vec[2],vec[1])) %>% t() %>% as.data.table() %>% `names<-`(c('start','end'))

plaque_mus.pep.seqinfo <- cbind(plaque_mus.pep.seqinfo, temp)
plaque_mus.pep.seqinfo[,`:=`(midpoint = (start + end)/2,
                             length = str_length(Sequence.pro))][,position := midpoint*100/length] 

plaque_mus.peplist <- split(plaque_mus.pep.seqinfo,by='Accession')

PosCor <- function(table){
  if(nrow(table) < 3){
    p1 <- 1
    p2 <- 1
    p3 <- 1
    p4 <- 1
  } else{
    fit1 <-lm(log2(AD_SPR1) ~ position + I(position ^ 2),
              weights = Abn_AD_SPR1,
              data = table)
    e1 <- summary(fit1)$coefficients
    if(nrow(e1) ==3){
      fit2 <-lm(log2(AD_SPR2) ~ position + I(position ^ 2),
                weights = AD_SPR2,
                data = table)
      fit3 <-lm(log2(AD_NPR1) ~ position + I(position ^ 2),
                weights = Abn_AD_NPR1,
                data = table)
      fit4 <-lm(log2(AD_NPR2) ~ position + I(position ^ 2),
                weights = Abn_AD_NPR2,
                data = table)
      e2 <- summary(fit2)$coefficients
      e3 <- summary(fit3)$coefficients
      e4 <- summary(fit4)$coefficients
      p1 <- min(e1[2:3,4])
      p2 <- min(e2[2:3,4])
      p3 <- min(e3[2:3,4])
      p4 <- min(e4[2:3,4])
    } else {
      p1 <- 1
      p2 <- 1
      p3 <- 1
      p4 <- 1
    }
  }
  data.table(
    Accession = table$Accession[1],
    Genename = table$Symbol.mus[1],
    AD_SPR1 = p1,
    AD_SPR2 = p2,
    AD_NPR1 = p3,
    AD_NPR2 = p4
  )
}
part.enriched <- lapply(plaque_mus.peplist,PosCor) %>% rbindlist() # 30 seconds.

part.enriched.adsp1.mus <- plaque_mus.pep.seqinfo[AD_SPR1 >=2,Accession] %>% 
  intersect(part.enriched[AD_SPR1 <= 0.05,Accession])
part.enriched.adsp2.mus <- plaque_mus.pep.seqinfo[AD_SPR2 >=2,Accession] %>% 
  intersect(part.enriched[AD_SPR2 <= 0.05,Accession])
part.enriched.adnp1.mus <- plaque_mus.pep.seqinfo[AD_NPR1 >=2,Accession] %>% 
  intersect(part.enriched[AD_NPR1 <= 0.05,Accession])
part.enriched.adnp2.mus <- plaque_mus.pep.seqinfo[AD_NPR2 >=2,Accession] %>% 
  intersect(part.enriched[AD_NPR2 <= 0.05,Accession])

AD_SPR1_outlier.mus <- plaque_mus.pep.seqinfo[ExtremOutlier(AD_SPR1) & AD_SPR1 >=2,.SD, by = Accession]
AD_SPR2_outlier.mus <- plaque_mus.pep.seqinfo[ExtremOutlier(AD_SPR2) & AD_SPR2 >=2,.SD, by = Accession]
AD_NPR1_outlier.mus <- plaque_mus.pep.seqinfo[ExtremOutlier(AD_NPR1) & AD_NPR1 >=2,.SD, by = Accession]
AD_NPR2_outlier.mus <- plaque_mus.pep.seqinfo[ExtremOutlier(AD_NPR2) & AD_NPR2 >=2,.SD, by = Accession]

AD_SPR1_partial <- AD_SPR1_outlier.mus %>% 
  rbind(plaque_mus.pep.seqinfo[Accession %in% part.enriched.adsp1.mus & AD_SPR1 >= 2]) %>% unique()
AD_SPR2_partial <- AD_SPR2_outlier.mus %>% 
  rbind(plaque_mus.pep.seqinfo[Accession %in% part.enriched.adsp2.mus & AD_SPR2 >= 2]) %>% unique()
AD_NPR1_partial <- AD_NPR1_outlier.mus %>% 
  rbind(plaque_mus.pep.seqinfo[Accession %in% part.enriched.adnp1.mus & AD_NPR1 >= 2]) %>% unique()
AD_NPR2_partial <- AD_NPR2_outlier.mus %>% 
  rbind(plaque_mus.pep.seqinfo[Accession %in% part.enriched.adnp2.mus & AD_NPR2 >= 2]) %>% unique()#Mouse

write.csv(AD_SPR1_partial,'./figures/AD_SPR1_partial.csv')
write.csv(AD_SPR2_partial,'./figures/AD_SPR2_partial.csv')
write.csv(AD_NPR1_partial,'./figures/AD_NPR1_partial.csv')
write.csv(AD_NPR2_partial,'./figures/AD_NPR2_partial.csv')

AD_SPR1_out_pro <- AD_SPR1_partial$Accession %>% unique()
AD_SPR2_out_pro <- AD_SPR2_partial$Accession %>% unique()
temp <- intersect(AD_SPR1_out_pro,AD_SPR2_out_pro)
AD_SPR_out_pro <- AD_SPR_partial$Accession %>% unique()
nonAD_SPR_out_pro <- nonAD_SPR_partial$Accession %>% unique()
temp1 <- intersect(AD_SPR_out_pro,nonAD_SPR_out_pro)
partial.core.pros <- plaque_hum.mus_orth[Accession.mus %in% temp & Accession.hum %in% temp1]
temp2 <- plaque_hum.mus_orth[Accession.hum %in% nonAD_SPR_out_pro & Accession.mus %in% temp]
write.csv(partial.core.pros,'./figures/Partial_Core_Pros.csv')
rm(AD_SPR1_out_pro,
   AD_SPR2_out_pro,
   AD_SPR_out_pro,
   nonAD_SPR_out_pro,
   partial.core.pros)

library(VennDiagram)
tiff(filename = './figures/Partial_core_venn.tiff', width = 4, height = 4, units = 'in', res = 300)
FigS4E_venn <- draw.triple.venn(
  area1 = 182,
  area2 = 75,
  area3 = 401,
  n12 = 41,
  n23 = 27,
  n13 = 42,
  n123 = 23,
  cex = rep(1.5, 7),
  fontfamily = rep('sans', 7),
  euler.d = T,
  scaled = T,
  margin = 0.05
)
grid.draw(FigS4E_venn)
dev.off()
rm(FigS4E_venn)
detach('package:VennDiagram', unload = T)


rm(
  part.enriched,
  part.enriched.adsp,
  part.enriched.nadsp,
  part.enriched.adbrain,
  ExtremOutlier,
  AD_SPR_outlier,
  nonAD_SPR_outlier,
  AD_NPR_outlier,
  PosCor,
  AD_SPR_partial,
  nonAD_SPR_partial,
  AD_NPR_partial,
  part.enriched.adsp1.mus,
  part.enriched.adsp2.mus, 
  part.enriched.adnp1.mus, 
  part.enriched.adnp2.mus, 
  AD_SPR1_outlier.mus,
  AD_SPR2_outlier.mus, 
  AD_NPR1_outlier.mus,
  AD_NPR2_outlier.mus, 
  AD_SPR1_partial, 
  AD_SPR2_partial,
  AD_NPR1_partial,
  AD_NPR2_partial
)

# Figure 3I plot. (Human)(Not updated.) -------------------------------------------
hum.pep.plot <- function(pepinfo.table,unipro.no){
  table <- filter(pepinfo.table,Master_Protein_Accessions == unipro.no)
  main <- paste0(table$Symbol.hum[1],'(',table$Master_Protein_Accessions[1],')')
  label1 <- 'aa1'
  label2 <- paste0('aa',table$length[1])
  label <- paste0('aa',table$start,'-',table$end)
  table$label <- label
  table <- melt(table,measure.vars = c('AD_SPR',
                                       'NonAD_SPR',
                                       'AD_PPR'),
                id.vars = c('Abn_AD_SPR',
                            'Abn_AD_NPR',
                            'Abn_nonAD_SPR',
                            'label',
                            'position'),
                variable.name = 'loci',
                value.name = 'ratio') %>% mutate(loci = as.character(loci))
  temp <- rep(0,nrow(table))
  for (i in 1:nrow(table)) {
    temp[i] <- switch (table$loci[i],
                       'AD_SPR' = log10(table$Abn_AD_SPR[i]),
                       'NonAD_SPR' = log10(table$Abn_nonAD_SPR[i]),
                       'AD_PPR' = log10(table$Abn_AD_NPR[i])
    )
  }
  table$Abn <- temp
  table$loci <- factor(table$loci,levels = c('AD_SPR','NonAD_SPR','AD_PPR'))
  label.table <- as.data.table(table) 
  label.table <- label.table[,.(ratio = max(ratio)), by = .(label,position)] %>% filter(ratio >=2 |ratio <= 0.5)
  plot <- ggplot(table, aes(x = position, y = log2(ratio)))
  plot + geom_line(aes(color = loci, linetype = loci),size=1) +
    geom_point(aes(size = Abn, color = loci),alpha = 0.5) +
    scale_size_continuous(range = c(1,4),guide = F) +
    scale_linetype_manual(values = c(1,1,2),guide=F)+
    scale_color_discrete(guide =F) +
    scale_x_continuous(breaks = c(1,100),
                       labels = c(label1,label2),
                       limits = c(0,101))+
    scale_y_continuous(limits = c(min(log2(table$ratio))-0.5,max(log2(table$ratio))+0.5),
                       breaks = c(-4:4)[between(-4:4, min(log2(table$ratio)) - 0.5, max(log2(table$ratio))+0.5)])+
    geom_text_repel(data = label.table, aes(label = label),
                    nudge_y = 0.2,
                    segment.size = 0.2)+
    ggtitle(main)+
    labs(x=NULL,y='Enrichment fold')+
    # geom_abline(slope = 0,intercept = 0.5)+
    theme_xf +
    theme(legend.position = c(0.15,0.7))
}
for (i in difpat.pro$Accession) {
  hum.pep.plot(hum.pep.seqinfo,unipro.no = i)
  ggsave(paste0(i,'.difpat','.tiff'),width = 5,height = 4,units = 'in',dpi = 150)
}
for (i in unique(hum.pep.outlier$Master_Protein_Accessions)) {
  hum.pep.plot(hum.pep.seqinfo,unipro.no = i)
  ggsave(paste0(i,'.pepoutlier','.tiff'),width = 5,height = 4,units = 'in',dpi = 150)
  
}
hum.pep.plot(hum.pep.seqinfo,unipro.no = 'P01024')
hum.pep.plot(hum.pep.seqinfo,unipro.no = 'Q07954')

# hum.pep.plot(hum.pep.seqinfo,unipro.no = 'P21741')
# Figure 3I plot. (Mouse)(Not updated.)-------------------------------------------------------
mus.pep.plot <- function(pepinfo.table,unipro.no){
  table <- filter(pepinfo.table,Master_Protein_Accessions == unipro.no,complete.cases(pepinfo.table))
  main <- paste0(table$Mouse.Associated.Gene.Name[1],'(',table$Master_Protein_Accessions[1],')')
  label1 <- 'aa1'
  label2 <- paste0('aa',table$length[1])
  label <- paste0('aa',table$start,'-',table$end)
  table$label <- label
  table <- melt(table,measure.vars = c('AD_SPR_1',
                                       'AD_SPR_2',
                                       'AD_PPR_1',
                                       'AD_PPR_2'),
                id.vars = c('Abn_Plaque_AD1',
                            'Abn_Plaque_AD2',
                            'Abn_Control_AD1',
                            'Abn_Control_AD2',
                            'label',
                            'position'),
                variable.name = 'loci',
                value.name = 'ratio') %>% mutate(loci = as.character(loci))
  temp <- rep(0,nrow(table))
  for (i in 1:nrow(table)) {
    temp[i] <- switch (table$loci[i],
                       'AD_SPR_1' = log10(table$Abn_Plaque_AD1[i]),
                       'AD_SPR_2' = log10(table$Abn_Plaque_AD2[i]),
                       'AD_PPR_1' = log10(table$Abn_Control_AD1[i]),
                       'AD_PPR_2' = log10(table$Abn_Control_AD2[i])
    )
  }
  table$Abn <- temp
  table$loci <- factor(table$loci,levels = c('AD_SPR_1','AD_SPR_2','AD_PPR_1','AD_PPR_2'))
  label.table <- as.data.table(table) 
  label.table <- label.table[,.(ratio = max(ratio)), by = .(label,position)] %>% filter(ratio >=2 |ratio <= 0.5)
  plot <- ggplot(table, aes(x = position, y = log2(ratio)))
  plot + geom_line(aes(color = loci, linetype = loci),size=1) +
    geom_point(aes(size = Abn, color = loci),alpha = 0.5) +
    scale_size_continuous(range = c(1,4),guide = F) +
    scale_linetype_manual(values = c(1,1,2,2),guide = F)+
    scale_color_discrete(guide = F) +
    labs(x=NULL,y='Enrichment fold')+
    scale_x_continuous(breaks = c(1,100),
                       labels = c(label1,label2),
                       limits = c(0,101))+
    scale_y_continuous(limits = c(min(log2(table$ratio))-0.5,max(log2(table$ratio))+0.5),
                       breaks = c(-4:4)[between(-4:4, min(log2(table$ratio)) - 0.5, max(log2(table$ratio))+0.5)])+
    geom_text_repel(data = label.table, aes(label = label),
                    nudge_y = 0.5,
                    segment.size = 0.6)+
    ggtitle(main)+
    # geom_abline(slope = 0,intercept = 0.5)+
    theme_xf +
    theme(legend.position = c(0.2,0.7))
}
for (i in mus.difpat.pro$Accession) {
  mus.pep.plot(mus.pep.seqinfo,unipro.no = i)
  ggsave(paste0(i,'.difpat','.mus','.tiff'),width = 5,height = 4,units = 'in',dpi = 150)
}
for (i in unique(mus.pep.outlier$Master_Protein_Accessions)) {
  mus.pep.plot(mus.pep.seqinfo,unipro.no = i)
  ggsave(paste0(i,'.pepoutlier','.mus','.tiff'),width = 5,height = 4,units = 'in',dpi = 150)
  
}
# mus.pep.plot(mus.pep.seqinfo,unipro.no = 'Q91ZX7')
mus.pep.plot(mus.pep.seqinfo,unipro.no = 'APPswe')
ggsave(paste0('APPswe','.mus','.tiff'),width = 5,height = 4,units = 'in',dpi = 150)
# Figure 3IAPPswe.pep.table (Specific,do not run,not updated.)--------------------------------------------------------
mus.ad.pep <- read.PD.csv('APPswe_Plaque_Mouse_Lumos_pep.csv') %>% 
  as.data.table() %>%
  left_join(mus_feature.new[,3:4],by = c('Master_Protein_Accessions' = 'Mouse.UniProt.SwissProt.Accession')) %>% 
  filter(Percolator_PEP_Sequest_HT < 0.01,Proteins ==1)
# mus.ad.pep <- mus.ad.pep[complete.cases(mus.ad.pep),]
mus.pep.seqinfo <- mus.ad.pep %>%
  group_by(Sequence,Master_Protein_Accessions,Mouse.Associated.Gene.Name) %>% 
  summarise(Abn_Plaque_AD1=sum(Abn__AD_SPR1),
            Abn_Plaque_AD2=sum(Abn__AD_SPR2),
            Abn_Control_AD1=sum(Abn__AD_NPR1),
            Abn_Control_AD2=sum(Abn__AD_NPR2),
            Abn_Control_WT = sum(Abn__WT_Control)) %>% 
  mutate(AD_SPR_1 = Abn_Plaque_AD1/Abn_Control_AD1,
         AD_SPR_2 = Abn_Plaque_AD2/Abn_Control_AD2,
         AD_PPR_1 = Abn_Control_AD1/Abn_Control_WT,
         AD_PPR_2 = Abn_Control_AD2/Abn_Control_WT)
mus.protein.seq <- read.csv('Mouse_Protein_sequence.txt',header = F,stringsAsFactors = F)
names(mus.protein.seq) <- c('Accession','ProSequence')
mus.pep.seqinfo <- left_join(mus.pep.seqinfo,mus.protein.seq,by = c('Master_Protein_Accessions'='Accession'))
result <- data.frame(start = NULL, end = NULL)
for (i in 1:nrow(mus.pep.seqinfo)) {
  temp <-
    str_locate(mus.pep.seqinfo$ProSequence[i] ,mus.pep.seqinfo$Sequence[i])%>% as.data.frame()
  result <- bind_rows(result, temp)
}
mus.pep.seqinfo <- cbind(as.data.frame(mus.pep.seqinfo),result) 
position <- function(vec){
  res <- rescale(vec)[2]
  res
}
mus.pep.seqinfo <- mus.pep.seqinfo %>% mutate(midpoint = (start + end)/2,
                                              length = str_length(ProSequence))
posi <- as.matrix(cbind(1,mus.pep.seqinfo[,16:17])) %>% apply(1,position)
mus.pep.seqinfo <- mus.pep.seqinfo %>% mutate(position= posi * 100)
rm(posi)
mus.pep.seqinfo$Mouse.Associated.Gene.Name[mus.pep.seqinfo$Master_Protein_Accessions == 'APPswe'] <- 'APPswe'
# Figure 4AB (See Figure 2D) ------------------------------------------
plaque.ratio <- plaque_hum.mus_orth[!duplicated(Accession.hum),.(Accession.hum,
                                                                 Rat_AD_SPR_AD_NPR.hum,
                                                                 Rat_nonAD_SPR_nonAD_NPR.hum,
                                                                 Rat_AD_SPR_AD_NPR.mus)]
temp <- plaque.ratio

plaque.ratio[,subset := rep(5,nrow(plaque.ratio))] 

plaque.ratio[Rat_AD_SPR_AD_NPR.hum > 1.149 &
               Rat_nonAD_SPR_nonAD_NPR.hum > 1.149 &
               Rat_AD_SPR_AD_NPR.mus > 1.149, subset := 1]
plaque.ratio[Rat_AD_SPR_AD_NPR.hum > 1.149 &
               Rat_nonAD_SPR_nonAD_NPR.hum <= 1.149 &
               Rat_AD_SPR_AD_NPR.mus > 1.149, subset := 2]
plaque.ratio[Rat_AD_SPR_AD_NPR.hum > 1.149 &
               Rat_AD_SPR_AD_NPR.mus <= 1.149, subset := 3]
plaque.ratio[Rat_AD_SPR_AD_NPR.hum <= 1.149 &
               Rat_AD_SPR_AD_NPR.mus > 1.149, subset := 4]
plaque.ratio[Rat_AD_SPR_AD_NPR.hum <= 1.149 &
               Rat_nonAD_SPR_nonAD_NPR.hum <= 1.149 &
               Rat_AD_SPR_AD_NPR.mus > 1.149, subset := 4]
plaque.ratio[Rat_AD_SPR_AD_NPR.hum <= 1.149 &
               Rat_nonAD_SPR_nonAD_NPR.hum <= 1.149 &
               Rat_AD_SPR_AD_NPR.mus <= 0.871, subset := 6]
plaque.ratio[Rat_AD_SPR_AD_NPR.hum < 0.871 &
               Rat_AD_SPR_AD_NPR.mus <= 0.871, subset := 7]


plaque_subs <- list(NULL)
for(i in 1:7){
  sub <- plaque.ratio[subset == i] %>% as.data.frame()
  sub <- sub[sample(1:nrow(sub),nrow(sub)),]
  plaque_subs[[i]] <- sub
} #Randomlization
plaque.ratio <- rbindlist(plaque_subs)
setorder(plaque.ratio,-subset)
plaque.ratio[,`:=`(Accession.hum=NULL, order=seq(nrow(plaque.ratio)), subset=NULL)]
GenOrderedTileplot(plaque.ratio,breaks = 0.5)
ggsave('./figures/Fig4_Plaqueheatmap.tiff',width = 3,height = 8,units = 'in',dpi = 300)

for(i in 1:length(plaque_subs)){
  plaque_subs_GOBP[[i]] <- gen.hum.table.GOBP.fis(plaque_subs[[i]],
                                                  background = 'table',
                                                  backgroundtable = temp)
} #GOBP analysis.

for(i in seq(length(plaque_subs_GOBP))){
  write.csv(plaque_subs_GOBP[[i]], file = paste0('./tables/','plaque_subs_GOBP',i,'.csv'),row.names = F)
}

library(topGO)
cc.table <- goterm[Name %in% c("cytoplasm",
                               'cytoskeleton',
                               "plasma membrane", 
                               "nucleus", 
                               "extracellular region",
                               "mitochondrion",
                               "endosome")]

plaque_subs.info <- lapply(plaque_subs,function(table)left_join(table,hum.info))
plaque_all <- rbindlist(plaque_subs.info)

hum.allgene <- plaque_all$Entrez.hum #Extract pros mapped to GO

CalcPerc <- function(vec){
  detect.genes <- factor(as.integer(hum.allgene %in% vec)) %>%
    `names<-`(hum.allgene) # Generate data for topGO.
  
  hum.GO.CC <- new('topGOdata',ontology = 'CC',
                   allGenes = detect.genes,
                   nodeSize = 10,
                   annot = annFUN.gene2GO,
                   gene2GO = hum.ent2go.list)
  res <- numeric(0)
  for(i in seq(nrow(cc.table))){
      test <- genesInTerm(hum.GO.CC,cc.table$ID[i]) %>%unlist()%>% as.integer()
      perc <- sum(test %in% vec)/length(vec)
      res <- c(res,perc)
  }
  names(res) <- cc.table$Name
  data.frame(as.list(res))
}
Res <- lapply(plaque_subs.info,function(table)unique(table$Entrez.hum)) %>% lapply(CalcPerc) %>% rbindlist(use.names = T) %>%mutate(subset = paste0('p',1:7))

GenBarPlot <- function(table){
  names(table) <- c('perc','subset')
  p <- ggplot(table,aes(x=subset, y=perc * 100))
  p + geom_bar(stat = 'identity',width = 0.8) + theme_xf
}
for(i in seq(nrow(cc.table))){
  go.name <- names(Res)[i]
  GenBarPlot(Res[,c(go.name,'subset')])
  ggsave(paste0('./figures/subsets.',go.name,'.tiff'),height = 2, width = 4,units = 'in', dpi = 150)
  
}

rm(Res,GenBarPlot,CalcPerc,plaque_all) #GO CC analysis
# Figure 4CD Wiki Pathway boxplot ------------------------------------------------
load('./database/Hum_wikipath_pros.Rdata')
load('./database/Mus_wikipath_pros.Rdata')
load('./database/Hum_path2Uni_pros.Rdata')
load('./database/Mus_Path2Uni_pros.Rdata')

RatioStat <- function(table) {
  do.call(rbind, lapply(table, function(vec)
    c(
      mean(vec),
      boxplot.stats(vec)$stats,
      length(vec),
      length(boxplot.stats(vec)$out)
    ))) %>% as.data.frame()
} #Define functions.

control <- abs(log2(plaque_mus$Abn_AD_NPR1/plaque_mus$Abn_AD_NPR2))
temp1 <- RatioStat(abs(log2(plaque_hum[,.(Rat_AD_SPR_AD_NPR,Rat_nonAD_SPR_nonAD_NPR)])))
temp2 <- RatioStat(abs(log2(plaque_mus[,.(Rat_AD_SPR_AD_NPR)])))
temp3 <- RatioStat(data.frame(a=control))
temp <- rbind(temp1,temp2) %>% rbind(temp3) %>% 
  mutate(loci = c('AD.SPs','nonAD.SPs','APP.SPs','Control'), pathname=c(rep('Total Proteome',3),'Control'))

names(temp)[1:8] <- c('mean','boxmin','boxlower','boxmedian','boxupper','boxmax','n','out')

plaque.hum <- plaque_hum %>% mutate(AD.brain = abs(log2(Rat_AD_NPR_nonAD_NPR)),
                                    AD.SPs = abs(log2(Rat_AD_SPR_AD_NPR)),
                                    nonAD.SPs = abs(log2(Rat_nonAD_SPR_nonAD_NPR))) %>%
  select(Accession,AD.brain,AD.SPs,nonAD.SPs) %>% as.data.table()

plaque.mus <- plaque_hum.mus_orth %>% mutate(APP.brain = abs(log2(Rat_AD_NPR_WT_Contol.mus)),
                                             APP.SPs = abs(log2(Rat_AD_SPR_AD_NPR.mus))) %>% as.data.table()

GenChangeReport.hum <- function(pros,table.hum= plaque.hum,table.mus=plaque.mus){
  path.hum <- table.hum[Accession %in% pros,.(AD.brain,AD.SPs,nonAD.SPs)]
  path.mus <- table.mus[Accession.hum %in% pros,.(APP.brain,APP.SPs)]
  hum.pros <- hum.info[Accession.hum %in% table.hum[Accession %in% pros, Accession],Symbol.hum]%>% paste(collapse = ',')
  mus.pros <- mus.info[Accession.mus %in% table.mus[Accession.hum %in% pros, Accession.mus], Symbol.mus] %>% paste(collapse = ',')
  
  if(nrow(path.hum) <5){
    return(NULL)
  }
  p1 <- t.test(path.hum$AD.brain, control)$p.value
  p2 <- t.test(path.hum$AD.SPs, control)$p.value
  p3 <- t.test(path.hum$nonAD.SPs, control)$p.value

  stat <- RatioStat(path.hum) %>% cbind(p.value = c(p1,p2,p3)) %>% as.data.frame()  
  stat <- stat %>% mutate(loci = rownames(stat),detected = hum.pros)
  if(nrow(path.mus) >=5){
    p4 <- t.test(path.mus$APP.brain, control)$p.value
    p5 <- t.test(path.mus$APP.SPs,control)$p.value
    stat.mus <- RatioStat(path.mus) %>% cbind(p.value = c(p4,p5))%>% as.data.frame 
    stat.mus <- stat.mus %>% mutate(loci = rownames(stat.mus),detected = mus.pros)
    stat <- rbind(stat,stat.mus)
  }
  stat <- as.data.frame(stat) %>% mutate(total = length(pros))
  names(stat)[1:8] <- c('mean','boxmin','boxlower','boxmedian','boxupper','boxmax','n','out')
  stat
} # Pathstat.hum
Path_analysis.hum <- function(path2uni.list,...){
  stat.list <- lapply(path2uni.list,GenChangeReport.hum,...)
  nullindex <- sapply(stat.list,is.null)
  stat.full.list <- stat.list[!nullindex]
  pathname <- names(stat.list)[!nullindex]
  for(i in seq(length(stat.full.list))){
    stat.full.list[[i]] <- stat.full.list[[i]] %>% mutate(pathname = pathname[i]) %>% as.data.table()
  }
  stat.table <- rbindlist(stat.full.list)
  return(stat.table)
} #Path_analysis.hum
  
wikipath.res.hum <- Path_analysis.hum(hum.wikipath.pros) %>% rbind(temp,fill = T)
sp.paths.wiki.hum <- wikipath.res.hum[loci =='AD.SPs' & p.value < 0.05 & n >=10 &  boxmedian > 0.0635]
sp.paths.wiki.hum <- wikipath.res.hum[pathname %in% sp.paths.wiki.hum$pathname]
# allpath.res.hum <- Path_analysis.hum(hum.path2uni.list) #Human result


plaque.hum <- plaque_hum.mus_orth %>% mutate(AD.brain = abs(log2(Rat_AD_NPR_nonAD_NPR.hum)),
                                    AD.SPs = abs(log2(Rat_AD_SPR_AD_NPR.hum)),
                                    nonAD.SPs = abs(log2(Rat_nonAD_SPR_nonAD_NPR.hum))) %>% as.data.table()
plaque.mus <- plaque_mus %>% mutate(APP.brain = abs(log2(Rat_AD_NPR_WT_Contol)),
                                             APP.SPs = abs(log2(Rat_AD_SPR_AD_NPR))) %>%
  select(Accession,APP.brain,APP.SPs) %>% as.data.table()

GenChangeReport.mus <- function(pros,table.hum= plaque.hum,table.mus=plaque.mus){
  path.hum <- table.hum[Accession.mus %in% pros,.(AD.brain,AD.SPs,nonAD.SPs)]
  path.mus <- table.mus[Accession %in% pros,.(APP.brain,APP.SPs)]
  hum.pros <- hum.info[Accession.hum %in% table.hum[Accession.mus %in% pros, Accession.hum],Symbol.hum]%>% paste(collapse = ',')
  mus.pros <- mus.info[Accession.mus %in% table.mus[Accession %in% pros, Accession], Symbol.mus] %>% paste(collapse = ',')

  
  if(nrow(path.mus) <5){
    return(NULL)
  }
  p4 <- t.test(path.mus$APP.brain, control)$p.value
  p5 <- t.test(path.mus$APP.SPs,control)$p.value
  stat <- RatioStat(path.mus) %>% cbind(p.value = c(p4,p5)) %>% as.data.frame()
  stat <- stat %>% mutate(loci = rownames(stat),detected = mus.pros)
  

  if(nrow(path.hum) >=5){
    
    p1 <- t.test(path.hum$AD.brain, control)$p.value
    p2 <- t.test(path.hum$AD.SPs, control)$p.value
    p3 <- t.test(path.hum$nonAD.SPs, control)$p.value
    stat.hum <- RatioStat(path.hum) %>% cbind(p.value = c(p1,p2,p3)) %>% as.data.frame 
    stat.hum <- stat.hum %>% mutate(loci = rownames(stat.hum),detected = hum.pros)
    stat <- rbind(stat.hum,stat)
  }
  stat <- as.data.frame(stat) %>% mutate(total = length(pros))
  names(stat)[1:8] <- c('mean','boxmin','boxlower','boxmedian','boxupper','boxmax','n','out')
  stat
} # Pathstat.mus
Path_analysis.mus <- function(path2uni.list,...){
  stat.list <- lapply(path2uni.list,GenChangeReport.mus,...)
  nullindex <- sapply(stat.list,is.null)
  stat.full.list <- stat.list[!nullindex]
  pathname <- names(stat.list)[!nullindex]
  for(i in seq(length(stat.full.list))){
    stat.full.list[[i]] <- stat.full.list[[i]] %>% mutate(pathname = pathname[i]) %>% as.data.table()
  }
  stat.table <- rbindlist(stat.full.list)
  return(stat.table)
} #Path_analysis.mus

wikipath.res.mus <- Path_analysis.mus(mus.wikipath.pros) %>% rbind(temp,fill= T)
sp.paths.wiki.mus <- wikipath.res.mus[loci == 'APP.SPs' & p.value <0.05 & n >= 10 & boxmedian > 0.0635]
sp.paths.wiki.mus <- wikipath.res.mus[pathname %in% sp.paths.wiki.mus$pathname]
# allpath.res.mus <- Path_analysis.mus(mus.path2uni.list) #Mouse result

rm(control,temp1,temp2,temp3,temp)

# Figure 4CD Plot changed pathways  (Human) (Not updated)------------------------------------------------

Top5.path <- wikipath.res.hum %>% filter(pathname %in% sp.paths.wiki.hum$pathname, loci =='AD.SPs') %>%
  arrange(desc(boxmedian)) %>% head(n=5)
plot.table <- wikipath.res.hum %>% 
  filter(pathname %in% c(Top5.path$pathname,'Total Proteome','Control'),
         loci %in% c('AD.SPs','nonAD.SPs','APP.SPs','Control')) %>% 
  mutate(pathname = factor(pathname, levels = c(Top5.path$pathname,'Total Proteome','Control')),
         loci = factor(loci, levels = c('AD.SPs','nonAD.SPs','APP.SPs','Control')))
p <- ggplot(plot.table, aes(x = pathname,fill = loci))
p + geom_boxplot(aes(lower = boxlower,
                     middle = boxmedian,
                     upper = boxupper,
                     ymin = boxmin,
                     ymax = boxmax),stat = 'identity') +
  coord_cartesian(ylim = c(0,2))+ #limits without removing data
  theme_xf 
  # theme(axis.text.x = element_text(angle =60))

ggsave('./figures/Figure4C_pathwaychange.hum.tiff', height = 5, width = 10, units = 'in',dpi = 150)

# Figure 4CD Plot changed pathways (Mouse) -------------------------------------------

Top5.path <- wikipath.res.mus %>% filter(n >=10, loci == 'APP.SPs') %>% arrange(desc(boxmedian)) %>% head(n=5)
plot.table <- wikipath.res.mus %>% 
  filter(pathname %in% c(Top5.path$pathname,'Total Proteome','Control'),
         loci %in% c('AD.SPs','nonAD.SPs','APP.SPs','Control')) %>% 
  mutate(pathname = factor(pathname, levels = c(Top5.path$pathname,'Total Proteome','Control')),
         loci = factor(loci, levels = c('AD.SPs','nonAD.SPs','APP.SPs','Control')))
p <- ggplot(plot.table, aes(x = pathname,fill = loci))
p + geom_boxplot(aes(lower = boxlower,
                     middle = boxmedian,
                     upper = boxupper,
                     ymin = boxmin,
                     ymax = boxmax),stat = 'identity') +
  coord_cartesian(ylim = c(0,2))+ #limits without removing data
  theme_xf 
# theme(axis.text.x = element_text(angle =60))

ggsave('./figures/Figure4D_pathwaychange.mus.tiff', height = 5, width = 10, units = 'in',dpi = 150)
# Figure 4EF pathway  Correlation analysis----------------------------------------------------
path.cor.table <- data.frame(
  pathname = character(0),
  AD_nonAD = numeric(0),
  AD_nonAD_dirc = numeric(0),
  AD_APP = numeric(0),
  AD_APP_dirc = numeric(0),
  N = numeric(0),
  orth.N = numeric(0)
)
Calcdirec <- function(table){
  direc <- sapply(log2(table),function(x) x>0)
  same <- sum(direc[,1]==direc[,2])
  perc <- same/nrow(table)*100
  perc
}
sp.paths.hum <- hum.wikipath.pros[names(hum.wikipath.pros)%in% sp.paths.wiki.hum$pathname]
sp.paths.mus <- mus.wikipath.pros[names(mus.wikipath.pros)%in% sp.paths.wiki.mus$pathname]
for(i in 1:length(sp.paths.hum)){
  pathname <- paste0(names(sp.paths.hum)[i],'.hum')
  pros.hum <- sp.paths.hum[[i]]
  path.hum <- plaque_hum[Accession %in% pros.hum]
  path.orth <- plaque_hum.mus_orth[Accession.hum %in% pros.hum] 
  # if (nrow(path.orth) <3){
  #   next()
  # }
  AD_nonAD <- cor(log2(path.hum$Rat_AD_SPR_AD_NPR), log2(path.hum$Rat_nonAD_SPR_nonAD_NPR))
  AD_nonAD_dirc <- Calcdirec(path.hum[,.(Rat_AD_SPR_AD_NPR,Rat_nonAD_SPR_nonAD_NPR)])
  AD_APP <- cor(log2(path.orth$Rat_AD_SPR_AD_NPR.hum),log2(path.orth$Rat_AD_SPR_AD_NPR.mus))
  AD_APP_dirc <- Calcdirec(path.orth[,.(Rat_AD_SPR_AD_NPR.hum,Rat_AD_SPR_AD_NPR.mus)])
  
  orth.N <- nrow(path.orth)
  N <- nrow(path.hum)
  path.cor.table <- rbindlist(list(path.cor.table, data.frame(pathname,AD_nonAD,AD_nonAD_dirc,AD_APP,AD_APP_dirc,N,orth.N)),use.names = F) 
}

for(i in 1:length(sp.paths.mus)){
  pathname <- paste0(names(sp.paths.mus)[i],'.mus')
  pros.mus <- sp.paths.mus[[i]]
  path.mus <- plaque_mus[Accession %in% pros.mus]
  path.orth <- plaque_hum.mus_orth[Accession.mus %in% pros.mus]
  AD_nonAD <- cor(log2(path.orth$Rat_AD_SPR_AD_NPR.hum), log2(path.orth$Rat_nonAD_SPR_nonAD_NPR.hum))
  AD_nonAD_dirc <- Calcdirec(path.orth[,.(Rat_AD_SPR_AD_NPR.hum,Rat_nonAD_SPR_nonAD_NPR.hum)])
  AD_APP <- cor(log2(path.orth$Rat_AD_SPR_AD_NPR.hum),path.orth$Rat_AD_SPR_AD_NPR.mus)
  AD_APP_dirc <- Calcdirec(path.orth[,.(Rat_AD_SPR_AD_NPR.hum,Rat_AD_SPR_AD_NPR.mus)])
  orth.N <- nrow(path.orth)
  N <- nrow(path.mus)
  path.cor.table <- rbindlist(list(path.cor.table, data.frame(pathname,AD_nonAD,AD_nonAD_dirc,AD_APP,AD_APP_dirc,N,orth.N)),use.names = F) 
}

path.cor.table <- filter(path.cor.table,orth.N>4) %>% filter(pathname != 'notch signaling pathway.hum' | orth.N >10) 
path.cor.table$pathname <- tolower(path.cor.table$pathname)

# Figure 4EF plot --------------------------------------------------------
interest.path <- c(
  "type ii interferon signaling (ifng).hum",
  'microglia pathogen phagocytosis pathway.mus',
  'alzheimers disease.hum',
  'tyrobp causal network.mus',
  'spinal cord injury.hum',
  'angiogenesis.hum',
  'mapk signaling pathway.hum',
  'toll-like receptor signaling pathway.hum',
  'glutathione metabolism.mus',
  "allograft rejection.hum",
  'apoptosis.hum',
  'complement activation.hum',
  'triacylglyceride synthesis.hum',
  'dna damage response.hum',
  'proteasome degradation.hum',
  'wnt signaling pathway.hum',
  'primary focal segmental glomerulosclerosis fsgs.hum',
  'notch signaling pathway.hum',
  'oxidative stress.hum',
  'inflammatory response pathway.hum'
)

p.table <-
  path.cor.table %>% 
  filter(pathname %in% interest.path)

p.table <- melt(p.table,measure.vars = c('AD_nonAD','AD_APP'),
                value.name = 'Correlation',
                variable.name = 'SPs') %>% mutate(pathname = factor(pathname, levels = rev(interest.path)),
                                                  SPs = factor(SPs, levels = c('AD_APP','AD_nonAD')))
p <- ggplot(p.table,aes(x=pathname,y= Correlation, fill = SPs))
p + geom_bar(color = 'black',stat = 'identity',position = 'dodge', width = 0.8)+
  coord_flip() +
  theme_xf +
  labs(x = NULL)

ggsave('./figures/Figure4E_SPs_pathway.tiff', width = 10, height = 10, units = 'in',dpi = 150)


# Figure 4E plot --------------------------------------------------------
interest.path <- c(
  'angiogenesis.hum',
  'microglia pathogen phagocytosis pathway.mus',
  'apoptosis.hum',
  "type ii interferon signaling (ifng).hum",
  'primary focal segmental glomerulosclerosis fsgs.hum',
  'spinal cord injury.hum',
  'tyrobp causal network.mus',
  'triacylglyceride synthesis.hum',
  'toll-like receptor signaling pathway.hum',
  'alzheimers disease.hum',
  'complement activation.hum',
  'mapk signaling pathway.hum',
  'glutathione metabolism.mus',
  "allograft rejection.hum",
  'dna damage response.hum',
  'proteasome degradation.hum',
  'wnt signaling pathway.hum',
  'notch signaling pathway.hum',
  'oxidative stress.hum',
  'inflammatory response pathway.hum'
)



p.table <-
  path.cor.table %>% 
  filter(pathname %in% interest.path) %>% arrange(desc(AD_APP_dirc))
levels <- p.table$pathname

p.table <- melt(p.table,measure.vars = c('AD_nonAD_dirc','AD_APP_dirc'),
                value.name = 'Direction',
                variable.name = 'SPs') %>% mutate(pathname = factor(pathname, levels = rev(levels)),
                                                  SPs = factor(SPs, levels = c('AD_APP_dirc','AD_nonAD_dirc'))) %>% 
  mutate(Direction = Direction-50)
p <- ggplot(p.table,aes(x=pathname,y= Direction, fill = SPs))
p + geom_bar(color = 'black',stat = 'identity',position = 'dodge', width = 0.8)+
  coord_flip() +
  theme_xf +
  ylim(-50,50)+
  labs(x = NULL)

ggsave('./figures/Figure4F_SPs_pathway_direct.tiff', width = 10, height = 10, units = 'in',dpi = 150)


# Start Aging analysis (Read data)---------------------------------------------------------
age_hip_hum <- ReadCleanMSProRes('Human_aging_HC.csv', bloodcontamine = T, commoncontamine = T)
age_tem_hum <- ReadCleanMSProRes('Human_aging_TL.csv', bloodcontamine = T, commoncontamine = T)

age_hip_hum[,`:=`(Rat_Old_Young = NULL, Abn_Old= NULL)]
setnames(age_hip_hum,c('Rat_Veryold_Young','Abn_Veryold'),c('Rat_Old_Young','Abn_Old'))
age_tem_hum[,`:=`(Rat_Old_Young = NULL, Abn_Old= NULL)]
setnames(age_tem_hum,c('Rat_Veryold_Young','Abn_Veryold'),c('Rat_Old_Young','Abn_Old')) #Human Aging



age_hip_mus <- ReadCleanMSProRes('Mouse_aging_HC.csv', bloodcontamine = T, commoncontamine = T)
age_tem_mus <- ReadCleanMSProRes('Mouse_aging_TL.csv', bloodcontamine = T, commoncontamine = T)

age_hip_mus <- age_hip_mus[,.(Accession,Description,Unique_Peptides,`Rat_50weeks,_WT_20weeks,_WT`,
                              `Rat_120weeks,_WT_20weeks,_WT`,`Abn_20weeks,_WT`,`Abn_50weeks,_WT`,`Abn_120weeks,_WT`,
                              Score_Sequest_HT,Genesym)]
setnames(age_hip_mus, 
         c('Rat_50weeks,_WT_20weeks,_WT','Rat_120weeks,_WT_20weeks,_WT','Abn_20weeks,_WT','Abn_50weeks,_WT','Abn_120weeks,_WT'),
         c('Rat_Middle_Young','Rat_Old_Young','Abn_Young','Abn_Middle','Abn_Old'))

age_tem_mus <- age_tem_mus[,.(Accession,Description,Unique_Peptides,`Rat_50weeks,_WT_20weeks,_WT`,
                              `Rat_120weeks,_WT_20weeks,_WT`,`Abn_20weeks,_WT`,`Abn_50weeks,_WT`,`Abn_120weeks,_WT`,
                              Score_Sequest_HT,Genesym)]
setnames(age_tem_mus, 
         c('Rat_50weeks,_WT_20weeks,_WT','Rat_120weeks,_WT_20weeks,_WT','Abn_20weeks,_WT','Abn_50weeks,_WT','Abn_120weeks,_WT'),
         c('Rat_Middle_Young','Rat_Old_Young','Abn_Young','Abn_Middle','Abn_Old')) #Mouse Aging


write.csv(age_hip_mus,'Mouse_aging_HC_filtered.csv')
write.csv(age_tem_mus,'Mouse_aging_TL_filtered.csv')
age_hh_mt <- genmatchtable(age_hip_hum,age_tem_mus)
age_ht_mh <- genmatchtable(age_tem_hum,age_hip_mus)
age_hip_hum.mus  <- genmatchtable(age_hip_hum,age_hip_mus)
age_tem_hum.mus  <- genmatchtable(age_tem_hum,age_tem_mus)


# Figure 5A corplot(Not updated, could be simplified) -------------------------------------------------------
library(corrplot)
age_hum <- merge(age_hip_hum,age_tem_hum,by = 'Accession')
age_mus <- merge(age_hip_mus,age_tem_mus,by = 'Accession')
age_min <- merge(age_hip_hum.mus,age_tem_hum.mus,by = 'Accession.hum')
names(age_min) <- names(age_min) %>%
  str_replace_all('\\.x','\\.hip') %>% 
  str_replace_all('\\.y','\\.tem')


abnmatrix.hum.hip <- select(age_hip_hum,starts_with('Abn'))
abnmatrix.hum.tem <- select(age_tem_hum,starts_with('Abn'))
abnmatrix.mus.hip <- select(age_hip_mus,starts_with('Abn'))
abnmatrix.mus.tem <- select(age_tem_mus,starts_with('Abn'))
abnmatrix.hum <- age_hum %>% filter(!duplicated(Accession)) %>% select(starts_with('Abn'))
abnmatrix.mus <- age_mus %>% filter(!duplicated(Accession)) %>% select(starts_with('Abn'))
abnmatrix.hip <- age_hip_hum.mus %>% filter(!duplicated(Accession.hum)) %>% select(starts_with('Abn'))
abnmatrix.tem <- age_tem_hum.mus %>% filter(!duplicated(Accession.hum)) %>% select(starts_with('Abn'))
abnmatrix.hhmt <- age_hh_mt %>% filter(!duplicated(Accession.hum)) %>% select(starts_with('Abn'))
abnmatrix.htmh <- age_ht_mh %>% filter(!duplicated(Accession.hum)) %>% select(starts_with('Abn'))
abnmatrix.min <- age_min %>% filter(!duplicated(Accession.hum)) %>% select(starts_with('Abn')) %>% select(contains('hum'),everything())
names(abnmatrix.min) <- c('HH.Y','HH.M','HH.O',
                          'HT.Y','HT.M','HT.O',
                          'MH.20A','MH.50A','MH.120A',
                          'MH.20W','MH.50W','MH.120W',
                          'MT.20A','MT.50A','MT.120A',
                          'MT.20W','MT.50W','MT.120W')
M <- cor(abnmatrix.min)
M.hum <- cor(abnmatrix.hum)
M.mus <- cor(abnmatrix.mus)
M.hip <- cor(abnmatrix.hip)
M.tem <- cor(abnmatrix.tem)
M.hhmt <- cor(abnmatrix.hhmt)
M.htmh <- cor(abnmatrix.htmh)
M.hh <- cor(abnmatrix.hum.hip)
M.ht <- cor(abnmatrix.hum.tem)
M.mh <- cor(abnmatrix.mus.hip)
M.mt <- cor(abnmatrix.mus.tem)
M[1:6,1:6] <- M.hum
M[7:18,7:18] <- M.mus
M[7:12,1:3] <- M.hip[4:9,1:3]
M[1:3,7:12] <- M.hip[1:3,4:9]
M[13:18,4:6] <- M.tem[4:9,1:3]
M[4:6,13:18] <- M.tem[1:3,4:9]
M[13:18,1:3] <- M.hhmt[4:9,1:3]
M[1:3,13:18] <- M.hhmt[1:3,4:9]
M[7:12,4:6] <- M.htmh[4:9,1:3]
M[4:6,7:12] <- M.htmh[1:3,4:9]
M[1:3,1:3] <- M.hh
M[4:6,4:6] <- M.ht
M[7:12,7:12] <- M.mh
M[13:18,13:18] <- M.mt

# 
# abnmatrix <- select(plaque_hum.mus, starts_with('Abn'))
# names(abnmatrix) <- c('AD_Plaque', 'AD_Brain', 'NonAD_Plaque', 'NonAD_Brain', 'AD_Mouse_Plaque1', 'AD_Mouse_Brain1',
#                       'AD_Mouse_Plaque2', 'AD_Mouse_Brain2', 'WT_Mouse_Brain')
# M <- cor(abnmatrix)
# abnmatrix.hum <- select(plaque_hum,starts_with('Abn'))
# M.hum <- cor(abnmatrix.hum)
# abnmatrix.mus <- select(plaque_mus,starts_with('Abn'))
# M.mus <- cor(abnmatrix.mus)
# M[1:4,1:4] <- M.hum
# M[5:9,5:9] <- M.mus
# 
tiff(filename = 'Fig6A_age_cor.tiff', width = 9, height = 9, units = 'in',res = 300)
corrplot.mixed(M, lower = 'ellipse',upper = 'number',
               tl.pos = 'lt',
               number.digits = 2,
               tl.col ='black'
)
dev.off()
# Figure 6E heat map
age_hip_hum <- ReadCleanMSProRes('Human_aging_HC.csv')

age_hip_hum <- age_hip_hum[,.(Accession,`Rat_(Middle)_(Young)`,`Rat_(Veryold)_(Young)`)] %>%
  `names<-`(c('Accession','Rat_Middle_Young','Rat_Old_Young'))
ad.age_hip_hum <- plaque_hum %>% merge(age_hip_hum,by = 'Accession')
ad.age_hip_hum <- ad.age_hip_hum[,.(Accession,Rat_AD_SPR_AD_NPR,Rat_Old_Young)][!duplicated(Accession)]

ad.age_hip_hum[,subset := 4] 

ad.age_hip_hum[Rat_AD_SPR_AD_NPR > 1.14 & Rat_Old_Young > 1.14 , subset := 1]
ad.age_hip_hum[Rat_AD_SPR_AD_NPR > 1.14 & Rat_Old_Young <= 1.14, subset := 2]
ad.age_hip_hum[Rat_AD_SPR_AD_NPR <= 1.14 & Rat_Old_Young > 1.14, subset := 3]
ad.age_hip_hum[Rat_AD_SPR_AD_NPR >= 0.88 & Rat_Old_Young < 0.88 & subset ==4 , subset := 5]
ad.age_hip_hum[Rat_AD_SPR_AD_NPR < 0.88 & Rat_Old_Young < 0.88 & subset ==4, subset := 6]



ad.age_hip_hum_subs <- list(NULL)
for(i in 1:6){
  sub <- ad.age_hip_hum[subset == i] %>% as.data.frame()
  sub <- sub[sample(1:nrow(sub),nrow(sub)),]
  ad.age_hip_hum_subs[[i]] <- sub
} #Randomlization
ad.age_hip_hum <- rbindlist(ad.age_hip_hum_subs)
ad.age_hip_hum_sub1pro <- ad.age_hip_hum[subset == 1, Accession]
setorder(ad.age_hip_hum,-subset)
ad.age_hip_hum[,`:=`(Accession=NULL, order=seq(nrow(ad.age_hip_hum)), subset=NULL)]
GenOrderedTileplot(ad.age_hip_hum,breaks = 0.5)
ggsave('./figures/Figure5E_Plaqueheatmap.tiff',width = 3,height = 8,units = 'in',dpi = 300)

ad.age_sub1_ppi <- hum.allppi.uni[Accession.hum.x %in% ad.age_hip_hum_sub1pro & Accession.hum.y %in% ad.age_hip_hum_sub1pro]

write.csv(ad.age_sub1_ppi,file = './figures/AD_Age_sub1_ppi.csv',quote = F,row.names = F)
temp <- plaque_hum_fea[,Description:=NULL]
write.csv(temp,file = './figures/hum.info.ad.csv',quote = F,row.names = F)


# Figure 5B corplot AD Age ------------------------------------------------
cor.table.humpart <- plaque_hum %>% merge(age_hip_hum, by = 'Accession', all = T) %>% 
  merge(age_tem_hum,by='Accession',suffix = c('.hip','.tem'),all = T) %>% 
  merge(ortholog.pairpro, by.x = 'Accession',by.y = 'Accession.hum',all.x = T)
cor.table.humpart <- cor.table.humpart[,.(Accession,Rat_AD_NPR_nonAD_NPR,Rat_AD_SPR_AD_NPR,Rat_nonAD_SPR_nonAD_NPR,
                                         Rat_Old_Young.hip,Rat_Old_Young.tem,Accession.mus)]
setnames(cor.table.humpart,
         c('Rat_AD_NPR_nonAD_NPR','Rat_AD_SPR_AD_NPR','Rat_nonAD_SPR_nonAD_NPR','Rat_Old_Young.hip','Rat_Old_Young.tem'),
         c('AD_brain','AD_SPs','nonAD_SPs','hum_HC_aging','hum_TL_aging')) #Human table for correlation


cor.table.muspart <- plaque_mus %>% merge(age_hip_mus, by = 'Accession', all = T) %>% 
  merge(age_tem_mus,by='Accession',suffix = c('.hip','.tem'),all = T) 
cor.table.muspart <- cor.table.muspart[,.(Accession,Rat_AD_NPR_WT_Contol,Rat_AD_SPR_AD_NPR,Rat_Old_Young.hip,Rat_Old_Young.tem)]

setnames(cor.table.muspart,
         c('Rat_AD_NPR_WT_Contol','Rat_AD_SPR_AD_NPR','Rat_Old_Young.hip','Rat_Old_Young.tem'),
         c('APPPS1_brain','APPPS1_SPs','mus_HC_aging','mus_TL_aging')) #Mouse table for correlation

cor.table <- merge(cor.table.humpart,cor.table.muspart,by.x = 'Accession.mus',by.y = 'Accession',all = T)
cor.table <- cor.table[,`:=`(Accession=NULL,Accession.mus=NULL)] %>% log2()


cor.result <- cor(cor.table, use='pairwise.complete.obs',method = 'spearman')


library(corrplot)

tiff(filename = './figures/AD_age_cor.tiff',width = 6, height = 6, units = 'in',res = 300)
corrplot.mixed(cor.result,lower = 'circle',upper = 'number',
               tl.pos = 'lt',
               number.digits = 2,
               tl.col ='black'
)
dev.off()

tiff(filename = './figures/AD_age_cor_pairs.tiff',width = 8, height = 8, units = 'in',res = 300)
library(car)
scatterplotMatrix(cor.table, smoother = NULL, diagonal = 'none',reg.line = F,use = "pairwise.complete.obs",ylim=c(-2,2),xlim=c(-2,2),pch=19,cex = 0.3,lwd=2,mar = c(0,0,0,0))
dev.off()

# Figure 5CD ---------------------------------------------------------------
plaque_sub.pro <- lapply(plaque_subs,function(table)table$Accession.hum)
sp.sub1.hum <- plaque_sub.pro[[1]]
sp.sub23.hum <- c(plaque_sub.pro[[2]],plaque_sub.pro[[3]])
sp.sub456.hum <- c(plaque_sub.pro[[4]],plaque_sub.pro[[5]],plaque_sub.pro[[6]])
sp.sub7.hum <- plaque_sub.pro[[7]]
sp.core.hum <- sp.core



age_hip_hum[Accession %in% sp.sub1.hum, class := 'sp.sub1'][Accession %in% sp.sub23.hum, class := 'sp.sub23'][Accession %in% sp.sub46.hum, class := 'sp.sub46'][Accession %in% sp.sub7.hum, class := 'sp.sub7']
# [Accession %in% sp.core.hum, class := 'sp.core']
age_hip_hum.sub <- age_hip_hum[complete.cases(age_hip_hum)]

age_tem_hum[Accession %in% sp.sub1.hum, class := 'sp.sub1'][Accession %in% sp.sub23.hum, class := 'sp.sub23'][Accession %in% sp.sub456.hum, class := 'sp.sub456'][Accession %in% sp.sub7.hum, class := 'sp.sub7']
# [Accession %in% sp.core.hum, class := 'sp.core']
age_tem_hum.sub <- age_tem_hum[complete.cases(age_tem_hum)] # classify


Hum2MusPro <- function(vec){
  res <- ortholog.pairpro[Accession.hum %in% vec, Accession.mus]
  unique(res)
}

sp.sub1.mus <- plaque_sub.pro[[1]] %>% Hum2MusPro()
sp.sub24.mus <- c(plaque_sub.pro[[2]],plaque_sub.pro[[4]]) %>% Hum2MusPro()
sp.sub35.mus <- c(plaque_sub.pro[[3]],plaque_sub.pro[[5]]) %>% Hum2MusPro()
sp.sub67.mus <- c(plaque_sub.pro[[6]],plaque_sub.pro[[7]]) %>% Hum2MusPro()
sp.core.mus <- sp.core %>% Hum2MusPro()


age_hip_mus[Accession %in% sp.sub1.mus, class := 'sp.sub1'][Accession %in% sp.sub24.mus, class := 'sp.sub24'][Accession %in% sp.sub35.mus, class := 'sp.sub35'][Accession %in% sp.sub67.mus, class := 'sp.sub67']
# [Accession %in% sp.core.mus, class := 'sp.core']
age_hip_mus.sub <- age_hip_mus[complete.cases(age_hip_mus)]

age_tem_mus[Accession %in% sp.sub1.mus, class := 'sp.sub1'][Accession %in% sp.sub24.mus, class := 'sp.sub24'][Accession %in% sp.sub35.mus, class := 'sp.sub35'][Accession %in% sp.sub67.mus, class := 'sp.sub67']
# [Accession %in% sp.core.mus, class := 'sp.core']
age_tem_mus.sub <- age_tem_mus[complete.cases(age_tem_mus)] # classify

# Figure 5CD plot ----------------------------------------------------------
Se <- function(vec){
  sd(vec)/sqrt(length(vec))
}

GenTimeLine <- function(table){
  table <- as.data.table(table)
  table <- table[,.(Rat_Middle_Young,Rat_Old_Young,class)][,Rat_Young_Young := 1] %>% melt(id.vars ='class')
  table <- table[,value := log2(value)][,.(Mean = mean(value), Se = Se(value)), by = .(variable,class)]
  table[,`:=`(variable = factor(variable,levels = c('Rat_Young_Young','Rat_Middle_Young','Rat_Old_Young')))]
  plot <- ggplot(table, aes(x=variable,y=Mean,group=class,color = class))
  plot + geom_line(size = 1)+
    geom_point(size =2)+
    geom_errorbar(aes(ymax = Mean + Se, ymin = Mean-Se),width = 0.2,size = 0.8)+
    theme_xf +
    ylim(c(-0.2,0.5))
}

GenTimeLine(age_hip_hum.sub)
ggsave('./figures/Figure5C_hip.hum.tiff', width = 5, height = 4, units = 'in',dpi = 150)
GenTimeLine(age_tem_hum.sub)
ggsave('./figures/Figure5C_tem.hum.tiff', width = 5, height = 4, units = 'in',dpi = 150)
GenTimeLine(age_hip_mus.sub)
ggsave('./figures/Figure5D_hip.mus.tiff', width = 5, height = 4, units = 'in',dpi = 150)
GenTimeLine(age_tem_mus.sub)
ggsave('./figures/Figure5D_tem.mus.tiff', width = 5, height = 4, units = 'in',dpi = 150)
# Figure 5E plot ----------------------------------------------------------
plaque_age_hip <- plaque_hum %>% merge(age_hip_hum,by='Accession')
plaque_age_hip <- plaque_age_hip[,.(Accession,Rat_AD_SPR_AD_NPR,Rat_Old_Young)]

plaque_age_hip[,subset := rep(4,nrow(plaque.ratio))] 

plaque_age_hip[Rat_AD_SPR_AD_NPR > 1.14 &
                 Rat_Old_Young > 1.14, subset := 1]
plaque_age_hip[Rat_AD_SPR_AD_NPR > 1.14 &
                 Rat_Old_Young <= 1.14, subset := 2]
plaque_age_hip[Rat_AD_SPR_AD_NPR <= 1.14 &
                 Rat_Old_Young > 1.14, subset := 3]
plaque_age_hip[Rat_AD_SPR_AD_NPR > 0.88 & subset == 4 &
                 Rat_Old_Young <= 0.88 , subset := 5]
plaque_age_hip[Rat_AD_SPR_AD_NPR <= 0.88 & subset == 4 &
                 Rat_Old_Young > 0.88 , subset := 6]
plaque_age_hip[Rat_AD_SPR_AD_NPR <= 0.88 & subset == 4 &
                 Rat_Old_Young <= 0.88 , subset := 7]

plaque_age_subs <- list(NULL)
for(i in 1:7){
  sub <- plaque_age_hip[subset == i] %>% as.data.frame()
  sub <- sub[sample(1:nrow(sub),nrow(sub)),]
  plaque_age_subs[[i]] <- sub
} #Randomlization
plaque_age_hip <- rbindlist(plaque_age_subs)
setorder(plaque_age_hip,-subset)
plaque_age_hip[,`:=`(Accession=NULL, order=seq(nrow(plaque_age_hip)), subset=NULL)]
GenOrderedTileplot(plaque_age_hip,breaks = 0.5)
ggsave('./figures/Figure5E_Plaqueheatmap.tiff',width = 3,height = 8,units = 'in',dpi = 300)
# Stoped here -------------------------------------------------------------
# Subset path (Result not good)-------------------------------------------------------
plaque_sub.pro <- lapply(plaque_subs,function(table)table$Accession.hum)
plaque_sub.ppi <- lapply(plaque_sub.pro,function(vec)hum.allppi.uni[Accession.hum.x %in% vec & Accession.hum.y %in% vec])
Pro2Path.hum <- function(pros,pathlist){
  hum.table<- plaque_hum[Accession %in% pros] %>% mutate(AD.brain = abs(log2(Rat_AD_NPR_nonAD_NPR)),
                                                           AD.SPs = abs(log2(Rat_AD_SPR_AD_NPR)),
                                                           nonAD.SPs = abs(log2(Rat_nonAD_SPR_nonAD_NPR))) %>%
    select(Accession,AD.brain,AD.SPs,nonAD.SPs) %>% as.data.table()
  
  mus.table <- plaque_hum.mus_orth[Accession.hum %in% pros] %>% mutate(APP.brain = abs(log2(Rat_AD_NPR_WT_Contol.mus)),
                                                                        APP.SPs = abs(log2(Rat_AD_SPR_AD_NPR.mus))) %>% as.data.table()
  Path_analysis(pathlist,table.hum = hum.table,table.mus= mus.table)
}
plaquesub.allpath <- lapply(plaque_sub.pro,Pro2Path.hum,pathlist = hum.path2uni.list)
plaquesub.wikipath <- lapply(plaque_sub.pro,Pro2Path.hum,pathlist = hum.wikipath.pros)

for(i in seq(length(plaque_sub.ppi))){
  temp <- plaque_sub.ppi[[i]] %>% as.data.table()
  temp <- temp[Interaction_type == 'in-complex-with'  | Interaction_type == 'interacts-with']
  write.csv(temp,file = paste0('./figures/','plaque_subset_',i,'.csv'),quote = F)
  temp <- hum.info[Accession.hum %in% plaque_sub.pro[[i]],Entrez.hum]
  write.table(temp,file = paste0('./figures/','plaque_subset_',i,'.txt'),quote = F,row.names = F)
}
write.csv(plaque_hum.mus_orth,file = './figures/ppi_node_info.csv',quote = F)

temp <- hum.info[Accession.hum %in% plaque_hum.mus_orth$Accession.hum]
write.csv(temp,file = './figures/hum_info.csv',quote = F,row.names = F)


temp <- rbindlist(plaque_subs) %>% merge(hum.info,by='Accession.hum',all.x = T)
write.csv(temp,file = './figures/orth_subsets_info.csv',row.names = F,quote = F)
allppi <- hum.allppi.uni[Accession.hum.x %in% plaque_hum.mus_orth$Accession.hum & 
                           Accession.hum.y %in% plaque_hum.mus_orth$Accession.hum]
write.csv(allppi,file = './figures/orth_ppi.csv',quote = F,row.names = F)




# Sup Figure 5A-D cluster -------------------------------------------------------
Gentimeline.hum <- function(cleanPD,clusternum,rearrange=c(1,2,3)){
  library(cluster)
  library(reshape2)
  cleanPD$Rat__Young <- 1
  matrix.rat <- select(cleanPD,starts_with('Rat'))
  clarax <- clara(matrix.rat[,2:1],clusternum)
  
  index <- clarax$clustering %>% factor(levels=rearrange,labels= c(1,2,3))
  matrix.rat$Cluster <- index 
  names(index) <- cleanPD$Accession
  assign('cluster.index',index,envir = parent.frame())
  grouped <- matrix.rat %>% 
    melt(id = 'Cluster', 
         value.name = 'Ratio',
         variable.name = 'Stage') %>% 
    group_by(Cluster,Stage) %>% 
    summarise(mean = mean(Ratio), se = sd(Ratio)/sqrt(n()),n=n())
  grouped$Stage <- factor(str_extract(grouped$Stage,'(?<=Rat__)\\w+'),
                          levels = c('Young','Middle','Veryold'),
                          labels = c('Young','Mid-age','Old'))
  grouped$Cluster <- as.factor(grouped$Cluster)
  plot <- ggplot(grouped, aes(x=Stage, y = mean, group = Cluster))
  plot + geom_line(aes(color= Cluster),size = 1) +
    geom_point(aes(color=Cluster),size = 2.5)+
    geom_errorbar(aes(ymin = (mean -se), ymax = (mean+se),color = Cluster),
                  width = 0.1,
                  size = 1) +
    geom_text(data = filter(grouped,Stage == tail(levels(grouped$Stage),1)),
              aes(label = paste0('N= ',n, '(Cluster ',Cluster,')'), 
                  color = Cluster),
              hjust = -0.3,
              fontface = 'bold') +
    scale_y_continuous(limits = c(0.8,1.5))+
    theme_xf +
    scale_x_discrete(expand = c(0,1.4))+
    scale_color_manual(values = c('firebrick2','dodgerblue2','forestgreen'),guide=F)+
    labs(x = NULL,y= expression('compared to Young group'))
  
}
Gentimeline.mus <- function(cleanPD,clusternum,rearrange=c(1,2,3),ADinclude = T){
  library(cluster)
  library(reshape2)
  cleanPD$Rat__20weeks_AD <- 1
  cleanPD$Rat__20weeks_WT <- 1
  matrix.rat <- select(cleanPD,starts_with('Rat'))[,4:9]
  clarax <- clara(matrix.rat[,1:4],clusternum)
  index <- clarax$clustering %>% factor(levels=rearrange,labels= c(1,2,3))
  matrix.rat$Cluster <- index 
  names(index) <- cleanPD$Accession
  assign('cluster.index',index,envir = parent.frame())
  grouped <- matrix.rat %>% 
    melt(id = 'Cluster', 
         value.name = 'Ratio',
         variable.name = 'Stage') %>% 
    group_by(Cluster,Stage) %>% 
    summarise(mean = mean(Ratio), se = sd(Ratio)/sqrt(n()),n=n())
  grouped$Sample <- 'WT'
  grouped$Sample[str_detect(grouped$Stage,'AD')] <- 'APP/PS1'
  grouped$Sample <- factor(grouped$Sample,levels = c('WT','APP/PS1'))
  grouped$Stage <- factor(str_extract(grouped$Stage,'(?<=Rat__)\\w+(?=_)'),
                          levels = c('20weeks','50weeks','120weeks'),
                          labels = c('Young','Mid-age','Old'))
  if(ADinclude == F){
    grouped <- grouped[grouped$Sample == 'WT',]
    grouped$Sample <- as.character(grouped$Sample)
  }
  grouped$Cluster <- as.factor(grouped$Cluster)
  plot <- ggplot(grouped, aes(x=Stage, y = mean, group = interaction(Cluster,Sample)))
  plot + geom_line(aes(color= Cluster,linetype = Sample),size = 1) +
    geom_point(aes(color=Cluster,shape = Sample),size = 2.5)+
    geom_errorbar(aes(ymin = (mean -se), ymax = (mean+se),color = Cluster),
                  width = 0.1,
                  size = 1) +
    geom_text(data = filter(grouped,Stage == tail(levels(grouped$Stage),1),Sample=='WT'),
              aes(label = paste0('N= ',n, '(Cluster ',Cluster,')'),  
                  color = Cluster),
              hjust = -0.3,
              fontface = 'bold') +
    scale_y_continuous(limits = c(0.9,1.2))+
    theme_xf +
    scale_x_discrete(expand = c(0,1.4))+
    scale_color_manual(values = c('firebrick2','dodgerblue2','forestgreen'),guide=F)+
    theme(legend.position = c(0.2,0.8))+
    labs(x = NULL,y= expression('compared to Young group'))
  
}

Gentimeline.mus(age_hip_mus,3,rearrange = c(1,3,2),ADinclude = F)
ggsave('Sup_Fig4C_mus_hip_time.tiff', width = 7, height = 4.5, units = 'in',dpi = 150)
cluster.age_hip_mus <- cluster.index

Gentimeline.mus(age_tem_mus,3,rearrange = c(1,3,2),ADinclude = F)
ggsave('Sup_Fig4D_mus_tem_time.tiff', width = 7, height = 4.5, units = 'in',dpi = 150)
cluster.age_tem_mus <- cluster.index

Gentimeline.hum(age_tem_hum,3,rearrange = c(1,2,3))
ggsave('Sup_Fig4B_hum_tem_time.tiff', width = 7, height = 4.5, units = 'in',dpi = 150)
cluster.age_tem_hum <- cluster.index

Gentimeline.hum(age_hip_hum,3,rearrange = c(2,1,3))
ggsave('Sup_Fig4A_hum_hip_time.tiff', width = 7, height = 4.5, units = 'in',dpi = 150)
cluster.age_hip_hum <- cluster.index
# clustered heatmap (Extremely slow, be careful!)----------------------
library(RColorBrewer)
cluster.age_hip_hum[cluster.index==2] <- 0
age_min_ratio <- age_min %>%
  mutate(Cluster = cluster.index[as.character(age_min$Accession.hum)]) %>% 
  filter(!duplicated(Genesym.hum.hip)) %>% 
  merge(hum_feature.new,
        by.x = 'Accession.hum', 
        by.y = 'Human.UniProt.SwissProt.Accession') %>% 
  arrange(Cluster)
rownames(age_min_ratio) <- age_min_ratio$Symbol.hum
age_min_rat <- select(age_min_ratio,starts_with('Rat'))[,c(1:2,10:11,6:9,15:18)]  %>% log2() %>% as.matrix()
colnames(age_min_rat) <- c('HH.M','HH.O',
                           'HT.M','HT.O',
                           'MH.50A','MH.120A','MH.50W','MH.120W',
                           'MT.50A','MT.120A','MT.50W','MT.120W')

my.palette <- colorRampPalette(c('forestgreen','forestgreen','grey90','firebrick4','firebrick4'))(n=399)
my.breaks <- c(seq(min(age_min_rat),-0.5,length = 100),
               seq(-0.49,0.49,length = 200),
               seq(0.5,max(age_min_rat),length =100))

tiff('Fig4C_AD_Heatmap.tiff',width = 10,height = 8,units = 'in',res = 150)

library(gplots)
heatmap.2(age_min_rat,
          density.info = 'none',
          margins = c(1.5,0),
          Colv = FALSE,
          Rowv = FALSE,
          col = my.palette,
          breaks = my.breaks,
          dendrogram = 'none',
          labRow = F,
          labCol = colnames(age_min_rat),
          cexCol = 1.5,
          trace = 'none',
          key.par = list(mar = c(0,0,0,0),
                         pin = c(0.8,0.2)),
          srtCol = 0,
          adjCol = c(0.5,0),
          key.title = F,
          keysize = 0.8)

dev.off()
# Sup Figure 5E-H Gen up/down proteins. -------------------------------------------------
Gentimeline.hum(age_hip_hum,3)
Gentimeline.hum(age_hip_hum,3,rearrange = c(2,1,3))
table(cluster.index)
age_hip_hum$Cluster <- cluster.index
age_hip_hum.up <- age_hip_hum %>% filter(Cluster==1)
age_hip_hum.down <- age_hip_hum %>% filter(Cluster==3)


Gentimeline.hum(age_tem_hum,3,rearrange = c(1,2,3))
table(cluster.index)
age_tem_hum$Cluster <- cluster.index
age_tem_hum.up <- age_tem_hum %>% filter(Cluster==1)
age_tem_hum.down <- age_tem_hum %>% filter(Cluster==3)

Gentimeline.mus(age_hip_mus,3,rearrange = c(1,3,2))
table(cluster.index)
age_hip_mus$Cluster <- cluster.index
age_hip_mus.up <- age_hip_mus %>% filter(Cluster==1)
age_hip_mus.down <- age_hip_mus %>% filter(Cluster==3)


Gentimeline.mus(age_tem_mus,3,rearrange = c(1,3,2))
table(cluster.index)
age_tem_mus$Cluster <- cluster.index
age_tem_mus.up <- age_tem_mus %>% filter(Cluster==1)
age_tem_mus.down <- age_tem_mus %>% filter(Cluster==3)
# Sup Figure 5E-H venn (Refer to Figure 1A)------------------------------------------------------
hum.up <- intersect(age_hip_hum.up$Accession,age_tem_hum.up$Accession)
age.hum.up <- age_hum %>% filter(Accession %in% hum.up)
mus.up <- intersect(age_hip_mus.up$Accession,age_tem_mus.up$Accession)
age.mus.up <- age_mus %>% filter(Accession %in% mus.up)
hip.up <- age_hip_hum.mus %>% filter(Accession.hum %in% age_hip_hum.up$Accession,Accession.mus %in% age_hip_mus.up$Accession) %>% 
  filter(!duplicated(Accession.hum))
tem.up <- age_tem_hum.mus %>% filter(Accession.hum %in% age_tem_hum.up$Accession,Accession.mus %in% age_tem_mus.up$Accession) %>% 
  filter(!duplicated(Accession.hum))
min.up <- age_min %>% filter(Accession.hum %in% hum.up, Accession.mus.hip %in% mus.up) %>% 
  filter(!duplicated(Accession.hum))

hum.down <- intersect(age_hip_hum.down$Accession,age_tem_hum.down$Accession)
age.hum.down <- age_hum %>% filter(Accession %in% hum.down)
mus.down <- intersect(age_hip_mus.down$Accession,age_tem_mus.down$Accession)
age.mus.down <- age_mus %>% filter(Accession %in% mus.down)
hip.down <- age_hip_hum.mus %>% filter(Accession.hum %in% age_hip_hum.down$Accession,Accession.mus %in% age_hip_mus.down$Accession) %>% 
  filter(!duplicated(Accession.hum))
tem.down <- age_tem_hum.mus %>% filter(Accession.hum %in% age_tem_hum.down$Accession,Accession.mus %in% age_tem_mus.down$Accession) %>% 
  filter(!duplicated(Accession.hum))
min.down <- age_min %>% filter(Accession.hum %in% hum.down, Accession.mus.hip %in% mus.down) %>% 
  filter(!duplicated(Accession.hum))
hum.cluster.alltable <- data.frame(Accession = union(age_hip_hum$Accession,age_tem_hum$Accession),
                                   Cluster =2,stringsAsFactors = F)
hum.cluster.alltable[hum.cluster.alltable$Accession %in% hum.up,2] <- 1
hum.cluster.alltable[hum.cluster.alltable$Accession %in% hum.down,2] <- 3


mus.cluster.alltable <- data.frame(Accession = union(age_hip_mus$Accession,age_tem_mus$Accession),
                                   Cluster =2,stringsAsFactors = F)
mus.cluster.alltable[mus.cluster.alltable$Accession %in% mus.up,2] <- 1
mus.cluster.alltable[mus.cluster.alltable$Accession %in% mus.down,2] <- 3 
# Sup Figure 5E-H plot -----------------------------------------------------------
Genvenn2 <- function(cleanPD1,cleanPD2,method = c('hh','mm','hm'),tiffname,cat1,cat2,angle = c(220, 140),
                     label.dist = c(0.05,0.05)){
  library(VennDiagram)
  t1 <-  cleanPD1$Accession[!duplicated(cleanPD1$Accession)]
  t2 <-  cleanPD2$Accession[!duplicated(cleanPD2$Accession)]
  a1 <- length(t1)
  a2 <- length(t2)
  if (method == 'hh' | method == 'mm') {
    tc <-  intersect(t1, t2)
    ac <- length(tc)
  } else if (method == 'hm') {
    tc.table <- genmatchtable(cleanPD1, cleanPD2,genclassified = F)
    ac <- tc.table %>% filter(!duplicated(Accession.hum)) %>% nrow()
  } else{
    stop('Not valid method')
  }
  tiff(
    filename = tiffname,
    width = 6,
    height = 6,
    units = 'in',
    res = 300
  )
  venn <- draw.pairwise.venn(
    area1 = a1,
    area2 = a2,
    cross.area = ac,
    category = c(paste(a1, '\n', cat1),
                 paste(a2, '\n', cat2)),
    col = c('#4682B4', '#771F1F'),
    alpha = c(0, 0),
    cex = rep(1.5, 3),
    fontfamily = rep('sans', 3),
    cat.pos = angle,
    cat.dist = label.dist,
    cat.prompts = T,
    cat.cex = rep(1.5, 2),
    cat.col = c('#4682B4', '#771F1F'),
    cat.just = c(list(c(0.8,-0.5), list(0.3,-0.5))),
    cat.fontfamily = rep('sans', 2),
    margin = 0.05
  )
  grid.draw(venn)
  dev.off()
}

Genvenn2(age_hip_hum.up,age_tem_hum.up,
         tiffname = 'Fig_S5E_humup_venn.tiff',
         method = 'hh',
         cat1 = 'hippocampus',
         cat2 = 'temporallobe')
Genvenn2(age_hip_mus.up,age_tem_mus.up,
         tiffname = 'Fig_S5G_musup_venn.tiff',
         method = 'mm',
         cat1 = 'hippocampus',
         cat2 = 'temporallobe')
Genvenn2(age_hip_hum.down,age_tem_hum.down,
         tiffname = 'Fig_S5F_humdown_venn.tiff',
         method = 'hh',
         angle = c(-30,30),
         label.dist = c(0.025,0.025),
         cat1 = 'hippocampus',
         cat2 = 'temporallobe')
Genvenn2(age_hip_mus.down,age_tem_mus.down,
         tiffname = 'Fig_S5H_musdown_venn.tiff',
         method = 'mm',
         angle = c(-30,30),
         label.dist = c(0.025,0.025),
         cat1 = 'hippocampus',
         cat2 = 'temporallobe')

Genvenn2(age.hum.up,age.mus.up,
         tiffname = 'Fig4D_up_venn.tiff',
         method = 'hm',
         angle = c(-40,40),
         label.dist = c(0.03,0.03),
         cat1 = 'Human',
         cat2 = 'Mouse')

Genvenn2(age.hum.down,age.mus.down,
         tiffname = 'Fig4D_down_venn.tiff',
         method = 'hm',
         label.dist = c(0.025,0.025),
         cat1 = 'Human',
         cat2 = 'Mouse')
# Figure 6B-C plot --------------------------------------------------------
GenPiePlot <- function(pros,species = c('hum','mus')){
  if(species == 'hum'){
    table <- left_join (data.frame(Accession = pros),hum.cluster.alltable)
  } else if(species == 'mus'){
    table <- left_join (data.frame(Accession = pros),mus.cluster.alltable)
  } else {
    stop('no information')
  }
  table[is.na(table$Cluster),2] <- 2
  x <- table$Cluster
  x1 <- sum(x==1)
  x2 <- sum(x==2)
  x3 <- sum(x==3)
  t1 <- round(x1*100/length(x),2)
  t2 <- round(x2*100/length(x),2)
  t3 <- round(x3*100/length(x),2)
  labels <- c(paste(c(x1,', ',t1,'%'),collapse = ''),
              paste(c(x3,', ',t3,'%'),collapse = ''),
              paste(c(x2,', ',t2,'%'),collapse = ''))
  breaks <- c((x1/2),
              (x1+x3/2),
              (x1+x3+x2/2))
  plot <- ggplot(table, aes(x = factor(1), fill = factor(Cluster, levels = c(1,3,2))))
  plot + geom_bar(width = 1) +
    scale_fill_manual(labels = c('aging-up','aging-down','unchanged'),
                      values = c('firebrick2','forestgreen','gray95'),drop = F) +
    scale_y_continuous(breaks = breaks, labels = labels)+
    coord_polar(theta = 'y') +
    theme_xf +
    labs(x = NULL, y=NULL)+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          # axis.text.x =element_blank(),
          axis.text.y = element_blank()
          )
}
GenPiePlot(hum.cluster.alltable$Accession,species = 'hum')
ggsave('Figure 6B.tiff', width = 6, height = 5, units = 'in',dpi = 150)
GenPiePlot(mus.cluster.alltable$Accession,species = 'mus')
ggsave('Figure 6C.tiff', width = 6, height = 5, units = 'in',dpi = 150)
# Figure 6D-G --------------------------------------------------------------
GenBpPoint2<- function(table.fis){
  plot <- ggplot(table.fis, aes(y = -log10(classicFisher), x = EnrichFold)) +
    geom_point(size=3.5,
               shape = 21,
               color = 'black') +
    # scale_fill_manual(values = c('white',muted('green'))) +
    # scale_x_continuous(limits = c(5,20),breaks = c(5,10,15,20))+
    labs(y = expression('P value (-log '[10]*')'),
         x = 'Enrichfold (Detected/Expected)') +
    expand_limits(y=1.5)+
    theme_xf
  plot
}

humageup <- gen.hum.table.GOBP.fis(data.frame(Accession = hum.up),background = 'table',backgroundtable = hum.cluster.alltable)
write.csv(humageup,file = 'Human_aging_up.csv')
GenBpPoint2(humageup[1:20,])
ggsave('Figure6E.tiff',height = 5,width = 6,units = 'in',dpi = 150)

humagedown <- gen.hum.table.GOBP.fis(data.frame(Accession = hum.down),background = 'table',backgroundtable = hum.cluster.alltable)
write.csv(humagedown,file = 'Human_aging_down.csv')
GenBpPoint2(humagedown[1:20,])
ggsave('Figure6F.tiff',height = 5,width = 6,units = 'in',dpi = 150)

musageup <- gen.mus.table.GOBP.fis(data.frame(Accession = mus.up),background = 'table',backgroundtable = mus.cluster.alltable)
write.csv(musageup,file = 'Mouse_aging_up.csv')
GenBpPoint2(musageup[1:20,])
ggsave('Figure6G.tiff',height = 5,width = 6,units = 'in',dpi = 150)

musagedown <- gen.mus.table.GOBP.fis(data.frame(Accession = mus.down),background = 'table',backgroundtable = mus.cluster.alltable)
write.csv(musagedown,file = 'Mouse_aging_down.csv')
GenBpPoint2(musagedown[1:20,])
ggsave('Figure6H.tiff',height = 5,width = 6,units = 'in',dpi = 150)
# Figure 7A-7D plot ------------------------------
Regentimeline.hum <- function(cleanPD,humup_pros,humdown_pros){
  library(reshape2)
  cleanPD$Rat__Young <- 1
  cleanPD$Cluster <- 2
  cleanPD$Cluster[cleanPD$Accession %in% humup_pros] <- 1
  cleanPD$Cluster[cleanPD$Accession %in% humdown_pros] <- 3
  matrix.rat <- select(cleanPD, starts_with('Rat'),Cluster)
  grouped <- matrix.rat %>%
    melt(id = 'Cluster',
         value.name = 'Ratio',
         variable.name = 'Stage') %>%
    group_by(Cluster, Stage) %>%
    summarise(mean = mean(Ratio),
              se = sd(Ratio) / sqrt(n()),
              n = n())
  grouped$Stage <- factor(
    str_extract(grouped$Stage, '(?<=Rat__)\\w+'),
    levels = c('Young', 'Middle', 'Veryold'),
    labels = c('Young', 'Mid-age', 'Old')
  )
  grouped$Cluster <- as.factor(grouped$Cluster)
  plot <- ggplot(grouped, aes(x = Stage, y = mean, group = Cluster))
  plot + geom_line(aes(color= Cluster),size = 1) +
    geom_point(aes(color=Cluster),size = 2.5)+
    geom_errorbar(aes(ymin = (mean -se), ymax = (mean+se),color = Cluster),
                  width = 0.1,
                  size = 1) +
    geom_text(data = filter(grouped,Stage == tail(levels(grouped$Stage),1)),
              aes(label = paste0('N= ',n, '(Cluster ',Cluster,')'), 
                  color = Cluster),
              hjust = -0.3,
              fontface = 'bold') +
    scale_y_continuous(limits = c(0.8,1.5))+
    theme_xf +
    scale_x_discrete(expand = c(0,1.4))+
    scale_color_manual(values = c('firebrick2','dodgerblue2','forestgreen'),guide=F)+
    labs(x = NULL,y= expression('compared to Young group (log '[2]*')')) 
  
  
}
Regentimeline.mus <- function(cleanPD,musup_pros,musdown_pros) {
  library(reshape2)
  cleanPD$Rat__20weeks_AD <- 1
  cleanPD$Rat__20weeks_WT <- 1
  
  cleanPD$Cluster <- 2
  cleanPD$Cluster[cleanPD$Accession %in% musup_pros] <- 1
  cleanPD$Cluster[cleanPD$Accession %in% musdown_pros] <- 3
  matrix.rat <- select(cleanPD, starts_with('Rat'))[, 4:9] %>% cbind(Cluster = cleanPD$Cluster)
  grouped <- matrix.rat %>%
    melt(id = 'Cluster',
         value.name = 'Ratio',
         variable.name = 'Stage') %>%
    group_by(Cluster, Stage) %>%
    summarise(
      mean = mean(Ratio),
      se = sd(Ratio) / sqrt(n()),
      n = n()
    )
  grouped$Sample <- 'WT'
  grouped$Sample[str_detect(grouped$Stage, 'AD')] <- 'PS1'
  grouped$Sample <- factor(grouped$Sample, levels = c('WT', 'PS1'))
  grouped$Stage <-
    factor(
      str_extract(grouped$Stage, '(?<=Rat__)\\w+(?=_)'),
      levels = c('20weeks', '50weeks', '120weeks'),
      labels = c('Young','Mid-age','Old')
    )
  grouped$Cluster <- as.factor(grouped$Cluster)
  plot <-
    ggplot(grouped, aes(
      x = Stage,
      y = mean,
      group = interaction(Cluster, Sample)
    ))
  plot + geom_line(aes(color= Cluster,linetype = Sample),size = 1) +
    geom_point(aes(color=Cluster,shape = Sample),size = 2.5)+
    geom_errorbar(aes(ymin = (mean -se), ymax = (mean+se),color = Cluster),
                  width = 0.1,
                  size = 1) +
    geom_text(data = filter(grouped,Stage == tail(levels(grouped$Stage),1),Sample=='WT'),
              aes(label = paste0('N= ',n, '(Cluster ',Cluster,')'),
                  color = Cluster),
              hjust = -0.3,
              fontface = 'bold') +
    scale_y_continuous(limits = c(0.8,1.5))+
    theme_xf +
    scale_x_discrete(expand = c(0,1.4))+
    scale_color_manual(values = c('firebrick2','dodgerblue2','forestgreen'),guide=F)+
    theme(legend.position = c(0.2,0.8))+
    labs(x = NULL,y= expression('compared to Young group (log '[2]*')'))
}

Regentimeline.hum(age_hip_hum,hum.up,hum.down)
ggsave('Fig7A_hum_hip_time.tiff', width = 7, height = 4.5, units = 'in',dpi = 150)

Regentimeline.hum(age_tem_hum,hum.up,hum.down)
ggsave('Fig7B_hum_tem_time.tiff', width = 7, height = 4.5, units = 'in',dpi = 150)

Regentimeline.mus(age_hip_mus,mus.up,mus.down)
ggsave('Fig7C_mus_hip_time.tiff', width = 7, height = 4.5, units = 'in',dpi = 150)

Regentimeline.mus(age_tem_mus,mus.up,mus.down)
ggsave('Fig7D_mus_tem_time.tiff', width = 7, height = 4.5, units = 'in',dpi = 150)
# Figure 7 adrelated gene in aging cluster -------------------------------
humad.all <- as.data.table(humad.all)
humadplaque.all <- as.data.table(humadplaque.all)
humnormplaque.all <- as.data.table(humnormplaque.all)
musad.all <- as.data.table(musad.all)
musadplaque1.all <- as.data.table(musadplaque1.all)
musadplaque2.all <- as.data.table(musadplaque2.all)

hum.adrelatedup.pro <- union(humad.all[humad.changed == 'up',Accession], 
                             humadplaque.all[humadplaque.changed == 'up',Accession]) %>% 
  union(humnormplaque.all[humnormplaque.changed == 'up', Accession])
mus.adrelatedup.pro <- union(musad.all[musad.changed == 'up',Accession],
                             musadplaque1.all[musadplaque1.changed == 'up',Accession]) %>% 
  union(musadplaque2.all[musadplaque2.changed == 'up',Accession])

hum.adrelateddown.pro <- union(humad.all[humad.changed == 'down',Accession], 
                               humadplaque.all[humadplaque.changed == 'down',Accession]) %>% 
  union(humnormplaque.all[humnormplaque.changed == 'down', Accession])
mus.adrelateddown.pro <- union(musad.all[musad.changed == 'down',Accession],
                               musadplaque1.all[musadplaque1.changed == 'down',Accession]) %>% 
  union(musadplaque2.all[musadplaque2.changed == 'down',Accession])
# Figure 7E-H plot ----------------------------------------------------------
GenPiePlot2 <- function(pros,species = c('hum','mus')){
  if(species == 'hum'){
    table <- left_join (data.frame(Accession = pros),hum.cluster.alltable)
  } else if(species == 'mus'){
    table <- left_join (data.frame(Accession = pros),mus.cluster.alltable)
  } else {
    stop('no information')
  }
  table[is.na(table$Cluster),2] <- 2
  x <- table$Cluster
  x1 <- sum(x==1)
  x2 <- sum(x==2)
  x3 <- sum(x==3)
  t1 <- round(x1*100/length(x),2)
  t2 <- round(x2*100/length(x),2)
  t3 <- round(x3*100/length(x),2)
  labels <- c(paste(c(x1,', ',t1,'%'),collapse = ''),
              paste(c(x3,', ',t3,'%'),collapse = ''),
              paste(c(x2,', ',t2,'%'),collapse = ''))
  breaks <- c((x1/2),
              (x1+x3/2),
              (x1+x3+x2/2))
  plot <- ggplot(table, aes(x = factor(1), fill = factor(Cluster, levels = c(1,3,2))))
  plot + geom_bar(width = 1) +
    scale_fill_manual(labels = c('aging-up','aging-down','unchanged'),
                      values = c('firebrick2','forestgreen','gray95'),drop = F) +
    scale_y_continuous(breaks = breaks, labels = labels)+
    coord_polar(theta = 'y') +
    theme_xf +
    labs(x = NULL, y=NULL)+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x =element_blank(),
          axis.text.y = element_blank()
    )
}
GenPiePlot(hum.adrelatedup.pro,species = 'hum')
GenPiePlot2(hum.adrelatedup.pro,species = 'hum')
ggsave('Fig7E_humadup_venn.tiff', width = 5, height = 4, units = 'in',dpi = 150)
GenPiePlot(hum.adrelateddown.pro,species = 'hum')
GenPiePlot2(hum.adrelateddown.pro,species = 'hum')
ggsave('Fig7F_humaddown_venn.tiff', width = 5, height = 4, units = 'in',dpi = 150)
GenPiePlot(mus.adrelatedup.pro,species = 'mus')
GenPiePlot2(mus.adrelatedup.pro,species = 'mus')
ggsave('Fig7G_musadup_venn.tiff', width = 5, height = 4, units = 'in',dpi = 150)
GenPiePlot(mus.adrelateddown.pro,species = 'mus')
GenPiePlot2(mus.adrelateddown.pro,species = 'mus')
ggsave('Fig7H_musaddown_venn.tiff', width = 5, height = 4, units = 'in',dpi = 150)
# Figure 7I Final tile plot (Human)-----------------------------------------------
musadplaque.changed.0.05 <- inner_join(musadplaque1.changed.0.05,musadplaque2.changed.0.05)
final.table <- left_join(adrelated.table,hum.cluster.alltable, by = c('Accession.hum'='Accession')) %>% 
  left_join(humad.changed.0.05[,c(1,ncol(humad.changed.0.05))],by =c('Accession.hum'='Accession')) %>% 
  left_join(humadplaque.changed.0.05[,c(1,ncol(humadplaque.changed.0.05))],by =c('Accession.hum'='Accession')) %>% 
  left_join(humnormplaque.changed.0.05[,c(1,ncol(humnormplaque.changed.0.05))],by =c('Accession.hum'='Accession')) %>% 
  left_join(mus.cluster.alltable, by = c('Accession.mus'='Accession')) %>% 
  left_join(musad.changed.0.05[,c(1,ncol(musad.changed.0.05))],by =c('Accession.mus'='Accession')) %>% 
  left_join(musadplaque.changed.0.05[,c(1,ncol(musadplaque.changed.0.05)-1)],by =c('Accession.mus'='Accession')) %>% 
  left_join(musadplaque.changed.0.05[,c(1,ncol(musadplaque.changed.0.05))],by =c('Accession.mus'='Accession'))
final.table$Cluster.x <- as.numeric(final.table$Cluster.x)
final.table$Cluster.y <- as.numeric(final.table$Cluster.y)
final.table$Cluster.x[final.table$Cluster.x ==3] <- 0.5
final.table$Cluster.x[final.table$Cluster.x ==1] <- 2.1
final.table$Cluster.x[final.table$Cluster.x ==2] <- 1
final.table$Cluster.y[final.table$Cluster.y ==1] <- 2.1
final.table$Cluster.y[final.table$Cluster.y ==3] <- 0.5
final.table$Cluster.y[final.table$Cluster.y ==2] <- 1

final.table.hum <- final.table[,c(2:5,11)] %>% 
  filter(Genesym.hum %in% adrelated.hum.table$Genesym.hum) %>% 
  filter((Cluster.x == 2.1)|(Cluster.x == 0.5))
TilePlot3 <- function(table){
  
  plot.t <- melt(table,
                 id.vars = 1,
                 value.name = 'Ratio',
                 variable.name = 'Site') %>% as.data.table()
  plot.t[,Ratio:=log2(Ratio)]
  plot.t$Genesym.hum <- factor(plot.t$Genesym.hum, levels = rev(table$Genesym.hum))
  
  plot <- ggplot(plot.t,aes(Site,Genesym.hum))
  plot + geom_tile(aes(fill = Ratio)) +
    scale_fill_gradientn(
      colours = c(
        'forestgreen',
        'forestgreen',
        'white',
        'firebrick4',
        'firebrick4'),
      breaks = c(trunc(min(plot.t$Ratio, na.rm = T)*100)/100,
                 -1,
                 0,
                 1,
                 trunc(max(plot.t$Ratio, na.rm = T)*100)/100),
      values = scales::rescale(c(min(plot.t$Ratio, na.rm = T),
                                 -1,
                                 0,
                                 1,
                                 max(plot.t$Ratio, na.rm = T))),
      labels = c(trunc(min(plot.t$Ratio, na.rm = T)*100)/100,
                 '-1',
                 '0',
                 '1',
                 trunc(max(plot.t$Ratio, na.rm = T)*100)/100)) +
    # scale_color_manual(values = c(up = 'firebrick',down = 'forestgreen',nc=NA),na.value = NA,guide =F)+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    labs(x=NULL,y=NULL)+
    theme(legend.position = "top", 
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = rel(0.8)),
          axis.text.x = element_text(size = rel(0.8), 
                                     angle = 330, 
                                     hjust = 0, 
                                     colour = "grey50"))
}
final.table.hum.index<-c("SPOCK2","HTRA1","GFAP","CSRP1","ACSF3","AQP1","CD44","ALDH2","LTF","DEFA1" )
final.table.hum <- final.table.hum[match(final.table.hum.index,final.table.hum$Genesym.hum),]
TilePlot3(final.table.hum)
ggsave('Figure7G.tiff',width = 3, height = 4, units = 'in',dpi = 150)

final.table.all <- final.table[,c(2:5,11,6:8,15)]
TilePlot3(final.table.all)
ggsave('Figure7G_all.tiff',width = 4, height = 25, units = 'in',dpi = 150)


final.table.hum <- data.frame(lapply(final.table.hum, na2nc))
names(final.table.hum) <- c('Genename.hum','Human AD brain','Human AD plaque','Human Non AD plaque','Human aging')
final.table.hum <- final.table.hum %>% arrange(desc(`Human AD plaque`),desc(`Human Non AD plaque`),desc(`Human AD brain`),desc(`Human aging`))



plot.table <- melt(final.table.hum,
                   id.vars = 'Genename.hum',
                   value.name = 'change',
                   variable.name = 'site')
plot.table$change <- factor(plot.table$change,levels = c('up','down','nc'))
plot.table$Genename.hum <- factor(plot.table$Genename.hum,levels = rev(final.table.hum$Genename.hum))
plot.table$site <- factor(plot.table$site,levels = names(final.table.hum[2:9]))
plot <- ggplot(plot.table,aes(site,Genename.hum))
plot + geom_tile(aes(fill = change),color = 'white',na.rm = F) +
  scale_fill_manual(values = c('red', 'green','grey95'),na.value = 'grey50',guide = F) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  # coord_flip()+
  labs(x=NULL,y=NULL)+
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = rel(1)),
        axis.text.x = element_text(size = rel(0.8),
                                   angle = 330,
                                   hjust = 0,
                                   colour = "grey50"),
        plot.margin = unit(c(0.1,0.5,0.1,0.1),'in'))
ggsave('Fig7G_final_pro.hum.tiff', width = 3, height = 7, units = 'in',dpi = 150)
# Regenerate HeatMap class (Extremely slow too....Final!)------------------------------------------------

age_min_ratio <- age_min %>% mutate(Cluster = 2) %>%  
  filter(!duplicated(Genesym.hum.hip)) %>% 
  merge(hum_feature.new,
        by.x = 'Accession.hum', 
        by.y = 'Human.UniProt.SwissProt.Accession') 
age_min_ratio$Cluster[age_min_ratio$Accession.hum %in% hum.up] <- 1
age_min_ratio$Cluster[age_min_ratio$Accession.hum %in% hum.down] <- 3
rownames(age_min_ratio) <- age_min_ratio$Symbol.hum

age_min_ratio <- age_min_ratio %>% arrange(Cluster)
age_min_rat <- select(age_min_ratio,starts_with('Rat'))[,c(1:2,10:11,6:9,15:18)]  %>% log2() %>% as.matrix()
colnames(age_min_rat) <- c('HH.M','HH.O',
                           'HT.M','HT.O',
                           'MH.50A','MH.120A','MH.50W','MH.120W',
                           'MT.50A','MT.120A','MT.50W','MT.120W')
library(RColorBrewer)
my.palette <- colorRampPalette(c('forestgreen','forestgreen','grey90','firebrick4','firebrick4'))(n=399)
my.breaks <- c(seq(min(age_min_rat),-0.5,length = 100),
               seq(-0.49,0.49,length = 200),
               seq(0.5,max(age_min_rat),length =100))


tiff('Fig4C_AD_Heatmap_regen.tiff',width = 10,height = 8,units = 'in',res = 150)

library(gplots)
heatmap.2(age_min_rat,
          density.info = 'none',
          margins = c(1.5,0),
          Colv = FALSE,
          Rowv = FALSE,
          col = my.palette,
          breaks = my.breaks,
          dendrogram = 'none',
          labRow = F,
          labCol = colnames(age_min_rat),
          cexCol = 1.5,
          trace = 'none',
          key.par = list(mar = c(0,0,0,0),
                         pin = c(0.8,0.2)),
          srtCol = 0,
          adjCol = c(0.5,0),
          key.title = F,
          keysize = 0.8)

dev.off()
# Sub Figure 5 I-J Density Plot --------------------------------------------------------
GenClusterDensity <- function(cleanPD,up_pro,down_pro){
  Cluster <- rep(2,nrow(cleanPD))
  Cluster[cleanPD$Accession %in% up_pro] <- 1
  Cluster[cleanPD$Accession %in% down_pro] <- 3
  abn_matrix <- select(cleanPD,starts_with('Abn'))
  ave_abn <- transmute(abn_matrix,ave = apply(abn_matrix, 1, mean)) %>% cbind(Cluster=Cluster)
  ave_abn$Cluster <- factor(ave_abn$Cluster, labels = c('Cluster1','Cluster2','Cluster3'))
  P <- ggplot(ave_abn,aes(x=log10(ave))) +
    geom_line(stat = 'density', aes(color=Cluster), size = 1.2) +
    scale_color_manual(values = c('firebrick2','dodgerblue2','forestgreen')) +
    labs(x= expression('Protein abundance (log'[2]*')'), y = 'Density')
  P + theme_xf + theme(legend.position = c(.8,.9))
  
}
GenClusterDensity(age_hum,hum.up, hum.down)
ggsave('Sup_Fig4I_hum_density.tiff', width = 6, height = 5, units = 'in',dpi = 150)
GenClusterDensity(age_mus,mus.up, mus.down)
ggsave('Sup_Fig4J_mus_density.tiff', width = 6, height = 5, units = 'in',dpi = 150)
# Generate Cluster heatmap (full)------------------------------------------------
GenTilePlot <- function(data,x = 'Loci',y = 'Accession',fill = 'cluster',
                        colorvalue =c('firebrick2','gray95','forestgreen'), na.value = 'gray5'){
  plot <- ggplot(Cluster.all, aes_string(x, y))
  plot + geom_tile(aes_string(fill = fill)) +
    scale_fill_manual(values = c('firebrick2','gray95','forestgreen'),na.value = 'gray50',guide=F) +
    scale_y_discrete(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    labs(x =NULL, y =NULL)+
    theme(
      legend.position = "top",
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(
        size = rel(0.8),
        angle = 330,
        hjust = 0,
        colour = "grey50"
      )
    )
}
ReformCluster <- function(cluster,species=c('hum','mus')){
  title <- deparse(substitute(cluster))
  Accession <- data.frame(Accession = as.character(names(cluster)))
  cluster <- data.frame(Accession,value = unname(cluster))
  if (species == 'hum') {
    result <- left_join(Accession,cluster)
    
  } else if(species == 'mus'){
    temp <- left_join(Accession,hum2mus.pair.protein[,c(4,8)], by = c('Accession' = 'Mouse.UniProt.SwissProt.Accession')) %>% 
      left_join(cluster)
    result <- data.frame(Accession = coalesce(temp[,2],temp[,1]),a = temp[,3])
  } else {
    stop('No species information')
  }
  names(result)[2] <- title
  result
}
Cluster.matrix <- full_join(ReformCluster(cluster.age_hip_hum,species = 'hum'),
                            ReformCluster(cluster.age_tem_hum,species = 'hum')) %>% 
  full_join(ReformCluster(cluster.age_hip_mus,species = 'mus')) %>% 
  full_join(ReformCluster(cluster.age_tem_mus,species = 'mus')) %>% 
  arrange(cluster.age_hip_hum, cluster.age_tem_hum,cluster.age_hip_mus,cluster.age_tem_mus)
Cluster.all <- melt(Cluster.matrix,id.vars = 'Accession', variable.name = 'Loci', value.name = 'cluster')
Cluster.all$Accession <- factor(Cluster.all$Accession, levels = rev(Cluster.matrix$Accession))
GenTilePlot(Cluster.all)
ggsave('Sup_FigNN_heatmap.tiff', width = 1.5, height = 5, units = 'in',dpi = 150)
# Figure 4F, Sup Figure 5K.Generate Cluster heatmap (hum vs. mus)  ---------------------------------
hum.cluster <- rep(2,nrow(age_hum))
names(hum.cluster) <- age_hum$Accession
hum.cluster[names(hum.cluster)%in% hum.up] <- 1
hum.cluster[names(hum.cluster)%in% hum.down] <- 3

mus.cluster <- rep(2,nrow(age_mus))
names(mus.cluster) <- age_mus$Accession
mus.cluster[names(mus.cluster)%in% mus.up] <- 1
mus.cluster[names(mus.cluster)%in% mus.down] <- 3

Cluster.matrix <- full_join(ReformCluster(hum.cluster,species = 'hum'),
                            ReformCluster(mus.cluster,species = 'mus')) %>% arrange(hum.cluster,mus.cluster)
Cluster.all <- melt(Cluster.matrix,id.vars = 'Accession', variable.name = 'Loci', value.name = 'cluster')
Cluster.all$Accession <- factor(Cluster.all$Accession, levels = rev(as.character(Cluster.matrix$Accession)))
Cluster.all$cluster <- as.factor(Cluster.all$cluster)
GenTilePlot(Cluster.all)
ggsave('Sup_Fig5K_heatmap_hummus.tiff', width = 1.5, height = 7.5, units = 'in',dpi = 150)

Cluster.matrix <-  Cluster.matrix %>% filter(Accession %in% age_min$Accession.hum)
Cluster.all <- melt(Cluster.matrix,id.vars = 'Accession', variable.name = 'Loci', value.name = 'cluster')
Cluster.all$Accession <- factor(Cluster.all$Accession, levels = rev(as.character(Cluster.matrix$Accession)))
Cluster.all$cluster <- as.factor(Cluster.all$cluster)
GenTilePlot(Cluster.all)
ggsave('Fig4F_heatmap_hummus_simple.tiff', width = 1, height = 5, units = 'in',dpi = 150)
# Sup Figure 6N (Final tile plot full) ---------------------------------------
hum.ad.info <- left_join(humad.all,humadplaque.all) %>% left_join(humnormplaque.all) %>%
  select(Accession,Symbol.hum,contains('changed')) 

mus.ad.info <- left_join(musad.all,musadplaque1.all) %>% 
  left_join(musadplaque2.all) %>% 
  select(Accession, Mouse.Associated.Gene.Name,contains('changed'))

age.sig.adpro <- mus.ad.info %>% 
  left_join(plaque_hum.mus_orth[,c('Accession.hum','Accession.mus')],by=c('Accession'='Accession.mus')) %>% 
  full_join(hum.ad.info,by = c('Accession.hum'='Accession')) %>% 
  filter(humad.changed != 'nochange' | humadplaque.changed != 'nochange' | humnormplaque.changed != 'nochange'|
           musad.changed != 'nochange' | musadplaque1.changed != 'nochange' | musadplaque2.changed != 'nochange') %>% 
  left_join(mus.cluster.alltable) %>% 
  left_join(hum.cluster.alltable,by=c('Accession.hum' = 'Accession'),suffix = c('.mus','.hum')) %>% 
  filter(Cluster.hum !=2 | Cluster.mus !=2)
age.sig.adpro[,sapply(age.sig.adpro, is.factor)] <- lapply(age.sig.adpro[,sapply(age.sig.adpro, is.factor)],as.character)
age.sig.adpro[age.sig.adpro == 1] <- 'up'
age.sig.adpro[age.sig.adpro == 2] <- 'nc'
age.sig.adpro[age.sig.adpro == 3] <- 'down'
age.sig.adpro[age.sig.adpro == 'nochange'] <- 'nc'

age.sig.adpro <- age.sig.adpro[,c(7:10,12,2:5,11)]
age.sig.adpro <- age.sig.adpro %>% arrange(desc(humadplaque.changed),desc(humnormplaque.changed),desc(humad.changed),desc(Cluster.hum),
                                           desc(musadplaque1.changed),desc(musadplaque2.changed),desc(musad.changed),desc(Cluster.mus))
names(age.sig.adpro) <- c('Genename.hum','Human AD brain','Human AD plaque','Human Non AD plaque','Human aging','Genename.mus',
                          'Mouse AD brain','Mouse AD plaque1','Mouse AD plaque2','Mouse aging')
age.sig.adpro$Genename <- coalesce(age.sig.adpro$Genename.hum,paste0(age.sig.adpro$Genename.mus,'.mus'))
age.sig.adpro <- age.sig.adpro[,c(-1,-6)]
plot.table <- melt(age.sig.adpro,
                   id.vars = 'Genename',
                   value.name = 'change',
                   variable.name = 'site')
plot.table$change <- factor(plot.table$change,levels = c('up','down','nc'))
plot.table$Genename.mus <- factor(plot.table$Genename,levels = rev(age.sig.adpro$Genename))

plot <- ggplot(plot.table,aes(site,Genename.mus))
plot + geom_tile(aes(fill = change),color = 'white',na.rm = F) +
  scale_fill_manual(values = c('red', 'green','grey95'),na.value = 'grey50',guide = F) +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  labs(x=NULL,y=NULL)+
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = rel(1.2)),
        axis.text.x = element_text(size = rel(1.2),
                                   angle = 320,
                                   hjust = 0,
                                   colour = "grey50"),
        plot.margin = unit(c(0.1,0.5,0.1,0.1),'in'))
ggsave('Fig5E_final_pro.mus.tiff', width = 4, height = 15, units = 'in',dpi = 150)





# pathvisio csv -----------------------------------------------------------
plaque.hum.full <- plaque_hum %>% left_join(hum2mus.pair.protein, by = c('Accession' = "Human.UniProt.SwissProt.Accession"))
plaque.mus.full <- plaque_mus %>% left_join(hum2mus.pair.protein, by = c('Accession' = "Mouse.UniProt.SwissProt.Accession"))

pathvisio.hum <- plaque.hum.full %>% mutate(AD.SPs = log2(`Rat_AD_SPR/AD_NPR`),
                                       nonAD.SPs = log2(`Rat_nonAD_SPR/nonAD_NPR`),
                                       AD.brain = log2(`Rat_AD_NPR/nonAD_NPR`)) %>% 
  select(Accession,
         Entrez.hum = Human.EntrezGene.ID,
         Entrez.mus = Mouse.EntrezGene.ID,
         AD.SPs,
         nonAD.SPs,
         AD.brain)
write.csv(pathvisio.hum,file = 'pathvisio.hum.csv',quote = F)

pathvisio.mus <- plaque.mus.full %>% mutate(APP.SPs = log2((Rat_AD_SPR1_AD_NPR1+Rat_AD_SPR2_AD_NPR2)/2),
                                       APP.brain = log2((`Rat_AD_NPR1/WT_Control` + `Rat_AD_NPR2/WT_Control`)/2)) %>% 
  select(Accession.mus = Accession,
         Entrez.hum = Human.EntrezGene.ID,
         Entrez.mus = Mouse.EntrezGene.ID,
         APP.SPs,
         APP.brain)
write.csv(pathvisio.mus,file = 'pathvisio.mus.csv',quote = F)

rm(plaque.hum.full,plaque.mus.full)


# Aging_pathway_analysis --------------------------------------------------
Na20 <- function(vec){
  vec[is.na(vec)] <- 0
  vec
}

hum.age <- full_join(age_hip_hum,age_tem_hum, by = 'Accession') %>% 
  select(Accession,
         hip.old = `Rat__Veryold/Young.x`,
         tem.old = `Rat__Veryold/Young.y`) %>% mutate(hum.old = log2((Na20(hip.old) + Na20(tem.old))/2)) %>% 
  left_join(hum2mus.pair.protein,by = c('Accession' = "Human.UniProt.SwissProt.Accession")) %>% 
  select(Accession.hum = Accession,
         Entrez.hum = Human.EntrezGene.ID,
         hum.old,
         Accession.mus = Mouse.UniProt.SwissProt.Accession,
         Entrez.mus = Mouse.EntrezGene.ID)%>% as.data.table()

mus.age <- full_join(age_hip_mus,age_tem_mus, by = 'Accession') %>% 
  select(Accession,
         hip.old = `Rat__120weeks_WT/20weeks_WT.x`,
         tem.old = `Rat__120weeks_WT/20weeks_WT.y`) %>% mutate(mus.old = log2((Na20(hip.old) + Na20(tem.old))/2)) %>% 
  left_join(hum2mus.pair.protein,by = c('Accession' = "Mouse.UniProt.SwissProt.Accession")) %>% 
  select(Accession.mus = Accession,
         Entrez.mus = Mouse.EntrezGene.ID,
         mus.old,
         Accession.hum = Human.UniProt.SwissProt.Accession,
         Entrez.hum = Human.EntrezGene.ID)%>% as.data.table()

write.csv(hum.age, file = 'pathvisio.humage.csv',quote = F)
write.csv(mus.age, file = 'pathvisio.musage.csv',quote = F)


orth.age <- inner_join(hum.age[,.(Accession.hum,hum.old,Accession.mus,Entrez.hum)],
                       mus.age[,.(Accession.mus,mus.old,Entrez.mus)], 
                       by = 'Accession.mus') 
write.csv(orth.age,file = 'pathvisio.orthage.csv',quote = F)

hum.path.stat <- list(NULL)

for(i in 1:length(hum.pathwaypros)){
  pathname <- names(hum.pathwaypros)[i]
  path.hum <- hum.age %>% filter(Accession.hum %in% hum.pathwaypros[[i]]) %>% select(hum.old)
  path.mus <- mus.age %>% filter(Accession.hum %in% hum.pathwaypros[[i]]) %>% select(mus.old)
  if(nrow(path.hum) <5) {next()}
  
  p1 <- t.test(path.hum$hum.old, hum.age$hum.old)$p.value
  stat <- RatioStat(path.hum) %>% cbind(p.value = p1)
  if(nrow(path.mus) >=5){
    p2 <- t.test(path.mus$mus.old, mus.age$mus.old)$p.value
    stat.mus <- RatioStat(path.mus) %>% cbind(p.value = p2)
    stat <- rbind(stat,stat.mus)
  }
  stat <- stat %>% as.data.frame() %>% mutate(pathname = pathname,
                                              loci = rownames(stat))
  names(stat)[1:8] <- c('mean','boxmin','boxlower','boxmedian','boxupper','boxmax','n','out')
  hum.path.stat <- c(hum.path.stat,list(pathname = stat))
}
hum.path.stat <- do.call(rbind,hum.path.stat) %>% `rownames<-`(NULL) #Human pathway stat

mus.path.stat <- list(NULL)

for(i in 1:length(mus.pathwaypros)){
  pathname <- names(mus.pathwaypros)[i]
  path.hum <- hum.age %>% filter(Accession.mus %in% mus.pathwaypros[[i]]) %>% select(hum.old)
  path.mus <- mus.age %>% filter(Accession.mus %in% mus.pathwaypros[[i]]) %>% select(mus.old)
  if(nrow(path.mus) <5) {next()}
  
  p1 <- t.test(path.mus$mus.old, mus.age$mus.old)$p.value
  stat <- RatioStat(path.mus) %>% cbind(p.value = p1)
  if(nrow(path.hum) >=5){
    p2 <- t.test(path.hum$hum.old, hum.age$hum.old)$p.value
    stat.hum <- RatioStat(path.hum) %>% cbind(p.value = p2)
    stat <- rbind(stat,stat.hum)
  }
  stat <- stat %>% as.data.frame() %>% mutate(pathname = pathname,
                                              loci = rownames(stat))
  names(stat)[1:8] <- c('mean','boxmin','boxlower','boxmedian','boxupper','boxmax','n','out')
  mus.path.stat <- c(mus.path.stat,list(pathname = stat))
}
mus.path.stat <- do.call(rbind,mus.path.stat) %>% `rownames<-`(NULL) #musan pathway stat



temp1 <- hum.age %>% 
  select(hum.old) %>% 
  RatioStat() %>%
  as.data.frame() %>% 
  mutate(p.value = 1, pathname = 'Total Proteome',loci= rownames(.))

temp2 <- mus.age %>% 
  select(mus.old) %>% 
  RatioStat() %>%
  as.data.frame() %>% 
  mutate(p.value = 1, pathname = 'Total Proteome',loci= rownames(.))

hum.path.stat <- rbindlist(list(hum.path.stat,temp1,temp2),use.names = F,fill = F)
mus.path.stat <- rbindlist(list(mus.path.stat,temp1,temp2),use.names = F,fill = F)
rm(hum.age,mus.age,temp1,temp2,i)


# Prepare data--------------------------------------------------------------------

hip_hum <- ReadCleanMSRes('Human_aging_HC.csv',bloodcontamine=T)
hip_hum <- hip_hum %>% mutate(Rat_Old_Young_ = NULL, Abn_F1_128_Sample_Old = NULL) %>%
  rename(Rat_Middle_Young = Rat_Middle_Young_,
         Rat_Old_Young = Rat_Veryold_Young_,
         Abn_Young = Abn_F1_126_Sample_Young,
         Abn_Middle = Abn_F1_127_Sample_Middle,
         Abn_Old = Abn_F1_129_Sample_Veryold)  %>% as.data.table()
tem_hum <- ReadCleanMSRes('Human_aging_TL.csv',bloodcontamine=T) 
tem_hum <- tem_hum %>% mutate(Rat_Old_Young_ = NULL, Abn_F1_130_Sample_Old = NULL) %>%
  rename(Rat_Middle_Young = Rat_Middle_Young_,
         Rat_Old_Young = Rat_Veryold_Young_,
         Abn_Young = Abn_F1_126_Sample_Young,
         Abn_Middle = Abn_F1_129_Sample_Middle,
         Abn_Old = Abn_F1_131_Sample_Veryold) %>% as.data.table()
hip_mus <- ReadCleanMSRes('Mouse_aging_HC.csv',bloodcontamine=T)
hip_mus <- hip_mus  %>%  mutate(Rat_20weeks_AD_20weeks_WT_ =NULL,
                             Rat_50weeks_AD_50weeks_WT_ =NULL,
                             Rat_120weeks_AD_120weeks_WT_ =NULL) %>% 
  rename(Rat_Middle_Young = Rat_50weeks_WT_20weeks_WT_,
         Rat_Old_Young = Rat_120weeks_WT_20weeks_WT_,
         Rat_Middle_Young_APP = Rat_50weeks_AD_20weeks_AD_,
         Rat_Old_Young_APP = Rat_120weeks_AD_20weeks_AD_,
         Abn_Young = Abn_F1_127_Sample_20weeks_WT,
         Abn_Middle = Abn_F1_129_Sample_50weeks_WT,
         Abn_Old = Abn_F1_131_Sample_120weeks_WT,
         Abn_Young_APP = Abn_F1_126_Sample_20weeks_AD,
         Abn_Middle_APP = Abn_F1_128_Sample_50weeks_AD,
         Abn_Old_APP = Abn_F1_130_Sample_120weeks_AD) %>% as.data.table()
  
tem_mus <- ReadCleanMSRes('Mouse_aging_TL.csv',bloodcontamine=T)
tem_mus <- tem_mus  %>%  mutate(Rat_20weeks_AD_20weeks_WT_ =NULL,
                                Rat_50weeks_AD_50weeks_WT_ =NULL,
                                Rat_120weeks_AD_120weeks_WT_ =NULL) %>% 
  rename(Rat_Middle_Young = Rat_50weeks_WT_20weeks_WT_,
         Rat_Old_Young = Rat_120weeks_WT_20weeks_WT_,
         Rat_Middle_Young_APP = Rat_50weeks_AD_20weeks_AD_,
         Rat_Old_Young_APP = Rat_120weeks_AD_20weeks_AD_,
         Abn_Young = Abn_F1_127_Sample_20weeks_WT,
         Abn_Middle = Abn_F1_129_Sample_50weeks_WT,
         Abn_Old = Abn_F1_131_Sample_120weeks_WT,
         Abn_Young_APP = Abn_F1_126_Sample_20weeks_AD,
         Abn_Middle_APP = Abn_F1_128_Sample_50weeks_AD,
         Abn_Old_APP = Abn_F1_130_Sample_120weeks_AD
  ) %>% as.data.table()

age_hum <- merge(hip_hum,tem_hum,by= c('Accession','Description','Genesym'),all = T,suffix = c('.hip','.tem'))
age_mus <- merge(hip_mus,tem_mus,by= c('Accession','Description','Genesym'),all = T, suffix = c('.hip','.tem'))

age_full <-
  merge(age_hum, ortholog.pairpro, all.x = T,by.x = 'Accession',by.y = 'Accession.hum') %>% 
  merge(age_mus,all = T, by.x = 'Accession.mus',by.y = 'Accession',sort = F,suffix = c('.hum','.mus'))

plaque_age_hum <- merge(age_hum,plaque_hum,by = c('Accession','Description'),all = T)
plaque_age_mus <- merge(age_mus,plaque_mus,by = c('Accession','Description'),all = T)



# Gen_hum_mus_heatmap -----------------------------------------------------
p.table <- plaque_age_hum[,.(Accession,Rat_AD_SPR_AD_NPR,Rat_AD_NPR_nonAD_NPR,Rat_Old_Young.hip,Rat_Old_Young.tem)] %>% filter(complete.cases(.))
GenClusterMap(p.table,4) #NOT SIGNIFICANT

p.table <- plaque_age_mus[,.(Accession,Rat_AD_SPR_AD_NPR,
                             Rat_AD_NPR_WT_Contol,
                             Rat_Old_Young.tem,
                             Rat_Old_Young.hip)] %>% filter(complete.cases(.))
GenClusterMap(p.table,5) # NOT SIGINIFCANT
CalcMeanSE <- function(vec){
  mean = mean(vec,na.rm =T)
  se = sd(vec,na.rm = T)/sqrt(length(na.omit(vec)))
  data.frame(mean = mean,se = se)
}

plaque_subs.hip.hum <- lapply(plaque_subs,function(table)hip_hum[Accession %in% table$Accession.hum]) %>% 
  lapply(function(table)CalcMeanSE(log2(table$Rat_Old_Young))) %>% rbindlist()


plaque_subs.tem.hum <- lapply(plaque_subs,function(table)tem_hum[Accession %in% table$Accession.hum]) %>% 
  lapply(function(table)CalcMeanSE(log2(table$Rat_Old_Young))) %>% rbindlist()


plaque_subs.hip.mus <- lapply(plaque_subs,function(table)left_join(table,ortholog.pairpro)) %>%
  lapply(function(table)hip_mus[Accession %in% table$Accession.mus]) %>% 
  lapply(function(table)CalcMeanSE(log2(table$Rat_Old_Young))) %>% rbindlist()

plaque_subs.tem.mus <- lapply(plaque_subs,function(table)left_join(table,ortholog.pairpro)) %>%
  lapply(function(table)tem_mus[Accession %in% table$Accession.mus]) %>% 
  lapply(function(table)CalcMeanSE(log2(table$Rat_Old_Young))) %>% rbindlist()

GenPointRange <- function(table){
  table <- cbind(table,subset = paste0('p',1:7))
  p <- ggplot(table, aes(x=subset, y = mean))
  p + geom_point(size =1.5) + geom_errorbar(aes(ymax = mean+se, ymin = mean-se),width = 0.5,size =0.8)+
    ylim(-0.3,0.3)+theme_xf
}

GenPointRange(plaque_subs.hip.hum)
ggsave('plaque_subs.hip.hum.tiff',height = 3,width = 4,units = 'in',dpi = 150)
GenPointRange(plaque_subs.tem.hum)
ggsave('plaque_subs.tem.hum.tiff',height = 3,width = 4,units = 'in',dpi = 150)
GenPointRange(plaque_subs.hip.mus)
ggsave('plaque_subs.hip.mus.tiff',height = 3,width = 4,units = 'in',dpi = 150)
GenPointRange(plaque_subs.tem.mus)
ggsave('plaque_subs.tem.mus.tiff',height = 3,width = 4,units = 'in',dpi = 150)

# Analyze the SPs and aging upregulated proteins in human -----------------

temp <- plaque_age_hum[,.(Accession,Rat_AD_SPR_AD_NPR,Rat_Old_Young.hip)] %>% filter(complete.cases(.))

GenClusterMap(temp,5,levels = c(3,5,2,1,4),reorder = T)
clustered.table[,cluster:= as.numeric(cluster)]
clustered.table[Rat_AD_SPR_AD_NPR < 1.14 & cluster == 5, cluster:= 3]
clustered.table[Rat_Old_Young.hip < 1.14 & cluster == 5, cluster:= 4]
clustered.table[Rat_AD_SPR_AD_NPR < 1.14 & cluster == 4, cluster:= 2]
clustered.table[Rat_Old_Young.hip < 0.88 & cluster == 4, cluster:= 1]
clustered.table[Rat_Old_Young.hip < 1.14 & cluster == 3, cluster:= 2]
clustered.table[Rat_Old_Young.hip > 1.14 & cluster == 2, cluster:= 4]
clustered.table[Rat_Old_Young.hip < 0.88 & cluster == 2, cluster:= 1]

setorder(clustered.table, cluster)

p.table <- list(NULL)
for(i in 1:5){
  sub <- clustered.table[cluster == i] %>% as.data.frame()
  sub <- sub[sample(1:nrow(sub),nrow(sub)),]
  p.table[[i]] <- sub
} #Randomlization
p.table <- rbindlist(p.table) %>% as.data.table()
p.table[,`:=` (order= seq(1:nrow(p.table)),Accession = NULL)]
GenOrderedTileplot(p.table,breaks = 0.5)

ggsave('plaque_age_hum_heatmap.tiff',width = 2.5,height = 7, units = 'in',dpi =150)

# GO analysis of both plaque and age up proteins --------------------------
sp_age_up <- clustered.table[,`:=` (Accession.hum = Accession, Accession = NULL)][Rat_AD_SPR_AD_NPR>1.14 & Rat_Old_Young.hip>1.14]
temp <- gen.hum.table.GOBP.fis(sp_age_up,background = 'table',backgroundtable = clustered.table)
write.table(sp_age_up$Accession.hum,file = 'sp_age_up.txt',quote = F,row.names = F) # String analysis
write.table(plaque_age_hum$Accession,file='sp_age.txt',quote = F,row.names = F)
# SPs proteins time line --------------------------------------------------

# SPs_age_up_interaction --------------------------------------------------
temp <- read.csv('SPs_age_up.tsv',sep = '\t',stringsAsFactors = F)
temp1 <- c(temp$X.node1,temp$node2) %>% table() %>% as.data.frame()
write.table(temp1,'nodeinfo.txt',quote = F,row.names = F,sep = '\t')

# aging_correlation -------------------------------------------------------
library(corrplot)
temp <- age_full %>% select(starts_with('Abn')) %>% select(-contains('APP'))
cor.result <- cor(temp, use='pairwise.complete.obs',method = 'spearman')
tiff(filename = 'Age_cor_full.tiff',width = 7, height = 7, units = 'in',res = 300)
corrplot.mixed(cor.result,lower = 'ellipse',upper = 'number',
               tl.pos = 'lt',
               number.digits = 2,
               tl.col ='black'
)
dev.off()



# aging_cluster -----------------------------------------------------------

Gentimeline <- function(table,clusternum,rearrange=c(1,2,3)){
  library(cluster)
  library(reshape2)
  table$Rat_Young_Young <- 1
  matrix.rat <- select(table,starts_with('Rat'))
  clarax <- clara(matrix.rat[,1:2],clusternum)
  
  index <- clarax$clustering %>% factor(levels=rearrange,labels= c(1,2,3))
  matrix.rat$Cluster <- index 
  index <- data.frame(Accession=table$Accession,index = index)
  assign('cluster.index',index,envir = parent.frame())
  grouped <- matrix.rat %>% 
    melt(id = 'Cluster', 
         value.name = 'Ratio',
         variable.name = 'Stage') %>% 
    group_by(Cluster,Stage) %>% 
    summarise(mean = mean(Ratio), se = sd(Ratio)/sqrt(n()),n=n())
  grouped$Stage <- factor(str_extract(grouped$Stage,'(?<=Rat_)\\w+'),
                          levels = c('Young_Young','Middle_Young','Old_Young'),
                          labels = c('Young','Middle age','Old'))
  grouped$Cluster <- as.factor(grouped$Cluster)
  plot <- ggplot(grouped, aes(x=Stage, y = mean, group = Cluster))
  plot + geom_line(aes(color= Cluster),size = 1) +
    geom_point(aes(color=Cluster),size = 2.5)+
    geom_errorbar(aes(ymin = (mean -se), ymax = (mean+se),color = Cluster),
                  width = 0.1,
                  size = 1) +
    geom_text(data = filter(grouped,Stage == tail(levels(grouped$Stage),1)),
              aes(label = paste0('N= ',n, '(Cluster ',Cluster,')'), 
                  color = Cluster),
              hjust = -0.3,
              fontface = 'bold') +
    scale_y_continuous(limits = c(0.8,1.5))+
    theme_xf +
    scale_x_discrete(expand = c(0,1.4))+
    scale_color_manual(values = c('firebrick2','dodgerblue2','forestgreen'),guide=F)+
    labs(x = NULL,y= expression('compared to Young group'))
  
}
Gentimeline(hip_hum,clusternum = 3,rearrange = c(2,1,3))
cluster.index.hip.hum <- cluster.index %>% as.data.table()
Gentimeline(tem_hum,clusternum = 3,rearrange = c(1,2,3))
cluster.index.tem.hum <- cluster.index %>% as.data.table()


Gentimeline.mus <- function(table,clusternum,rearrange=c(1,2,3),ADinclude = T){
  library(cluster)
  library(reshape2)
  table$Rat_Young_Young <- 1
  table$Rat_Young_Young_APP <- 1
  matrix.rat <- select(table,starts_with('Rat'))
  clarax <- clara(matrix.rat[,1:4],clusternum)
  index <- clarax$clustering %>% factor(levels=rearrange,labels= c(1,2,3))
  matrix.rat$Cluster <- index 
  index <- data.frame(Accession=table$Accession,index = index)
  assign('cluster.index',index,envir = parent.frame())
  grouped <- matrix.rat %>% 
    melt(id = 'Cluster', 
         value.name = 'Ratio',
         variable.name = 'Stage') %>% 
    group_by(Cluster,Stage) %>% 
    summarise(mean = mean(Ratio), se = sd(Ratio)/sqrt(n()),n=n())
  grouped$Sample <- 'WT'
  grouped$Sample[str_detect(grouped$Stage,'APP')] <- 'APP/PS1'
  grouped$Sample[!str_detect(grouped$Stage,'APP')] <- 'WT'
  grouped$Sample <- factor(grouped$Sample,levels = c('WT','APP/PS1'))
  grouped$Stage <- grouped$Stage %>% str_replace_all('Rat_','') %>% 
    str_replace_all('_APP','') %>% factor(levels = c('Young_Young','Middle_Young','Old_Young'),
                                          labels = c('Young','Middle age','Old'))
  if(ADinclude == F){
    grouped <- grouped[grouped$Sample == 'WT',]
    grouped$Sample <- as.character(grouped$Sample)
  }
  grouped$Cluster <- as.factor(grouped$Cluster)
  plot <- ggplot(grouped, aes(x=Stage, y = mean, group = interaction(Cluster,Sample)))
  plot + geom_line(aes(color= Cluster,linetype = Sample),size = 1) +
    geom_point(aes(color=Cluster,shape = Sample),size = 2.5)+
    geom_errorbar(aes(ymin = (mean -se), ymax = (mean+se),color = Cluster),
                  width = 0.1,
                  size = 1) +
    geom_text(data = filter(grouped,Stage == tail(levels(grouped$Stage),1),Sample=='WT'),
              aes(label = paste0('N= ',n, '(Cluster ',Cluster,')'),  
                  color = Cluster),
              hjust = -0.3,
              fontface = 'bold') +
    scale_y_continuous(limits = c(0.8,1.4))+
    theme_xf +
    scale_x_discrete(expand = c(0,1.4))+
    scale_color_manual(values = c('firebrick2','dodgerblue2','forestgreen'),guide=F)+
    # theme(legend.position = c(0.3,0.8))+
    labs(x = NULL,y= expression('compared to Young group'))
  
}
Gentimeline.mus(hip_mus,3,rearrange = c(1,3,2),ADinclude = F)
cluster.index.hip.mus <- cluster.index %>% as.data.table()
Gentimeline.mus(tem_mus,3,rearrange = c(1,3,2),ADinclude = F)
cluster.index.tem.mus <- cluster.index %>% as.data.table()

cluster.index.hum <- merge(cluster.index.hip.hum,cluster.index.tem.hum,by = 'Accession',all = T)
cluster.index.hum[, index := 'nochange'][index.x == 1 &
                                           index.y == 1, index := 'age_up'][index.x == 3 & index.y == 3, index := 'age_down']

cluster.index.mus <- merge(cluster.index.hip.mus,cluster.index.tem.mus,by = 'Accession',all = T)
cluster.index.mus[, index := 'nochange'][index.x == 1 &
                                           index.y == 1, index := 'age_up'][index.x == 3 & index.y == 3, index := 'age_down']
test <- cluster.index.hum[Accession %in% core.pro$Accession.hum]
age_cluster_table <- merge(cluster.index.hum[,.(Accession,index)],
                           ortholog.pairpro,
                           all.x = T,
                           by.x = 'Accession',
                           by.y = 'Accession.hum') %>% merge(cluster.index.mus[,.(Accession,index)],
                                                             by.x = 'Accession.mus',
                                                             by.y = 'Accession',
                                                             all = T)
p.table <- age_cluster_table[complete.cases(age_cluster_table)][,Accession.mus := NULL]
p.table[,`:=`(index.x = factor(index.x,levels = c('age_up','nochange','age_down')),
              index.y = factor(index.y,levels = c('age_up','nochange','age_down')))]
setorder(p.table,index.x,index.y)
levels <- p.table$Accession
p.table <- melt(p.table, id.vars = 'Accession',value.name = 'class',variable.name = 'loci')
p.table[,class:= as.character(class)]
test <- p.table[Accession %in% core.pro$Accession.hum]
p.table[Accession %in% core.pro$Accession.hum, class := 'core'][,Accession:= factor(Accession,levels = rev(levels))][,class:=factor(class,levels = c('age_up','core','nochange','age_down'))]

plot <- ggplot(p.table,aes(x=loci,y=Accession))
plot + geom_tile(aes(fill = class)) +
  scale_fill_manual(values = c('firebrick2','steelblue','gray95','forestgreen'),
                    na.value = 'gray50',guide=F) +
  scale_y_discrete(expand = c(0,0))+
  scale_x_discrete(expand = c(0,0))+
  labs(x =NULL, y =NULL)+
  theme(legend.position = "top",
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      size = rel(0.8),
      angle = 330,
      hjust = 0,
      colour = "grey50"
    ))
ggsave('test.tiff')






