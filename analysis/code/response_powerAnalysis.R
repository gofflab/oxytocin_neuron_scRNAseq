#code for response including power analysis 
---
title: "Dolen_OT_neurons_response"
author: "Genevieve Stein-O'Brien"
date: "03/08/2020"
output: html_document
---

# Initialize Environment
```{r init}
source('code/init.R')
library(reticulate)
use_condaenv('OT-reticulate')
set.seed(1234)
library("AnnotationDbi")
library("org.Mm.eg.db")
library("biomaRt")
library("RColorBrewer")
library(corrplot)
library(monocle)
library('ncells')
library(reshape2)
library(tidyverse)
library(knitr)
library(Biobase)
library(edgeR)
```

```{r load_data}
dat<-readRDS(file="../data/filtered10kOxy.rds")

#call cell types by clusters
pData(dat)$CellType<-"Magnocellular"
pData(dat)$CellType[pData(dat)$Cluster==1]<-"Parvocellular"
PMcolors<-c("orangered","darkcyan")

nCellCutoff<-5
expressed_genes<-rownames(fData(dat))[fData(dat)$num_cells_expressed>=nCellCutoff]

PV_vs_Magno_res<-readRDS("PV_vs_Magno_res.rds")
qval_cutoff<-0.001
PV_vs_Magno_sigGeneIDs<-PV_vs_Magno_res$gene_id[PV_vs_Magno_res$qval<=qval_cutoff]
length(PV_vs_Magno_sigGeneIDs)

qval_cutoff<-0.001
sigGeneIDs<-PV_vs_Magno_res$gene_id[PV_vs_Magno_res$qval<=qval_cutoff]
length(sigGeneIDs)

magExp<-apply(exprs(dat)[sigGeneIDs,pData(dat)$CellType=="Magnocellular"],1,mean)
pavExp<-apply(exprs(dat)[sigGeneIDs,pData(dat)$CellType=="Parvocellular"],1,mean)

("Kcnmb4","Calb1","Cnr1","Reln")

fc<-magExp/pavExp
str(fData(dat))
tmp<-fData(dat)[names(fc),"gene_short_name"]
names(fc)<-tmp

fc[c("Kcnmb4","Calb1","Cnr1","Reln")]

fc[c("Calb1", "Kcnmb4", "Reln", "Cnr1")]
delta<-magExp-pavExp

CT<-ifelse(delta<0,"Parvo","Magno")

fData(dat)[fData(dat)$gene_short_name=="Oxt",]

length(which(exprs(dat)[fData(dat)$gene_short_name=="Oxt",]>0))
dim(exprs(dat))
```

```{r response}

fData(dat)$isSignificant<-fData(dat)$gene_short_name %in% lookupGeneName(dat,PV_vs_Magno_sigGeneIDs)


#IUPHAR.genes<-read.delim("IUPHAR_genes.txt",header=F,stringsAsFactors=F)$V1
#IUPHAR.genes<-read.delim("IUPHAR_genes.txt",header=F,stringsAsFactors=F)$V1
All_channels<-c("Asic1","Asic2","Asic3","Asic4","Asic5","Ano1","Ano2","Ano3","Ano4","Ano5","Ano6","Ano7","Ano8","Ano9","Ano10","Mip","Aqp1","Aqp2","Aqp3","Aqp4","Aqp5","Aqp6","Aqp7","Aqp8","Aqp9","Aqp10","Aqp11","Aqp12a","Aqp12b","Best1","Best2","Best3","Best4","Cacna1a","Cacna1b","Cacna1c","Cacna1d","Cacna1e","Cacna1f","Cacna1g","Cacna1h","Cacna1i","Cacna1s","Cacna2d1","Cacna2d2","Cacna2d3","Cacna2d4","Cacnb1","Cacnb2","Cacnb3","Cacnb4","Cacng1","Cacng2","Cacng3","Cacng4","Cacng5","Cacng6","Cacng7","Cacng8","Catsperb","Catsperd","Catsperg","Catsper1","Catsper2","Catsper3","Catsper4","Cftr","Clic1","Clic2","Clic3","Clic4","Clic5","Clic6","Clcnka","Clcnkb","Bsnd","Clcn1","Clcn2","Clcn3","Clcn4","Clcn5","Clcn6","Clcn7","Gja1","Gja3","Gja4","Gja5","Gja6p","Gja8","Gja9","Gja10","Gjb1","Gjb2","Gjb3","Gjb4","Gjb5","Gjb6","Gjb7","Gjc1","Gjc2","Gjc3","Gjd2","Gjd3","Gjd4","Gje1","Itpr1","Itpr2","Itpr3","Kcnma1","Kcnn1","Kcnn2","Kcnn3","Kcnn4","Kcnu1","Kcnt1","Kcnt2","Kcnk1","Kcnk2","Kcnk3","Kcnk4","Kcnk5","Kcnk6","Kcnk7","Kcnk9","Kcnk10","Kcnk12","Kcnk13","Kcnk15","Kcnk16","Kcnk17","Kcnk18","Kcnj1","Kcnj2","Kcnj3","Kcnj4","Kcnj5","Kcnj6","Kcnj8","Kcnj9","Kcnj10","Kcnj11","Kcnj12","Kcnj13","Kcnj14","Kcnj15","Kcnj16","Kcnj18","Kcna1","Kcna2","Kcna3","Kcna4","Kcna5","Kcna6","Kcna7","Kcna10","Kcnb1","Kcnb2","Kcnc1","Kcnc2","Kcnc3","Kcnc4","Kcnd1","Kcnd2","Kcnd3","Kcnf1","Kcng1","Kcng2","Kcng3","Kcng4","Kcnh1","Kcnh2","Kcnh3","Kcnh4","Kcnh5","Kcnh6","Kcnh7","Kcnh8","Kcnq1","Kcnq2","Kcnq3","Kcnq4","Kcnq5","Kcns1","Kcns2","Kcns3","Kcnv1","Kcnv2","Ryr1","Ryr2","Ryr3","Scnn1a","Scnn1b","Scnn1d","Scnn1g","Nalcn","Scn1a","Scn2a","Scn3a","Scn4a","Scn5a","Scn8a","Scn9a","Scn10a","Scn11a","Scn1b","Scn2b","Scn3b","Scn4b","Tpcn1","Tpcn2","Vdac1","Vdac2","Vdac3","Dpp6","Dpp10","Kcnab1","Kcnab2","Kcnab3","Kcne1","Kcne1b","Kcne2","Kcne3","Kcne4","Kcne5","Kcnip1","Kcnip2","Kcnip3","Kcnip4","Kcnmb1","Kcnmb2","Kcnmb3","Kcnmb4")

Synaptic_transmission_genes<-c("2610042L04Rik","9530002B09Rik","Abat","Abhd6","Adcyap1","Adgrl1","Adipoq","Adnp","Adora1","Adora2a","Adra1a","Adrb1","Adrb2","Agrn","Agt","Als2","Anapc2","Apba1","Apba2","Apba3","Arc","Arf1","Arrb2","Asic1","Atad1","Atp2a2","Atp2b2","Atxn1","Baiap2","Bche","Bdnf","Braf","Brsk1","Btbd9","Cacna1a","Cacna1b","Cacna1c","Cacna1g","Cacna2d2","Cacnb1","Cacnb2","Cacnb3","Cacnb4","Cacng2","Cacng3","Cacng4","Cacng5","Cacng7","Cacng8","Cadps","Cadps2","Calb1","Camk2a","Camk2b","Camk4","Car2","Car7","Cartpt","Cav2","Cckbr","Ccl2","Cd24a","Cd38","Cdc20","Cdh8","Cdk5","Cel","Celf4","Celsr1","Chat","Chrm1","Chrm2","Chrm3","Chrm4","Chrm5","Chrna1","Chrna2","Chrna3","Chrna4","Chrna5","Chrna6","Chrna7","Chrna9","Chrna10","Chrnb1","Chrnb2","Chrnb3","Chrnb4","Chrnd","Chrne","Chrng","Clcn3","Clstn1","Clstn2","Clstn3","Cnih2","Cnih3","Cnr1","Cnr2","Cntn2","Cntnap4","Cpeb1","Cpeb3","Cplx1","Cplx2","Cplx3","Cplx4","Creb1","Crh","Crhbp","Crhr1","Crhr2","Crtc1","Cspg5","Ctnnb1","Ctnnd2","Cux2","Cx3cl1","Cx3cr1","D730048I06Rik","Dapk1","Dbi","Dbn1","Dgki","Dkk1","Dlg1","Dlg2","Dlg3","Dlg4","Dlgap1","Dlgap2","Dmpk","Dnm1","Doc2a","Doc2b","Doc2g","Drd1","Drd2","Drd3","Drd4","Drd5","Dtnbp1","Dvl1","Edn1",
  "Edn3","Egfr","Egr1","Egr3","Eif2ak4","Eif4a3","Eif4ebp2","Ephb2","Erbb4","Etv5","Exoc4","Fgf12","Fgf14","Flot1","Fmr1","Gabbr1","Gabra1","Gabra3","Gabra5","Gabra6","Gabrb2","Gabrb3","Gabrd","Gabrg1","Gabrg2","Gdnf","Gfap","Ghrl","Gip","Gipc1","Gjc1","Gjd2","Glra1","Glra2","Glra3","Glra4","Glrb","Gls","Glul","Gm12886","Gm12887","Gm12888","Gnai1","Gnai2","Gpm6b","Gpr1","Gpr21","Gpr52","Gpr149","Gria1","Gria2","Gria4","Grid2","Grid2ip","Grik1","Grik2","Grik3","Grik5","Grin1","Grin2a","Grin2b","Grin2c","Grin2d","Grm1","Grm2","Grm3","Grm4","Grm5","Grm6","Grm7","Grm8","Gsg1l","Gsk3b","Hap1","Hcrt","Hcrtr1","Hcrtr2","Hnrnpk","Hras","Hrh1","Hrh2","Hrh3","Hrh4","Htr1b","Htr1d","Htr2a","Htr2c","Htr4","Htr6","Htr7","Htt","Ica1","Ifng","Igsf9b","Iqsec2","Itpka","Itpr3","Jph3","Jph4","Kalrn","Kcnb1","Kcnc3","Kcnc4","Kcnj10","Kcnma1","Kcnmb1","Kcnmb2","Kcnmb3","Kcnmb4","Kcnn2","Kcnq4","Kdr","Kif1b","Kif5b","Kiss1","Kiss1r","Kit","Klhl24","Kmt2a","Kras","Lama2","Lgi1","Lin7a","Lin7b","Lin7c","Lrp6","Lrp8","Lrrk2","Lrrtm1","Lrrtm2","Lynx1","Lypd1","Lzts1","Mapk1","Mapk8ip2","Mapt","Mecp2","Mef2c","Met","Mgll","Mink1","Mmp9","Mtmr2","Musk","Mylk2","Myo5a","Myo6","Napa","Napb","Nat8l","Ncam1","Ncdn","Ncs1","Neto1","Neto2","Neurl1a","Neurod2","Nf1","Nfatc4","Ngf","Ngfr","Nlgn1","Nlgn2","Nlgn3","Nlgn4l","Nmu","Nos1","Npas4","Npbwr1","Npff","Npffr1","Npffr2","Nps","Nptn","Npy2r","Npy5r","Nr2e1",
  "Nr3c1","Nrgn","Nrxn1","Nrxn2","Nrxn3","Nsmf","Ntf3","Ntrk1","Ntrk2","Ntsr1","Omp","Ophn1","Oprd1","Oprk1","Oprl1","Oprm1","Otof","Oxt","Oxtr","P2rx1","P2rx2","P2rx3","P2rx4","P2rx7","Pafah1b1","Paip2","Pak1","Park2","Park7","Pate4","Pcdh8","Pcdh17","Pcdhb16","Pclo","Pdyn","Pdzd11","Pebp1","Penk","Pfn2","Pick1","Pink1","Pla2g6","Plat","Plcb3","Plcl1","Plcl2","Plk2","Pmch","Pnkd","Pnoc","Ppargc1a","Ppfia3","Ppp1r9a","Ppp1r9b","Ppp3ca","Ppt1","Prkaca","Prkce","Prkcg","Prkcz","Psen1","Psen2","Pten","Ptgdr2","Ptger4","Ptgs2","Ptk2","Ptk2b","Ptpn5","Ptprn2","Rab3a","Rab3b","Rab3gap1","Rab5a","Rab8a","Rab11a","Rac1","Rac3","Rap1a","Rap1b","Rapgef2","Rapgef4","Rapsn","Rara","Rasd2","Rasgrf1","Rasgrf2","Reln","Retn","Rgs14","Ric3","Rims1","Rims2","Rims3","Rims4","Rin1","Rph3a","S1pr2","S100b","Scrib","Sdcbp","Serpine2","Sez6","Shank1","Shank2","Shank3","Shc3","Shisa6","Shisa7","Shisa8","Shisa9","Shisa9","Sipa1l1","Slc1a3","Slc5a7","Slc6a1","Slc6a3","Slc6a4","Slc6a5","Slc6a9","Slc6a11","Slc6a12","Slc6a13","Slc6a20a","Slc6a20b","Slc8a2","Slc8a3","Slc12a4","Slc12a5","Slc12a6","Slc12a7","Slc17a7","Slc18a1","Slc24a2","Slc29a1","Slc30a1","Slitrk5","Snap23","Snap25","Snap29","Snap47","Snap91",
  "Snapin","Snca","Sncaip","Sncb","Sncg","Sorcs3","Spg11","Sphk1","Sptbn2","Srf","Sstr1","Sstr2","Sstr3","Sstr4","Sstr5","Stac3","Star","Stau1","Stau2","Ston2","Stx1a","Stx1b","Stx2","Stx3","Stx4a","Stx11","Stx19","Stxbp1","Sv2a","Sv2b","Sv2c","Syde1","Syn1","Syn2","Syn3","Syngap1","Syngr1","Synj1","Syp","Syt1","Syt2","Syt3","Syt4","Syt5","Syt6","Syt7","Syt8","Syt9","Syt10","Syt11","Syt12","Syt13",
  "Syt15","Syt17","Sytl1","Sytl3","Sytl4","Sytl5","Tac1","Tacr1","Tacr2","Th","Tmod2","Tnf","Tnr","Tor1a","Tpgs1","Trim9","Ucn","Unc13a","Unc13b","Unc13c","Usp14","Usp46","Uts2","Vamp2","Vdac1","Vdac3","Vgf","Wnt7a","Xbp1","Ywhag","Zmynd8")

fData(dat)$Synaptic_transmission_genes<-fData(dat)$gene_short_name %in% Synaptic_transmission_genes
indxSTG<-fData(dat)[fData(dat)$stg & fData(dat)$isSignificant,"gene_short_name"]
length(indxSTG) #0


Amine_receptors<-c("Htr1a","Htr1b","Htr1d","Htr1e","Htr1f","Htr2a","Htr2b","Htr2c","Htr4","Htr5a","Htr5bp","Htr6","Htr7","Adra1a","Adra1b","Adra1d","Adra2a","Adra2b","Adra2c","Adrb1","Adrb2","Adrb3","Chrm1","Chrm2","Chrm3","Chrm4","Chrm5","Drd1","Drd2","Drd3","Drd4","Drd5","Hrh1","Hrh2","Hrh3","Hrh4","Taar1","Taar2","Taar3","Taar4p","Taar5","Taar6","Taar7p","Taar8","Taar9")
Peptide_receptors<-c("Aplnr","Ghsr","Kiss1r","Mlnr","Prlhr","Qrfpr","Trhr","Agtr1","Agtr2","Avpr1a","Avpr1b","Avpr2","Oxtr","Bdkrb1","Bdkrb2","Cckar","Cckbr","Ednra","Ednrb","Galr1","Galr2","Galr3","Gnrhr","Gnrhr2","Hcrtr1","Hcrtr2","Mchr1","Mchr2","Mc1r","Mc2r","Mc3r","Mc4r","Mc5r","Nmur1","Nmur2","Npbwr1","Npbwr2","Npffr1","Npffr2","Npsr1","Npy1r","Npy2r","Npy4r","Npy5r","Npy6r","Ntsr1","Ntsr2","Oprd1","Oprk1","Oprl1","Oprm1","Sstr1","Sstr2","Sstr3","Sstr4","Sstr5","Tacr1","Tacr2","Tacr3")
VIP_calcitonin<-c("Adcyap1r1","Vipr1","Vipr2","Calcr","Calcrl","Crhr1","Crhr2","Casr","Gprc6a")
Other_receptors<-c(Amine_receptors,Peptide_receptors,VIP_calcitonin,"Pgrmc1","Pgrmc2","Nenf","Cyb5d2","Glp2r")

fData(dat)$Other_receptors<-fData(dat)$gene_short_name %in% Other_receptors
indxOR<-fData(dat)[fData(dat)$GPCRs & fData(dat)$isSignificant,"gene_short_name"]
length(indxOR) #6


GPCRs<-c("Htr1a","Htr1b","Htr1d","Htr1e","Htr1f","Htr2a","Htr2b","Htr2c","Htr4","Htr5a","Htr5bp","Htr6","Htr7","Gpr107","Gpr108","Gpr137","Gpr143","Gpr157","Tpra1","Adora1","Adora2a","Adora2b","Adora3","Adgra1","Adgra2","Adgra3","Adgrb1","Adgrb2","Adgrb3","Celsr1","Celsr2","Celsr3","Adgrd1","Adgrd2","Adgre1","Adgre2","Adgre3","Adgre4p","Adgre5","Adgrf1","Adgrf2","Adgrf3","Adgrf4","Adgrf5","Adgrg1","Adgrg2","Adgrg3","Adgrg4","Adgrg5","Adgrg6","Adgrg7","Adgrl1","Adgrl2","Adgrl3","Adgrl4","Adgrv1","Adra1a","Adra1b","Adra1d","Adra2a","Adra2b","Adra2c","Adrb1","Adrb2","Adrb3","Agtr1","Agtr2","Avpr1a","Avpr1b","Avpr2","Oxtr","Ackr1","Ackr2","Ackr3","Ackr4","Ccrl2","Pitpnm3","Bdkrb1","Bdkrb2","Ccr1","Ccr10","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccr9","Cx3cr1","Cxcr1","Cxcr2","Cxcr3"
         ,"Cxcr4","Cxcr5","Cxcr6","Calcr","Calcrl","Casr","Gprc6a","Cnr1","Cnr2","Cmklr1","Cckar","Cckbr","Chrm1","Chrm2","Chrm3","Chrm4","Chrm5","C3ar1","C5ar1","C5ar2","Crhr1","Crhr2","Drd1","Drd2","Drd3","Drd4","Drd5","Ednra","Ednrb","F2r","F2rl1","F2rl2","F2rl3","Fpr1","Fpr2","Fpr3","Ffar1","Ffar2","Ffar3","Ffar4","Gpbar1","Gper1","Gpr1","Gpr101","Gpr119","Gpr12","Gpr132","Gpr135","Gpr139","Gpr141","Gpr142","Gpr146","Gpr148","Gpr149","Gpr15","Gpr150","Gpr151","Gpr152","Gpr153","Gpr160","Gpr161","Gpr162","Gpr17","Gpr171","Gpr173","Gpr174","Gpr176","Gpr18","Gpr182","Gpr183","Gpr19","Gpr20","Gpr21","Gpr22","Gpr25","Gpr26","Gpr27","Gpr3","Gpr31","Gpr32","Gpr33","Gpr34","Gpr35","Gpr37","Gpr37l1","Gpr39","Gpr4","Gpr42","Gpr45","Gpr50","Gpr52","Gpr55","Gpr6","Gpr61","Gpr62","Gpr63","Gpr65","Gpr68","Gpr75","Gpr78","Gpr79","Gpr82","Gpr83","Gpr84","Gpr85","Gpr87","Gpr88","Lgr4","Lgr5","Lgr6","Mas1","Mas1l","Mrgprd","Mrgpre","Mrgprf","Mrgprg","Mrgprx1","Mrgprx2","Mrgprx3","Mrgprx4","Gpr156","Gpr158","Gpr179","Gprc5a","Gprc5b","Gprc5c","Gprc5d","Fzd1","Fzd10","Smo","Fzd2","Fzd3","Fzd4","Fzd5","Fzd6","Fzd7","Fzd8","Fzd9","Galr1","Galr2","Galr3","Gabbr1","Gabbr2","Gcgr","Ghrhr","Gipr","Glp1r","Glp2r","Sctr","Grm1","Grm2","Grm3","Grm4","Grm5","Grm6","Grm7","Grm8","Fshr","Lhcgr","Tshr","Gnrhr","Gnrhr2","Hrh1","Hrh2","Hrh3","Hrh4","Hcar1","Hcar2","Hcar3","Hcrtr1","Hcrtr2","Cysltr1","Cysltr2","Fpr2","Ltb4r","Ltb4r2","Oxer1","Lpar1","Lpar2",
         "Lpar3","Lpar4","Lpar5","Lpar6","Mchr1","Mchr2","Mc1r","Mc2r","Mc3r","Mc4r","Mc5r","Mtnr1a","Mtnr1b","Nmur1","Nmur2","Npbwr1","Npbwr2","Npffr1","Npffr2","Npsr1","Npy1r","Npy2r","Npy4r","Npy5r","Npy6r","Ntsr1","Ntsr2","Or1a1","Or1a2","Or1aa1p","Or1ab1p","Or1ac1p","Or1b1","Or1c1","Or1d2","Or1d3p","Or1d4","Or1d5","Or1e1","Or1e2","Or1e3","Or1f1","Or1f12","Or1f2p","Or1g1","Or1h1p","Or1i1","Or1j1","Or1j2","Or1j4","Or1k1","Or1l1","Or1l3","Or1l4","Or1l6","Or1l8","Or1m1","Or1m4p","Or1n1","Or1n2","Or1p1","Or1q1","Or1r1p","Or1s1","Or1s2","Or1x1p","Or1x5p","Or10a2","Or10a3","Or10a4","Or10a5","Or10a6","Or10a7","Or10aa1p","Or10ab1p","Or10ac1","Or10ad1","Or10ae1p","Or10ae3p","Or10af1p","Or10ag1","Or10ah1p","Or10ak1p","Or10b1p","Or10c1","Or10d1p","Or10d3","Or10d4p","Or10d5p","Or10g1p","Or10g2","Or10g3","Or10g4","Or10g5p","Or10g6","Or10g7","Or10g8","Or10g9","Or10h1","Or10h2","Or10h3","Or10h4","Or10h5","Or10j1","Or10j2p","Or10j3","Or10j4","Or10j5","Or10j6p","Or10j7p","Or10j8p","Or10j9p","Or10k1","Or10k2","Or10n1p","Or10p1","Or10q1","Or10q2p","Or10r1p","Or10r2","Or10r3p","Or10s1","Or10t1p","Or10t2","Or10u1p","Or10v1","Or10v2p","Or10v3p","Or10v7p","Or10w1","Or10x1","Or10y1p","Or10z1","Or11a1","Or11g1p","Or11g2","Or11h1","Or11h12","Or11h13p","Or11h2","Or11h3p","Or11h4","Or11h5p","Or11h6","Or11h7","Or11i1p","Or11j1p","Or11j2p","Or11j5p","Or11k1p","Or11k2p","Or11l1","Or11m1p","Or11n1p","Or11p1p","Or11q1p","Or12d1","Or12d2","Or12d3","Or13a1"
         ,"Or13c1p","Or13c2","Or13c3","Or13c4","Or13c5","Or13c6p","Or13c7","Or13c8","Or13c9","Or13d1","Or13d2p","Or13d3p","Or13e1p","Or13f1","Or13g1","Or13h1","Or13i1p","Or13j1","Or13k1p","Or13z1p","Or13z2p","Or13z3p","Or14a16","Or14a2","Or14c36","Or14i1","Or14j1","Or14k1","Or14l1p","Or2a1","Or2a12","Or2a13p","Or2a14","Or2a15p","Or2a2","Or2a20p","Or2a25","Or2a3p","Or2a4","Or2a41p","Or2a42","Or2a5","Or2a7","Or2a9p","Or2ad1p","Or2ae1","Or2af1p","Or2ag1","Or2ag2","Or2ah1p","Or2ai1p","Or2aj1","Or2ak2","Or2al1p","Or2am1p","Or2ao1p","Or2ap1","Or2aq1p","Or2as1p","Or2as2p","Or2at1p","Or2at2p","Or2at4","Or2b11","Or2b2","Or2b3","Or2b4p","Or2b6","Or2b7p","Or2b8p","Or2bh1p","Or2c1","Or2c3","Or2d2","Or2d3","Or2e1p","Or2f1","Or2f2","Or2g1p","Or2g2","Or2g3","Or2g6","Or2h1","Or2h2","Or2h4p","Or2h5p","Or2i1p","Or2j1","Or2j2","Or2j3","Or2j4p","Or2k2","Or2l13","Or2l1p","Or2l2","Or2l3","Or2l5","Or2l6p","Or2l8","Or2l9p","Or2m1p","Or2m2","Or2m3","Or2m4","Or2m5","Or2m7","Or2n1p","Or2p1p","Or2q1p","Or2r1p","Or2s1p","Or2s2","Or2t1","Or2t10","Or2t11","Or2t12","Or2t2","Or2t27","Or2t29","Or2t3","Or2t32p","Or2t33","Or2t34","Or2t35","Or2t4","Or2t5","Or2t6","Or2t7","Or2t8","Or2u1p","Or2u2p","Or2v1","Or2v2","Or2w1","Or2w2p","Or2w3","Or2w4p","Or2w5","Or2w6p","Or2x1p","Or2y1","Or2z1","Or3a1","Or3a2","Or3a3","Or3a4p","Or3b1p","Or3d1p","Or4a10p","Or4a11p","Or4a12p","Or4a13p","Or4a14p","Or4a15","Or4a16","Or4a17p","Or4a18p","Or4a19p","Or4a1p","Or4a21p","Or4a2p","Or4a3p","Or4a40p","Or4a41p","Or4a42p","Or4a43p","Or4a44p","Or4a45p","Or4a46p","Or4a47","Or4a48p","Or4a49p","Or4a4p","Or4a5","Or4a50p","Or4a6p","Or4a7p","Or4a8",
         "Or4a9p","Or4b1","Or4b2p","Or4c10p","Or4c11","Or4c12","Or4c13","Or4c14p","Or4c15","Or4c16","Or4c1p","Or4c2p","Or4c3","Or4c45","Or4c46","Or4c48p","Or4c49p","Or4c4p","Or4c5","Or4c50p","Or4c6","Or4c7p","Or4c9p","Or4d1","Or4d10","Or4d11","Or4d12p","Or4d2","Or4d5","Or4d6","Or4d7p","Or4d8p","Or4d9","Or4e1","Or4e2","Or4f13p","Or4f14p","Or4f15","Or4f16","Or4f17","Or4f1p","Or4f21","Or4f28p","Or4f29","Or4f2p","Or4f3","Or4f4","Or4f5","Or4f6","Or4f7p","Or4f8p","Or4g11p","Or4g1p","Or4g2p","Or4g3p","Or4g4p","Or4g6p","Or4h12p","Or4h6p","Or4k1","Or4k11p","Or4k12p","Or4k13","Or4k14","Or4k15","Or4k16p","Or4k17","Or4k2","Or4k3","Or4k4p","Or4k5","Or4k6p","Or4k7p","Or4k8p","Or4l1","Or4m1","Or4m2","Or4n1p","Or4n2","Or4n3p","Or4n4","Or4n5","Or4p1p","Or4p4","Or4q1p","Or4q2","Or4q3","Or4r1p","Or4r2p","Or4r3p","Or4s1","Or4s2","Or4t1p","Or4u1p","Or4v1p","Or4w1p","Or4x1","Or4x2","Or4x7p","Or5a1","Or5a2","Or5ac1","Or5ac2","Or5ac4p","Or5ah1p","Or5ak1p","Or5ak2","Or5ak3p","Or5ak4p","Or5al1","Or5al2p","Or5am1p","Or5an1","Or5an2p","Or5ao1p","Or5ap1p","Or5ap2","Or5aq1p","Or5ar1","Or5as1","Or5au1","Or5aw1p","Or5az1p","Or5b10p","Or5b12","Or5b15p","Or5b17","Or5b19p","Or5b1p","Or5b2","Or5b21","Or5b3","Or5ba1p","Or5bb1p","Or5bc1p","Or5bd1p","Or5be1p","Or5bh1p","Or5bj1p","Or5bk1p","Or5bl1p","Or5bm1p"
         ,"Or5bn1p","Or5bn2p","Or5bp1p","Or5bq1p","Or5br1p","Or5bs1p","Or5bt1p","Or5c1","Or5d13","Or5d14","Or5d15p","Or5d16","Or5d17p","Or5d18","Or5d2p","Or5d3p","Or5e1p","Or5f1","Or5f2p","Or5g1p","Or5g3","Or5g4p","Or5g5p","Or5h1","Or5h14","Or5h15","Or5h2","Or5h3p","Or5h4p","Or5h5p","Or5h6","Or5h7p","Or5h8","Or5i1","Or5j1p","Or5j2","Or5j7p","Or5k1","Or5k2","Or5k3","Or5k4","Or5l1","Or5l2","Or5m1","Or5m10","Or5m11","Or5m12p","Or5m13p","Or5m14p","Or5m2p","Or5m3","Or5m4p","Or5m5p","Or5m6p","Or5m7p","Or5m8","Or5m9","Or5p1p","Or5p2","Or5p3","Or5p4p","Or5r1","Or5s1p","Or5t1","Or5t2","Or5t3","Or5v1","Or5w1p","Or5w2","Or51a10p","Or51a1p","Or51a2","Or51a3p","Or51a4","Or51a5p","Or51a6p","Or51a7","Or51a8p","Or51a9p","Or51ab1p","Or51b2","Or51b3p","Or51b4","Or51b5","Or51b6","Or51b8p","Or51c1p","Or51c4p","Or51d1","Or51e1","Or51e2","Or51f1","Or51f2","Or51f3p","Or51f4p","Or51f5p","Or51g1","Or51g2","Or51h1","Or51h2p","Or51i1","Or51i2","Or51j1","Or51k1p","Or51l1","Or51m1","Or51n1p","Or51p1p","Or51q1","Or51r1p","Or51s1","Or51t1","Or51v1","Or52a1","Or52a4p","Or52a5","Or52b1p","Or52b2","Or52b3p","Or52b4","Or52b5p","Or52b6","Or52d1","Or52e1","Or52e2","Or52e3p","Or52e4","Or52e5","Or52e6","Or52e7p","Or52e8","Or52h1","Or52h2p","Or52i1","Or52i2","Or52j1p","Or52j2p","Or52j3","Or52k1","Or52k2","Or52k3p","Or52l1","Or52l2p","Or52m1","Or52m2p","Or52n1","Or52n2","Or52n3p","Or52n4","Or52n5","Or52p1p","Or52p2p","Or52q1p","Or52r1","Or52s1p","Or52t1p","Or52u1p",
         "Or52v1p","Or52w1","Or52x1p","Or52y1p","Or52z1","Or55b1p","Or56a1","Or56a3","Or56a4","Or56a5","Or56a7p","Or56b1","Or56b2p","Or56b3p","Or56b4","Or6a2","Or6b1","Or6b2","Or6b3","Or6c1","Or6c2","Or6c3","Or6c4","Or6c5p","Or6c6","Or6c64p","Or6c65","Or6c66p","Or6c68","Or6c69p","Or6c70","Or6c71p","Or6c72p","Or6c73p","Or6c74","Or6c75","Or6c76","Or6c7p","Or6d1p","Or6e1p","Or6f1","Or6j1","Or6k1p","Or6k2","Or6k3","Or6k4p","Or6k5p","Or6k6","Or6l1p","Or6l2p","Or6m1","Or6m2p","Or6m3p","Or6n1","Or6n2","Or6p1","Or6q1","Or6r1p","Or6r2p","Or6s1","Or6t1","Or6u2p","Or6v1","Or6w1p","Or6x1","Or6y1","Or7a10","Or7a11p","Or7a15p","Or7a17","Or7a18p","Or7a19p","Or7a1p","Or7a2p","Or7a3p","Or7a5","Or7a8p","Or7c1","Or7c2","Or7d11p","Or7d1p","Or7d2","Or7d4","Or7e100p","Or7e101p","Or7e102p","Or7e104p","Or7e105p","Or7e106p","Or7e108p","Or7e109p","Or7e10p","Or7e110p","Or7e111p","Or7e115p","Or7e116p","Or7e117p","Or7e11p","Or7e121p","Or7e122p","Or7e125p","Or7e126p","Or7e128p","Or7e129p","Or7e12p","Or7e130p","Or7e136p","Or7e13p","Or7e140p","Or7e145p","Or7e148p","Or7e149p","Or7e14p","Or7e154p","Or7e155p","Or7e156p","Or7e157p","Or7e158p","Or7e159p","Or7e15p","Or7e160p","Or7e161p","Or7e162p","Or7e163p","Or7e16p","Or7e18p","Or7e19p","Or7e1p","Or7e21p","Or7e22p","Or7e23p","Or7e24","Or7e25p","Or7e26p","Or7e28p","Or7e29p","Or7e2p","Or7e31p","Or7e33p","Or7e35p","Or7e36p","Or7e37p","Or7e38p","Or7e39p","Or7e41p","Or7e43p","Or7e46p","Or7e47p","Or7e4p","Or7e53p","Or7e55p","Or7e59p","Or7e5p","Or7e62p","Or7e66p","Or7e7p","Or7e83p","Or7e84p","Or7e85p","Or7e86p","Or7e87p","Or7e89p","Or7e8p","Or7e90p","Or7e91p","Or7e93p","Or7e94p","Or7e96p","Or7e97p","Or7e99p","Or7g1","Or7g15p","Or7g2","Or7g3","Or7h1p","Or7h2p","Or7k1p","Or7l1p","Or7m1p","Or8a1","Or8a2p","Or8a3p","Or8b10p","Or8b12","Or8b1p","Or8b2","Or8b3","Or8b4","Or8b5p","Or8b6p","Or8b7p","Or8b8","Or8b9p","Or8c1p","Or8d1","Or8d2","Or8d4","Or8f1p","Or8g1","Or8g2","Or8g3p","Or8g5","Or8g7p","Or8h1","Or8h2","Or8h3","Or8i1p","Or8i2","Or8i4p","Or8j1","Or8j2","Or8j3","Or8k1","Or8k2p","Or8k3","Or8k4p","Or8k5","Or8l1p","Or8q1p","Or8r1p","Or8s1","Or8s21p","Or8t1p","Or8u1","Or8u8","Or8u9","Or8v1p","Or8x1p","Or9a1p","Or9a2","Or9a3p","Or9a4","Or9g1","Or9g2p","Or9g3p","Or9g4","Or9g9","Or9h1p","Or9i1","Or9i2p","Or9i3p"
         ,"Or9k1p","Or9k2","Or9l1p","Or9m1p","Or9n1p","Or9p1p","Or9q1","Or9q2","Or9r1p","Or9s24p","Oprd1","Oprk1","Oprl1","Oprm1","Opn1lw","Opn1mw","Opn1mw2","Opn1mw3","Opn1sw","Rho","Opn3","Opn4","Opn5","Rgr","Rrh","Oxgr1","Pth1r","Pth2r","Aplnr","Ghsr","Kiss1r","Mlnr","Prlhr","Qrfpr","Trhr","Ptafr","Prokr1","Prokr2","Ptgdr","Ptgdr2","Ptger1","Ptger2","Ptger3","Ptger4","Ptgfr","Ptgir","Tbxa2r","P2ry1","P2ry10","P2ry11","P2ry12","P2ry13","P2ry14","P2ry2","P2ry4","P2ry6","P2ry8","Rxfp1","Rxfp2","Rxfp3","Rxfp4","Sstr1","Sstr2","Sstr3","Sstr4","Sstr5","S1pr1","S1pr2","S1pr3","S1pr4","S1pr5","Sucnr1","Tacr1","Tacr2","Tacr3","Tas1r1","Tas1r2","Tas1r3","Tas2r1","Tas2r10","Tas2r12p","Tas2r13","Tas2r14","Tas2r15p","Tas2r16","Tas2r18p","Tas2r19","Tas2r20","Tas2r22","Tas2r2p","Tas2r3","Tas2r30","Tas2r31","Tas2r33","Tas2r36","Tas2r37","Tas2r38","Tas2r39","Tas2r4","Tas2r40","Tas2r41","Tas2r42","Tas2r43","Tas2r45","Tas2r46","Tas2r5","Tas2r50","Tas2r60","Tas2r62p","Tas2r63p","Tas2r64p","Tas2r67p","Tas2r68p","Tas2r6p","Tas2r7","Tas2r8","Tas2r9","Taar1","Taar2","Taar3","Taar4p","Taar5","Taar6","Taar7p","Taar8","Taar9","Adcyap1r1","Vipr1","Vipr2","Vn1r1","Vn1r100p","Vn1r101p","Vn1r102p","Vn1r103p","Vn1r104p","Vn1r105p","Vn1r106p","Vn1r107p","Vn1r108p","Vn1r109p","Vn1r10p","Vn1r110p","Vn1r111p","Vn1r112p","Vn1r11p","Vn1r12p","Vn1r13p","Vn1r14p","Vn1r15p","Vn1r16p","Vn1r17p","Vn1r18p","Vn1r19p","Vn1r2","Vn1r20p","Vn1r21p","Vn1r22p","Vn1r23p","Vn1r24p","Vn1r25p","Vn1r26p","Vn1r27p","Vn1r28p","Vn1r29p","Vn1r3","Vn1r30p","Vn1r31p","Vn1r32p","Vn1r33p","Vn1r34p","Vn1r35p","Vn1r36p","Vn1r37p","Vn1r38p","Vn1r39p","Vn1r4","Vn1r40p","Vn1r41p","Vn1r42p","Vn1r43p","Vn1r44p","Vn1r45p","Vn1r46p","Vn1r47p","Vn1r48p","Vn1r49p","Vn1r5","Vn1r51p","Vn1r52p","Vn1r53p","Vn1r54p"
         ,"Vn1r55p","Vn1r56p","Vn1r57p","Vn1r58p","Vn1r59p","Vn1r60p","Vn1r61p","Vn1r62p","Vn1r63p","Vn1r64p","Vn1r65p","Vn1r66p","Vn1r67p","Vn1r68p","Vn1r69p","Vn1r6p","Vn1r70p","Vn1r71p","Vn1r72p","Vn1r73p","Vn1r74p","Vn1r75p","Vn1r76p","Vn1r77p","Vn1r78p","Vn1r79p","Vn1r7p","Vn1r80p","Vn1r81p","Vn1r82p","Vn1r83p","Vn1r84p","Vn1r85p","Vn1r86p","Vn1r87p","Vn1r88p","Vn1r89p","Vn1r8p","Vn1r90p","Vn1r91p","Vn1r92p","Vn1r93p","Vn1r94p","Vn1r95p","Vn1r96p","Vn1r97p","Vn1r98p","Vn1r99p","Vn1r9p","Vn2r10p","Vn2r11p","Vn2r12p"
         ,"Vn2r13p","Vn2r14p","Vn2r15p","Vn2r16p","Vn2r17p","Vn2r18p","Vn2r19p","Vn2r1p","Vn2r20p","Vn2r21p","Vn2r2p","Vn2r3p","Vn2r6p","Vn2r7p","Vn2r9p","Xcr1")


fData(dat)$GPCRs<-fData(dat)$gene_short_name %in% GPCRs
indxGPCRs<-fData(dat)[fData(dat)$GPCRs & fData(dat)$isSignificant,"gene_short_name"]
length(indxGPCRs) #6 

Neuropeptides<-c("Penk","Pomc","Pdyn","Pnoc","Avp","Oxt","Gast","Cck","Sst","Cort","Npvf","Npff","Npy","Ppy","Pyy","Prlh","Calca","Calcb","Iapp","Adm","Adm2","Nppa","Nppb","Nppc","Grp","Nmb","Edn1","Edn2","Edn3","Cgc","Sct","Vip","Adcyap1","Ghrh","Gip","Crh","Ucn","Ucn2","Ucn3","Uts2","Uts2b","Tac1","Tac2","Nms","Nmu","Kng1","Agt","Nts","Chga","Chgb","Scg2","Scg3","Scg5","Vgf","Mln","Ghrl","Gal","Galp","Gnrh1","Npb","Npw","Nps","Nxph1","Nxph2","Nxph3","Nxph4","Ins","Ins1","Ins2","Igf1","Igf2","Rln1","Rln3","Trh","Pthlh","Pmch","Hcrt","Cartpt","Agrp","Prl","Apln","Kiss1","Dbi","Cbln1","Cbln2","Cbln3","Cbln4","Lep","Adipoq","Nampt","Retn","Retnla","Retnlb","Retnlg","Nucb2","Ubl5")
Na_channels_Mouse<-c("Asic1","Asic2","Asic3","Asic4","Asic5","Scnn1a","Scnn1b","Scnn1d","Scnn1g","Nalcn","Scn1a","Scn2a1","Scn3a","Scn4a","Scn5a","Scn8a","Scn9a","Scn10a","Scn11a","Scn1b","Scn2b","Scn3b","Scn4b","Scnm1","Scn7a")
K_channels_Mouse<-c("Hcn1","Hcn2","Hcn3","Hcn4","Kcnma1","Kcnn1","Kcnn2","Kcnn3","Kcnn4","Kcnu1","Kcnt1","Kcnt2","Kcnk1","Kcnk2","Kcnk3","Kcnk4","Kcnk5","Kcnk6","Kcnk7","Kcnk9","Kcnk10","Kcnk12","Kcnk13","Kcnk15","Kcnk16","Kcnk17","Kcnk18","Kcnj1","Kcnj2","Kcnj3","Kcnj4","Kcnj5","Kcnj6","Kcnj8","Kcnj9","Kcnj10","Kcnj11","Kcnj12","Kcnj13","Kcnj14","Kcnj15","Kcnj16","Kcnj18","Kcna1","Kcna2","Kcna3","Kcna4","Kcna5","Kcna6","Kcna7","Kcna10","Kcnb1","Kcnb2","Kcnc1","Kcnc2","Kcnc3","Kcnc4","Kcnd1","Kcnd2","Kcnd3","Kcnf1","Kcng1","Kcng2","Kcng3","Kcng4","Kcnh1","Kcnh2","Kcnh3","Kcnh4","Kcnh5","Kcnh6","Kcnh7","Kcnh8","Kcnq1","Kcnq2","Kcnq3","Kcnq4","Kcnq5","Kcns1","Kcns2","Kcns3","Kcnv1","Kcnv2","Kcnab1","Kcnab2","Kcnab3","Kcne1","Kcne2","Kcne3","Kcne4","Kcne1l","Kcnip1","Kcnip2","Kcnip3","Kcnip4","Kcnmb1","Kcnmb2","Kcnmb3","Kcnmb4","Dpp6","Dpp10","Kcnrg","Kcne1l")
Voltage_gated_channels_HUGO<-c("Cacna1a","Cacna1b","Cacna1c","Cacna1d","Cacna1e","Cacna1f","Cacna1g","Cacna1h","Cacna1i","Cacna1s","Cacna2d1","Cacna2d2","Cacna2d3","Cacna2d4","Cacnb1","Cacnb2","Cacnb3","Cacnb4","Cacng1","Cacng2","Cacng3","Cacng4","Cacng5","Cacng6","Cacng7","Cacng8","Catsperb","Catsperd","Catsperg","Catsper1","Catsper2","Catsper3","Catsper4","Cnga1","Cnga2","Cnga3","Cnga4","Cngb1","Cngb3","Hcn1","Hcn2","Hcn3","Hcn4","Hvcn1","Kcnma1","Kcnn1","Kcnn2","Kcnn3","Kcnn4","Kcnu1","Kcnk1","Kcnk2","Kcnk3","Kcnk4","Kcnk5","Kcnk6","Kcnk7","Kcnk9","Kcnk10","Kcnk12","Kcnk13","Kcnk15","Kcnk16","Kcnk17","Kcnk18","Kcnj1","Kcnj2","Kcnj3","Kcnj4","Kcnj5","Kcnj6","Kcnj8","Kcnj9","Kcnj10","Kcnj11","Kcnj12","Kcnj13","Kcnj14","Kcnj15","Kcnj16","Kcnj18","Kcna1","Kcna2","Kcna3","Kcna4","Kcna5","Kcna6","Kcna7","Kcna10","Kcnb1","Kcnb2","Kcnc1","Kcnc2","Kcnc3","Kcnc4","Kcnd1","Kcnd2","Kcnd3","Kcnf1","Kcng1","Kcng2","Kcng3","Kcng4","Kcnh1","Kcnh2","Kcnh3","Kcnh4","Kcnh5","Kcnh6","Kcnh7","Kcnh8","Kcnq1","Kcnq2","Kcnq3","Kcnq4","Kcnq5","Kcns1","Kcns2","Kcns3","Kcnv1","Kcnv2","Scn1a","Scn2a","Scn3a","Scn4a","Scn5a","Scn8a","Scn9a","Scn10a","Scn11a","Scn1b","Scn2b","Scn3b","Scn4b","Trpa1","Trpc1","Trpc2","Trpc3","Trpc4","Trpc5","Trpc6","Trpc7","Mcoln1","Mcoln2","Mcoln3","Trpm1","Trpm2","Trpm3","Trpm4","Trpm5","Trpm6","Trpm7","Trpm8","Pkd2","Pkd2l1","Pkd2l2","Trpv1","Trpv2","Trpv3","Trpv4","Trpv5","Trpv6","Tpcn1","Tpcn2")
Ligand_gated_ion_channels_HUGO<-c("Htr3a","Htr3b","Htr3c","Htr3d","Htr3e","Cftr","Chrna1","Chrna2","Chrna3","Chrna4","Chrna5","Chrna6","Chrna7","Chrna9","Chrna10","Chrnb1","Chrnb2","Chrnb3","Chrnb4","Chrnd","Chrne","Chrng","Gabra1","Gabra2","Gabra3","Gabra4","Gabra5","Gabra6","Gabrb1","Gabrb2","Gabrb3","Gabrd","Gabre","Gabrg1","Gabrg2","Gabrg3","Gabrp","Gabrq","Gabrr1","Gabrr2","Gabrr3","Gria1","Gria2","Gria3","Gria4","Grid1","Grid2","Grik1","Grik2","Grik3","Grik4","Grik5","Grin1","Grin2a","Grin2b","Grin2c","Grin2d","Grin3a","Grin3b","Glra1","Glra2","Glra3","Glra4","Glrb","Itpr1","Itpr2","Itpr3","P2rx1","P2rx2","P2rx3","P2rx4","P2rx5","P2rx6","P2rx7","Ryr1","Ryr2","Ryr3","Zacn")

NNKVL<-c(Neuropeptides,Na_channels_Mouse,K_channels_Mouse,Voltage_gated_channels_HUGO,Ligand_gated_ion_channels_HUGO)

fData(dat)$NNKVL<-fData(dat)$gene_short_name %in% NNKVL
indxNNKVL<-fData(dat)[fData(dat)$NNKVL & fData(dat)$isSignificant,"gene_short_name"]
length(indxNNKVL)#0

### Ion Channels 
KChannels<-c("Kcnma1","Kcnn1","Kcnn2","Kcnn3","Kcnn4","Kcnu1","Kcnt1","Kcnt2","Kcnk1","Kcnk2","Kcnk3","Kcnk4","Kcnk5","Kcnk6","Kcnk7","Kcnk9","Kcnk10","Kcnk12","Kcnk13","Kcnk15","Kcnk16","Kcnk17","Kcnk18","Kcnj1","Kcnj2","Kcnj3","Kcnj4","Kcnj5","Kcnj6","Kcnj8","Kcnj9","Kcnj10","Kcnj11","Kcnj12","Kcnj13","Kcnj14","Kcnj15","Kcnj16","Kcnj18","Kcna1","Kcna2","Kcna3","Kcna4","Kcna5","Kcna6","Kcna7","Kcna10","Kcnb1","Kcnb2","Kcnc1","Kcnc2","Kcnc3","Kcnc4","Kcnd1","Kcnd2","Kcnd3","Kcnf1","Kcng1","Kcng2","Kcng3","Kcng4","Kcnh1","Kcnh2","Kcnh3","Kcnh4","Kcnh5","Kcnh6","Kcnh7","Kcnh8","Kcnq1","Kcnq2","Kcnq3","Kcnq4","Kcnq5","Kcns1","Kcns2","Kcns3","Kcnv1","Kcnv2")
KChannel_reg_subunits<-c("Dpp6","Dpp10","Kcnab1","Kcnab2","Kcnab3","Kcne1","Kcne1b","Kcne2","Kcne3","Kcne4","Kcne5","Kcnip1","Kcnip2","Kcnip3","Kcnip4","Kcnmb1","Kcnmb2","Kcnmb3","Kcnmb4", "Hcn1","Hcn2","Hcn3","Hcn4")
NaChannels<-c("Asic1","Asic2","Asic3","Asic4","Asic5","Scnn1a","Scnn1b","Scnn1d","Scnn1g","Nalcn","Scn1a","Scn2a","Scn3a","Scn4a","Scn5a","Scn8a","Scn9a","Scn10a","Scn11a","Scn1b","Scn2b","Scn3b","Scn4b","Scnm1")
CaChannels<-c("Cacna1a","Cacna1b","Cacna1c","Cacna1d","Cacna1e","Cacna1f","Cacna1g","Cacna1h","Cacna1i","Cacna1s","Cacna2d1","Cacna2d2","Cacna2d3","Cacna2d4","Cacnb1","Cacnb2","Cacnb3","Cacnb4","Cacng1","Cacng2","Cacng3","Cacng4","Cacng5","Cacng6","Cacng7","Cacng8","Catsperb","Catsperd","Catsperg","Catsper1","Catsper2","Catsper3","Catsper4","Itpr1","Itpr2","Itpr3","Ryr1","Ryr2","Ryr3","Tpcn1","Tpcn2")
ClChannels<-c("Ano1","Ano2","Ano3","Ano4","Ano5","Ano6","Ano7","Ano8","Ano9","Ano10","Best1","Best2","Best3","Best4","Cftr","Clic1","Clic2","Clic3","Clic4","Clic5","Clic6","Clcnka","Clcnkb","Bsnd","Clcn1","Clcn2","Clcn3","Clcn4","Clcn5","Clcn6","Clcn7")
# all of the above
Ion_channels<-c(KChannels,KChannel_reg_subunits,NaChannels,CaChannels,ClChannels)

fData(dat)$IonChannels<-fData(dat)$gene_short_name %in% Ion_channels
indxIC<-fData(dat)[fData(dat)$IonChannels & fData(dat)$isSignificant,"gene_short_name"]
length(indxIC)#0

paste(indxIC,indxNNKVL,indxGPCRs,indxOR,indxSTG,sep=",")
tmp<-union(indxIC,indxNNKVL)
tmp<-union(tmp,indxGPCRs)
tmp<-union(tmp,indxOR)
tmp<-union(tmp,indxSTG)
indx<-which(fData(dat)$gene_short_name %in% tmp)
length(indx)
#indx<-c(indxIC,indxNNKVL,indxGPCRs,indxOR,indxSTG)
#indx<-indx[-duplicated(indx)]


  sub<-dat[indx,]
  mat<-as.matrix(exprs(sub))
  mat<-log10(mat+1)
  #mat_df <- mat[apply(mat, MARGIN = 1, FUN = function(x) sd(x) != 0),]

pdf("plots/efiz_genes.pdf")
  pheatmap(mat=mat,
            scale="row",
            labels_row=fData(sub)$gene_short_name,
            annotation_col=pData(sub)[,c("source_plate","Total_mRNAs","num_genes_expressed","Cluster","Fluorogold")],
            labelsCol=FALSE,  
            #annotation_colors=c(label_colors,PMcolors),
            #color = colorRampPalette(piratepal(palette="brave"))(100), 
            #color = rev(brewer.pal(n = 11, name ="RdBu")))(100), # magma(100),
            clustering_distance_cols = "correlation",
            clustering_distance_rows = "correlation",
            show_rownames = T, show_colnames = F,
            cutree_rows=2, cutree_cols=2, breaks=seq(-3,3,length=101))
dev.off()

pdf("plots/efiz_genes_bycluster.pdf",height=20,width=20)
plot_cell_clusters(dat,markers=indx) + scale_color_viridis(option="viridis")
dev.off()


magExp<-apply(exprs(dat)[indx,pData(dat)$CellType=="Magnocellular"],1,mean)
pavExp<-apply(exprs(dat)[indx,pData(dat)$CellType=="Parvocellular"],1,mean)

mExp<-data.frame("magnocellular"=magExp,"parvocellular"=pavExp)
rownames(mExp)<-fData(dat)$gene_short_name[indx]
mExp

pdf("plots/efiz_genes_mean.pdf")
  pheatmap(mat=log10(mExp+1),
            #scale="row",
            clustering_distance_cols = "correlation",
            clustering_distance_rows = "correlation",
            show_rownames = T, show_colnames = T)
dev.off()

pdf("plots/efiz_genes_violin.pdf")
	monocle::plot_genes_violin(sub, grouping="CellType", ncol=3, min_expr=0.1,
		color_by="CellType",panel_order=tmp) +
		scale_fill_manual(values=PMcolors) + theme(legend.position = "top")
dev.off()

```

 power calculations use a permutation test with `r nPerm` simulations to estimate
  the power of a t-test with Benjamini Hotchberg adjusted p-values to detect a 
  specified number of significant guides from the total number of guides considered 
  for analysis at a given false discovery rate, with the specified number of 
  replicates per group. True positives are assumed to have a log fold change of 2 
  with specified standard deviation.


```{r e_code}
library('printr')
library('limma')

# m = number of genes that are differentially expressed
# g = total number of genes, sd = sd 
# mn = log fold change, n = number of samples per group

.markPwr=function(mn,sd=mn,n=50,m=100,g=25000,fdr=0.1){
  pT=2*pt(abs(rnorm(m,mean=mn,sd=sd)/sqrt(2/n)),lower.tail=F,df=n-1)
  pN=runif(g)
  padj=p.adjust(c(pT,pN),method="BH")
  ansout <- rbind(100*sum(padj[1:m]<=fdr)/m, 
                  sum(padj[seq(from=m+1, to=g)]<=fdr))
  row.names(ansout) <- c('TP','FP')
  return(ansout)
}

#given total number of cells 
pwr <- matrix(nrow=18,ncol=2, dimnames=list(paste0('g',seq(from=10,to=180,by=10)),c('TP','FP')))
nPerm <- 1000
for (midx in row.names(pwr)) {
  message(midx)
  m <- as.numeric(sub('g','',midx))
  ans=sapply(rep(2,nPerm),.markPwr,n=dim(pData(dat))[1],sd=0.5,
  	m=length(PV_vs_Magno_sigGeneIDs),g=length(expressed_genes),fdr=0.001)
  pwr[midx,] <- apply(ans,1,mean)
}

# using number of parvocellular cells
pwr <- matrix(nrow=18,ncol=2, dimnames=list(paste0('g',seq(from=10,to=180,by=10)),c('TP','FP')))
nPerm <- 1000
for (midx in row.names(pwr)) {
  message(midx)
  m <- as.numeric(sub('g','',midx))
  ans=sapply(rep(2,nPerm),.markPwr,n=sum(pData(dat)$CellType=="Parvocellular")*2,sd=0.5,
  	m=length(PV_vs_Magno_sigGeneIDs),g=length(expressed_genes),fdr=0.001)
  pwr[midx,] <- apply(ans,1,mean)
}

pwr

```

```{r }
param <- vector("list", 2)
names(param) <- c("Magnocellular", "Parvocellular")

diffanal <- function(eset, alpha) {
    y <- DGEList(counts=exprs(eset), group=pData(eset)$population)
    keep <- rowSums(cpm(y) > 1) >= 2
    y <- y[keep, keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y)

    tab <- exactTest(y)$table
    tab$FDR <- with(tab, p.adjust(PValue, "BH"))
    sigtab <- with(tab, tab[FDR < alpha & logFC > 0,])

    p <- nrow(eset)
    p1 <- nrow(sigtab)/p*100
    dropout <- (1 - sum(exprs(eset) > 0)/prod(dim(eset)))*100
    return(list(p=p, p1=p1, dropout=dropout, tab=sigtab))
}

estimateFC <- function(eset, siggenes) {
    emat <- log2(exprs(eset[siggenes,]) + 1)
    x0 <- with(pData(eset), emat[,population == 0])
    x1 <- with(pData(eset), emat[,population == 1])
    n0 <- ncol(x0)
    n1 <- ncol(x1)
    s0 <- apply(x0,1,sd)
    s1 <- apply(x1,1,sd)
    sp <- sqrt(((n0-1)*s0 + (n1-1)*s1)/(n0+n1-2))

    y0 <- apply(x0,2,function(x) x/sp)
    y1 <- apply(x1,2,function(x) x/sp)
    m0 <- rowMeans(y0)
    m1 <- rowMeans(y1)

    return(list(fc=2^mean(m1-m0),mu=mean(m0), sigma=sd(m0)))


ans <- diffanal(eset_sel, 0.05)
ans$n0 <- with(pData(eset_sel), sum(population == 0))
ans$n1 <- with(pData(eset_sel), sum(population == 1))
est <- estimateFC(eset_sel, rownames(ans$tab))
ans$fc <- est$fc
ans$mu <- est$mu
ans$sigma <- est$sigma

param[[1]] <- ans


```



```{r power_analysis_ncells}

#Expressed genes with expression in >=5 cells
numCellThreshold<-5
expressed_genes<-row.names(subset(fData(dat),num_cells_expressed >= numCellThreshold))
length(expressed_genes)

#from manuscript given Expected: in the PVN, 66% (952 ± 69 neurons) of OTergic neurons are magnocellular,
# while 34% (488 ± 25 neurons) of OTergic neurons are parvocellular 

pData(dat)$Cluster==1

sum(pData(dat)$CellType=="Parvocellular")/sum(pData(dat)$CellType=="Magnocellular")*100

pcalc<-ncells(m1=length(PV_vs_Magno_sigGeneIDs)/length(expressed_genes)*100,
	 ncore=4, ,dropout=.6, #dfactor, type1=qval_cutoff
	pi1=sum(pData(dat)$CellType=="Parvocellular")/sum(pData(dat)$CellType=="Magnocellular")*100, 
	p=length(expressed_genes),foldchange=2, 
	mu=1, sigma=5)


```
