### Install and library Packages----
# install.packages('meta')
# install.packages('metafor')
library(meta)
library(metafor)

### Read data----
df_miRPoly <- read.csv(file="miRPoly_data.csv", header = T)
Polymorphisms <- c("DICER A>G rs3742330", "DROSHA T>C rs10719", "RAN C>T rs14035", 
                   "miR-10a A>T rs3809783", "miR-124 C>G rs531564", "miR-125a C>T rs12976445", 
                   "miR-125a G>A rs41275794", "miR-146a C>G rs2910164", "miR-149 T>C rs2292832", 
                   "miR-196a2 C>T rs11614913", "miR-27a A>G rs895819", "miR-323b T>A rs56103835", 
                   "miR-423 C>A rs6505162", "miR-499a A>G rs3746444")
Genetic_models <- c("Allele model", "Dominant model", "Recessive model", "Homozygotic model", "Heterozygotic model")
df_miRPoly$polymorphism <- factor(df_miRPoly$polymorphism, levels = Polymorphisms)
df_miRPoly$model <- factor(df_miRPoly$model, levels = Genetic_models)

### Meta-analyses----
for (i in 1:length(Polymorphisms)){
  meta <- metabin(case1, case.totol, control1, control.total, data=df_miRPoly, sm="OR", studlab=paste(Author, Year), 
                    subgroup = model, allstudies = T, overall = F, overall.hetstat = F, label.e = "RPL", label.c = "Control", 
                    incr=1, test.subgroup=F, fixed=F, subset = polymorphism==Polymorphisms[i])
  df.result <- data.frame(polymorphism=Polymorphisms[i], model=meta[["bylevs"]], n=meta[["k.w"]], 
                          OR.random=exp(meta[["TE.random.w"]]), lower.OR.random=exp(meta[["lower.random.w"]]), upper.OR.random=exp(meta[["upper.random.w"]]), p.random=meta[["pval.random.w"]],
                          OR.fixed=exp(meta[["TE.fixed.w"]]), lower.OR.fixed=exp(meta[["lower.fixed.w"]]), upper.OR.fixed=exp(meta[["upper.fixed.w"]]), p.fixed=meta[["pval.fixed.w"]],
                          I2=meta[["I2.w"]], p.heterogeneity=meta[["pval.Q.w"]])
  write.csv(df.result, file = paste("Result.", Polymorphisms[i], ".csv", sep = ""), row.names=F)
}

### Forest plots for meta-analyses
for (i in 1:length(Polymorphisms)){
  meta <- metabin(case1, case.totol, control1, control.total, data=df_miRPoly, sm="OR", studlab=paste(Author, Year), 
                  subgroup = model, allstudies = T, overall = F, overall.hetstat = F, label.e = "RPL", label.c = "Control", 
                  incr=1, test.subgroup=F, subset = polymorphism==Polymorphisms[i])
  pdf(paste("forest.", Polymorphisms[i], ".pdf", sep = ""), width = 18, height = 20, title = paste("Forest plot showing meta-analysis of", Polymorphisms[i], sep = " "))
  forest(meta, col.diamond.fixed="blue", 
         col.diamond.random="red", digits=2, print.subgroup.labels = T)
  dev.off()
}

### sensitivity analysis----
for (i in 1:length(Polymorphisms)) {
  for (j in 1:length(Genetic_models)) {
    meta_sensitivity <- metabin(case1, case.totol, control1, control.total, data=df_miRPoly, sm="OR", studlab=paste(Author, Year), 
                                allstudies = T, label.e = "RPL", label.c = "Control", 
                                incr=1, fixed=F, subset = polymorphism==Polymorphisms[i] & model==Genetic_models[j])
    sensitivity <- metainf(meta_sensitivity, pooled = "random")
    pdf(paste("Sensitivity.", Polymorphisms[i], Genetic_models[j], ".pdf", sep = ""), width = 12, height = 8, title = paste("Sensitivity analysis of", Polymorphisms[i], "under", Genetic_models[j], sep = " "))
    forest(sensitivity)
    dev.off()
  }
}

### Subgroup meta-analysis----
# Subgroup based on ethnicity
for (i in 1:length(Polymorphisms)) {
  for (j in 1:length(Genetic_models)) {
    meta_subgroup <- metabin(case1, case.totol, control1, control.total, data=df_miRPoly, sm="OR", studlab=paste(Author, Year), 
                             subgroup = ethnicity, allstudies = T, overall = T, overall.hetstat = T, label.e = "RPL", label.c = "Control", 
                             incr=1, test.subgroup=F, subset = polymorphism==Polymorphisms[i] & model==Genetic_models[j])
    df.subgroup.result <- data.frame(polymorphism=Polymorphisms[i], model=Genetic_models[j], ethnicity.1=meta_subgroup[["bylevs"]][1], n.1=meta_subgroup[["k.w"]][1], 
                            OR.random.1=exp(meta_subgroup[["TE.random.w"]][1]), lower.OR.random.1=exp(meta_subgroup[["lower.random.w"]][1]), upper.OR.random.1=exp(meta_subgroup[["upper.random.w"]][1]), p.random.1=meta_subgroup[["pval.random.w"]][1],
                            OR.fixed.1=exp(meta_subgroup[["TE.fixed.w"]][1]), lower.OR.fixed.1=exp(meta_subgroup[["lower.fixed.w"]][1]), upper.OR.fixed.1=exp(meta_subgroup[["upper.fixed.w"]][1]), p.fixed.1=meta[["pval.fixed.w"]][1],
                            I2.1=meta_subgroup[["I2.w"]][1], p.heterogeneity.1=meta_subgroup[["pval.Q.w"]][1],
                            ethnicity.2=meta_subgroup[["bylevs"]][2], n.2=meta_subgroup[["k.w"]][2], 
                            OR.random.2=exp(meta_subgroup[["TE.random.w"]][2]), lower.OR.random.2=exp(meta_subgroup[["lower.random.w"]][2]), upper.OR.random.2=exp(meta_subgroup[["upper.random.w"]][2]), p.random.2=meta_subgroup[["pval.random.w"]][2],
                            OR.fixed.2=exp(meta_subgroup[["TE.fixed.w"]][2]), lower.OR.fixed.1=exp(meta_subgroup[["lower.fixed.w"]][2]), upper.OR.fixed.2=exp(meta_subgroup[["upper.fixed.w"]][2]), p.fixed.2=meta[["pval.fixed.w"]][2],
                            I2.2=meta_subgroup[["I2.w"]][2], p.heterogeneity.2=meta_subgroup[["pval.Q.w"]][2])
    write.csv(df.subgroup.result, file = paste("Result.", Polymorphisms[i], Genetic_models[j],  ".csv", sep = ""), row.names=F)
    pdf(paste("Subgroup", Polymorphisms[i], Genetic_models[j], ".pdf", sep = ""), width = 12, height = 10, title = paste("Forest plot showing subgroup analysis of", Polymorphisms[i], "under", Genetic_models[j], sep = " "))
    forest(meta_subgroup, col.diamond.fixed="blue",
           col.diamond.random="red", digits=2, print.subgroup.labels = T)
    dev.off()
  }
}


