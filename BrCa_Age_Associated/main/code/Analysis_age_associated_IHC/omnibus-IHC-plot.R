# setwd("~/Shared/manuscripts/BrCa_Age_Associated-GIT/main/code/Analysis_NormalBreast13_TMA_IHC_Trends/Tables")

normal.ihc.corr <- 
  read.csv("./main/code/Analysis_NormalBreast13_TMA_IHC_Trends/Tables/NormalBreast13_TMA_Correlations.csv", 
            header = T, stringsAsFactors = F)

# setwd("~/Shared/manuscripts/BrCa_Age_Associated-GIT/main/code/Analysis_MB09_TMA_IHC_Trends/Tables")

mb.ihc.corr <- 
  read.csv("./main/code/Analysis_MB09_TMA_IHC_Trends/Tables/MB09_TMA_Correlations.csv", 
           header = T, stringsAsFactors = F)

# setwd("~/Shared/manuscripts/BrCa_Age_Associated-GIT/main/code/Analysis_BigSeries_TMA_IHC_Trends/Tables")

big.ihc.corr <- 
  read.csv("./main/code/Analysis_BigSeries_TMA_IHC_Trends/Tables/BigSeries_TMA_Correlations.csv", 
           header = T, stringsAsFactors = F)

# setwd("~/Shared/manuscripts/BrCa_Age_Associated-GIT/main/code/Analysis_AgeAssociated_GeneCounts")

##TMA_correlations_omnibus.csv from the above

# omnibus.ihc.corr <- read.csv("../Analysis_AgeAssociated_GeneCounts/TMA_correlations_omnibus.csv",
#                              header = T, stringsAsFactors = F)
omnibus.ihc.corr <- read.csv("./main/code/Analysis_AgeAssociated_GeneCounts/TMA_correlations_omnibus.csv",
                             header = T, stringsAsFactors = F)

omnibus.ihc.corr$Gene.f <- as.factor(gsub("ER-PR", "ER-IHC", omnibus.ihc.corr$Canonical.Gene))

require("reshape")
require("ggplot2")
require("forcats")
omnibus.ihc.corr.melt <-
  melt(omnibus.ihc.corr, id = c("Gene.f"),
###        measure.vars = c("BIGSERIES.WHOLE", "BIGSERIES.TRAINING", "BIGSERIES.VALIDATION", "METABRIC", "NORMAL.BREAST"),
       measure.vars = c("BIGSERIES.WHOLE", "METABRIC", "NORMAL.BREAST"),
###       value_name = "correlation", variable_name = "dataset", na.rm = T)
              variable_name = "Dataset", na.rm = T)
names(omnibus.ihc.corr.melt)[which(names(omnibus.ihc.corr.melt) == "value")] <- "correlation"
head(omnibus.ihc.corr.melt)
write.csv(omnibus.ihc.corr.melt, file="main/code/Analysis_age_associated_IHC/omnibus-protein-express.csv")

p <- ggplot(omnibus.ihc.corr.melt, aes(x = fct_reorder(Gene.f, -correlation, median), y = correlation))

png(filename = "./main/code/Analysis_age_associated_IHC/omnibus-protein-express.png",
    width = 6, height = 3, units = "in", pointsize = 10,
    bg = "white",  res = 300)

print(p + geom_point(aes(color = Dataset), pch = 19,  size = 5, alpha = 0.5) +
        scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#999999")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        theme(legend.position = c(0.6, 0.8))+labs(x = "Gene", y = "IHC dynamic range association")
)
dev.off()

pdf(file = "./main/code/Analysis_age_associated_IHC/omnibus-protein-express.pdf",
    width = 6, height = 3.5, bg = "white")

print(p + geom_point(aes(color = Dataset), pch = 19,  size = 5, alpha = 0.5) +
        scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#999999")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        theme(legend.position = c(0.6, 0.8))+labs(x = "Gene", y = "IHC dynamic range association")
)
dev.off()
