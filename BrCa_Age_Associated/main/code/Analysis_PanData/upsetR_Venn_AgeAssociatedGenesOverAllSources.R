### upsetR venn graph for age associated genes
library("UpSetR")
##movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
##                   header = T, sep = ";")
#gsets <- read.table(file = "main/data/Data_PanData/AgeAssociatedGenes_AllSources_FC1.25.csv", sep = ",",
#                    stringsAsFactors = FALSE, header = TRUE)
gsets <- read.table(file = "main/data/Data_PanData/AgeAssociatedGenes_AllSources_FC1.25.csv", sep = ",",
                    header = TRUE)

gsetsl <- tapply(as.integer(gsets$Gene), gsets$Source, FUN = function(x) as.integer(unlist(list(x))), simplify = FALSE)


#png(filename="main/figures/upsetVennAllSourcesFC1p25.png", res = 300, width=7, height=5, units="in")
pdf(file="main/figures/upsetVennAllSourcesFC1p25v3.pdf",  width=12, height=5, useDingbats = FALSE)
print(upset(fromList(gsetsl), order.by = "freq", nsets = length(gsetsl), mb.ratio = c(0.5, 0.5)))
dev.off()

sort(gsetsl[["Breast_MB"]][[1]])

sort(intersect(gsetsl[["Breast_TCGA"]][[1]], gsetsl[["Breast_MB"]][[1]]))

sort(gsetsl[["Breast_TCGA"]][[1]])

sort(gsetsl[["Lung"]][[1]])

sort(gsetsl[["Thyroid"]][[1]])

sort(intersect(intersect(gsetsl[["Lung"]][[1]], gsetsl[["Breast_TCGA"]][[1]]), gsetsl[["Breast_MB"]][[1]]))

sort(intersect(gsetsl[["Lung"]][[1]], gsetsl[["Breast_MB"]][[1]]))

sort(intersect(gsetsl[["Lung"]][[1]], gsetsl[["Breast_TCGA"]][[1]]))

sort(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Breast_TCGA"]][[1]]))

sort(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Breast_MB"]][[1]]))

sort(intersect(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Breast_MB"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]))

sort(gsetsl[["Prostate"]][[1]])

sort(gsetsl[["Kidney"]][[1]])

sort(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Lung"]][[1]]))

sort(intersect(intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Breast_MB"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]))

sort(intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Breast_MB"]][[1]]))

sort(intersect(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Lung"]][[1]]), gsetsl[["Breast_MB"]][[1]]))

sort(intersect(intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Breast_MB"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]))

sort(intersect(intersect(intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Thyroid"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]), gsetsl[["Breast_MB"]][[1]]))

sort(intersect(intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Thyroid"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]))

sort(intersect(intersect(intersect(gsetsl[["Lung"]][[1]], gsetsl[["Thyroid"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]), gsetsl[["Breast_MB"]][[1]]))

paste(sort(intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Lung"]][[1]])), collapse = ",")

  sort(intersect(intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Breast_MB"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]))

sort(intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Breast_TCGA"]][[1]]))

sort(intersect(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Lung"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]))

sort(intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Breast_MB"]][[1]]))

sort(intersect(intersect(intersect(gsetsl[["Lung"]][[1]], gsetsl[["Prostate"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]), gsetsl[["Breast_MB"]][[1]]))


write.table(paste(c("Breast_MB",
                    sort(gsetsl[["Breast_MB"]][[1]])), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = FALSE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Breast_TCGA and Breast_MB",
                    sort(intersect(gsetsl[["Breast_TCGA"]][[1]], gsetsl[["Breast_MB"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Breast_TCGA",
                    sort(gsetsl[["Breast_TCGA"]][[1]])), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung",
                    sort(gsetsl[["Lung"]][[1]])), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid",
                    sort(gsetsl[["Thyroid"]][[1]])), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Breast_TCGA and Breast_MB",
                    sort(intersect(intersect(gsetsl[["Lung"]][[1]], gsetsl[["Breast_TCGA"]][[1]]), gsetsl[["Breast_MB"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Breast_MB",
                    sort(intersect(gsetsl[["Lung"]][[1]], gsetsl[["Breast_MB"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Breast_TCGA",
                    sort(intersect(gsetsl[["Lung"]][[1]], gsetsl[["Breast_TCGA"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid and Breast_TCGA",
                    sort(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Breast_TCGA"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid and Breast_MB",
                    sort(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Breast_MB"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid and Breast_TCGA and Breast_MB",
                    sort(intersect(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Breast_MB"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Prostate",
                    sort(gsetsl[["Prostate"]][[1]])), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Kidney",
                    sort(gsetsl[["Kidney"]][[1]])), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid and Lung",
                    sort(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Lung"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Kidney and Breast_TCGA and Breast_MB",
                    sort(intersect(intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Breast_MB"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Prostate and Breast_MB",
                    sort(intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Breast_MB"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid and Lung and Breast_MB",
                    sort(intersect(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Lung"]][[1]]), gsetsl[["Breast_MB"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Prostate and Breast_TCGA and Breast_MB",
                    sort(intersect(intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Breast_MB"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Kidney and Thyroid and Breast_TCGA and Breast_MB",
                    sort(intersect(intersect(intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Thyroid"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]), gsetsl[["Breast_MB"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Kidney and Thyroid and Breast_TCGA",
                    sort(intersect(intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Thyroid"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Thyroid and Breast_TCGA and Breast_MB",
                    sort(intersect(intersect(intersect(gsetsl[["Lung"]][[1]], gsetsl[["Thyroid"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]), gsetsl[["Breast_MB"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Prostate and Lung",
                    sort(intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Lung"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Prostate and Breast_TCGA",
                    sort(intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Breast_TCGA"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Thyroid and Breast_TCGA",
                    sort(intersect(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Lung"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Kidney and Breast_MB",
                    sort(intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Breast_MB"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Prostate and Breast_TCGA and Breast_MB",
                    sort(intersect(intersect(intersect(gsetsl[["Lung"]][[1]], gsetsl[["Prostate"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]), gsetsl[["Breast_MB"]][[1]]))), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)


###

sbmb <- sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(gsetsl[["Breast_MB"]][[1]], gsetsl[["Breast_TCGA"]][[1]]), 
          gsetsl[["Kidney"]][[1]]),
        gsetsl[["Lung"]][[1]]),
      gsetsl[["Prostate"]][[1]]),
    gsetsl[["Thyroid"]][[1]])
)
sbtcmb <- sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          intersect(gsetsl[["Breast_TCGA"]][[1]], gsetsl[["Breast_MB"]][[1]]), 
          gsetsl[["Kidney"]][[1]]),
        gsetsl[["Lung"]][[1]]),
      gsetsl[["Prostate"]][[1]]),
    gsetsl[["Thyroid"]][[1]])
)
sbtc <- sort(
    setdiff(
      setdiff(
        setdiff(
          setdiff(
            setdiff(gsetsl[["Breast_TCGA"]][[1]], gsetsl[["Breast_MB"]][[1]]), 
            gsetsl[["Kidney"]][[1]]),
          gsetsl[["Lung"]][[1]]),
        gsetsl[["Prostate"]][[1]]),
      gsetsl[["Thyroid"]][[1]])
  )
sl <- sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(gsetsl[["Lung"]][[1]], gsetsl[["Breast_MB"]][[1]]),
          gsetsl[["Breast_TCGA"]][[1]]), 
        gsetsl[["Kidney"]][[1]]),
      gsetsl[["Prostate"]][[1]]),
    gsetsl[["Thyroid"]][[1]])
)

st <- sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(gsetsl[["Thyroid"]][[1]], gsetsl[["Breast_MB"]][[1]]),
          gsetsl[["Breast_TCGA"]][[1]]), 
        gsetsl[["Kidney"]][[1]]),
      gsetsl[["Prostate"]][[1]]),
    gsetsl[["Lung"]][[1]])
)

slbtcmb <- sort(
  setdiff(
    setdiff(
      setdiff(
        intersect(
          intersect(gsetsl[["Lung"]][[1]], gsetsl[["Breast_TCGA"]][[1]]), 
          gsetsl[["Breast_MB"]][[1]]), 
        gsetsl[["Kidney"]][[1]]),
      gsetsl[["Prostate"]][[1]]),
    gsetsl[["Thyroid"]][[1]])
)
## length(slbtcmb) ## 119

slbmb <- sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          intersect(gsetsl[["Lung"]][[1]], gsetsl[["Breast_MB"]][[1]]), 
          gsetsl[["Breast_TCGA"]][[1]]), 
        gsetsl[["Kidney"]][[1]]),
      gsetsl[["Prostate"]][[1]]),
    gsetsl[["Thyroid"]][[1]])
)
## length(slbmb) ## 113

slbtc <- sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          intersect(gsetsl[["Lung"]][[1]], gsetsl[["Breast_TCGA"]][[1]]), 
          gsetsl[["Breast_MB"]][[1]]), 
        gsetsl[["Kidney"]][[1]]),
      gsetsl[["Prostate"]][[1]]),
    gsetsl[["Thyroid"]][[1]])
)
## length(slbtc) ## 92

stbtc <- sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Breast_TCGA"]][[1]]), 
          gsetsl[["Breast_MB"]][[1]]), 
        gsetsl[["Kidney"]][[1]]),
      gsetsl[["Prostate"]][[1]]),
    gsetsl[["Lung"]][[1]])
)
## length(stbtc) ## 55


stbmb <- sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Breast_MB"]][[1]]), 
          gsetsl[["Breast_TCGA"]][[1]]), 
        gsetsl[["Kidney"]][[1]]),
      gsetsl[["Prostate"]][[1]]),
    gsetsl[["Lung"]][[1]])
)
## length(stbmb) ## 55

stbmbtc <- sort(
  setdiff(
    setdiff(
      setdiff(
        sort(intersect(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Breast_MB"]][[1]]), gsetsl[["Breast_TCGA"]][[1]])), 
        gsetsl[["Kidney"]][[1]]),
      gsetsl[["Prostate"]][[1]]),
    gsetsl[["Lung"]][[1]])
)
## length(stbmbtc) ## 52

sp <- sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(gsetsl[["Prostate"]][[1]], gsetsl[["Breast_MB"]][[1]]),
          gsetsl[["Breast_TCGA"]][[1]]), 
        gsetsl[["Kidney"]][[1]]),
      gsetsl[["Thyroid"]][[1]]),
    gsetsl[["Lung"]][[1]])
)
## length(sp) ## 20

sk <-sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(gsetsl[["Kidney"]][[1]], gsetsl[["Breast_MB"]][[1]]),
          gsetsl[["Breast_TCGA"]][[1]]), 
        gsetsl[["Prostate"]][[1]]),
      gsetsl[["Thyroid"]][[1]]),
    gsetsl[["Lung"]][[1]])
)
## length(sk) ## 12

stl <-sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Lung"]][[1]]), 
          gsetsl[["Breast_MB"]][[1]]),
        gsetsl[["Breast_TCGA"]][[1]]), 
      gsetsl[["Prostate"]][[1]]), 
    gsetsl[["Kidney"]][[1]])
)
## length(stl) ## 12

skbtcmb <-sort(
  setdiff(
    setdiff(
      setdiff(
        intersect(
          intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Breast_TCGA"]][[1]]), 
          gsetsl[["Breast_MB"]][[1]]),
        gsetsl[["Lung"]][[1]]), 
      gsetsl[["Prostate"]][[1]]), 
    gsetsl[["Thyroid"]][[1]])
)
## length(skbtcmb) ## 5

spbmb <- sort(
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Breast_MB"]][[1]]), 
          gsetsl[["Breast_TCGA"]][[1]]),
        gsetsl[["Lung"]][[1]]), 
      gsetsl[["Kidney"]][[1]]), 
    gsetsl[["Thyroid"]][[1]])
)
## length(spbmb) ## 3

stlbmb <- sort(
  setdiff(
    setdiff(
      setdiff(
        intersect(intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Lung"]][[1]]), gsetsl[["Breast_MB"]][[1]]), 
        gsetsl[["Breast_TCGA"]][[1]]),
      gsetsl[["Kidney"]][[1]]), 
    gsetsl[["Prostate"]][[1]])
)
## length(stlbmb) ## 3

spbmbtc <- sort(
  setdiff(
    setdiff(
      setdiff(
        intersect(intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Breast_MB"]][[1]]), gsetsl[["Breast_TCGA"]][[1]]), 
        gsetsl[["Kidney"]][[1]]),
      gsetsl[["Lung"]][[1]]), 
    gsetsl[["Thyroid"]][[1]])
)
## length(spbmbtc) ## 3

sktbmbtc <- sort(
  setdiff(
    setdiff(
      intersect(
        intersect(
          intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Thyroid"]][[1]]), 
          gsetsl[["Breast_TCGA"]][[1]]), 
        gsetsl[["Breast_MB"]][[1]]), 
      gsetsl[["Lung"]][[1]]), 
    gsetsl[["Prostate"]][[1]])
)
## length(sktbmbtc) ## 3

sktbtc <- sort(
  setdiff(
    setdiff(
      setdiff(
        intersect(
          intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Thyroid"]][[1]]), 
          gsetsl[["Breast_TCGA"]][[1]]), 
        gsetsl[["Breast_MB"]][[1]]),
      gsetsl[["Lung"]][[1]]), 
    gsetsl[["Prostate"]][[1]])
)
## length(sktbtc) ## 2

stlbmbtc <- sort(
  setdiff(
    setdiff(
      intersect(
        intersect(
          intersect(gsetsl[["Lung"]][[1]], gsetsl[["Thyroid"]][[1]]), 
          gsetsl[["Breast_TCGA"]][[1]]), 
        gsetsl[["Breast_MB"]][[1]]), 
      gsetsl[["Kidney"]][[1]]), 
    gsetsl[["Prostate"]][[1]])
)
## length(stlbmbtc) ## 2





write.table(paste(c("Breast_MB", sbmb), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = FALSE, row.names = FALSE, col.names = FALSE, quote = FALSE) ## 1808 entries

write.table(paste(c("Breast_TCGA and Breast_MB", sbtcmb), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Breast_TCGA", sbtc), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)  ## 1354 entries

write.table(paste(c("Lung", sl), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid", st), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Breast_TCGA and Breast_MB", slbtcmb), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Breast_MB", slbmb), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Breast_TCGA", slbtc), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid and Breast_TCGA", stbtc), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid and Breast_MB", stbmb), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid and Breast_TCGA and Breast_MB", stbmbtc), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Prostate", sp), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Kidney", sk), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid and Lung", stl), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Kidney and Breast_TCGA and Breast_MB", skbtcmb), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Prostate and Breast_MB", spbmb), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Thyroid and Lung and Breast_MB", stlbmb), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Prostate and Breast_TCGA and Breast_MB", spbmbtc), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Kidney and Thyroid and Breast_TCGA and Breast_MB", sktbmbtc), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Kidney and Thyroid and Breast_TCGA", sktbtc), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Thyroid and Breast_TCGA and Breast_MB", stlbmbtc), collapse = ","),
            file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Prostate and Lung",
                    sort(
                      setdiff(
                        setdiff(
                          setdiff(
                            setdiff(
                              intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Lung"]][[1]]), 
                            gsetsl[["Breast_MB"]][[1]]),
                          gsetsl[["Breast_TCGA"]][[1]]), 
                        gsetsl[["Thyroid"]][[1]]), 
                      gsetsl[["Kidney"]][[1]]))
), collapse = ","),
file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Prostate and Breast_TCGA",
                    sort(
                      setdiff(
                        setdiff(
                          setdiff(
                            setdiff(
                              intersect(gsetsl[["Prostate"]][[1]], gsetsl[["Breast_TCGA"]][[1]]), 
                              gsetsl[["Breast_MB"]][[1]]),
                            gsetsl[["Lung"]][[1]]), 
                          gsetsl[["Thyroid"]][[1]]), 
                        gsetsl[["Kidney"]][[1]]))
), collapse = ","),
file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Thyroid and Breast_TCGA",
                    sort(
                      setdiff(
                        setdiff(
                          setdiff(
                            intersect(
                              intersect(gsetsl[["Thyroid"]][[1]], gsetsl[["Lung"]][[1]]), 
                              gsetsl[["Breast_TCGA"]][[1]]),
                            gsetsl[["Breast_MB"]][[1]]), 
                          gsetsl[["Prostate"]][[1]]), 
                        gsetsl[["Kidney"]][[1]]))
), collapse = ","),
file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Kidney and Breast_MB",
                    sort(
                      setdiff(
                        setdiff(
                          setdiff(
                            setdiff(
                              intersect(gsetsl[["Kidney"]][[1]], gsetsl[["Breast_MB"]][[1]]), 
                              gsetsl[["Breast_TCGA"]][[1]]),
                            gsetsl[["Lung"]][[1]]), 
                          gsetsl[["Prostate"]][[1]]), 
                        gsetsl[["Thyroid"]][[1]]))
), collapse = ","),
file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(paste(c("Lung and Prostate and Breast_TCGA and Breast_MB",
                    sort(
                      setdiff(
                        setdiff(
                          intersect(
                            intersect(
                              intersect(gsetsl[["Lung"]][[1]], gsetsl[["Prostate"]][[1]]), 
                              gsetsl[["Breast_TCGA"]][[1]]), 
                            gsetsl[["Breast_MB"]][[1]]), 
                          gsetsl[["Thyroid"]][[1]]), 
                        gsetsl[["Kidney"]][[1]]))
), collapse = ","),
file = "UpSetR_MB_TCGA_GeneVennSets.csv", append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)


gsr <- read.table(file = "UpSetR_MB_TCGA_GeneVennSets.csv", fill = TRUE, sep = ",", header = FALSE, stringsAsFactors = FALSE)
gsrmat <- t(as.matrix(gsr))
write.table(gsrmat, 
            file = paste("./main/code/Analysis_PanData/",
                         "Supplemental_UpSetR_MB_TCGA_GeneVennSetCols.csv", 
                         sep = ""),
            row.names = FALSE, col.names = FALSE, sep = ",")
file.remove(file = "UpSetR_MB_TCGA_GeneVennSets.csv")
file.remove(file = "UpSetR_MB_TCGA_GeneSets.csv")
file.remove(file = "UpSetR_MB_TCGA_GeneVennSetCols.csv")

gsrvec <- as.vector(gsrmat[-1, ])
str(gsrvec)
gsrlab <- gsub(" ", ".", rep(gsrmat[1, ], each = nrow(gsrmat[-1, ])))
head(gsrlab)
tail(gsrlab)
length(gsrlab)
gsrlmat <- cbind(gsrvec, gsrlab)
head(gsrlmat)
dimnames(gsrlmat)[[2]] <- c("Hugo_Symbol", "AgeAssociatedVennSet")
if(any(gsrlmat[, "Hugo_Symbol"] == "")) { gsrlmat <- gsrlmat[which(gsrlmat[, "Hugo_Symbol"] != ""), ] }
gsrlmat[, "AgeAssociatedVennSet"] <- gsub("\\.", " ", gsrlmat[, "AgeAssociatedVennSet"] )
tail(gsrlmat, 30)
write.table(gsrlmat, file = paste("./main/code/Analysis_PanData/",
                                  "UpSetR_MB_TCGA_Gene_VennSet.csv", 
                                  sep = ""),
            sep = ",", row.names = FALSE
            )
