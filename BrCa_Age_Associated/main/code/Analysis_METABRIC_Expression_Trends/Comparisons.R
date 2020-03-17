### Compare METABRIC and TCGA expression findings

TCGAintClustN <- c(75, 38, 181, 165, 84, 60, 100, 145, 74, 157)
MBintClustN <- c(140, 72, 294, 344, 191, 86, 193, 300, 146, 226)
cor(TCGAintClustN, MBintClustN, method = "spearman")
cor.test(TCGAintClustN, MBintClustN, method = "spearman")
### > cor.test(TCGAintClustN, MBintClustN, method = "spearman")
### 
### 	Spearman's rank correlation rho
### 
### data:  TCGAintClustN and MBintClustN
### S = 12, p-value = 0.0001302
### alternative hypothesis: true rho is not equal to 0
### sample estimates:
###       rho 
### 0.9272727 


TCGAintClustAgeAssocN <- c(0, 0, 177, 1, 0, 0, 1, 258, 0, 0)
MBintClustAgeAssocN <- c(8, 0, 878, 536, 4, 0, 12, 636, 1, 0)
cor.test(TCGAintClustAgeAssocN, MBintClustAgeAssocN, method = "spearman")
### > cor.test(TCGAintClustAgeAssocN, MBintClustAgeAssocN, method = "spearman")
### 
### 	Spearman's rank correlation rho
### 
### data:  TCGAintClustAgeAssocN and MBintClustAgeAssocN
### S = 19.595, p-value = 0.0007522
### alternative hypothesis: true rho is not equal to 0
### sample estimates:
###       rho 
### 0.8812435 
### 
### Warning message:
### In cor.test.default(TCGAintClustAgeAssocN, MBintClustAgeAssocN,  :
###   Cannot compute exact p-value with ties



