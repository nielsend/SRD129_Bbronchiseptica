library(dplyr)
dayneg0 <- read.csv(file = 'day0andNeg.csv')

dayneg0list <- unique(dayneg0$X.NAME.)
print(dayneg0list)
