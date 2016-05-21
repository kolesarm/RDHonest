## 2. Lee data

## Dataset from Mostly Harmless Econometrics website
lee <- readstata13::read.dta13("~/teaching/Datasets/Lee2008/table_two_final.dta")
lee <- subset(lee, difdemshare!=0 & use ==1)
s <- sort(lee$difdemshare, index.return=TRUE)       # sort

lee08 <- data.frame(voteshare=100*lee$demsharenext[s$ix], margin=100*s$x)
devtools::use_data(lee08, overwrite=TRUE, internal=FALSE)
