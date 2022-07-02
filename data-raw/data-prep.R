## 1. Lee data from Mostly Harmless Econometrics website
dir1 <- "~/teaching/Datasets/Lee2008/table_two_final.dta"
lee <- readstata13::read.dta13(dir1)
lee <- subset(lee, difdemshare!=0 & use ==1)
s <- sort(lee$difdemshare, index.return=TRUE)       # sort

lee08 <- data.frame(voteshare=100*lee$demsharenext[s$ix], margin=100*s$x)
devtools::use_data(lee08, overwrite=TRUE, internal=FALSE)

## 2. Oreopoulos data from AER website
dir2 <- "~/teaching/Datasets/Oreopoulos2006/uk/combined general household survey.dta"
cghs <- readstata13::read.dta13(dir2, generate.factors=TRUE,
                                nonint.factors=TRUE)

cghs$yearat14 <- cghs$yobirth+14
d <- within(cghs,
            keep <- yearat14 >= 35 & yearat14 <= 65 & age <= 64 & agelfted >=
                10 & !is.na(agelfted) & !is.na(earn) & nireland==0)
d <- d[d$keep, ]
cghs <- data.frame(earnings=d$earn, yearat14=1900+d$yearat14)
devtools::use_data(cghs, overwrite=TRUE, internal=FALSE)

## 3. Lalive data from Rafael Lalive's website,
## https://sites.google.com/site/rafaellalive/research
dir3 <- "~/teaching/Datasets/Lalive2008/releaseData.dta"
rebp <- readstata13::read.dta13(dir3)
## Only treated region
rebp <- rebp[rebp$tr==1, c("age", "period", "female", "unemployment_duration")]
names(rebp)[4] <- "duration"
rebp$period <- rebp$period==1
rebp$female <- rebp$female==1
## Fix rounding in stata
rebp$age <- round(rebp$age*12)/12
devtools::use_data(rebp, overwrite=TRUE, internal=FALSE)

## 4. LM data from Douglas Miller's website

dir4 <- "~/teaching/Datasets/LudwigMiller2007/analysis data/"
## Table 3 and 4
d3 <- foreign::read.dta(paste0(dir4, "census3.dta"))
d4 <- foreign::read.dta(paste0(dir4, "census_1990.dta"))

## observation 3133 is full of NA's
d4 <- d4[!is.na(d4$oldcode), ]
## now both d3 and d4 have the same number of obs, 3138

## 27056 Yellowstone County, MT is there twice, remove it
dup <- d3$oldcode[duplicated(d3$oldcode)]
dup == d4$oldcode[duplicated(d4$oldcode)]
d3 <- d3[-which.max(d3$oldcode==dup), ]
d4 <- d4[-which.max(d4$oldcode==dup), ]

## left with 3137 observations, merge them
d3$mortHS <- d3$age5_9_sum2
d3$mortInj <- d3$age5_9_injury_rate
d4$highSchool <- d4$hsplus18_24
d <- merge(d3[, c("oldcode", "state", "povrate60", "mortHS", "mortInj")],
           d4[, c("oldcode", "state", "povrate60", "highSchool")])
## remove 6 with all three outcomes missing
d$oldcode[(is.na(d$mortHS) & is.na(d$highSchool))]
## 12003 13045 42065 47061 47102 51024, most of these have most data missing
d <- d[!(is.na(d$mortHS) & is.na(d$highSchool)), ]

## remove 3 with poverty rate missing
d <- d[!is.na(d$povrate60), ]

## left with 3128 observations

## Add county names to dataset
url <- "https://www2.census.gov/geo/docs/reference/codes/files/national_county.txt"
if (!file.exists("fips.txt")) download.file(url, destfile="fips.txt")
fips <- read.csv("fips.txt", header=FALSE)
names(fips) <- c("statepc", "statefp", "countyfp", "county", "classfp")
fips$classfp <- NULL
fips$cfips <- 1000*fips$statefp+fips$countyfp
## These correspond to
## statefp         State FIPS Code         12
## countyfp        County FIPS Code        011
## county      County Name      Broward County
## classfp         FIPS Class Code         H1
## drop non-us

## Add deleted counties, see https://www.census.gov/geo/reference/county-changes.html
fips <- rbind(fips,
              data.frame(statepc=c("AK", "MT", "VA", "FL", "AK", "AL"),
                         statefp=c(2, 30, 51, 12, 2, 2),
                         countyfp=c(280, 113, 780, 25, 231, 201),
                         county=c("Wrangell-Petersburg",
                             "Yellowstone National Park",
                             "South Boston", "Date County",
                             "Skagway-Yakutat-Angoon",
                             "Prince of Wales-Outer Ketchikan"),
                         cfips=c(2280, 30113, 51780, 12025, 2231, 2201)))
fips <- fips[fips$statefp<57, ]


## 2.1 Add county names to dateset
lmfips <- foreign::read.dta(paste0(dir, "oldcode_fips.dta"))
names(lmfips) <- c("oldcode", "statefp", "countyfp", "cfips")
lmfips <- lmfips[!is.na(lmfips$oldcode), ]
## Drop first batch of duplicates and non-matching
c(46135, 42068) %in%  d$oldcode
lmfips <- lmfips[-c(313, 548, 2427, 2429, 2825, 2432, 1652), ] # drop duplicates
lmfips[!complete.cases(lmfips), ]

## clean this
lmfips["312", c("statefp", "countyfp")] <- c(10, 1)
lmfips["549", c("statefp", "countyfp")] <- c(16, 1)
lmfips["2428", c("statefp", "countyfp")] <- c(46, 131)
lmfips["2430", c("statefp", "countyfp")] <- c(46, 135)
lmfips[!complete.cases(lmfips), ]

## still some duplicates
lmfipsnd <- lmfips[!duplicated(lmfips$oldcode), ]
lmfipsd <- lmfips[duplicated(lmfips$oldcode), ]

d2 <- merge(x=d, y=lmfipsnd, by="oldcode", all.x=TRUE)
## merged everything

dd <- merge(d2, fips, all.x=TRUE)
dd[!complete.cases(dd), ]
## merged everything

cutoff <- sort(dd$povrate)[nrow(dd)-299]
dd$povrate60 <- dd$povrate60-cutoff
## order by povrate
headst <- dd[order(dd[, "povrate60"]), ]
headst$cfips <- NULL
headst$state <- NULL

## Check that we match Table III before removing the Yellowstone outlier
headst <- headst[headst$oldcode!=27057, ] # Yellowstone National Park

devtools::use_data(headst, overwrite=TRUE, internal=FALSE)

## 4. Battistin et al data from AER website
rcp <- readstata13::read.dta13("~/teaching/Datasets/BattistinEtAl2009/datapaper_ab.dta", generate.factors=TRUE, nonint.factors=TRUE)
rcp <- rcp[, c(2, 29, 27, 8, 4, 6)]
## Survey year
names(rcp) <- c("survey_year", "elig_year", "retired", "food", "c", "cn")
rcp$retired <- rcp$retired == "retired"
## drop if at
rcp <- rcp[rcp$elig_year!=0, ]
usethis::use_data(rcp, overwrite=TRUE, internal=FALSE)

## // profiles by S: collapse (mean) mc=lnc mcn=lncn mf=lncf ret=retired
## hh=hh_size age=age minors=minors children=children mcfp=lncfp, by(time anno)

## generate elig = time>=0
## keep if abs(time)<=10 & time!=0

## IV regressions with year dummies: Table 5
## ivregress 2sls mcn (ret=elig) time c.time#c.time i.anno, robust first
## ivregress 2sls mf (ret=elig) time c.time#c.time i.anno, robust
