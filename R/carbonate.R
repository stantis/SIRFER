library(readxl)
library(openxlsx)
source("R/helpers.R")

d = read_xls("data/220510_COND.xls")

#remove on-off
d = d[d$Method != "CO2 On Off",]

#re-register run sequece #s
d$Row = d$Row - min(d$Row) + 1

#check for the right # of peaks
for(i in unique(d$Row)){
  r = max(d$`Peak Nr`[d$Row == i])
  if(r != 9){message(paste("Missing peak(s) for seq #", i))}
}

#get sample peaks
d = d[d$`Peak Nr` %in% c(4:8),]

#detect and remove outliers
d$Outliler = NA
for(i in unique(d$Row)){
  d$Outliler[d$Row == i] = dout(d[d$Row == i,])
}
d.no = d[!d$Outliler,]

#calculate stats per sample
l = length(unique(d.no$Row))
d.sam = data.frame(
  "Row" = unique(d.no$Row),
  "Identifier_1" = character(l),
  "Weight" = numeric(l),
  "Area.max" = numeric(l),
  "Area.sd" = numeric(l),
  "d13C.mean" = numeric(l),
  "d13C.sd" = numeric(l),
  "d18O.mean" = numeric(l),
  "d18O.sd" = numeric(l),
  "Ignore" = rep(FALSE)
)
for(i in d.sam$Row){
  d.sam$Identifier_1[i] = d.no$`Identifier 1`[match(i, d.no$Row)]
  d.sam$Weight[i] = d.no$`Identifier 2`[match(i, d.no$Row)]
  #areas include outliers for standardization over all peaks
  d.sam$Area.max[i] = max(d$`Area All`[d$Row == i])
  d.sam$Area.sd[i] = sd(d$`Area All`[d$Row == i])
  #isotopes exclude outliers
  d.sam$d13C.mean[i] = mean(d.no$`d 13C/12C`[d.no$Row == i])
  d.sam$d13C.sd[i] = sd(d.no$`d 13C/12C`[d.no$Row == i])
  d.sam$d18O.mean[i] = mean(d.no$`d 18O/16O`[d.no$Row == i])
  d.sam$d18O.sd[i] = sd(d.no$`d 18O/16O`[d.no$Row == i])
}

#screen for high SDs
d.sam$Ignore = d.sam$d13C.sd > 0.25 | d.sam$d18O.sd > 0.3
d.good = d.sam[!d.sam$Ignore,]

#some values to use
plrm1 = list("ID" = "CARRARA", "d13C" = 2.1, "d18O" = -1.8, "pCO3" = 0.6)
plrm2 = list("ID" = "LSVEC", "d13C" = -46.6, "d18O" = -26.7, "pCO3" = 0.6)
slrm = list("ID" = "MARBLE", "d13C" = 1.9, "d18O" = -11.3, "pCO3" = 0.6)

#drift correction
##fits and plots the spine
dfit = drift(d.good, plrm1, plrm2, slrm)

##apply the drift correction 
d.good$d13C.dc = d.good$d13C.mean - predict(dfit[[1]], d.good$Row)$y
d.good$d18O.dc = d.good$d18O.mean - predict(dfit[[2]], d.good$Row)$y

#calibration
##calibration fit
cfit = cal(d.good, plrm1, plrm2, slrm)

##apply the calibration
d.good$d13C.cal = predict(cfit$c.cal, d.good)
d.good$d13C.cal.se = predict(cfit$c.cal, d.good, se.fit = TRUE)$se.fit
d.good$d18O.cal = predict(cfit$o.cal, d.good)
d.good$d18O.cal.se = predict(cfit$o.cal, d.good, se.fit = TRUE)$se.fit
d.good$pCO3 = predict(cfit$k.cal, d.good) / d.good$Weight
d.good$pCO3.se = predict(cfit$k.cal, d.good, se.fit = TRUE)$se.fit / 
  d.good$Weight

#write the results
##parse the results
plrm1.data = d.good[d.good$Identifier_1 == plrm1$ID,
                    c("Row", "Weight", "Area.max", "d13C.cal",
                      "d13C.cal.se", "d18O.cal", "d18O.cal.se",
                      "pCO3", "pCO3.se")]
plrm2.data = d.good[d.good$Identifier_1 == plrm2$ID,
                    c("Row", "Weight", "Area.max", "d13C.cal",
                      "d13C.cal.se", "d18O.cal", "d18O.cal.se",
                      "pCO3", "pCO3.se")]
slrm.data = d.good[d.good$Identifier_1 == slrm$ID,
                   c("Row", "Weight", "Area.max", "d13C.cal",
                     "d13C.cal.se", "d18O.cal", "d18O.cal.se",
                     "pCO3", "pCO3.se")]
sam.data = d.good[!(d.good$Identifier_1 %in% c(plrm1$ID, plrm2$ID, slrm$ID)),
                  c("Row", "Identifier_1", "Weight", "Area.max", 
                    "d13C.cal", "d13C.cal.se", "d18O.cal", 
                    "d18O.cal.se", "pCO3", "pCO3.se")]

##excel object
wb = createWorkbook()
options("openxlsx.numFmt" = "0.00")
addWorksheet(wb, "Samples")
addWorksheet(wb, "Standards")
addWorksheet(wb, "All")
addWorksheet(wb, "Raw")

##rounding
int = createStyle(numFmt = "0")

##sample results
writeData(wb, "Samples", sam.data)
addStyle(wb, "Samples", int, rows = 1:100, cols = 1)

##standard results
y = 1
writeData(wb, "Standards", plrm1$ID)
y = y + 1
writeData(wb, "Standards", plrm1.data, startRow = y)
y = y + nrow(plrm1.data) + 1
writeData(wb, "Standards", "Known", startRow = y)
writeData(wb, "Standards", plrm1$d13C, startRow = y, startCol = 4)
writeData(wb, "Standards", plrm1$d18O, startRow = y, startCol = 6)
writeData(wb, "Standards", plrm1$pCO3, startRow = y, startCol = 8)
y = y + 2

writeData(wb, "Standards", plrm2$ID, startRow = y)
y = y + 1
writeData(wb, "Standards", plrm2.data, startRow = y)
y = y + nrow(plrm2.data) + 1
writeData(wb, "Standards", "Known", startRow = y)
writeData(wb, "Standards", plrm2$d13C, startRow = y, startCol = 4)
writeData(wb, "Standards", plrm2$d18O, startRow = y, startCol = 6)
writeData(wb, "Standards", plrm2$pCO3, startRow = y, startCol = 8)
y = y + 2

writeData(wb, "Standards", slrm$ID, startRow = y)
y = y + 1
writeData(wb, "Standards", slrm.data, startRow = y)
y = y + nrow(slrm.data) + 1
writeData(wb, "Standards", "Known", startRow = y)
writeData(wb, "Standards", slrm$d13C, startRow = y, startCol = 4)
writeData(wb, "Standards", slrm$d18O, startRow = y, startCol = 6)
writeData(wb, "Standards", slrm$pCO3, startRow = y, startCol = 8)

addStyle(wb, "Standards", int, rows = 1:100, cols = 1)


##all processed
writeData(wb, "All", d.good)
addStyle(wb, "All", int, rows = 1:100, cols = 1)

##all raw
writeData(wb, "Raw", d)
addStyle(wb, "Raw", int, rows = 1:100, cols = c(1, 6, 7),
         gridExpand = TRUE)

##save it
saveWorkbook(wb, "out/testing.xlsx", overwrite = TRUE)
