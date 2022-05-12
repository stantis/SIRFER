library(readxl)
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
  "Identifier_2" = numeric(l),
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
  d.sam$Identifier_2[i] = d.no$`Identifier 2`[match(i, d.no$Row)]
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
cfit = cal(d.good, plrm1, plrm2)





plot(plrm1.data$`Area All`, plrm1.data$`d 18O/16O`)
plot(plrm2.data$`Area All`, plrm2.data$`d 13C/12C`)
plot(slrm.data$`Area All`, slrm.data$`d 13C/12C`)

for(i in unique(d$Row)){
  dsub = d[d$Row == i,]
  plot(dsub$Start, dsub$`d 13C/12C`, main = i)
  points(dsub$Start[dsub$Outliler], dsub$`d 13C/12C`[dsub$Outliler], col = "red")
}


