library(readxl); library(openxlsx); library(stringr); library(dplyr); library(tidyr)

# Setting up functions ----------------------------------------------------

##outlier function
dout = function(dsub, thold = 3.5){
  ol = rep(NA, nrow(dsub))
  for(i in 1:nrow(dsub)){
    c.sd = sd(dsub$`d 13C/12C`[-i])
    o.sd = sd(dsub$`d 18O/16O`[-i])
    c.z = (dsub$`d 13C/12C`[i] - mean(dsub$`d 13C/12C`[-i])) / c.sd
    o.z = (dsub$`d 18O/16O`[i] - mean(dsub$`d 18O/16O`[-i])) / o.sd
    ol[i] = abs(c.z) > thold | abs(o.z) > thold
  }
  return(ol)
}

##offset function
off = function(d){
  return(d - mean(d, na.rm = TRUE))
}

##drift function
drift = function(d.good, plrm1, plrm2, slrm){
  #parse primary and secondary RMs
  dl = list(d.good[d.good$Identifier_1 == plrm1$ID,],
            d.good[d.good$Identifier_1 == plrm2$ID,])
  sl = d.good[d.good$Identifier_1 == slrm$ID,]
  
  #offsets from mean
  for(i in 1:length(dl)){
    dl[[i]]$d13C.off = off(dl[[i]]$d13C.mean)
    dl[[i]]$d18O.off = off(dl[[i]]$d18O.mean)
  }
  sl$d13C.off = off(sl$d13C.mean)
  sl$d18O.off = off(sl$d18O.mean)
  
  #bundle primary RMs
  r.all = c.all = o.all = numeric()
  for(i in dl){
    r.all = c(r.all, i$Row)
    c.all = c(c.all, i$d13C.off)
    o.all = c(o.all, i$d18O.off)
  }
  
  #ranges for plotting
  drr = range(c(r.all, sl$Row))
  drc = range(c(c.all, sl$d13C.off))
  dro = range(c(o.all, sl$d18O.off))
  
  #plot d13C drift
  plot(dl[[1]]$Row, dl[[1]]$d13C.off, xlim = drr, ylim = drc, pch = 21,
       bg = 1, main = expression(delta^{13}*"C drift"),
       xlab = "Sequence #", ylab = expression(delta^{13}*"C offset"))
  points(dl[[2]]$Row, dl[[2]]$d13C.off, bg = 2, pch = 21)
  points(sl$Row, sl$d13C.off, col = 3)
  text(dl[[1]]$Row, dl[[1]]$d13C.off, dl[[1]]$Row)
  text(dl[[2]]$Row, dl[[2]]$d13C.off, dl[[2]]$Row, col = 2)
  text(sl$Row, sl$d13C.off, sl$Row, col = 3)
  
  #spline fit
  ssc = smooth.spline(r.all, c.all, df = floor(length(r.all)/3))
  ssc.p = predict(ssc, seq(drr[1]:drr[2]))
  lines(ssc.p)
  points(sl$Row, sl$d13C.off - predict(ssc, sl$Row)$y, pch = 21, bg = 3)
  
  #plot d18O drift
  plot(dl[[1]]$Row, dl[[1]]$d18O.off, xlim = drr, ylim = dro, pch = 21,
       bg = 1, main = expression(delta^{18}*"O drift"),
       xlab = "Sequence #", ylab = expression(delta^{18}*"O offset"))
  points(dl[[2]]$Row, dl[[2]]$d18O.off, bg = 2, pch = 21)
  points(sl$Row, sl$d18O.off, col = 3)
  text(dl[[1]]$Row, dl[[1]]$d18O.off, dl[[1]]$Row)
  text(dl[[2]]$Row, dl[[2]]$d18O.off, dl[[2]]$Row, col = 2)
  text(sl$Row, sl$d18O.off, sl$Row, col = 3)
  
  #spline fit
  sso = smooth.spline(r.all, o.all, df = floor(length(r.all)/3))
  sso.p = predict(sso, seq(drr[1]:drr[2]))
  lines(sso.p)
  points(sl$Row, sl$d18O.off - predict(sso, sl$Row)$y, pch = 21, bg = 3)
  
  #report correction results for slrm
  cat("SLRM SD: C =", round(sd(sl$d13C.off), 2), 
      "O =", round(sd(sl$d18O.off), 2), "\n")
  cat("After drift: C =", 
      round(sd(sl$d13C.off - predict(ssc, sl$Row)$y), 2),
      "O =", round(sd(sl$d18O.off - predict(sso, sl$Row)$y), 2))
  
  return(list("ssc" = ssc, "sso" = sso))
}

##calibration function
cal = function(dl, plrm1, plrm2, slrm){
  #parse plrms
  dl$d13C.known[dl$Identifier_1 == plrm1$ID] = plrm1$d13C
  dl$d13C.known[dl$Identifier_1 == plrm2$ID] = plrm2$d13C
  
  dl$d18O.known[dl$Identifier_1 == plrm1$ID] = plrm1$d18O
  dl$d18O.known[dl$Identifier_1 == plrm2$ID] = plrm2$d18O
  
  dl$CO3.known[dl$Identifier_1 == plrm1$ID] = plrm1$pCO3 * 
    dl$Weight[dl$Identifier_1 == plrm1$ID]
  dl$CO3.known[dl$Identifier_1 == plrm2$ID] = plrm2$pCO3 * 
    dl$Weight[dl$Identifier_1 == plrm2$ID]
  
  #plot d13C and calibration fit
  plot(dl$d13C.dc, dl$d13C.known, pch = 21, bg = 1, 
       main = expression(delta^{13}*"C calibration"), 
       xlab = expression(delta^{13}*"C measured"),
       ylab = expression(delta^{13}*"C known"))
  c.cal = lm(d13C.known ~ d13C.dc, data = dl)
  abline(c.cal)
  pb = par("usr")
  cvals = round(coef(c.cal), 2)
  text(pb[1] + 0.05 * diff(pb[c(1,2)]),
       pb[4] - 0.15 * diff(pb[c(3,4)]),
       paste("y =", cvals[2], "* x +", cvals[1], 
             "\nRMSE =", round(sqrt(summary(c.cal)$sigma), 2)), 
       adj = 0)
  xs = data.frame("d13C.dc" = 
                    c(floor(range(dl$d13C.dc)[1]):
                        ceiling(range(dl$d13C.dc)[2])))
  cci = predict(c.cal, xs, interval = "confidence", level = 0.95)
  lines(xs$d13C.dc, cci[,2], lty = 3)
  lines(xs$d13C.dc, cci[,3], lty = 3)
  
  #plot d18O and calibration fit
  plot(dl$d18O.dc, dl$d18O.known, pch = 21, bg = 1, 
       main = expression(delta^{18}*"O calibration"), 
       xlab = expression(delta^{18}*"O measured"),
       ylab = expression(delta^{18}*"O known"))
  o.cal = lm(d18O.known ~ d18O.dc, data = dl)
  abline(o.cal)
  pb = par("usr")
  cvals = round(coef(o.cal), 2)
  text(pb[1] + 0.05 * diff(pb[c(1,2)]),
       pb[4] - 0.15 * diff(pb[c(3,4)]),
       paste("y =", cvals[2], "* x +", cvals[1], 
             "\nRMSE =", round(sqrt(summary(o.cal)$sigma), 2)), 
       adj = 0)
  xs = data.frame("d18O.dc" = 
                    c(floor(range(dl$d18O.dc)[1]):
                        ceiling(range(dl$d18O.dc)[2])))
  oci = predict(o.cal, xs, interval = "confidence", level = 0.95)
  lines(xs$d18O.dc, oci[,2], lty = 3)
  lines(xs$d18O.dc, oci[,3], lty = 3)
  
  #plot CO3 and calibration fit
  plot(dl$Area.max, dl$CO3.known, 
       xlim = c(0, max(dl$Area.max, na.rm = TRUE)),
       ylim = c(0, max(dl$CO3.known, na.rm = TRUE)),
       xlab = "Peak area (Vs)", 
       ylab = expression("CO"[3]*" (mg)"))
  kfit.free = lm(CO3.known ~ Area.max, data = dl)
  kfit.fix = lm(CO3.known ~ 0 + Area.max, data = dl)
  abline(kfit.free)
  abline(kfit.fix, lty = 2)
  
  #report results for SLRM
  d.s = dl[dl$Identifier_1 == slrm$ID,]
  d.s$d13C.cal = predict(c.cal, d.s)
  d.s$d18O.cal = predict(o.cal, d.s)
  d.s$pCO3 = predict(kfit.fix, d.s) / d.s$Weight
  
  cat("SLRM d13C: mean:", round(mean(d.s$d13C.cal), 2),
      "sd:", round(sd(d.s$d13C.cal), 2), "known:", 
      slrm$d13C, "\n")
  cat("SLRM d18O: mean:", round(mean(d.s$d18O.cal), 2),
      "sd:", round(sd(d.s$d18O.cal), 2), "known:", 
      slrm$d18O, "\n")
  cat("SLRM %CO3: mean:", round(mean(d.s$pCO3), 2),
      "sd:", round(sd(d.s$pCO3), 2), "known:", 
      slrm$pCO3, "\n")
  
  return(list("c.cal" = c.cal, "o.cal" = o.cal, 
              "k.cal" = kfit.fix))  
}

# Cleaning raw data -------------------------------------------------------

path <- file.path(getwd(), 'input')
files <- list.files(path, full.names = F)

for (file in files) {
  d =  read_excel(paste0('input/',file), 
                  col_types = c("numeric", "text", "numeric", 
                                "text", "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "text", "text", 
                                "text", "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "text", "text"))
  
  #remove on-off
  d = d[d$Method != "CO2 On Off",]
  
  #re-register run sequece #s
  d$Row = d$Row - min(d$Row) + 1
  
  #remove conditioners
  d = d[d$`Identifier 1` != "COND",]
  
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
  d.no = d.no[!is.na(d.no$Outliler),]
  
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
  for(i in seq_along(d.sam$Row)){
    d.sam$Identifier_1[i] = d.no$`Identifier 1`[match(d.sam$Row[i], 
                                                      d.no$Row)]
    d.sam$Weight[i] = d.no$`Identifier 2`[match(d.sam$Row[i], 
                                                d.no$Row)]
    #areas include outliers for standardization over all peaks
    d.sam$Area.max[i] = max(d$`Area All`[d$Row == d.sam$Row[i]])
    d.sam$Area.sd[i] = sd(d$`Area All`[d$Row == d.sam$Row[i]])
    #isotopes exclude outliers
    d.sam$d13C.mean[i] = mean(d.no$`d 13C/12C`[d.no$Row == 
                                                 d.sam$Row[i]])
    d.sam$d13C.sd[i] = sd(d.no$`d 13C/12C`[d.no$Row == d.sam$Row[i]])
    d.sam$d18O.mean[i] = mean(d.no$`d 18O/16O`[d.no$Row == 
                                                 d.sam$Row[i]])
    d.sam$d18O.sd[i] = sd(d.no$`d 18O/16O`[d.no$Row == d.sam$Row[i]])
  }
  
  #screen for high SDs
  d.sam$Ignore = d.sam$d13C.sd > 0.25 | d.sam$d18O.sd > 0.3
  d.good = d.sam[!d.sam$Ignore,]
  
  #some values to use
  plrm1 = list("ID" = "Carrara", "d13C" = 2.1, "d18O" = -1.8, "pCO3" = 0.6)
  #plrm2 = list("ID" = "LSVEC", "d13C" = -46.6, "d18O" = -26.7, "pCO3" = 0.6)
  plrm2 = list("ID" = "CO-8", "d13C" = -5.764, "d18O" = -22.7, "pCO3" = 0.6)
  slrm = list("ID" = "Marble", "d13C" = 1.9, "d18O" = -11.3, "pCO3" = 0.6)
  
  #drift correction
  ##fits and plots the spine
  dfit = drift(d.good, plrm1, plrm2, slrm)
  
  ##apply the drift correction 
  cd = TRUE
  if(cd){
    d.good$d13C.dc = d.good$d13C.mean - predict(dfit[[1]], d.good$Row)$y
    d.good$d18O.dc = d.good$d18O.mean - predict(dfit[[2]], d.good$Row)$y
  } else{
    d.good$d13C.dc = d.good$d13C.mean
    d.good$d18O.dc = d.good$d18O.mean
  }
  
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
  saveWorkbook(wb, file.path("xlsxOutput", paste0(file, "")), overwrite = TRUE)
  
}

# Cleaning Output ---------------------------------------------------------

path <- file.path(getwd(), 'xlsxOutput')
files <- list.files(path, full.names = F)

df <- data.frame() 
for (file in files) {
  r <- read_excel(paste0('xlsxOutput/',file))
  r$batch_id <- basename(file)
  df <- rbind(df, r)
}

df$batch_id <- df$batch_id %>% str_replace(".xlsx", "") #remove .xlsx from end of batch_id

df$participant_id <- df$sample_id %>% substr(0,11) #this will only work when FINDEM sample_ids are properly set up, so take care

#renaming some variables and removing some unneccessary ones
df <-  df %>% 
  rename(d13C = d13C.cal, # can't have .s in SQL files it will cause all sorts of trouble
         d18O = d18O.cal, 
         se_O = d18O.cal.se, 
         se_C = d13C.cal.se, 
         se_CO3 = pCO3.se, 
         weight = Weight, 
         sample_id = Identifier_1, 
         area_max = Area.max) %>% 
  select(-c(Row)) %>% 
  filter(str_detect(sample_id, "LAR|FIND-EM|TSU-|CMU")) #get rid of other samples run in batches

write.csv(df, "SQLFiles/isotopes.csv")
