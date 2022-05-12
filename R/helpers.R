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
  dl = list(d.good[d.good$Identifier_1 == plrm1$ID,],
            d.good[d.good$Identifier_1 == plrm2$ID,])
  sl = d.good[d.good$Identifier_1 == slrm$ID,]
  
  for(i in 1:length(dl)){
    dl[[i]]$d13C.off = off(dl[[i]]$d13C.mean)
    dl[[i]]$d18O.off = off(dl[[i]]$d18O.mean)
  }
  sl$d13C.off = off(sl$d13C.mean)
  sl$d18O.off = off(sl$d18O.mean)
  
  r.all = c.all = o.all = numeric()
  for(i in dl){
    r.all = c(r.all, i$Row)
    c.all = c(c.all, i$d13C.off)
    o.all = c(o.all, i$d18O.off)
  }
  drr = range(c(r.all, sl$Row))
  drc = range(c(c.all, sl$d13C.off))
  dro = range(c(o.all, sl$d18O.off))

  plot(dl[[1]]$Row, dl[[1]]$d13C.off, xlim = drr, ylim = drc, pch = 21,
       bg = 1, main = expression(delta^{13}*"C drift"),
       xlab = "Sequence #", ylab = expression(delta^{13}*"C offset"))
  points(dl[[2]]$Row, dl[[2]]$d13C.off, bg = 2, pch = 21)
  points(sl$Row, sl$d13C.off, col = 3)

  ssc = smooth.spline(r.all, c.all, df = floor(length(r.all)/3))
  ssc.p = predict(ssc, seq(drr[1]:drr[2]))
  lines(ssc.p)
  points(sl$Row, sl$d13C.off - predict(ssc, sl$Row)$y, pch = 21, bg = 3)
  
  plot(dl[[1]]$Row, dl[[1]]$d18O.off, xlim = drr, ylim = dro, pch = 21,
       bg = 1, main = expression(delta^{18}*"O drift"),
       xlab = "Sequence #", ylab = expression(delta^{18}*"O offset"))
  points(dl[[2]]$Row, dl[[2]]$d18O.off, bg = 2, pch = 21)
  points(sl$Row, sl$d18O.off, col = 3)
  
  sso = smooth.spline(r.all, o.all, df = floor(length(r.all)/3))
  sso.p = predict(sso, seq(drr[1]:drr[2]))
  lines(sso.p)
  points(sl$Row, sl$d18O.off - predict(sso, sl$Row)$y, pch = 21, bg = 3)
  
  return(list(ssc, sso))
}

##calibration function
cal = function(dl, plrm1, plrm2){
  dl$d13C.known[dl$Identifier_1 == plrm1$ID] = plrm1$d13C
  dl$d13C.known[dl$Identifier_1 == plrm2$ID] = plrm2$d13C
  
  dl$d18O.known[dl$Identifier_1 == plrm1$ID] = plrm1$d18O
  dl$d18O.known[dl$Identifier_1 == plrm2$ID] = plrm2$d18O
  
  dl$CO3.known[dl$Identifier_1 == plrm1$ID] = plrm1$pCO3 * 
    dl$Identifier_2[dl$Identifier_1 == plrm1$ID]
  dl$CO3.known[dl$Identifier_1 == plrm2$ID] = plrm2$pCO3 * 
    dl$Identifier_2[dl$Identifier_1 == plrm2$ID]
  
  plot(dl$d13C.dc, dl$d13C.known, pch = 21, bg = 1, 
       main = expression(delta^{13}*"C calibration"), 
       xlab = expression(delta^{13}*"C measured"),
       ylab = expression(delta^{13}*"C known"))
  c.cal = lm(dl$d13C.known ~ dl$d13C.dc)
  abline(c.cal)
  pb = par("usr")
  cvals = round(coef(c.cal), 2)
  text(pb[1] + 0.05 * diff(pb[c(1,2)]),
       pb[4] - 0.1 * diff(pb[c(3,4)]),
       paste("y =", cvals[2], "* x +", cvals[1]), adj = 0)
  
  plot(dl$d18O.dc, dl$d18O.known, pch = 21, bg = 1, 
       main = expression(delta^{18}*"O calibration"), 
       xlab = expression(delta^{18}*"O measured"),
       ylab = expression(delta^{18}*"O known"))
  o.cal = lm(dl$d18O.known ~ dl$d18O.dc)
  abline(o.cal)
  pb = par("usr")
  cvals = round(coef(o.cal), 2)
  text(pb[1] + 0.05 * diff(pb[c(1,2)]),
       pb[4] - 0.1 * diff(pb[c(3,4)]),
       paste("y =", cvals[2], "* x +", cvals[1]), adj = 0)
  
  

  return(list(c.cal, o.cal))  
}