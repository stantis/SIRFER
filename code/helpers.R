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
