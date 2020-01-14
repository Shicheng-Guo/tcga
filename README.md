

Timeline: 

*  es<-escalc(m1i=m1i, sd1i=sd1i, n1i=n1i, m2i=m2i, sd2i=sd2i, n2i=n2i,measure="MD",data=output)
*  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output,digits=5,control=list(stepadj=0.5))
