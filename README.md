

Timeline: 
* 2020/01/17: send Dr. Xiong drug responseincluding [chemotherapy](./extdata/clinical/Pancancer.radiationresponse2020.txt), [radiationtherapy](./extdata/clinical/Pancancer.radiationresponse2020.txt), immune and targeted therapy.
*  2020/01/14: Updated script, command and code can be found in my Rscript deposition, see and download [here](https://github.com/Shicheng-Guo/GscRbasement/blob/master/TcGaOverallDGEmeta.R)
*  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output,digits=5,control=list(stepadj=0.5))
*  2020/01/14: update md/rma command `digits` and `control` so that the model fitting can be converged as above
