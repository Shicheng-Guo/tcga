

Timeline: 

* [1470](/extdata/clinical/chemo-radio/) unique patients with responses to radiation-therapy and [1538](/extdata/clinical/chemo-radio/) patients with the response to chemotherapy. 
* When we merge both of them, there are [703](/extdata/clinical/chemo-radio/) unique patients with clear response to both chemo- and radio-therapy. 
* 2020/01/19: identify all the individual have both [chemotherapy]() and [radiotherapy]() treatment
* 2020/01/19:  [immune](./extdata/clinical/immunotherapy.txt) and [targeted therapy response](./extdata/clinical/TargetedMoleculartherapy.txt) phenotypes were collected with [script](./R/immunotherapy.tcga.R) and saved [here](./extdata/clinical/)
* 2020/01/17:  Dr. Xiong let me check do we have immune and targeted therapy response phenotypes. 
* 2020/01/17: send Dr. Xiong drug responseincluding [chemotherapy](./extdata/clinical/pancancer.chemotherapy.response2020.txt), [radiationtherapy](./extdata/clinical/Pancancer.radiationresponse2020.txt), immune and targeted therapy.
*  2020/01/14: Updated script, command and code can be found in my Rscript deposition, see and download [here](https://github.com/Shicheng-Guo/GscRbasement/blob/master/TcGaOverallDGEmeta.R)
*  md <- rma(es,slab=source,method = "REML", measure = "SMD",data=output,digits=5,control=list(stepadj=0.5))
*  2020/01/14: update md/rma command `digits` and `control` so that the model fitting can be converged as above
