upsilon <- function(x, mediator, dv, conf.level = 0.95, bootstrap.lavaan = TRUE, bootstrap.lavaan.type = "ordinary", bootstrap.boot = FALSE, bootstrap.boot.type = "perc", B = 1000, boot.data.out=FALSE){
  if("lavaan" %in% rownames(installed.packages())){
    require(lavaan)
    data <- data.frame(x=x,mediator=mediator,dv=dv)
    med.model <- paste0(paste0(colnames(data)[2],'~',colnames(data)[1]),' \n ',paste0(colnames(data)[3],'~',colnames(data)[1],'+',colnames(data)[2]))
    med.res <- sem(med.model,data)
    upsilon <- coef(med.res)[1]**2*coef(med.res)[3]**2*inspectSampleCov(med.model,data)$cov['x','x']/inspectSampleCov(med.model,data)$cov['dv','dv']
    adj.upsilon <- (coef(med.res)[1]**2-vcov(med.res)[1,1])*(coef(med.res)[3]**2-vcov(med.res)[3,3])*inspectSampleCov(med.model,data)$cov['x','x']/inspectSampleCov(med.model,data)$cov['dv','dv']
    
    if(bootstrap.lavaan==TRUE){
      med.boot.fun <- function(out){
        require(lavaan)
        data.boot <- lavInspect(out,what = "data")
        colnames(data.boot) <- lavNames(out)
        upsilon.boot <- coef(out)[1]**2*coef(out)[3]**2*inspectSampleCov(med.model,data.boot)$cov['x','x']/inspectSampleCov(med.model,data.boot)$cov['dv','dv']
        adj.upsilon.boot <- (coef(out)[1]**2-vcov(out)[1,1])*(coef(out)[3]**2-vcov(out)[3,3])*inspectSampleCov(med.model,data.boot)$cov['x','x']/inspectSampleCov(med.model,data.boot)$cov['dv','dv']
        as.vector(c(upsilon.boot,adj.upsilon.boot))
      }
      cat('Bootstrapping may take several minutes \n \n')
      boot.med.res <- bootstrapLavaan(med.res,R=B,type=bootstrap.lavaan.type,FUN=med.boot.fun)
      
      upsES.out <- data.frame(est=c(upsilon,adj.upsilon),
                              lcl=c(ups.LCL=quantile(boot.med.res[,1],probs=(1-conf.level)/2),
                                    adj.ups.UCL=quantile(boot.med.res[,2],probs=(1-conf.level)/2)),
                              ucl=c(ups.LCL=quantile(boot.med.res[,1],probs=1-(1-conf.level)/2),
                                    adj.ups.UCL=quantile(boot.med.res[,2],probs=1-(1-conf.level)/2)),
                              row.names=c('Upsilon','Adj Upsilon'))
      colnames(upsES.out) <- c('Estimate',paste0(conf.level*100,'% ',bootstrap.lavaan.type,' LCL'),paste0(conf.level*100,'% ',bootstrap.lavaan.type,' UCL'))
      return(upsES.out)
    }
    if(bootstrap.boot==TRUE){
      if("boot" %in% rownames(installed.packages())){
        med.boot.fun <- function(out,i){
          require(boot)
          require(lavaan)
          data.boot <- out[i,]
          colnames(data.boot) <- colnames(out)
          out.boot <- sem(med.model,data.boot)
          upsilon.boot <- coef(out.boot)[1]**2*coef(out.boot)[3]**2*inspectSampleCov(med.model,data.boot)$cov['x','x']/inspectSampleCov(med.model,data.boot)$cov['dv','dv']
          adj.upsilon.boot <- (coef(out.boot)[1]**2-vcov(out.boot)[1,1])*(coef(out.boot)[3]**2-vcov(out.boot)[3,3])*inspectSampleCov(med.model,data.boot)$cov['x','x']/inspectSampleCov(med.model,data.boot)$cov['dv','dv']
          as.vector(c(upsilon.boot,adj.upsilon.boot))
        }
        cat('Bootstrapping may take several minutes \n \n')
        boot.med.res <- boot(data,med.boot.fun,R=B)
        upsES.out <- data.frame(est=c(upsilon,adj.upsilon),
                                lcl=c(ups.LCL=boot.ci(boot.med.res,conf=conf.level,type=bootstrap.boot.type,index=1)[[4]][4],
                                    adj.ups.UCL=boot.ci(boot.med.res,conf=conf.level,type=bootstrap.boot.type,index=2)[[4]][4]),
                                ucl=c(ups.LCL=boot.ci(boot.med.res,conf=conf.level,type=bootstrap.boot.type,index=1)[[4]][5],
                                    adj.ups.UCL=boot.ci(boot.med.res,conf=conf.level,type=bootstrap.boot.type,index=2)[[4]][5]),
                                row.names=c('Upsilon','Adj Upsilon'))
        colnames(upsES.out) <- c('Estimate',paste0(conf.level*100,'% ',bootstrap.boot.type,' LCL'),paste0(conf.level*100,'% ',bootstrap.boot.type,' UCL'))
        if(boot.data.out==TRUE){
          return(list(upsilon=upsES.out,bootstrap.data=boot.med.res))
                 } else return(upsES.out)
        } else stop('boot CIs requires package boot')
      
      } else {
        upsES.out <- data.frame(Est=c(upsilon,adj.upsilon),row.names=c('Upsilon','Adj Upsilon'))
        colnames(upsES.out) <- 'Estimate'
        return(upsES.out)
        }
    } else stop('upsilon requires package lavaan')
}
