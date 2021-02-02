FitInlaAll <- function (forms, dat, fams = "zinb", logdisp = c(0, 0.01), precerr = c(1, 10^(-5)), curvedispfun = NULL, 
                        logitp0 = c(0, 0.01), ncpus = 2, 
            effoutput = TRUE, keepmargrand = FALSE, keepmarghyper = TRUE, 
            setthreads1 = TRUE, showupdate = FALSE, silentINLA = 2L, 
            updateby = 5000, ndigits = 5, addpackage = NULL, safemode = TRUE, 
            cf = NULL, designlist=NULL, ...) 
  {
  # forms <- formb;dat=datsim[1:3,];fams="gaussian";logdisp = c(0, 0.01); precerr = c(1, 10^(-5)); curvedispfun = NULL;
  # logitp0 = c(0, 0.01); ncpus = 1;
  # effoutput = TRUE; keepmargrand = FALSE; keepmarghyper = TRUE;
  # setthreads1 = TRUE; showupdate = FALSE; silentINLA = 2L;
  # updateby = 5000; ndigits = 5; addpackage = NULL; safemode = TRUE; cf = NULL
  # 
  # dat <- lapply(1:nrow(datsim[1:3,]),function(i) datsim[i,])
  # designlist <- lapply(1:nrow(datsim[1:3,]),function(i) group)
  
  
  
    print("NOTE: Warnings from INLA (eigenvalues, convergence, abort) can currently not be surpressed. Please ignore (generally)")
    cat("\n")
    if (setthreads1) 
      inla.setOption("num.threads", 1)
    if(class(dat)[1] =="list")  ngene <- length(dat) else ngene <- nrow(dat) #NEW
    
    if (mode(forms) == "call") 
      forms <- rep(list(forms), ngene)
    if (length(fams == 1)) 
      fams <- rep(fams, ngene)
    if (is.null(cf)) 
      cf <- list(prec.intercept = 0.001) else cf <- c(cf, prec.intercept = 0.001)
    form <- forms[[1]]
    
    
  #creates design from formula if designlist is unknown 
  if(is.null(designlist)){ 
    frmchr <- as.character(form)[[3]]
    sp <- strsplit(frmchr, "\\+")
    sp <- as.vector(sapply(sp, function(tt) {
      gsub(" ", "", tt)
    }))
    wht <- which(!is.na(match(sp, objects(envir = .GlobalEnv))))
    if (length(wht) > 0) {
      sp2 <- sp[wht]
      dfr <- data.frame(get(sp2[1]))
      if (length(sp2) > 1) {
        for (j in 2:length(sp2)) {
          dfr <- cbind(dfr, get(sp2[j]))
        }
      }
      names(dfr) <- sp2
    } else dfr <- data.frame()
    if (dim(dfr)[1] > 0) {
      elmiss <- length(which(is.na(dfr)))
      if (elmiss > 0) {
        print("Design contains missing fixed effects. Please use function 'ReDefMiss' first.")
        return(NULL)
      }
    }
  }
    fitinlaseqi <- function(i, ...) {
      #i<-1
      print(i)
      form <- forms[[i]]
      fam = fams[i]
      if(class(dat)[1] == "list") di <- as.numeric(dat[[i]]) else di <- as.numeric(dat[i, ])
      
      if(is.null(designlist)) {
        if (dim(dfr)[1] == 0) dattag <- data.frame(y = di) else dattag <- cbind(data.frame(y = di), dfr)
        } else { #design is different per tag
        dfr <- designlist[[i]]
        dattag <- cbind(data.frame(y = di), dfr)
      }
      if (!is.null(curvedispfun)) {
        mulogdisp <- curvedispfun(log(sum(di)))
        logdispi <- logdisp + c(mulogdisp, 0)
      } else {
        logdispi <- logdisp
      }
      if (fam == "zinb") {
        faminla = "zeroinflatednbinomial1"
        cd <- list(prior = c("gaussian", "gaussian"), param = c(logdispi, 
                                                                logitp0))
      }
      if (fam == "zip") {
        faminla = "zeroinflatedpoisson1"
        cd <- list(prior = "gaussian", param = logitp0)
      }
      if (fam == "nb") {
        faminla = "nbinomial"
        cd <- list(prior = "gaussian", param = logdispi)
      }
      if (fam == "poisson") {
        faminla = "poisson"
        cd <- list()
      }
      if (fam == "gaussian") {
        faminla = "gaussian"
        cd = list(hyper = list(prec = list(prior = "loggamma", 
                                           param = precerr)))
      }
      INLA:::inla.dynload.workaround() #NEW 11-1-2019
      result <- try(inla(formula = form, family = faminla, 
                         data = dattag, control.family = cd, silent = silentINLA, 
                         control.fixed = cf, ...))
      #mydfr <- data.frame(y=di,groupfac,batch,pers)
      #result <- try(inla(formula = form, family = faminla, data = mydfr, control.family = cd, silent = silentINLA, control.fixed = cf))
      if (class(result) == "try-error") 
        result <- NULL
      if (is.null(result$mlik)) {
        if ((faminla == "nbinomial" | faminla == "zeroinflatednbinomial1") & 
            max(di, na.rm = T) > 10^5) {
          maxval <- max(di, na.rm = T)
          if (maxval > 10^7) 
            di <- round(di/1000)
          if (maxval > 10^6 & maxval <= 10^7) 
            di <- round(di/100)
          if (maxval > 10^5 & maxval <= 10^6) 
            di <- round(di/10)
          if(is.null(designlist)) {
            if (dim(dfr)[1] == 0) dattag <- data.frame(y = di) else dattag <- cbind(data.frame(y = di), dfr)
          } else { #design is different per tag
            dfr <- designlist[[i]]
            dattag <- cbind(data.frame(y = di), dfr)
          }
          INLA:::inla.dynload.workaround() #NEW 11-1-2019
          result <- try(inla(formula = form, family = faminla, 
                             data = dattag, control.family = cd, silent = silentINLA, 
                             control.fixed = cf, ...))
        }
      }
      if (class(result) == "try-error") 
        result <- NULL else {
        if (effoutput) {
          nm = names(result)
          todel <- c(".control.defaults", "control.inla", 
                     "model.matrix", "lincomb", "control.expert", 
                     "control.mode", "control.results", "control.lincomb", 
                     "control.predictor", "joint.hyper", "misc", 
                     "logfile")
          if (!keepmargrand) 
            todel <- c(todel, "marginals.random")
          if (!keepmarghyper) 
            todel <- c(todel, "marginals.hyperpar", "internal.marginals.hyperpar")
          wh <- match(todel, nm)
          wh <- wh[!is.na(wh)]
          result <- result[-wh]
        }
        if (is.list(result)) 
          result <- rapply(result, signif, digits = ndigits, 
                           how = "replace", classes = "matrix")
      }
      return(list(result))
    } #END fitinlaseqi
    #f1 <- fitinlaseqi(1)
    
    
    if (ncpus == 1 | ngene == 1) {
      results <- sapply(1:ngene, fitinlaseqi, ...)
    }
    else {
      sfInit(parallel = TRUE, cpus = ncpus)
      sfLibrary(INLA)
      if (!is.null(addpackage)) {
        for (i in 1:length(addpackage)) {
          api <- addpackage[i]
          sfLibrary(api, character.only = TRUE)
        }
      }
      print("Exporting to slaves")
      mysfExport(forceexport = c("dat"))
      print("Exporting done")
      print("Started fitting")
      if (showupdate & updateby < ngene) {
        results <- as.list(rep(NA, ngene))
        sec <- seq(1, ngene, by = updateby)
        secl <- length(sec)
        if (sec[secl] > ngene - 2) 
          sec[secl] <- ngene + 1
        else sec <- c(sec, ngene + 1)
        secl <- length(sec)
        for (k in 1:(secl - 1)) {
          pmt <- proc.time()
          resultk <- sfSapply((sec[k]):(sec[k + 1] - 1), 
                              fitinlaseqi, ...)
          print(paste(sec[k + 1] - 1, "data rows done"))
          print(proc.time() - pmt)
          results[(sec[k]):(sec[k + 1] - 1)] <- resultk
        }
      }
      else {
        results <- sfSapply(1:ngene, fitinlaseqi, ...)
      }
    }
    if ((fams[1] == "zinb" | fams[1] == "nb") & safemode) {
      whichna <- which(is.na(mliks(results)))
      if (length(whichna) > 0) {
        logdisp <- c(logdisp[1], 100)
        if (length(whichna) == 1 | ncpus == 1) {
          newresults <- sapply(whichna, fitinlaseqi, ...)
        }
        else {
          sfExport("logdisp")
          newresults <- sfSapply(whichna, fitinlaseqi, 
                                 ...)
        }
        results[whichna] <- newresults
      }
    }
    if (ncpus > 1) {
      sfRemoveAll()
      sfStop()
    }
    return(results)
    if (setthreads1) 
      rm("inla.options")
  }