
library(lme4)
library(optimx)

# main function
scebe = function(nullmodel = NULL,     # fitted null model from lme4.
                 Phenofile,            # a file path to read in the phenotype.
                 IDcolname,            # a character to specifie the column name of the subject ID.
                 Phenocolname,         # a character to specifie the column name of the phenotype.
                 Covacolname,          # a character (vector) to specifie the column name of the covariates (except for age).
                 Timecolname,          # a character to specifie the column name of the age.
                 Genofile,             # a file path to read in the genotype (without file suffix like ".bim", "bed" or "fam").
                 GenofileIndex = NULL, # c("prefix.bim", "prefix.fam").
                 SNPIDfile = NULL,     # a file path to read in the SNP to be tested.
                 MinMAF = 0.01,        # a cutoff for variants with MAF > MinMAF to be tested.
                 Outputfile = NULL)    # a file path to store the output (with file suffix like ".txt").
{
  cat("Read in the phenotype...\n")
  if(is.data.frame(Phenofile))
  {
    Pheno_data = Phenofile
  }else
  {
    Pheno_data = data.table::fread(Phenofile)
  }
  
  cmd = paste0("Pheno_data = Pheno_data %>% arrange(", IDcolname, ")"); eval(parse(text = cmd))
  
  if(is.null(nullmodel))
  {
    cat("Model formula is...\n")
    model_formula = as.formula(paste(Phenocolname, "~", paste(Covacolname, collapse = "+"), "+", Timecolname, "+ (", Timecolname, "|", IDcolname, ")"))
    print(model_formula)
    
    cat("Fit the null model...\n")
    null_model = lmer(model_formula, data = Pheno_data, REML = FALSE,
                      control = lmerControl(optimizer ='optimx', optCtrl = list(method = 'nlminb')))
  }else
  {
    null_model = nullmodel
  }
  
  cat("Process the null model product...\n")
  
  varcor = VarCorr(null_model)
  cmd = paste0("varcor$", IDcolname); R = matrix(as.numeric(eval(parse(text = cmd))), ncol = 2)
  Ip=diag(1,2)
  
  #allS=allW=list()
  #allM=allISMIS=list()	  
  alls11=alls12=alls21=alls22=NULL
  allw11=allw12=allw21=allw22=NULL
  allsms11=allsms12=allsms21=allsms22=NULL
  
  Time.lst = split(Pheno_data[[Timecolname]], Pheno_data[[IDcolname]])
  N = length(Time.lst); n_measure = unlist(lapply(Time.lst,length))
  cmd = paste0("unique(Pheno_data$", IDcolname, ")"); SubjID = eval(parse(text = cmd)); SubjID = as.character(SubjID)
  
  for(i in 1:N)
  {
    if(i %% 10000 == 0) cat("Process", i, "/", N, "\n")
    
    nt = n_measure[i]
    time <- Time.lst[[i]] #Pheno_data[Pheno_data$ID ==i, Time]##
    Zi <- cbind(1, time)
    Gi = diag(sigma(null_model) ^ 2, nt)
    ZinvGZ <- crossprod(Zi, solve(Gi, Zi)) ##when n measure==1, return error!
    Mi=R+MASS::ginv(ZinvGZ)
    
    SIGMAi=Zi%*%R%*%t(Zi)+Gi
    Wi <- crossprod(Zi, solve(SIGMAi, Zi))
    #allW[[i]]=Wi
    allw11=c(allw11,Wi[1,1])
    allw12=c(allw12,Wi[1,2])
    allw21=c(allw21,Wi[2,1])
    allw22=c(allw22,Wi[2,2])
    
    invSi <- R %*% ZinvGZ + Ip
    Si=solve(invSi)
    alls11<-c(alls11,Si[1,1])
    alls12<-c(alls12,Si[1,2])
    alls21<-c(alls21,Si[2,1])
    alls22<-c(alls22,Si[2,2])
    
    Ipi=diag(1,2)
    SMS=(Ipi-Si)%*%Mi%*%t(Ipi-Si)
    
    allsms11<-c(allsms11,SMS[1,1])
    allsms12<-c(allsms12,SMS[1,2])
    allsms21<-c(allsms21,SMS[2,1])
    allsms22<-c(allsms22,SMS[2,2])
  }
  
  C = matrix(c(sum(allw11),sum(allw12),sum(allw21),sum(allw22)),2,byrow=T)
  
  inv_C = solve(C)
  
  invc11<-inv_C[1,1]
  invc12<-inv_C[1,2]
  invc21<-inv_C[2,1]
  invc22<-inv_C[2,2]
  
  cat("Start screen the whole genome...\n")
  
  if(is.null(GenofileIndex))
  {
    bedfile = paste0(Genofile, ".bed")
    bimfile = paste0(Genofile, ".bim")
    
    if(is.null(SNPIDfile))
    {
      GenoMatInfo = GRAB::GRAB.ReadGeno(GenoFile = bedfile,
                                        SampleIDs = SubjID,
                                        control = list(AllMarkers = TRUE,
                                                       ImputeMethod = "mean"))
    }else if(length(SNPIDfile) == 1)
    {
      GenoMatInfo = GRAB::GRAB.ReadGeno(GenoFile = bedfile,
                                        SampleIDs = SubjID,
                                        control = list(IDsToIncludeFile = SNPIDfile,
                                                       ImputeMethod = "mean"))
    }else
    {
      cmd = paste0("ranef(null_model)$", IDcolname); rf = eval(parse(text = cmd))
      names(rf)[c(1:2)] = c("Inter","Slope")
      
      results = c()
      for(i in 1:length(SNPIDfile))
      {
        cat(i, "\n")
        
        GenoMatInfo = GRAB::GRAB.ReadGeno(GenoFile = bedfile,
                                          SampleIDs = SubjID,
                                          control = list(IDsToIncludeFile = SNPIDfile[i],
                                                         ImputeMethod = "mean"))
        
        means = colMeans(GenoMatInfo$GenoMat)/2
        
        GenoMatInfo$GenoMat = GenoMatInfo$GenoMat[, which(means > MinMAF)]; gc()
        GenoMatInfo$markerInfo = GenoMatInfo$markerInfo[which(means > MinMAF), ]; gc()
        
        spr = slmm(y1 = rf$Inter, y2 = rf$Slope, X = GenoMatInfo$GenoMat)
        ebe = spr[,c(1,4)]
        
        Xc<-apply(GenoMatInfo$GenoMat, 2, function(x){scale(x,scale=F)})
        
        XcX<-Xc*GenoMatInfo$GenoMat
        Xc2<-Xc^2
        Sxc2<-colSums(Xc2)
        
        xstar11<-colSums(allw11*Xc)
        xstar12<-colSums(allw12*Xc)
        xstar21<-colSums(allw21*Xc)
        xstar22<-colSums(allw22*Xc)
        
        A11<-colSums(Xc2*alls11)/Sxc2
        A12<-colSums(Xc2*alls12)/Sxc2
        A21<-colSums(Xc2*alls21)/Sxc2
        A22<-colSums(Xc2*alls22)/Sxc2
        
        B11<-colSums(Xc*alls11)/Sxc2
        B12<-colSums(Xc*alls12)/Sxc2
        B21<-colSums(Xc*alls21)/Sxc2
        B22<-colSums(Xc*alls22)/Sxc2
        
        Sxc22<-(Sxc2)^2
        
        D11=colSums(Xc2*allsms11)/Sxc22
        D12=colSums(Xc2*allsms12)/Sxc22
        D21=colSums(Xc2*allsms21)/Sxc22
        D22=colSums(Xc2*allsms22)/Sxc22
        
        ssh211<-B11*invc11*xstar11+B12*invc21*xstar11+B11*invc12*xstar21+B12*invc22*xstar21
        ssh212<-B11*invc11*xstar12+B12*invc21*xstar12+B11*invc12*xstar22+B12*invc22*xstar22
        ssh221<-B21*invc11*xstar11+B22*invc21*xstar11+B21*invc12*xstar21+B22*invc22*xstar21
        ssh222<-B21*invc11*xstar12+B22*invc21*xstar12+B21*invc12*xstar22+B22*invc22*xstar22
        
        ssh11=1-A11+ssh211
        ssh12=0-A12+ssh212
        ssh21=0-A21+ssh221
        ssh22=1-A22+ssh222
        
        detssh=ssh11*ssh22-ssh12*ssh21
        
        invssh11<-ssh22/detssh
        invssh12<-(-ssh12/detssh)
        invssh21<-(-ssh21/detssh)
        invssh22<-ssh11/detssh
        
        Vgamma211<-B11*invc11*B11+B12*invc21*B11+B11*invc12*B12+B12*invc22*B12
        #Vgamma211<-Vgamma211/SXc2
        Vgamma212<-B11*invc11*B21+B12*invc21*B21+B11*invc12*B22+B12*invc22*B22
        #Vgamma212<-Vgamma212/SXc2
        Vgamma221<-B21*invc11*B11+B22*invc21*B11+B21*invc12*B12+B22*invc22*B12
        #Vgamma221<-Vgamma221/SXc2
        Vgamma222<-B21*invc11*B21+B22*invc21*B21+B21*invc12*B22+B22*invc22*B22
        #Vgamma222<-Vgamma222/SXc2
        
        Vgamma11=D11-Vgamma211
        Vgamma12=D12-Vgamma212
        Vgamma21=D21-Vgamma221
        Vgamma22=D22-Vgamma222 ###VAR of NEBE
        
        Vssh11<-invssh11*Vgamma11*invssh11+invssh12*Vgamma21*invssh11+invssh11*Vgamma12*invssh12+invssh12*Vgamma22*invssh12
        Vssh12<-invssh11*Vgamma11*invssh21+invssh12*Vgamma21*invssh21+invssh11*Vgamma12*invssh22+invssh12*Vgamma22*invssh22
        Vssh21<-invssh21*Vgamma11*invssh11+invssh22*Vgamma21*invssh11+invssh21*Vgamma12*invssh12+invssh22*Vgamma22*invssh12
        Vssh22<-invssh21*Vgamma11*invssh21+invssh22*Vgamma21*invssh21+invssh21*Vgamma12*invssh22+invssh22*Vgamma22*invssh22
        
        est.ssh.inter<-invssh11*ebe[,1]+invssh12*ebe[,2]
        est.ssh.slope<-invssh21*ebe[,1]+invssh22*ebe[,2]
        
        t.ssh.inter<-est.ssh.inter/sqrt(Vssh11)
        t.ssh.slope<-est.ssh.slope/sqrt(Vssh22)
        
        p.ssh.inter<-2*pnorm(abs(t.ssh.inter),lower.tail=F)
        p.ssh.slope<-2*pnorm(abs(t.ssh.slope),lower.tail=F)
        
        Output = cbind(GenoMatInfo$markerInfo, 
                       scebeCSpval = p.ssh.inter, ebeCSpval = spr[,3], 
                       scebeLTpval = p.ssh.slope, ebeLTpval = spr[,6])
        
        results = rbind(results, Output); gc()
      }
      
      data.table::fwrite(results, file = Outputfile,
                         row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
      
      return("Complete the analysis!")
    }
  }else
  {
    if(is.null(SNPIDfile))
    {
      GenoMatInfo = GRAB::GRAB.ReadGeno(GenoFile = Genofile,
                                        GenoFileIndex = GenofileIndex,
                                        SampleIDs = SubjID,
                                        control = list(AllMarkers = TRUE,
                                                       ImputeMethod = "mean"))
    }else
    {
      SNPIDs = data.table::fread(SNPIDfile, header = FALSE); SNPIDs = SNPIDs$V1
      
      if(length(SNPIDs) > 500)
      {
        nGroups = length(SNPIDs) %/% 500
        SNPIDsGroups = split(SNPIDs, cut(seq_along(SNPIDs), nGroups, labels = FALSE))
        
        cmd = paste0("ranef(null_model)$", IDcolname); rf = eval(parse(text = cmd))
        names(rf)[c(1:2)] = c("Inter","Slope")
        
        results = c()
        for(i in 1:length(SNPIDsGroups))
        {
          cat(i, "\n")
          
          temp_SNPIDfile = paste0(Outputfile, "_", i, ".txt")
          
          data.table::fwrite(data.frame(SNPIDsGroups[[i]]), file = temp_SNPIDfile, 
                             row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
          
          GenoMatInfo = GRAB::GRAB.ReadGeno(GenoFile = Genofile,
                                            GenoFileIndex = GenofileIndex,
                                            SampleIDs = SubjID,
                                            control = list(IDsToIncludeFile = temp_SNPIDfile,
                                                           ImputeMethod = "mean"))
          
          means = colMeans(GenoMatInfo$GenoMat)/2
          
          GenoMatInfo$GenoMat = GenoMatInfo$GenoMat[, which(means > MinMAF)]; gc()
          GenoMatInfo$markerInfo = GenoMatInfo$markerInfo[which(means > MinMAF), ]; gc()
          
          spr = slmm(y1 = rf$Inter, y2 = rf$Slope, X = GenoMatInfo$GenoMat)
          ebe = spr[,c(1,4)]
          
          Xc<-apply(GenoMatInfo$GenoMat, 2, function(x){scale(x,scale=F)})
          
          XcX<-Xc*GenoMatInfo$GenoMat
          Xc2<-Xc^2
          Sxc2<-colSums(Xc2)
          
          xstar11<-colSums(allw11*Xc)
          xstar12<-colSums(allw12*Xc)
          xstar21<-colSums(allw21*Xc)
          xstar22<-colSums(allw22*Xc)
          
          A11<-colSums(Xc2*alls11)/Sxc2
          A12<-colSums(Xc2*alls12)/Sxc2
          A21<-colSums(Xc2*alls21)/Sxc2
          A22<-colSums(Xc2*alls22)/Sxc2
          
          B11<-colSums(Xc*alls11)/Sxc2
          B12<-colSums(Xc*alls12)/Sxc2
          B21<-colSums(Xc*alls21)/Sxc2
          B22<-colSums(Xc*alls22)/Sxc2
          
          Sxc22<-(Sxc2)^2
          
          D11=colSums(Xc2*allsms11)/Sxc22
          D12=colSums(Xc2*allsms12)/Sxc22
          D21=colSums(Xc2*allsms21)/Sxc22
          D22=colSums(Xc2*allsms22)/Sxc22
          
          ssh211<-B11*invc11*xstar11+B12*invc21*xstar11+B11*invc12*xstar21+B12*invc22*xstar21
          ssh212<-B11*invc11*xstar12+B12*invc21*xstar12+B11*invc12*xstar22+B12*invc22*xstar22
          ssh221<-B21*invc11*xstar11+B22*invc21*xstar11+B21*invc12*xstar21+B22*invc22*xstar21
          ssh222<-B21*invc11*xstar12+B22*invc21*xstar12+B21*invc12*xstar22+B22*invc22*xstar22
          
          ssh11=1-A11+ssh211
          ssh12=0-A12+ssh212
          ssh21=0-A21+ssh221
          ssh22=1-A22+ssh222
          
          detssh=ssh11*ssh22-ssh12*ssh21
          
          invssh11<-ssh22/detssh
          invssh12<-(-ssh12/detssh)
          invssh21<-(-ssh21/detssh)
          invssh22<-ssh11/detssh
          
          Vgamma211<-B11*invc11*B11+B12*invc21*B11+B11*invc12*B12+B12*invc22*B12
          #Vgamma211<-Vgamma211/SXc2
          Vgamma212<-B11*invc11*B21+B12*invc21*B21+B11*invc12*B22+B12*invc22*B22
          #Vgamma212<-Vgamma212/SXc2
          Vgamma221<-B21*invc11*B11+B22*invc21*B11+B21*invc12*B12+B22*invc22*B12
          #Vgamma221<-Vgamma221/SXc2
          Vgamma222<-B21*invc11*B21+B22*invc21*B21+B21*invc12*B22+B22*invc22*B22
          #Vgamma222<-Vgamma222/SXc2
          
          Vgamma11=D11-Vgamma211
          Vgamma12=D12-Vgamma212
          Vgamma21=D21-Vgamma221
          Vgamma22=D22-Vgamma222 ###VAR of NEBE
          
          Vssh11<-invssh11*Vgamma11*invssh11+invssh12*Vgamma21*invssh11+invssh11*Vgamma12*invssh12+invssh12*Vgamma22*invssh12
          Vssh12<-invssh11*Vgamma11*invssh21+invssh12*Vgamma21*invssh21+invssh11*Vgamma12*invssh22+invssh12*Vgamma22*invssh22
          Vssh21<-invssh21*Vgamma11*invssh11+invssh22*Vgamma21*invssh11+invssh21*Vgamma12*invssh12+invssh22*Vgamma22*invssh12
          Vssh22<-invssh21*Vgamma11*invssh21+invssh22*Vgamma21*invssh21+invssh21*Vgamma12*invssh22+invssh22*Vgamma22*invssh22
          
          est.ssh.inter<-invssh11*ebe[,1]+invssh12*ebe[,2]
          est.ssh.slope<-invssh21*ebe[,1]+invssh22*ebe[,2]
          
          t.ssh.inter<-est.ssh.inter/sqrt(Vssh11)
          t.ssh.slope<-est.ssh.slope/sqrt(Vssh22)
          
          p.ssh.inter<-2*pnorm(abs(t.ssh.inter),lower.tail=F)
          p.ssh.slope<-2*pnorm(abs(t.ssh.slope),lower.tail=F)
          
          Output = cbind(GenoMatInfo$markerInfo, 
                         scebeCSpval = p.ssh.inter, ebeCSpval = spr[,3], 
                         scebeLTpval = p.ssh.slope, ebeLTpval = spr[,6])
          
          results = rbind(results, Output)
          
          file.remove(temp_SNPIDfile); rm(GenoMatInfo); gc()
        }
        
        if(is.null(Outputfile))
        {
          return(results)
        }else
        {
          data.table::fwrite(results, file = Outputfile,
                             row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        }
        
        return("Complete the analysis!")
        
      }else
      {
        GenoMatInfo = GRAB::GRAB.ReadGeno(GenoFile = Genofile,
                                          GenoFileIndex = GenofileIndex,
                                          SampleIDs = SubjID,
                                          control = list(IDsToIncludeFile = SNPIDfile,
                                                         ImputeMethod = "mean"))
      }
    }
  }
  
  means = colMeans(GenoMatInfo$GenoMat)/2
  
  GenoMatInfo$GenoMat = GenoMatInfo$GenoMat[, which(means > MinMAF)]; gc()
  GenoMatInfo$markerInfo = GenoMatInfo$markerInfo[which(means > MinMAF), ]; gc()
  
  cmd = paste0("ranef(null_model)$", IDcolname); rf = eval(parse(text = cmd))
  names(rf)[c(1:2)] = c("Inter","Slope")
  spr = slmm(y1 = rf$Inter, y2 = rf$Slope, X = GenoMatInfo$GenoMat)
  ebe = spr[,c(1,4)]
  
  Xc<-apply(GenoMatInfo$GenoMat, 2, function(x){scale(x,scale=F)})
  
  XcX<-Xc*GenoMatInfo$GenoMat
  Xc2<-Xc^2
  Sxc2<-colSums(Xc2)
  
  xstar11<-colSums(allw11*Xc)
  xstar12<-colSums(allw12*Xc)
  xstar21<-colSums(allw21*Xc)
  xstar22<-colSums(allw22*Xc)
  
  A11<-colSums(Xc2*alls11)/Sxc2
  A12<-colSums(Xc2*alls12)/Sxc2
  A21<-colSums(Xc2*alls21)/Sxc2
  A22<-colSums(Xc2*alls22)/Sxc2
  
  B11<-colSums(Xc*alls11)/Sxc2
  B12<-colSums(Xc*alls12)/Sxc2
  B21<-colSums(Xc*alls21)/Sxc2
  B22<-colSums(Xc*alls22)/Sxc2
  
  Sxc22<-(Sxc2)^2
  
  D11=colSums(Xc2*allsms11)/Sxc22
  D12=colSums(Xc2*allsms12)/Sxc22
  D21=colSums(Xc2*allsms21)/Sxc22
  D22=colSums(Xc2*allsms22)/Sxc22
  
  ssh211<-B11*invc11*xstar11+B12*invc21*xstar11+B11*invc12*xstar21+B12*invc22*xstar21
  ssh212<-B11*invc11*xstar12+B12*invc21*xstar12+B11*invc12*xstar22+B12*invc22*xstar22
  ssh221<-B21*invc11*xstar11+B22*invc21*xstar11+B21*invc12*xstar21+B22*invc22*xstar21
  ssh222<-B21*invc11*xstar12+B22*invc21*xstar12+B21*invc12*xstar22+B22*invc22*xstar22
  
  ssh11=1-A11+ssh211
  ssh12=0-A12+ssh212
  ssh21=0-A21+ssh221
  ssh22=1-A22+ssh222
  
  detssh=ssh11*ssh22-ssh12*ssh21
  
  invssh11<-ssh22/detssh
  invssh12<-(-ssh12/detssh)
  invssh21<-(-ssh21/detssh)
  invssh22<-ssh11/detssh
  
  Vgamma211<-B11*invc11*B11+B12*invc21*B11+B11*invc12*B12+B12*invc22*B12
  #Vgamma211<-Vgamma211/SXc2
  Vgamma212<-B11*invc11*B21+B12*invc21*B21+B11*invc12*B22+B12*invc22*B22
  #Vgamma212<-Vgamma212/SXc2
  Vgamma221<-B21*invc11*B11+B22*invc21*B11+B21*invc12*B12+B22*invc22*B12
  #Vgamma221<-Vgamma221/SXc2
  Vgamma222<-B21*invc11*B21+B22*invc21*B21+B21*invc12*B22+B22*invc22*B22
  #Vgamma222<-Vgamma222/SXc2
  
  Vgamma11=D11-Vgamma211
  Vgamma12=D12-Vgamma212
  Vgamma21=D21-Vgamma221
  Vgamma22=D22-Vgamma222 ###VAR of NEBE
  
  Vssh11<-invssh11*Vgamma11*invssh11+invssh12*Vgamma21*invssh11+invssh11*Vgamma12*invssh12+invssh12*Vgamma22*invssh12
  Vssh12<-invssh11*Vgamma11*invssh21+invssh12*Vgamma21*invssh21+invssh11*Vgamma12*invssh22+invssh12*Vgamma22*invssh22
  Vssh21<-invssh21*Vgamma11*invssh11+invssh22*Vgamma21*invssh11+invssh21*Vgamma12*invssh12+invssh22*Vgamma22*invssh12
  Vssh22<-invssh21*Vgamma11*invssh21+invssh22*Vgamma21*invssh21+invssh21*Vgamma12*invssh22+invssh22*Vgamma22*invssh22
  
  est.ssh.inter<-invssh11*ebe[,1]+invssh12*ebe[,2]
  est.ssh.slope<-invssh21*ebe[,1]+invssh22*ebe[,2]
  
  t.ssh.inter<-est.ssh.inter/sqrt(Vssh11)
  t.ssh.slope<-est.ssh.slope/sqrt(Vssh22)
  
  p.ssh.inter<-2*pnorm(abs(t.ssh.inter),lower.tail=F)
  p.ssh.slope<-2*pnorm(abs(t.ssh.slope),lower.tail=F)
  
  Output = cbind(GenoMatInfo$markerInfo, 
                 scebeCSpval = p.ssh.inter, ebeCSpval = spr[,3], 
                 scebeLTpval = p.ssh.slope, ebeLTpval = spr[,6])
  
  if(is.null(Outputfile))
  {
    return(Output)
  }else
  {
    data.table::fwrite(Output, file = Outputfile,
                       row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  }
  
  print("Complete the analysis!")
}


slmm = function(y1,y2, X)
{
  n1 = length(y1)
  n = n2 = n1
  
  y1bar = mean(y1)
  y2bar = mean(y2)
  
  y1c = y1 - y1bar
  y2c = y2 - y2bar
  s1 = colSums(X)
  s2 = colSums(X ^ 2)
  sxx = (s2 - (s1 ^ 2) / n)
  
  b1 = as.vector(colSums(y1c*X) / sxx)
  b2 = as.vector(colSums(y2c*X) / sxx)
  y1_sum = sum(y1)
  y2_sum = sum(y2)
  y1X = colSums(y1 * X)
  y2X = colSums(y2 * X)
  
  #Residuals1 <- (y1 - X %*% diag(b1) - t(replicate(n1, a1))) ## still potential for improvement
  #b1_sd <- sqrt((colSums(Residuals1 ^ 2) / (n1 - 2)) / sxx)
  #Residuals2 <- (y2 - X %*% diag(b2) - t(replicate(n2, a2))) ## still potential for improvement
  #b2_sd <- sqrt((colSums(Residuals2 ^ 2) / (n2 - 2)) / sxx)
  
  sigma11 = sum(y1 ^ 2) - (y1_sum ^ 2 * s2 - 2 * s1 * y1_sum * y1X + 
                             n * (y1X) ^ 2) / (n * sxx)
  var.gamma11 = sigma11 / (sxx * (n - 2))
  
  #sigma12 = sum(y1 * y2) - (y1_sum * y2_sum * s2 - y2_sum * s1 * y1X - 
  #                           y1_sum * s1 * y2X + n * y1X * y2X) / (n * sxx)
  #var.gamma12 = sigma12 / (sxx * (n - 2))
  
  sigma22 = sum(y2 ^ 2) - (y2_sum ^ 2 * s2 - 2 * s1 * y2_sum * y2X + 
                             n * (y2X) ^ 2) / (n * sxx)
  var.gamma22 = sigma22 / (sxx * (n - 2))
  
  b1_sd = sqrt(var.gamma11)
  b2_sd = sqrt(var.gamma22)
  
  t1_value = b1 / b1_sd
  p1_value = pt(q = -abs(t1_value), df =  n1 - 2) * 2 # Pr(>|t|)
  
  t2_value = b2 / b2_sd
  p2_value = pt(q = -abs(t2_value), df =  n2 - 2) * 2 # Pr(>|t|)
  
  result = cbind(b1, b1_sd, p1_value,b2, b2_sd,  p2_value)
  colnames(result) = c("b1_Est", "b1_Std",  "Pr1(>|t|)","b2_Est", "b2_Std",  "Pr2(>|t|)")
  rownames(result) = colnames(X)
  
  return(result)
}