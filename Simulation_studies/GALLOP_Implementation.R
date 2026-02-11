
library(lme4)
library(glmmTMB)
library(optimx)

# main function
gallop = function(nullmodel = NULL,     # fitted null model from lme4.
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
  
  # Extract variance components and compute penalty matrix (P)
  if(inherits(null_model, "merMod"))
  {
    varcor = VarCorr(null_model)
    cmd = paste0("varcor$", IDcolname); G = eval(parse(text = cmd))
    sig = attr(varcor, "sc")
    P = solve(G / sig ^ 2)
  }else if(inherits(null_model, "glmmTMB"))
  {
    varcor = VarCorr(null_model)$cond
    cmd = paste0("varcor$", IDcolname); G = eval(parse(text = cmd))
    sig = attr(varcor, "sc")
    P = solve(G / sig ^ 2)
  }else
  {
    stop("Currently we only support fitted models fitted by 'LME4' and 'glmmTMB'.")
  }
  
  # Put data in convenient arrays and vector
  cmd = paste0("unique(Pheno_data$", IDcolname, ")"); SubjID = eval(parse(text = cmd)); SubjID = as.character(SubjID)
  cmd = paste0("table(Pheno_data$", IDcolname, ")"); k = eval(parse(text = cmd))
  n = length(SubjID)
  k = k[SubjID]
  XY = bind_cols(intercept = 1, Pheno_data %>% select(all_of(Timecolname)), 
             Pheno_data %>% select(all_of(Covacolname)), Pheno_data %>% select(all_of(Phenocolname)))
  XY = as.matrix(XY)
  TT = XY[, 1:2]
  nxy = ncol(XY)
  X = XY[, -nxy]
  nx = ncol(X)
  y = XY[, nxy]
  ncov = nx - ncol(TT)
  
  # Compute components of block-diagonal system with covariates  (without SNP)
  # Additionally compute and store object SS = Rot %*% Si, which is used later on
  A21 = matrix(NA, 2*n, 2 + ncov)
  q2 = rep(NA, 2*n)
  SS = matrix(NA, 2*n, 2)
  
  uk = 0
  for(i in 1:n)
  {
    if(i %% 1000 == 0)
      cat(i, "\n")

    ki = k[i]
    uk = max(uk) + 1:ki
    u2 = (i - 1)*2 + 1:2
    Ti = TT[uk,]
    if(ki == 1) Ti = t(as.matrix(Ti))
    Si = crossprod(Ti, Ti)
    sv = svd(Si + P)
    Rot = sqrt(1/sv$d) * sv$u
    Q = Rot %*% t(Ti)
    SS[u2,] = Rot %*% Si
    A21[u2,] = Q %*% X[uk, ]
    q2[u2] = Q %*% y[uk]
  }
  
  q1 = crossprod(X, y)
  A11 = crossprod(X)
  Q = A11 - crossprod(A21)
  q = q1 - crossprod(A21, q2)
  sol = solve(Q, q)
  blups = q2 - A21 %*% sol
  
  # Compute sums of products per subject involved in the crossprod(X, G), crossprod(G) and crossprod(G, y).
  ex = matrix(1, 1, ncov + 2)
  et = matrix(1, 1, 2)
  XTk = kronecker(et, X) * kronecker(TT, ex)
  TTk = kronecker(et, TT) * kronecker(TT, et)
  Tyk = y * TT
  XTs = matrix(0, n, ncol(XTk))
  TTs = matrix(0, n, ncol(TTk))
  Tys = matrix(0, n, 2)
  AtS = matrix(0, n, 2 * nx)
  
  uk = 0
  for (i in 1:n) 
  {
    if(i %% 1000 == 0)
      cat(i, "\n")
    
    ki = k[i]
    uk = max(uk) + 1:ki
    if(ki == 1)
    {
      XTs[i, ] = XTk[uk, ]
      TTs[i, ] = TTk[uk, ]
      Tys[i, ] = Tyk[uk, ]
    }else
    {
      XTs[i, ] = apply(XTk[uk, ], 2, sum)
      TTs[i, ] = apply(TTk[uk, ], 2, sum)
      Tys[i, ] = apply(Tyk[uk, ], 2, sum)
    }
    u2 = (i - 1) * 2 + (1 : 2)
    AtS[i, ] = c(crossprod(A21[u2, ], SS[u2, ]))
  }
  
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
        MAF = colMeans(GenoMatInfo$GenoMat)/2
        
        Theta = D = matrix(NA, nrow(GenoMatInfo$markerInfo), 2)
        for (i in 1:nrow(D)) 
        {
          if(i %% 1000 == 0) cat("Process", i, "variants in this chunk.\n")
          
          si = GenoMatInfo$GenoMat[, i]
          snp2 = rep(si, each = 2)
          H1 = matrix(crossprod(si, XTs), nx, 2)
          H2 = snp2 * SS
          AtH = matrix(crossprod(si, AtS), nx, 2)
          R = H1 - AtH
          Cfix = solve(Q, R)
          Cran = H2 - A21 %*% Cfix
          GtG = matrix(crossprod(si ^ 2, TTs), 2, 2)
          Gty = matrix(crossprod(si, Tys), 2, 1)
          V = GtG - crossprod(H1, Cfix) - crossprod(H2, Cran)
          v = Gty - crossprod(H1, sol) - crossprod(H2, blups)
          Theta[i, ] = solve(V, v)
          D[i, ] = diag(solve(V))
        }
        
        SE = sig * sqrt(D)
        Pval = 2 * pnorm(-abs(Theta/SE))
        
        Output = cbind(as.matrix(GenoMatInfo$markerInfo), MAF, Theta, Pval)
        colnames(Output) = c("Chr", "Pos", "SNPID", "Ref", "Alt", "MAF", "CSeffect", "LTeffect", "CSpval", "LTpval")
        
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
          MAF = colMeans(GenoMatInfo$GenoMat)/2
          
          Theta = D = matrix(NA, nrow(GenoMatInfo$markerInfo), 2)
          for (i in 1:nrow(D)) 
          {
            if(i %% 1000 == 0) cat("Process", i, "variants in this chunk.\n")
            
            si = GenoMatInfo$GenoMat[, i]
            snp2 = rep(si, each = 2)
            H1 = matrix(crossprod(si, XTs), nx, 2)
            H2 = snp2 * SS
            AtH = matrix(crossprod(si, AtS), nx, 2)
            R = H1 - AtH
            Cfix = solve(Q, R)
            Cran = H2 - A21 %*% Cfix
            GtG = matrix(crossprod(si ^ 2, TTs), 2, 2)
            Gty = matrix(crossprod(si, Tys), 2, 1)
            V = GtG - crossprod(H1, Cfix) - crossprod(H2, Cran)
            v = Gty - crossprod(H1, sol) - crossprod(H2, blups)
            Theta[i, ] = solve(V, v)
            D[i, ] = diag(solve(V))
          }
          
          SE = sig * sqrt(D)
          Pval = 2 * pnorm(-abs(Theta/SE))
          
          Output = cbind(as.matrix(GenoMatInfo$markerInfo), MAF, Theta, Pval)
          colnames(Output) = c("Chr", "Pos", "SNPID", "Ref", "Alt", "MAF", "CSeffect", "LTeffect", "CSpval", "LTpval")
          
          results = rbind(results, Output)
          
          file.remove(temp_SNPIDfile); rm(GenoMatInfo); gc()
        }
        
        data.table::fwrite(results, file = Outputfile,
                           row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        
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
  MAF = colMeans(GenoMatInfo$GenoMat)/2
  
  Theta = D = matrix(NA, nrow(GenoMatInfo$markerInfo), 2)
  for (i in 1:nrow(D)) 
  {
    if(i %% 1000 == 0) cat("Process", i, "variants in this chunk.\n")
    
    si = GenoMatInfo$GenoMat[, i]
    snp2 = rep(si, each = 2)
    H1 = matrix(crossprod(si, XTs), nx, 2)
    H2 = snp2 * SS
    AtH = matrix(crossprod(si, AtS), nx, 2)
    R = H1 - AtH
    Cfix = solve(Q, R)
    Cran = H2 - A21 %*% Cfix
    GtG = matrix(crossprod(si ^ 2, TTs), 2, 2)
    Gty = matrix(crossprod(si, Tys), 2, 1)
    V = GtG - crossprod(H1, Cfix) - crossprod(H2, Cran)
    v = Gty - crossprod(H1, sol) - crossprod(H2, blups)
    Theta[i, ] = solve(V, v)
    D[i, ] = diag(solve(V))
  }
  
  SE = sig * sqrt(D)
  Pval = 2 * pnorm(-abs(Theta/SE))
  
  Output = cbind(as.matrix(GenoMatInfo$markerInfo), MAF, Theta, Pval)
  colnames(Output) = c("Chr", "Pos", "SNPID", "Ref", "Alt", "MAF", "CSeffect", "LTeffect", "CSpval", "LTpval")
  
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