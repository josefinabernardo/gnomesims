# Ignore warnings for variables without global binding
utils::globalVariables(c("expMeandz", "Av", "Amz", "Cv", "Ev", "E", "Adz", "Fdz",
                         "Fmz", "G", "Pdz", "Ydz", "Pmz", "Ymz", "Smz1",
                         "GDPmz", "Sdz1", "GDPdz", "b0", "pred", "bs1",
                         "bs2", "SD", "Rdz", "B0", "g", "b"))

# Define defaults outside of function
default_a <- sqrt(c(.4))
default_c <- sqrt(c(.3))
default_e <- sqrt(c(.3))
default_ct <- sqrt(c(0, .0025, .01))
default_si <- sqrt(c(0, .0025, .01))
default_x <- 0
default_assortm <- 0

#' OpenMx Simulation Function
#'
#' Function to run OpenMx models estimating cultural transmission and sibling interaction on simulated data
#'
#' @import OpenMx
#' @importFrom MASS mvrnorm
#' @importFrom dplyr mutate distinct
#' @importFrom stats D cor var
#' @importFrom magrittr %>%
#'
#' @param alpha Alpha used to calculate power
#' @param seed Set a seed if desired
#' @param standPGS Boolean to standardize the PGS
#' @param nmz Sample size monozygotic twins
#' @param ndz Sample size dizygotic twins
#' @param a Additive genetic path coefficient
#' @param c Shared environmental path coefficient
#' @param e Unique environmental path coefficient
#' @param ct Cultural Transmission - Parent genotype to child phenotype
#' @param si Sibling Interaction - Sibling 1 genotype to sibling 2 phenotype
#' @param x Sibling Interaction at the phenotypic level
#' @param nloci Number of diallelic loci
#' @param npgsloci Number of loci comprising the PGS
#' @param assortm Assortative mating - genetic correlation between the parents
#'
#' @return List with data frame of power estimates and data frame of parameter estimates
#' @export
#'
#' @examples
#' gnome_mx_simulation(ct = .01, si = .025, npgsloci = 10)
gnome_mx_simulation <- function(
    alpha = .05, # Alpha for power
    seed = NA, # Set a seed if desired
    standPGS = FALSE, # Standardize the PGS
    nmz = 4000, # Sample size monozygotic twins
    ndz = 4000, # Sample size dizygotic twins
    a = default_a, # Additive genetic path coefficient
    c = default_c, # Shared environmental path coefficient
    e = default_e, # Unique environmental path coefficient
    ct = default_ct, # Cultural Transmission - Parent genotype to child phenotype
    si = default_si, # Sibling Interaction - Sibling 1 genotype to sibling 2 phenotype
    x = default_x, # Sibling interaction at the phenotypic level
    nloci = 100, # Number of diallelic loci
    npgsloci = c(2, 5, 10, 15), # Number of loci comprising the PGS
    assortm = default_assortm # Assortative mating - genetic correlation
){
  gc()
  # Set option OpenMx
  # mxOption(model = NULL, key = "Default optimizer", "NPSOL")
  # mxOption(model = NULL, key = "Default optimizer", "CSOLNP")
  mxOption(model = NULL, key = "Default optimizer", "SLSQP")
  mxOption(NULL, "Number of Threads", parallel::detectCores())

  # Logical for assortative mating
  assortm_logical = all(assortm == 0)

  # Create all possible parameter combinations
  ifelse(assortm_logical,
         param_combinations <- expand.grid(a = a, c = c, e = e, x = x, ct = ct, si = si),
         param_combinations <- expand.grid(a = a, c = c, e = e, x = x, ct = ct, si = si, assortm = assortm))

  # Number of settings we iterate through
  n_set <- nrow(param_combinations)

  # R2 of the polygenic scores
  global_ppgs <- npgsloci/nloci

  # Initiate counters
  counter_within <- 0  # counts sets within PGS setting
  counter_overall <- 0 # counts sets overall

  # Determine number of rows for data frames
  ifelse(assortm_logical,
         n_cols <- 18,
         n_cols <- 19)
  n_rows <- n_set * length(global_ppgs)

  # Pre-allocate data frames with the appropriate dimensions
  final_mx_estimates <- data.frame(matrix(NA, nrow = n_rows, ncol = n_cols))
  final_mx_power <- data.frame(matrix(NA, nrow = n_rows, ncol = n_cols))

  for (ngp_i in seq_along(npgsloci)) {

    ngp <- nloci[ngp_i]
    p_pgs <- global_ppgs[ngp_i]  # percentage of genetic variance explained by pgs

    print(paste('Running simulation proportion of genetic variance explained by the PGS is:', p_pgs, "."))

    p_A <- 1 - p_pgs # not explained = A without pgs effect
    # e.g., if par_as^2 = .4, then this A variance is due to ng genes
    #                         of .4*(npg/ng) is due to the PGS
    #                         Given var(PH) = 1 (assuming no covAC), the PGS explained {.4*(npg/ng)}/1 of the phenotypic variance

    # Create
    ifelse(assortm_logical,
          setkeep <- matrix(NA, n_set, 10),   # to keep settings
          setkeep <- matrix(NA, n_set, 11))
    mxkeep <- matrix(NA, n_set, 16) # openmx results

    # Print number of settings to the user
    print(paste('The factorial design has', n_set, 'setting(s).'))

    for (i in 1:n_set) {
      par_a <- param_combinations$a[i]
      par_c <- param_combinations$c[i]
      par_e <- param_combinations$e[i]
      par_x <- param_combinations$x[i]
      par_g <- param_combinations$ct[i]
      par_b <- param_combinations$si[i]
      if(assortm_logical == FALSE) {
        par_assortm <- param_combinations$assortm[i]
      }

      counter_within <- counter_within + 1 # count sets in factorial design
      counter_overall <- counter_overall + 1 # count sets overall
      #
      print(c(counter_overall))
      #
      ifelse(assortm_logical,
            setkeep[counter_within,1:10] <- c(nmz, ndz, par_a, par_c, par_e, par_g, par_b, par_x, p_pgs, p_A),
            setkeep[counter_within,1:11] <- c(nmz, ndz, par_a, par_c, par_e, par_g, par_b, par_x, p_pgs, p_A, par_assortm))

      #
      VA1=p_A; VP=p_pgs;VC=1; VE=1 # .... VA1+VP = par_as^2
      ##
      # simulate data exactly:
      # the sample MZ and DZ phenotypic covariance = the population matrices
      #                 m    m   v    v   t1    t1    t2    t2
      #                 1   2    3    4   5     6     7      8   9  10   11
      SLmz=SLdz=diag(c(VA1, VP, VA1, VP, VA1/2, VP/2,VA1/2, VP/2,VC, VE, VE))
      #
      SLmz[5,7]=SLmz[7,5]=VA1/2
      SLmz[6,8]=SLmz[8,6]=VP/2

      # Assortative mating
      if(assortm_logical == FALSE) {
        SLmz[1,3]=SLmz[3,1]=SLdz[1,3]=SLdz[3,1]=par_assortm*sqrt(SLmz[1,1]*SLmz[3,3])
        SLmz[2,4]=SLmz[4,2]=SLdz[2,4]=SLdz[4,2]=par_assortm*sqrt(SLmz[2,2]*SLmz[4,4])
      }

      #
      # simulate the latent variables exactly
      #
      Ldz=MASS::mvrnorm(ndz, rep(0,11), Sigma=SLdz, emp=TRUE)  # emp=TRUE means exact data simulation, cov(Ldz) = SLdz
      Lmz=MASS::mvrnorm(nmz, rep(0,11), Sigma=SLmz, emp=TRUE)  # emp=TRUE means exact data simulation, cov(Lmz) = SLmz  #
      # build the exact simulated data DZ
      Am_=Ldz[,1]
      Pm=Ldz[,2] # PGS
      Af_=Ldz[,3]
      Pf=Ldz[,4] # PGS
      At1r=Ldz[,5]
      Pt1r=Ldz[,6]
      At2r=Ldz[,7]
      Pt2r=Ldz[,8]
      #
      At1_=.5*Am_+.5*Af_+At1r # A res t1
      Pt1=.5*Pm+.5*Pf+Pt1r # PGS t1
      At2_=.5*Am_+.5*Af_+At2r #t2
      Pt2=.5*Pm+.5*Pf+Pt2r# t2
      C=Ldz[,9]
      E1=Ldz[,10]
      E2=Ldz[,11]
      #
      Am=Am_+Pm
      Af=Af_+Pf
      At1=At1_+Pt1
      At2=At2_+Pt2
      #
      Pht1=par_a*At1 + par_e*E1 + par_c*C +  par_g*(Am+Af) + par_b*(At2)
      Pht2=par_a*At2 + par_e*E2 + par_c*C +  par_g*(Am+Af) + par_b*(At1)
      Pht1=Pht1 + par_x*Pht2
      Pht2=Pht2 + par_x*Pht1
      # standardize PGS
      if (standPGS) {
        Pm=scale(Pm)
        Pf=scale(Pf)
        Pt1=scale(Pt1)
        Pt2=scale(Pt2)
      }
      #
      exdatdz=as.data.frame(cbind(Pm, Pf, Pt1, Pt2, Pht1, Pht2))  # exact data sim
      exdatdzA=as.data.frame(cbind(Am, Af, At1, At2, Pht1, Pht2))
      #
      #
      # build the exact simulated data MZ
      Am_=Lmz[,1]
      Pm=Lmz[,2]
      Af_=Lmz[,3]
      Pf=Lmz[,4]
      At1r=Lmz[,5]
      Pt1r=Lmz[,6]
      At2r=Lmz[,5] # MZ MZ MZ MZ MZ Lmz[,7]
      Pt2r=Lmz[,6] # MZ MZ MZ MZ MZ Lmz[,8]
      #
      At1_=.5*Am_+.5*Af_+At1r # A res
      Pt1=.5*Pm+.5*Pf+Pt1r # PGS
      At2_=At1_ #.5*Am_+.5*Af_+At2r
      Pt2 = Pt1 # MZMZMZMZMZ .5*Pm+.5*Pf+Pt2r #
      C=Lmz[,9]
      E1=Lmz[,10]
      E2=Lmz[,11]
      #
      Am=Am_+Pm
      Af=Af_+Pf
      At1=At1_+Pt1
      At2=At1 # At2_+Pt2
      #
      Pht1=par_a*At1 + par_e*E1 + par_c*C +  par_g*(Am+Af) + par_b*(At2)
      Pht2=par_a*At2 + par_e*E2 + par_c*C +  par_g*(Am+Af) + par_b*(At1)
      Pht1=Pht1 + par_x*Pht2
      Pht2=Pht2 + par_x*Pht1
      #
      ## standardize PGS
      #
      if (standPGS) {
        Pm=scale(Pm)
        Pf=scale(Pf)
        Pt1=scale(Pt1)
        Pt2=scale(Pt2)
      }
      #
      exdatmz=as.data.frame(cbind(Pm, Pf, Pt1, Pt2, Pht1, Pht2))  #
      exdatmzA=as.data.frame(cbind(Am, Af, At1, At2, Pht1, Pht2))
      #
      colnames(exdatdz) =colnames(exdatmz) = c('pgsm','pgsf','pgst1','pgst2','pht1','pht2')
      # keep following
      round(cor(exdatmz),4) -> phrmz_
      round(cor(exdatdz),4)      -> phrdz_
      round(apply(exdatmz,2,var),4)[6] -> phvarmz_
      round(apply(exdatdz,2,var),4)[6] -> phvardz_
      #
      # --------------------------------- end exact data simualtion
      # exdatdz and exdatmz are dataframes, exact data simulation
      # when analysing these data, we should "get out" what we "put in"
      #
      # exact simulated data in the wide (horizontal) format
      #
      phdatmz_e = as.data.frame(exdatmz)
      phdatdz_e = as.data.frame(exdatdz)
      colnames(phdatmz_e)=colnames(phdatdz_e) =c('pgsm','pgsf','pgst1','pgst2','pht1','pht2')
      #
      # add sum of mother and father PGS and add mean of twins PGS ... we need these additional variables to
      #                                                                detect cov(AC)
      #
      addmz_e=cbind(phdatmz_e$pgsm+phdatmz_e$pgsf, (phdatmz_e$pgst1+phdatmz_e$pgst2)/2)
      colnames(addmz_e) = c('pgsmf','mpgst')
      phdatmz_e = cbind(phdatmz_e, addmz_e)
      adddz_e=cbind(phdatdz_e$pgsm+phdatdz_e$pgsf, (phdatdz_e$pgst1+phdatdz_e$pgst2)/2)
      colnames(adddz_e) = c('pgsmf','mpgst')
      phdatdz_e = cbind(phdatdz_e, adddz_e)
      # add sum and mean
      c(1,2,3,4,5,7,9,10) -> i1
      c(5,6,1,2,3,4,7,8) -> i2
      #
      apply(phdatmz_e,2,var)[i2]
      round(cor(phdatmz_e),3)[i2,i2]
      #
      apply(phdatdz_e,2,var)[i2]
      round(cor(phdatdz_e),3)[i2,i2]
      #
      #
      #
      # [1] "pht1"   "pht2"   "pgsm"   "pgsf"   "pgst1"  "pgsnt1" "pgst2"  "pgsnt2" "pgsmf"  "mpgst"
      #
      # Organize data in long format simulated data long format or vertical format
      #
      #
      # Long format exact simulated data
      phdatmzL_e = matrix(1,nmz*2,8)
      phdatmzL_e[,2]=nmz + c(c(1:nmz),c(1:nmz))
      phdatmzL_e[,3]=c(phdatmz_e$pht1,phdatmz_e$pht2)
      phdatmzL_e[,4]=c(phdatmz_e$pgsm,phdatmz_e$pgsm)
      phdatmzL_e[,5]=c(phdatmz_e$pgsf,phdatmz_e$pgsf)
      phdatmzL_e[,6]=c(phdatmz_e$pgst1,phdatmz_e$pgst2)
      phdatmzL_e[,7]=c(phdatmz_e$pgsmf,phdatmz_e$pgsmf)
      phdatmzL_e[,8]=c(phdatmz_e$mpgst,phdatmz_e$mpgst)
      #
      ix_ = sort.int(phdatmzL_e[,2], index.return=TRUE)
      phdatmzL_e = phdatmzL_e[ix_$ix,]
      colnames(phdatmzL_e)=c("zyg","famnr","ph","pgsm","pgsf","pgst","pgsmf","mpgst")
      phdatmzL_e=as.data.frame(phdatmzL_e)
      #
      #
      # Long format exact simulated data
      #
      phdatdzL_e = matrix(1,ndz*2,8)
      phdatdzL_e[,2]=c(c(1:ndz),c(1:ndz))
      phdatdzL_e[,3]=c(phdatdz_e$pht1,phdatdz_e$pht2)
      phdatdzL_e[,4]=c(phdatdz_e$pgsm,phdatdz_e$pgsm)
      phdatdzL_e[,5]=c(phdatdz_e$pgsf,phdatdz_e$pgsf)
      phdatdzL_e[,6]=c(phdatdz_e$pgst1,phdatdz_e$pgst2)
      phdatdzL_e[,7]=c(phdatdz_e$pgsmf,phdatdz_e$pgsmf)
      phdatdzL_e[,8]=c(phdatdz_e$mpgst,phdatdz_e$mpgst)
      #
      ix_ = sort.int(phdatdzL_e[,2], index.return=TRUE)
      phdatdzL_e = phdatdzL_e[ix_$ix,]
      colnames(phdatdzL_e)=c("zyg","famnr","ph","pgsm","pgsf","pgst","pgsmf","mpgst")
      phdatdzL_e=as.data.frame(phdatdzL_e)
      #
      phdatL_e=rbind(phdatmzL_e, phdatdzL_e)

      # OpenMx covariance structure modeling
      vnamesdz=c('pht1','pht2','pgst1','pgst2','pgsm','pgsf')
      vnamesmz=c('pht1','pht2','pgst1','pgsm','pgsf')
      #
      ##         1       2    3       4      5      6      7       8  .... 7 and 8 are correlated 1 in MZs
      w_mzdat=phdatmz_e[,vnamesmz]
      w_dzdat=phdatdz_e[,vnamesdz]
      # ==============================================================================================
      #
      # model covariance matrix:
      #  GDFD'G' + Y |  GDFD'
      #  FD'G'D      |  DFD'
      #
      # DFD' is the covariance matrix 2x2 of the parental PGS
      # D 2x2 diagonal contains the stds of the parental PGS
      # F is the 2x2 correlation matrix of the parental PGS.
      #
      # G contains the regressions of offspring phenotype on parental and offspring PGS
      # Y is the residual phenotypic covariance matrix modeled using ACE model
      #
      mz1out=4 # remove 8th variable PGS of mz twin 2.
      Filter=diag(6)[-mz1out,]
      #
      RAdz=matrix(c(.5),4,4)
      diag(RAdz)=1
      RAdz[3,4] <- RAdz[4,3] <- 0# m f assortative mating
      RAmz=RAdz
      RAmz[1,2]=RAmz[2,1]=1 # MZ twins
      RAfree=matrix(FALSE,4,4)
      RAlabels=matrix(NA,4,4)
      #
      M1 <- OpenMx::mxModel("M1",
                    mxMatrix(type='Full', nrow=5, ncol=6, free=FALSE, values=Filter, labels=c(NA), name='Filter'),
                    # mean
                    mxMatrix(type="Full", nrow=1, ncol=6,
                             #
                             free=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
                             values=c(0,0,0,0,0,0),
                             labels=c("mph","mph","mpgs","mpgs","mpgs","mpgs"),
                             name="expMeandz"),
                    mxAlgebra(expression=(expMeandz%*%t(Filter)), name="expMeanmz"),
                    # pgs stdevs matrix D
                    mxMatrix(type="Diag", nrow=4, ncol=4,
                             free=c(TRUE,TRUE,TRUE,TRUE),
                             values=c(.7,.7,.7,.7),
                             labels=c("sdp","sdp","sdp","sdp"), lbound=.01,
                             name="D"),
                    # pgs correlation matris F
                    mxMatrix(type="Symm", nrow=4, ncol=4,
                             free=RAfree, values=RAdz, labels=RAlabels, name="Fdz"),
                    mxMatrix(type="Symm", nrow=4, ncol=4,
                             free=RAfree, values=RAmz, labels=RAlabels, name="Fmz"),
                    # A matrix for P = A + C + E
                    mxMatrix(type="Symm", nrow=2, ncol=2,
                             free=c(FALSE), values=RAdz[1:2,1:2], labels=c(NA), name="Adz"),
                    mxMatrix(type="Symm", nrow=2, ncol=2,
                             free=c(FALSE), values=RAmz[1:2,1:2], labels=c(NA), name="Amz"),
                    # C matrix for P = A + C + E
                    mxMatrix(type="Symm", nrow=2, ncol=2,
                             free=c(FALSE), values=c(1), labels=c(NA), name="C"),
                    # E
                    mxMatrix(type="Iden", nrow=2, ncol=2, name="E"),
                    #
                    # A matrix for P = A + C + E
                    mxMatrix(type="Symm", nrow=1, ncol=1,
                             free=c(TRUE), values=c(.33), labels=c("sa2"), name="Av"),
                    # C matrix for P = A + C + E
                    mxMatrix(type="Symm", nrow=1, ncol=1,
                             free=c(TRUE), values=c(.33), labels=c("sc2"), name="Cv"),
                    # E matrix for P = A + C + E
                    mxMatrix(type="Symm", nrow=1, ncol=1,
                             free=c(TRUE), values=c(.33), labels=c("se2"), name="Ev"),
                    #
                    mxMatrix(type="Full", nrow=2, ncol=4,
                             free=matrix(c(
                               TRUE,TRUE,TRUE,TRUE,
                               TRUE,TRUE,TRUE,TRUE),2,4,byrow=TRUE),
                             labels=matrix(c(
                               'a1','b1','g1','g1',
                               'b1','a1','g1','g1'),2,4,byrow=TRUE),
                             values=matrix(c(
                               .1, 0, .01,.01,
                               0,.1, .01,.01),2,4,byrow=TRUE), name='G'),
                    #
                    mxAlgebra(expression=Av%x%Amz + Cv%x%C + Ev%x%E, name="Ymz"),
                    mxAlgebra(expression=Av%x%Adz + Cv%x%C + Ev%x%E, name="Ydz"),
                    mxAlgebra(expression=D%*%Fdz%*%t(D), name="Pdz"),
                    mxAlgebra(expression=D%*%Fmz%*%t(D), name="Pmz"),
                    mxAlgebra(expression=G%*%Pdz%*%t(G)+Ydz, name="Sdz1"),
                    mxAlgebra(expression=G%*%Pmz%*%t(G)+Ymz, name="Smz1"),
                    #
                    #   mxAlgebra(expression=G%*%D%*%Pmz, name="GDPmz"),
                    #   mxAlgebra(expression=G%*%D%*%Pdz, name="GDPdz"),
                    mxAlgebra(expression=G%*%Pmz, name="GDPmz"),
                    mxAlgebra(expression=G%*%Pdz, name="GDPdz"),
                    mxAlgebra(expression=Filter%*%rbind(cbind(Smz1,GDPmz),cbind(t(GDPmz),Pmz))%*%t(Filter), name="Smz"),
                    mxAlgebra(expression=rbind(cbind(Sdz1,GDPdz),cbind(t(GDPdz),Pdz)), name="Sdz"),
                    #
                    mxCI(c('b1','a1','g1'))
                    #
      )

      DZ <-  OpenMx::mxModel('DZ',
                     mxData( observed=w_dzdat, type="raw"),
                     mxExpectationNormal( covariance="M1.Sdz", means="M1.expMeandz",
                                          dimnames=vnamesdz),  # the fit function
                     mxFitFunctionML()
      )
      MZ <-  OpenMx::mxModel('MZ',
                     mxData( observed=w_mzdat, type="raw"),
                     mxExpectationNormal( covariance="M1.Smz", means="M1.expMeanmz",
                                          dimnames=vnamesmz),  # the fit function
                     mxFitFunctionML()
      )

      #
      Model_1 <-  OpenMx::mxModel("MZDZModel", M1, MZ, DZ,
                          mxFitFunctionMultigroup( c("MZ","DZ"))
      )
      #
      #
      # fit the model
      Model_1out <- mxRun(Model_1, intervals=FALSE, silent = TRUE)
      #
      #summary(Model_1out)
      sat_1out <- mxRefModels(Model_1out, run=TRUE)
      mxCompare(sat_1out, Model_1out)
      #
      #
      Model_1b <- omxSetParameters(Model_1out, labels=c('b1'), free=FALSE, values=c(0))
      Model_1b_out <- mxRun(Model_1b, intervals=TRUE, silent = TRUE)
      mxCompare(Model_1out, Model_1b_out)

      Model_1g <- omxSetParameters(Model_1out, labels=c('g1'), free=FALSE, values=c(0))
      Model_1g_out <- mxRun(Model_1g, intervals=TRUE, silent = TRUE)
      mxCompare(Model_1out, Model_1g_out)

      Model_1bg <- omxSetParameters(Model_1out, labels=c('g1','b1'), free=FALSE, values=c(0))
      Model_1bg_out <- mxRun(Model_1bg, intervals=TRUE, silent = TRUE)
      mxCompare(Model_1out, Model_1bg_out)
      #
      mxRefModels(Model_1out, run=TRUE) -> sat_1out
      mxCompare(sat_1out, Model_1out)
      #
      #
      # -----------------------------------------------   vars
      varnames=c('pht1','pht2')#
      #
      # the model to calculate expected summary statistics
      # this is the twin model phenotypic
      # [1] "pht1"   "pht2"   "pgsm"   "pgsf"   "pgst1"  "pgsnt1" "pgst2"  "pgsnt2" "pgsmf"  "mpgst"
      # [1] "pgsm"  "pgsf"  "pgst1" "pgst2" "pht1"  "pht2"  "pgsmf" "mpgst"

      # a model the data, the fit function (MZ)
      MZmodel <- OpenMx::mxModel("MZ",
                         #
                         # Matrix expMean for expected mean vector for MZ and DZ twins
                         mxMatrix(type="Full", nrow=1, ncol=4, free=FALSE, labels=c("data.pgst1","data.pgst2","data.pgsm","data.pgsf"), name="pred"),
                         mxMatrix(type="Full", nrow=1, ncol=4, free=c(TRUE,TRUE,TRUE,TRUE), values=c(0,0,0,0),
                                  labels=c("bpgst","bpgsb","bpgsg","bpgsg"), name="bs1"),
                         mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(0,0,0,0),
                                  labels=c("bpgsb","bpgst","bpgsg","bpgsg"), name="bs2"),
                         mxMatrix(type="Full", nrow=1, ncol=1,
                                  free=c(TRUE),values=c(0),labels=c("b0"),
                                  name="Int"),
                         mxAlgebra(expression=cbind(b0+pred%*%t(bs1), b0+pred%*%t(bs2)), name='expMean'),
                         mxData(observed=phdatmz_e, type="raw"),
                         mxExpectationNormal(covariance="ACE.expCovMZ",
                                             means = "expMean", varnames),
                         mxFitFunctionML()
      )
      # a model the data, the fit function (DZ)
      DZmodel <- OpenMx::mxModel("DZ",
                         #
                         # Matrix expMean for expected mean vector for MZ and DZ twins
                         mxMatrix(type="Full", nrow=1, ncol=4, free=FALSE, labels=c("data.pgst1","data.pgst2","data.pgsm","data.pgsf"), name="pred"),
                         mxMatrix(type="Full", nrow=1, ncol=4, free=c(TRUE,TRUE,TRUE,TRUE), values=c(0,0,0,0),
                                  labels=c("bpgst","bpgsb","bpgsg","bpgsg"), name="bs1"),
                         mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(0,0,0,0),
                                  labels=c("bpgsb","bpgst","bpgsg","bpgsg"), name="bs2"),
                         mxMatrix(type="Full", nrow=1, ncol=1,
                                  free=c(TRUE),values=c(0),labels=c("b0"),
                                  name="Int"),
                         mxAlgebra(expression=cbind(b0+pred%*%t(bs1), b0+pred%*%t(bs2)), name='expMean'),
                         mxData(observed=phdatdz_e, type="raw"),
                         mxExpectationNormal(covariance="ACE.expCovDZ",
                                             means = "expMean", varnames),
                         mxFitFunctionML()
      )

      varnames=c('pht1','pht2')#
      #
      # the model to calculate expected summary statistics
      # this is the twin model phenotypic
      # [1] "pht1"   "pht2"   "pgsm"   "pgsf"   "pgst1"  "pgsnt1" "pgst2"  "pgsnt2" "pgsmf"  "mpgst"
      # [1] "pgsm"  "pgsf"  "pgst1" "pgst2" "pht1"  "pht2"  "pgsmf" "mpgst"
      # a model the data, the fit function (MZ)
      MZmodel <- OpenMx::mxModel("MZ",
                         #
                         # Matrix expMean for expected mean vector for MZ and DZ twins
                         mxMatrix(type="Full", nrow=1, ncol=4, free=FALSE, labels=c("data.pgst1","data.pgst2","data.pgsm","data.pgsf"), name="pred"),
                         mxMatrix(type="Full", nrow=1, ncol=4, free=c(TRUE,TRUE,TRUE,TRUE), values=c(0,0,0,0),
                                  labels=c("bpgst","bpgsb","bpgsg","bpgsg"), name="bs1"),
                         mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(0,0,0,0),
                                  labels=c("bpgsb","bpgst","bpgsg","bpgsg"), name="bs2"),
                         mxMatrix(type="Full", nrow=1, ncol=1,
                                  free=c(TRUE),values=c(0),labels=c("b0"),
                                  name="Int"),
                         mxAlgebra(expression=cbind(b0+pred%*%t(bs1), b0+pred%*%t(bs2)), name='expMean'),
                         mxData(observed=phdatmz_e, type="raw"),
                         mxExpectationNormal(covariance="SAT.expCovMZ",
                                             means = "expMean", varnames),
                         mxFitFunctionML()
      )
      # a model the data, the fit function (DZ)
      DZmodel=OpenMx::mxModel("DZ",
                      #
                      # Matrix expMean for expected mean vector for MZ and DZ twins
                      mxMatrix(type="Full", nrow=1, ncol=4, free=FALSE, labels=c("data.pgst1","data.pgst2","data.pgsm","data.pgsf"), name="pred"),
                      mxMatrix(type="Full", nrow=1, ncol=4, free=c(TRUE,TRUE,TRUE,TRUE), values=c(0,0,0,0),
                               labels=c("bpgst","bpgsb","bpgsg","bpgsg"), name="bs1"),
                      mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(0,0,0,0),
                               labels=c("bpgsb","bpgst","bpgsg","bpgsg"), name="bs2"),
                      mxMatrix(type="Full", nrow=1, ncol=1,
                               free=c(TRUE),values=c(0),labels=c("b0"),
                               name="Int"),
                      mxAlgebra(expression=cbind(b0+pred%*%t(bs1), b0+pred%*%t(bs2)), name='expMean'),
                      mxData(observed=phdatdz_e, type="raw"),
                      mxExpectationNormal(covariance="SAT.expCovDZ",
                                          means = "expMean", varnames),
                      mxFitFunctionML()
      )

      varnames <- c('pht1','pht2')#
      #
      # the model to calculate expected summary statistics
      # this is the twin model phenotypic
      # [1] "pht1"   "pht2"   "pgsm"   "pgsf"   "pgst1"  "pgsnt1" "pgst2"  "pgsnt2" "pgsmf"  "mpgst"
      # [1] "pgsm"  "pgsf"  "pgst1" "pgst2" "pht1"  "pht2"  "pgsmf" "mpgst"
      nphen1 <- 1
      nphen2 <- 2
      DZModel  <-  OpenMx::mxModel("DZonly",
                           #
                           # Matrices a, c, and e to store the a, c, and e path coefficients
                           mxMatrix(type="Stand", nrow=nphen2, ncol=nphen2,
                                    free=c(TRUE), values=c(.25),
                                    labels=c("rdz"),name="Rdz"),
                           mxMatrix(type="Diag", nrow=nphen2, ncol=nphen2,
                                    free=c(TRUE), values=c(.7),
                                    labels=c("sd","sd"),name="SD"),
                           #
                           #
                           # Matrix expCovMZ for expected covariance matrix for DZ twins
                           #
                           mxAlgebra( expression=
                                        SD%*%Rdz%*%SD,,name="expCovDZ"),
                           #
                           # Matrix expMean for expected mean vector for DZ twins
                           #
                           mxMatrix(type="Full", nrow=1, ncol=4, free=FALSE, labels=c("data.pgst1","data.pgst2","data.pgsm","data.pgsf"), name="pred"),
                           mxMatrix(type="Full", nrow=1, ncol=4, free=c(TRUE,TRUE,TRUE,TRUE), values=c(0,0,0,0),
                                    labels=c("bpgst","bpgsb","bpgsg","bpgsg"), name="bs1"),
                           mxMatrix(type="Full", nrow=1, ncol=4, free=TRUE, values=c(0,0,0,0),
                                    labels=c("bpgsb","bpgst","bpgsg","bpgsg"), name="bs2"),
                           mxMatrix(type="Full", nrow=1, ncol=1,
                                    free=c(TRUE),values=c(0),labels=c("b0"),
                                    name="B0"),
                           mxAlgebra(expression=cbind(B0+pred%*%t(bs1), B0+pred%*%t(bs2)), name='expMean'),
                           mxData(observed=phdatdz_e, type="raw"),
                           mxExpectationNormal(covariance="expCovDZ",
                                               means ="expMean", varnames),
                           mxFitFunctionML()
      )
      # Model_4 <-  OpenMx::mxModel(name="DZ1SAT", DZModel)
      Model_4 <-  OpenMx::mxModel(DZModel)
      # fit the model
      Model_4out <- mxRun(Model_4, silent = TRUE)

      Model_4g <- omxSetParameters(Model_4out, labels='bpgsg', values=0, free=FALSE)
      Model_4g_out = mxRun(Model_4g, silent = TRUE)
      #
      Model_4b=omxSetParameters(Model_4out, labels='bpgsb', values=0, free=FALSE)
      Model_4b_out <- mxRun(Model_4b, silent = TRUE)
      # "bpgsb","bpgsg"
      Model_4bg <- omxSetParameters(Model_4out, labels=c('bpgsb','bpgsg'), free=FALSE, values=c(0))
      Model_4bg_out <- mxRun(Model_4bg, silent = TRUE)

      # Power
      ncp_tmp <- c(
        mxCompare(Model_1b_out, Model_1bg_out)[2,7],
        mxCompare(Model_1g_out, Model_1bg_out)[2,7],
        mxCompare(Model_1out, Model_1g_out)[2,7],
        mxCompare(Model_1out, Model_1b_out)[2,7],

        mxCompare(Model_4b_out, Model_4bg_out)[2,7],
        mxCompare(Model_4g_out, Model_4bg_out)[2,7],
        mxCompare(Model_4out, Model_4g_out)[2,7],
        mxCompare(Model_4out, Model_4b_out)[2,7]
      )

      ncp_power_tmp <- sapply(ncp_tmp, function(ncp) {
        gnome_power(alpha, 1, ncp)
      })

      # Estimates
      estimates_tmp <- c(
        summary(Model_1b_out)$parameters[8, "Estimate"], # g1 - CT Model 1 CT only
        summary(Model_1g_out)$parameters[8,'Estimate'], # b1 - SI Model 1 SI only
        summary(Model_1out)$parameters[9,'Estimate'], # g1 - CT Model 1 both
        summary(Model_1out)$parameters[8,'Estimate'], # b1 - SI Model 1 both

        summary(Model_4b_out)$parameters[4, 'Estimate'], # bpgsg - CT Model 4 CT only
        summary(Model_4g_out)$parameters[4, 'Estimate'], # bpgsb - SI Model 4 SI only
        summary(Model_4out)$parameters[5, 'Estimate'], # bpgsg - CT Model 4
        summary(Model_4out)$parameters[4, 'Estimate'] # bpgsb - SI Model 4
      )

      mxkeep[counter_within, ] <- c(ncp_power_tmp, estimates_tmp)

    }

    # For mxkeep
    jpow <- 1:8
    jest <- 9:16

    row_index <- counter_overall - n_set + 1

    if (assortm_logical == TRUE) {
      # Execute these lines when assortm_logical is TRUE
      final_mx_estimates[row_index : counter_overall, 1:10] <- setkeep[, 1:10]
      final_mx_power[row_index : counter_overall, 1:10] <- setkeep[, 1:10]
      final_mx_estimates[row_index : counter_overall, 11:18] <- round(mxkeep[, jest], 3)
      final_mx_power[row_index : counter_overall, 11:18] <- round(mxkeep[, jpow], 3)
    } else {
      # Execute these lines when assortm_logical is FALSE
      final_mx_estimates[row_index : counter_overall, 1:11] <- setkeep[, 1:11]
      final_mx_power[row_index : counter_overall, 1:11] <- setkeep[, 1:11]
      final_mx_estimates[row_index : counter_overall, 12:19] <- round(mxkeep[, jest], 3)
      final_mx_power[row_index : counter_overall, 12:19] <- round(mxkeep[, jpow], 3)
    }

    counter_within = 0 # reset set counter for each PGS setting
  }
  # Re-name columns
  ifelse(assortm_logical,
         setnames <- c('nmz','ndz','a','c','e','g','b','x','PGS','A'),
         setnames <- c('nmz','ndz','a','c','e','g','b','x','PGS','A','assortm'))
  colnames(final_mx_estimates) <- c(setnames, paste0("e", 1:8))
  colnames(final_mx_power) <- c(setnames, paste0("p", 1:8))

  # Use effect size function on the data sets
  final_mx_estimates <- final_mx_estimates %>%
    mutate(Smz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$mz,
           Sdz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$dz)
  final_mx_power <- final_mx_power %>%
    mutate(Smz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$mz,
           Sdz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$dz)

  return(list(power = final_mx_power, params = final_mx_estimates))
}
