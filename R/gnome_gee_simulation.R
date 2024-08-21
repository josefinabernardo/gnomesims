# Ignore warnings for variables without global binding
utils::globalVariables(c("famnr", "g", "b"))

# Define defaults outside of function
default_a <- sqrt(c(.4))
default_c <- sqrt(c(.3))
default_e <- sqrt(c(.3))
default_ct <- sqrt(c(0, .0025, .01))
default_si <- sqrt(c(0, .0025, .01))
default_x <- 0

#' Gee Simulation Function
#'
#' Function to run generalized estimating equations (gee) models estimating cultural transmission and sibling interaction on simulated data
#'
#' @importFrom MASS mvrnorm
#' @importFrom dplyr mutate
#' @importFrom stats D cor var
#' @importFrom magrittr %>%
#' @importFrom stats lm
#' @importFrom geepack geeglm
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
#' @param cmethod Gee error covariance structure
#'
#' @return List with data frame of power estimates and data frame of parameter estimates
#' @export
#'
#' @examples
#' gnome_gee_simulation(ct = .01, si = .025, npgsloci = 10)
gnome_gee_simulation <- function(
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
    cmethod = 'independence' # Gee error covariance structure
){

  print("Updated via copy-paste")

  #Create all possible parameter combinations
  param_combinations <- expand.grid(a = a, c = c, e = e, x = x, ct = ct, si = si)

  # Number of settings we iterate through
  n_set <- nrow(param_combinations)

  # R2 of the polygenic scores
  global_ppgs <- npgsloci/nloci

  # Initiate counters
  counter_within <- 0  # counts sets within PGS setting
  counter_overall <- 0 # counts sets overall

  # Determine number of rows for data frames
  n_rows <- n_set * length(global_ppgs)

  # Pre-allocate data frames with the appropriate dimensions
  final_gee_estimates <- data.frame(matrix(NA, nrow = n_rows, ncol = 18))
  final_gee_power <- data.frame(matrix(NA, nrow = n_rows, ncol = 18))

  for (ngp_i in seq_along(npgsloci)) {

    ngp <- nloci[ngp_i]
    p_pgs <- global_ppgs[ngp_i]  # percentage of genetic variance explained by pgs

    print(paste('Running simulation proportion of genetic variance explained by the PGS is:', p_pgs, "."))

    p_A <- 1-p_pgs # not explained = A without pgs effect
    # e.g., if par_as^2 = .4, then this A variance is due to ng genes
    #                         of .4*(npg/ng) is due to the PGS
    #                         Given var(PH) = 1 (assuming no covAC), the PGS explained {.4*(npg/ng)}/1 of the phenotypic variance

    # Create
    setkeep <- matrix(NA, n_set, 10)   # to keep settings
    geekeep <- matrix(NA, n_set, 16) # gee results

    # Print number of settings to the user
    print(paste('The factorial design has', n_set, 'setting(s).'))

    for (i in 1:n_set) {
      par_a <- param_combinations$a[i]
      par_c <- param_combinations$c[i]
      par_e <- param_combinations$e[i]
      par_x <- param_combinations$x[i]
      par_g <- param_combinations$ct[i]
      par_b <- param_combinations$si[i]

      counter_within <- counter_within + 1 # count sets in factorial design
      counter_overall <- counter_overall + 1 # count sets overall
      #
      print(c(counter_overall))
      #
      setkeep[counter_within,1:10] <- c(nmz, ndz, par_a, par_c, par_e, par_g, par_b, par_x, p_pgs, p_A)
      #
      VA1=p_A; VP=p_pgs;VC=1; VE=1 # .... VA1+VP = par_as^2
      ##
      # simulate data exactly:
      # the sample MZ and DZ phenotypic covariance = the population matrices
      #                 m    m   v    v   t1    t1    t2    t2
      #                 1   2    3    4   5     6     7      8   9  10   11
      SLmz=SLdz=diag(c(VA1, VP, VA1, VP, VA1/2, VP/2,VA1/2, VP/2,VC, VE, VE))
      #
      #SLmz[5,7]=SLmz[7,5]=VA1/2   #
      #SLmz[6,8]=SLmz[8,6]=VP/2    #
      #
      # simulate the latent variables exactly
      #
      Ldz=mvrnorm(ndz, rep(0,11), Sigma=SLdz, emp=T)  # emp=T means exact data simulation, cov(Ldz) = SLdz
      Lmz=mvrnorm(nmz, rep(0,11), Sigma=SLmz, emp=T)  # emp=T means exact data simulation, cov(Lmz) = SLmz  #
      # build the exact simulated data DZ
      Am_=Ldz[,1]
      Pm=Ldz[,2] # prs
      Af_=Ldz[,3]
      Pf=Ldz[,4] # prs
      At1r=Ldz[,5]
      Pt1r=Ldz[,6]
      At2r=Ldz[,7]
      Pt2r=Ldz[,8]
      #
      At1_=.5*Am_+.5*Af_+At1r # A res t1
      Pt1=.5*Pm+.5*Pf+Pt1r # Prs t1
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
      # standardize prs
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
      Pt1=.5*Pm+.5*Pf+Pt1r # Prs
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
      ## standardize prs
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
      # add sum of mother and father prs and add mean of twins prs ... we need these additional variables to
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
      ix_ = sort.int(phdatmzL_e[,2], index.return=T)
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
      ix_ = sort.int(phdatdzL_e[,2], index.return=T)
      phdatdzL_e = phdatdzL_e[ix_$ix,]
      colnames(phdatdzL_e)=c("zyg","famnr","ph","pgsm","pgsf","pgst","pgsmf","mpgst")
      phdatdzL_e=as.data.frame(phdatdzL_e)
      #
      phdatL_e=rbind(phdatmzL_e, phdatdzL_e)
      #
      #
      # simulated stochastically
      #                wide        wide     long      long      long mz+dz
      # data sets are  phdatdz and phdatdz, phdatdzL, phdatmzL, phdatL
      #                wide          wide        long       long        long mz+dz
      # simulated exactly
      # data sets are  phdatdz_e and phdatdz_e, phdatdzL_e, phdatmzL_e, phdatL_e
      #              pheno t1 pheno t2 mother    father   twin 1  nt twin1  twin2   nt twin2  m+f prs  mean twin prs
      # colnames [1] "pht1"   "pht2"   "pgsm"   "pgsf"   "pgst1"  "pgsnt1" "pgst2"  "pgsnt2" "pgsmf"  "mpgst"
      #
      # --------------------------------------------------- end data sim
      # start the analyses
      #
      # regression analyses. based on simulated data (not exact) ...
      #
      #
      # regression analyses. based on simulated data exact ... with power
      # DZ twin 1 only ...
      #
      eM0dz=lm(pht1~pgst1, data=phdatdz_e) #)$coefficients 			#   just regression of pheno on prs
      eM1dz=lm(pht1~pgsmf+pgst1, data=phdatdz_e) #)$coefficients 		#   with pgsmf ... this is equivalent to the transmitted / non-transmitted design
      eM2dz=lm(pht1~mpgst+pgst1, data=phdatdz_e) #)$coefficients 		#   with mpgst ...
      eM3dz=lm(pht1~pgsmf+mpgst+pgst1, data=phdatdz_e) # )$coefficients 	#   both detects two sources of cov(AC): twin interaction and cult transmission
      #round(summary(eM3dz)$coefficients,5)
      #
      # DZ twin 1 and twin 2 ... switch to geeglm
      #
      egeeM0dzL=geeglm(ph~pgst, id=famnr, corstr=cmethod,data=phdatdzL_e)#)$coefficients #
      egeeM1dzL=geeglm(ph~pgsmf+pgst, id=famnr, corstr=cmethod,data=phdatdzL_e)#)$coefficients #
      egeeM2dzL=geeglm(ph~mpgst+pgst, id=famnr, corstr=cmethod,data=phdatdzL_e)#)$coefficients  #
      egeeM3dzL=geeglm(ph~pgsmf+mpgst+pgst, id=famnr, corstr=cmethod,data=phdatdzL_e)#)$coefficients  #
      #
      # exact twins mz + dz gee
      # DZ and MZ twins
      egeeM0mzdzL=geeglm(ph~pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients #
      egeeM1mzdzL=geeglm(ph~pgsmf+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients #
      egeeM2mzdzL=geeglm(ph~mpgst+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients  #
      egeeM3mzdzL=geeglm(ph~pgsmf+mpgst+pgst, id=famnr, corstr=cmethod,data=phdatL_e)#)$coefficients  #

      # Power
      test_tmp <- c(
        summary(egeeM1mzdzL)$coefficients[2,3], ### ERROR
        summary(egeeM2mzdzL)$coefficients[2,3],
        summary(egeeM3mzdzL)$coefficients[2,3],
        summary(egeeM3mzdzL)$coefficients[3,3],

        summary(egeeM1dzL)$coefficients[2,3],
        summary(egeeM2dzL)$coefficients[2,3],
        summary(egeeM3dzL)$coefficients[2,3],
        summary(egeeM3dzL)$coefficients[3,3]
      )

      test_power_tmp <- sapply(test_tmp, function(ncp) {
        gnome_power(alpha, 1, ncp)
      })

      # Estimates
      estimates_tmp <- c(
        summary(egeeM1mzdzL)$coefficients[2,1], # CT Model 1 CT only
        summary(egeeM2mzdzL)$coefficients[2,1], # SI Model 1 SI only
        summary(egeeM3mzdzL)$coefficients[2,1], # CT Model 1 both
        summary(egeeM3mzdzL)$coefficients[3,1], # SI Model 1 both

        summary(egeeM1dzL)$coefficients[2,1], # CT Model 4 CT only
        summary(egeeM2dzL)$coefficients[2,1], # SI Model 1 SI only
        summary(egeeM3dzL)$coefficients[2,1], # CT Model 4
        summary(egeeM3dzL)$coefficients[3,1] # SI Model 4
      )

      geekeep[counter_within, ] <- c(test_power_tmp, estimates_tmp)

    }

    # For geekeep
    jpow <- 1:8
    jest <- 9:16

    row_index <- counter_overall - n_set + 1

    final_gee_estimates[row_index : counter_overall,1:10] <- setkeep[,1:10]
    final_gee_estimates[row_index : counter_overall,11:18] <- round(geekeep[,jest],3)
    final_gee_power[row_index : counter_overall,1:10] <- setkeep[,1:10]
    final_gee_power[row_index : counter_overall,11:18] <- round(geekeep[,jpow],3)

    counter_within = 0 # reset set counter for each PGS setting
  }
  # Re-name columns
  setnames <- c('nmz','ndz','a','c','e','g','b','x','PGS','A')
  colnames(final_gee_estimates) <- c(setnames, paste0("e", 1:8))
  colnames(final_gee_power) <- c(setnames, paste0("p", 1:8))

  # Use effect size function on the data sets
  final_gee_estimates <- final_gee_estimates %>%
    mutate(Smz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$mz,
           Sdz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$dz)
  final_gee_power <- final_gee_power %>%
    mutate(Smz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$mz,
           Sdz = gnome_effect(a = a, c = c, e = e, g = g, b = b)$dz)

  return(list(power = final_gee_power, params = final_gee_estimates))
}
