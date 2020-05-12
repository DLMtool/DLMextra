#' Length-based Integrated Mixed Effects (LIME) assessment model
#'
#' A wrapper function for the LIME model, a stock assessment model designed for data-limited scenarios with
#' length composition data and biological parameters. Optionally, catch or index abundance data can
#' also be included. Catch data should be used in order to estimate stock biomass for output-based harvest
#' control rules. A steepness value of one is used in typical applications. The reported FMSY reference point is a
#' proxy using either spawning potential ratio or Fmax (from yield-per-recruit). The BMSY reference point is
#' calculated using estimated mean recruitment and fishing at FMSY.
#'
#' @aliases lime
#' @param x An index for the objects in \code{Data} when running in \link[DLMtool]{runMSE}.
#' Otherwise, equals to 1 when running an assessment interactively.
#' @param Data An object of class Data.
#' @param add_catch Logical, whether to include catch data.
#' @param add_index Logical, whether to include index data.
#' @param ESS The maximum annual sample size of the length composition data.
#' @param CAL_dist Whether the model uses a multinomial \code{"mult"} or Dirichlet-multinomial \code{"Dirmult"}
#' likelihood for the length composition.
#' @param SigmaC The standard deviation of the catch in the likelihood.
#' @param SigmaI The standard deviation of the index in the likelihood.
#' @param SigmaR The standard deviation of the recruitment deviates in the likelihood.
#' @param SigmaF The standard deviation of the F random walk in the likelihood.
#' @param nseas The number of season in the model. In high F situations, multi-seasons may be needed to model smooth
#' length distributions.
#' @param Fproxy How the FMSY proxy reference point is calculated. The default is F40\%.
#' @param yind Optional, vector of years for the model. A subset of the data will be taken instead of the full
#' time series in the Data object.
#' @param integrate Logical, whether recruitment deviates are random effects in the model (TRUE) or penalized effects (FALSE).
#' FALSE by default.
#' @param run_LIME_args A named list of additional arguments for \link[LIME]{run_LIME}. Only arguments that are passed to \link[LIME]{format_input} are used.
#' @param control A named list of agruments for optimization to be passed to \code{\link[stats]{nlminb}}.
#' @param inner.control A named list of arguments for optimization of the random effects, which
#' is passed on to \code{\link[TMB]{newton}}.
#' @param ... Additional arguments (not currently used).
#' @return An object of \code{\linkS4class{Assessment}} containing assessment output.
#' @note
#' This is a wrapper function intended for a DLMtool \linkS4class{Data} object. The \code{LIME} package can be
#' downloaded from Github with \code{devtools::install_github("merrillrudd/LIME")}.
#'
#' @author Q. Huynh
#' @references
#' Rudd, M.B., and Thorson, J.T. 2017. Accounting for variable recruitment and fishing mortality in length-based
#' stock assessments for data-limited fisheries. Canadian Journal of Fisheries and Aquatic Sciences 75:1019-1035.
#' \url{https://doi.org/10.1139/cjfas-2017-0143}
#'
#' @section Required Data:
#' \itemize{
#' \item \code{LIME}: Mort, L50, L95, CAL, CAL_bins, vbK, vbLinf, vbt0, wla, wlb, LenCV, sigmaR
#' }
#' @section Optional Data:
#' \itemize{
#' \item \code{LIME}: Cat, Ind
#' }
#' @examples
#' \donttest{
#' library(MSEtool)
#' data(SimulatedData)
#'
#' res <- LIME(Data = SimulatedData)
#' plot(res)
#' summary(res)
#'
#' ## Use additional LIME functions
#' output <- res@@info
#' LIME::plot_LCfits(Inputs = output$Inputs, Report = output$Report)
#' LIME::plot_output(Inputs = output$Inputs, Report = output$Report, Sdreport = output$Sdreport)
#'
#' ## Create an MP that uses a F30% proxy and catch data
#' LIME_MP <- make_MP(DLMextra::LIME, HCR40_10, Fproxy = "F30", add_catch = TRUE)
#' }
#' @import MSEtool
#' @import dplyr
#' @importClassesFrom MSEtool Assessment
#' @importMethodsFrom MSEtool summary plot
#' @seealso \link[MSEtool]{plot.Assessment} \link[MSEtool]{summary.Assessment} \link[MSEtool]{make_MP}
#' @export
LIME <- function(x = 1, Data, add_catch = FALSE, add_index = FALSE, ESS = 50, CAL_dist = c("mult", "Dirmult"),
                 SigmaC = 0.2, SigmaI = 0.2, SigmaR = Data@sigmaR[x], SigmaF = 0.3, nseas = 1L,
                 Fproxy = c("F40", "F30", "Fmax"), yind = expression(1:ncol(Data@Cat)),
                 integrate = FALSE, run_LIME_args, control = list(iter.max = 2e+05, eval.max = 4e+05),
                 inner.control = list(), ...) {
  check_dependencies("LIME", "LIME")
  dots <- list(...)
  dependencies <- "Data@Mort, Data@L50, Data@L95, Data@CAL, Data@CAL_bins, Data@vbK, Data@vbLinf, Data@vbt0, Data@wla, Data@wlb, Data@LenCV, Data@sigmaR"
  CAL_dist <- match.arg(CAL_dist)
  Fproxy <- match.arg(Fproxy)
  yind <- eval(yind)

  LH_args <- list(vbk = Data@vbK[x], linf = Data@vbLinf[x], t0 = Data@vbt0[x],
                  M = Data@Mort[x], lwa = Data@wla[x], lwb = Data@wlb[x],
                  M50 = Data@L50[x], M95 = Data@L95[x], maturity_input = "length",
                  S50 = Data@L50[x], S95 = Data@L95[x], selex_input = "length", #placeholders; sel is estimated
                  nseasons = nseas, nfleets = 1L, CVlen = Data@LenCV[x],
                  SigmaC = SigmaC, SigmaI = SigmaI, SigmaR = SigmaR, SigmaF = SigmaF,
                  binwidth = unique(diff(Data@CAL_bins)))
  LH <- do.call(LIME::create_lh_list, LH_args)

  data_all <- list(years = yind, LF = structure(Data@CAL[x, yind, ], dimnames = list(yind, Data@CAL_bins[-1])),
                   neff_ft = pmin(ESS, rowSums(Data@CAL[x, yind, ]), na.rm = TRUE) %>% matrix(nrow = 1))
  data_avail <- paste0("LC", length(yind))

  if(add_index) {
    data_all$I_ft <- Data@Ind[x, yind, drop = FALSE]
    data_avail <- paste0(data_avail, "_Index")
  }
  if(add_catch) {
    data_all$C_ft <- Data@Cat[x, yind, drop = FALSE]
    data_avail <- paste0(data_avail, "_Catch")
  }

  LIME_args <- formals(LIME::run_LIME)
  input <- LIME::create_inputs(LH, data_all)
  format_input_args <- list(input = input, data_avail = data_avail, Fpen = 1L,
                            SigRpen = 1L, SigRprior = eval(LIME_args$SigRprior),
                            LFdist = as.integer(CAL_dist == "Dirmult"), C_type = ifelse(add_catch, 2, 0),
                            est_more = FALSE, fix_more = FALSE, est_F_ft = TRUE,
                            f_startval_ft = NULL, rdev_startval_t = NULL, est_selex_f = TRUE,
                            vals_selex_ft = matrix(-1, input$nfleets, length(input$highs)),
                            est_rdev_t = TRUE, mirror = NULL, est_totalF = FALSE,
                            prop_f = rep(1/input$nfleets, input$nfleets))

  if(!missing(run_LIME_args)) {
    args_ind <- match(names(run_LIME_args), names(format_input_args), nomatch = 0)
    format_input_args[args_ind] <- run_LIME_args[as.logical(args_ind)]
  }

  if(!is.null(dots$LIME_int) && dots$LIME_int) {
    format_input_args$modpath <- NULL
    out <- do.call(LIME::run_LIME, c(format_input_args, list(modpath = NULL)))
    return(out)
  }

  TmbList <- do.call(LIME::format_input, format_input_args)
  if(!integrate) TmbList$Random <- NULL

  obj <- TMB::MakeADFun(data = TmbList[["Data"]], parameters = TmbList[["Parameters"]],
                        random = TmbList[["Random"]], map = TmbList[["Map"]],
                        inner.control = inner.control, DLL = "LIME", silent = TRUE)

  F_up <- 3 # Maximum F

  Upr = rep(Inf, length(obj$par))
  Upr[match("log_sigma_R",names(obj$par))] = log(2)
  #if(is.null(S50_up)==FALSE) Upr[which(names(obj$par)=="log_S50_f")] <- log(S50_up)
  #if(is.null(S50_up)) Upr[which(names(obj$par)=="log_S50_f")] <- log(input$linf)
  Upr[which(names(obj$par)=="log_S50_f")] <- log(input$linf)
  Upr[which(names(obj$par)=="log_F_ft")] = log(F_up)
  Upr[match("log_sigma_F", names(obj$par))] <- log(2)
  # Upr[which(names(obj$par)=="log_theta")] <- log(10)

  Lwr <- rep(-Inf, length(obj$par))
  Lwr[match("log_CV_L",names(obj$par))] = log(0.001)
  Lwr[match("log_sigma_C",names(obj$par))] = log(0.001)
  Lwr[match("log_sigma_I",names(obj$par))] = log(0.001)
  Lwr[which(names(obj$par)=="log_S50_f")] = log(1)

  opt <- try(nlminb(obj$par, obj$fn, obj$gr, control = control, lower = Lwr, upper = Upr), silent = TRUE)
  SD <- TMB::sdreport(obj, bias.correct = length(TmbList[["Random"]]) > 0)
  Report <- obj$report(obj$env$last.par.best)

  # This is Fmax, not FMSY
  Derived <- list()
  Derived$Fmsy <- optimize(LIME::calc_msy, ages=input$ages, M=input$M, R0=exp(Report$beta), W_a=input$W_a,
                           S_fa=Report$S_fa, lower=0, upper=10, maximum=TRUE)$maximum
  Derived$FFmsy <- Report$F_t/Derived$Fmsy
  Derived$msy <- LIME::calc_msy(F=Derived$Fmsy, ages=input$ages, M=input$M, R0=exp(Report$beta), W_a=input$W_a, S_fa=Report$S_fa)
  Derived$Bmsy <- sum(LIME::calc_equil_abund(ages=input$ages, M=input$M, F=Derived$Fmsy, S_fa=Report$S_fa, R0=exp(Report$beta)) * input$W_a)
  Derived$BBmsy <- Report$TB_t/Derived$Bmsy
  Derived$SBmsy <- sum(LIME::calc_equil_abund(ages=input$ages, M=input$M, F=Derived$Fmsy, S_fa=Report$S_fa, R0=exp(Report$beta)) * input$W_a * input$Mat_a)
  Derived$SBBmsy <- Report$SB_t/Derived$SBmsy

  F30 <- tryCatch(uniroot(LIME::calc_ref, lower=0, upper=50, ages=input$ages, Mat_a=input$Mat_a, W_a=input$W_a, M=input$M, S_fa=Report$S_fa, ref=0.3)$root, error=function(e) NA) * input$nseasons
  F40 <- tryCatch(uniroot(LIME::calc_ref, lower=0, upper=50, ages=input$ages, Mat_a=input$Mat_a, W_a=input$W_a, M=input$M, S_fa=Report$S_fa, ref=0.4)$root, error=function(e) NA) * input$nseasons
  FF30 <- FF40 <- NULL
  if(is.na(F30)==FALSE) FF30 <- Report$F_t/F30
  if(is.na(F40)==FALSE) FF40 <- Report$F_t/F40

  Derived$F30 <- F30
  Derived$F40 <- F40
  Derived$FF30 <- FF30
  Derived$FF40 <- FF40

  df <- data.frame("gradient" = as.vector(obj$gr(opt$par)), "parameter" = names(obj$par),
                   "estimate"= opt$par, "transformed" = exp(opt$par))

  output <- list(input = input, data_avail = data_avail, Inputs = TmbList, Report = Report, Sdreport = SD,
                 obj = obj, opt = opt, Derived = Derived, df = df)

  # Assessment
  Year <- yind
  nll_report <- ifelse(is.character(opt), Report$jnll, opt$objective)
  sel_at_age <- Report$S_fa[1, ]
  VB <- colSums(sel_at_age * t(Report$N_ta) * Report$W_a)
  N0 <- sum(LIME::calc_equil_abund(ages=input$ages, M=input$M, F=0, S_fa=Report$S_fa, R0=exp(Report$beta)))
  B0 <- sum(LIME::calc_equil_abund(ages=input$ages, M=input$M, F=0, S_fa=Report$S_fa, R0=exp(Report$beta)) *
              input$W_a)
  VB0 <- sum(LIME::calc_equil_abund(ages=input$ages, M=input$M, F=0, S_fa=Report$S_fa, R0=exp(Report$beta)) *
               input$W_a * sel_at_age)

  Assessment <- new("Assessment", Model = "LIME", Name = Data@Name, conv = !is.character(SD) && SD$pdHess,
                    N0 = N0, R0 = exp(Report$beta), SSB0 = Report$SB0, B0 = B0,
                    h = Report$h, FMort = structure(Report$F_t, names = Year),
                    B = structure(Report$TB_t, names = Year),
                    B_B0 = structure(Report$TB_t/B0, names = Year),
                    SSB = structure(Report$SB_t, names = Year),
                    SSB_SSB0 = structure(Report$SB_t/Report$SB0, names = Year),
                    VB = structure(VB, names = Year),
                    VB_VB0 = structure(VB/VB0, names = Year),
                    R = structure(Report$R_t, names = Year),
                    N = structure(Report$N_t, names = Year),
                    N_at_age = Report$N_ta,
                    Selectivity = t(sel_at_age),
                    Catch = structure(Report$Cw_t_hat, names = Year),
                    Index = structure(Report$I_ft_hat[1, ], names = Year),
                    Obs_C_at_age = Report$Cn_ta,
                    Dev = Report$Nu_input, Dev_type = "log-Recruitment deviations",
                    NLL = structure(nll_report, names = "Total"),
                    info = output, obj = obj, opt = opt, SD = SD, TMB_report = Report,
                    dependencies = dependencies)

  if(add_catch) Assessment@Obs_Catch <- structure(data_all$C_ft[1, ], names = Year)
  if(add_index) Assessment@Obs_Index <- structure(data_all$I_ft[1, ], names = Year)

  if(Fproxy == "F30") {
    Assessment@FMSY <- F30
  } else if(Fproxy == "F40") {
    Assessment@FMSY <- F40
  } else {
    Assessment@FMSY <- Fmsy
  }

  Assessment@SSBMSY <- sum(LIME::calc_equil_abund(ages=input$ages, M=input$M, F=Assessment@FMSY, S_fa=Report$S_fa, R0=exp(Report$beta)) *
                             input$W_a * input$Mat_a)
  Assessment@BMSY <- sum(LIME::calc_equil_abund(ages=input$ages, M=input$M, F=Assessment@FMSY, S_fa=Report$S_fa, R0=exp(Report$beta)) *
                           input$W_a)
  Assessment@VBMSY <- sum(LIME::calc_equil_abund(ages=input$ages, M=input$M, F=Assessment@FMSY, S_fa=Report$S_fa, R0=exp(Report$beta)) *
                            input$W_a * sel_at_age)

  Assessment@F_FMSY <- Assessment@FMort/Assessment@FMSY
  Assessment@SSB_SSBMSY <- Assessment@SSB/Assessment@SSBMSY
  Assessment@B_BMSY <- Assessment@B/Assessment@BMSY
  Assessment@VB_VBMSY <- Assessment@VB/Assessment@VBMSY

  return(Assessment)
}
class(LIME) <- "Assess"

