#' @export
mbd_alpha_loglik = function(s, path = "F://Dropbox//University//R//MBD//results//[2017-05-29]",
                            alpha.min = 1, alpha.max = 15){
  if(.Platform$OS.type == "unix")
  {
    home_dir <- substring(getwd(), 1, 13)
    mle_file <- paste0(home_dir,'/mbd_like/mbd_MLE.txt')
  }else
  {
    mle_file <- paste0(path,'/mbd_MLE.txt')
  }
  res <- unname(as.matrix(read.csv(mle_file, stringsAsFactors = F)))

  load(file = "data/general_settings")
  load(file = "data/sim_data")

  #ALPHA ANALYSIS
  LL.alpha = rep(NA, times = alpha.max)
  if ( any(res[,dim(res)[2]]==s) )
  {
    for (alpha in alpha.min:alpha.max)
    {
      if (.Platform$OS.type == "windows"){print(alpha)}
      best_pars <- res[which(res[,dim(res)[2]] == s), 1:3]
      LL.alpha[alpha] <- try( MBD:::mbd_loglik0(brts = sim_data[[s]], pars = best_pars, soc = soc,
                                                cond = cond, tips_interval = tips_interval, alpha = alpha) )
    }
    out <- c(LL.alpha,s)
    if(.Platform$OS.type == "unix"){file_name <- "mbd_alpha_loglik.txt"}else{file_name <- paste0(path,'/mbd_alpha_loglik.txt')}
    write.table(matrix(out,ncol = length(out)),file = file_name,append = T,row.names = F,col.names = F, sep = ",")
  }else
  {
    out <- NULL
  }
  return(out)
}


#' @export
mbd_alpha_ML = function(s, path = "F://Dropbox//University//R//MBD//results//[2017-05-29]",
                        alpha.min = 1, alpha.max = 15){
  initparsopt <- c(1.8, 0.3, 0.15); idparsopt <- 1:3; parsfix <- NULL;
  load(file = "data/general_settings")
  load(file = "data/sim_data")

  #ALPHA ANALYSIS
  MLE.alpha = vector("list",alpha.max);
  time.alpha = rep(NA,alpha.max);
    for (alpha in alpha.min:alpha.max)
      {
      if(.Platform$OS.type=="windows"){print(alpha)}
      time.alpha[alpha] = system.time( MLE.alpha[[alpha]] <- tryCatch( MBD:::mbd_ML0(brts=sim_data[[s]], alpha = alpha,
                                               initparsopt=initparsopt,
                                               idparsopt=idparsopt,
                                               idparsfix = (1:3)[-idparsopt],
                                               parsfix = parsfix,
                                               missnumspec=0,
                                               cond= cond,
                                               soc = soc,
                                               tips_interval=tips_interval,
                                               res = 10*(1+length(brts)+missnumspec),
                                               tol = c(1E-3, 1E-4, 1E-6),
                                               maxiter = 1000 * round((1.25)^length(idparsopt)),
                                               changeloglikifnoconv = FALSE,
                                               optimmethod = 'subplex')
                                     ))
    }
    out <- matrix(unlist(MLE.alpha),nrow=length(time.alpha),ncol=6,byrow = T);
    out[,5] <- time.alpha
    out[,6] <- alpha.min:alpha.max
    out <- cbind(out, rep(s, length(time.alpha)))

    if(.Platform$OS.type=="unix"){file_name <- "mbd_alpha_MLE.txt"}else{file_name <- paste0(path,'/mbd_alpha_MLE.txt')}
    write.table(matrix(out, ncol = dim(out)[2]), file = file_name,
                append = T, row.names = F, col.names = F, sep = ",")
  return(out)
}
