mbd_setup=function(){
  #home_dir=NULL

  #loading previously saved data and settings
  load("data/general_settings", envir=globalenv())
  load("data/sim_data", envir=globalenv())

  #updating DDD to latest version (not yet on CRAN)
  if (packageVersion("DDD") != 3.4){
    if (.Platform$OS.type == "windows")
    {
      # install.packages(paste(home_dir,'libs/DDD_3.4.zip',sep = ''), repos = NULL, type = "source")
      install.packages("libs/DDD_3.4.zip", repos = NULL, type = "source")
    }
    if (.Platform$OS.type == "unix")
    {
      # install.packages(paste(home_dir,'libs/DDD_3.4.tar.gz',sep = ''), repos = NULL, type = "source")
      install.packages("libs/DDD_3.4.tar.gz", repos = NULL, type = "source")
    }
  }
}
