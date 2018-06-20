########################################################################################################################
#set variables

install.packages("RCurl", dependencies = TRUE); library("RCurl")
install.packages("pbdRPC", dependencies = TRUE); library("pbdRPC")

# local_path  = "F:/Dropbox/University/Progress/RQ1 - Multiple Births/Results + Reports"
# local_path2 = "F:/Dropbox/University/Progress/"
# local_path3 = "~/"

# system("dir F:/Dropbox/University/R/MBD")
# setwd(dir = "F:/Dropbox/University/R/MBD")

host_address = "peregrine.hpc.rug.nl"
remote_path = "/home/p274829/mbd_like2/sims/0.2-0.15-2.5-0.1/"
key_path = "C:\\Users\\P274829\\Documents\\keys\\"
key_path2 = "C:/Users/P274829/Documents/keys/"
file_name = "mbd_MLE1.txt"
remote_file = paste(remote_path,file_name,sep = '')
my.private.key = paste(key_path,"id_rsa",sep = '')
my.private.key2 = paste(key_path2,"id_rsa",sep = '')

########################################################################################################################
#test pbdRPC

rpc(cmd = "type /home/p274829/mbd_like2/zzz_ML.R",
    user = "p274829",
    hostname = "peregrine.hpc.rug.nl",
    priv.key = paste(key_path2,"id_rsa",sep = '')
    )

rpc("ls")
rpcopt_set( user = "p274829" , hostname = host_address , pport ="22", priv.key = my.private.key );.pbd_env$RPC.LI
# rpcopt_set( user = "p274829" , hostname = host_address , pport ="22",priv.key = my.private.key2 );.pbd_env$RPC.LI
.pbd_env$RPC.CT$use.shell.exec <- FALSE
rpc(machine(host_address),"whoami")
pbdRPC::
plink(args = "p274829@peregrine.hpc.rug.nl whoami",
      use.shell.exec = .pbd_env$RPC.CT$use.shell.exec,
      pause = .pbd_env$RPC.CT$pause, intern = .pbd_env$RPC.CT$intern,
      wait = .pbd_env$RPC.CT$wait)
check_exec("ping www.google.nl")

plink("peregrine.hpc.rug.nl 'whoami' ")
rpc("ping www.google.nl")
system('dir "C:\\Users\\P274829\\Documents\\keys"')

args <- "snoweye@192.168.56.101 Rscript -e 'whoami'"
ssh(args)      # Note ssh uses "-p" for server port.
plink(args)    # Note plink uses "-P" for server port.


########################################################################################################################
#test shell
shell(cmd = "pscp -scp ~/zzz_ML.R p274829@peregrine.hpc.rug.nl:/home/p274829/mbd_like2/zzz_ML.R"
      # shell = NULL
      # translate = T
      )
shell(cmd = "ping www.google.nl")
# key = sprintf(c(
#   ssh.public.keyfile  = paste(key_path,"id_rsa.pub",sep = ''),
#   ssh.private.keyfile = paste(key_path,"id_rsa",sep = ''),
#   "rsa"
# ))

shell.exec("p274829@peregrine.hpc.rug.nl:/home/p274829/mbd_like2/zzz_ML.R")

########################################################################################################################
#test scp

x=RCurl:::scp(host = host_address,path = remote_file,user = "p274829",verbose = T,
              key = sprintf(c(
                ssh.public.keyfile  = paste(key_path,"id_rsa.pub",sep = ''),
                ssh.private.keyfile = paste(key_path,"id_rsa",sep = ''),
                "rsa"
              ))
); rawToChar(x) #does something!!!!

# x=RCurl:::scp(host = host_address,path = remote_file,user = "p274829",keypasswd = get("G.chiave",envir=.pkgglobalenv),verbose = T,
#             key = sprintf(c(
#               ssh.public.keyfile  = paste(key_path,"id_rsa.pub",sep = ''),
#               ssh.private.keyfile = paste(key_path,"id_rsa",sep = ''),
#               "rsa"
#             )
#             )
# ); rawToChar(x) #does something!!!!

# remoteLogin(
#   host_address,
#   session = TRUE,
#   diff = TRUE,
#   commandline = TRUE,
#   username = p274829,
#   password = get("G.chiave",envir=.pkgglobalenv)
# )

########################################################################################################################
#alpha tests
#setup #rm(list=ls())
base_path=dirname(dirname(getwd()))
if( !exists("path") && interactive() ){
  path=choose.dir(paste(base_path,"/Progress/RQ1 - Multiple Births/Results + Reports/4_parameters/",sep = ''), "Choose a suitable folder")
}
load( as.character(paste(path,"/data/general_settings",sep = '')) );load( as.character(paste(path,"/data/sim_data",sep = '')) )
i=2;brts=sim_data[[i]];pars=sim_pars;brts;
pars;abstol=1e-16; reltol=1e-10;methode = "expo";missnumspec=0;
lambda=pars[1]; mu = pars[2]; nu=pars[3]; q=pars[4]; min_tips=tips_interval[1]; max_tips=tips_interval[2];
soc=2;cond=1 #watch out!

# i=1;sim_data[[i]]
# system.time(test1<-MBD:::pmb_loglik(pars=sim_pars,brts=sim_data[[i]]));test1
system.time(test2<-MBD:::mbd_loglik(pars=sim_pars,brts=sim_data[[i]],cond=0));test2
system.time(test3<-MBD:::mbd_loglik(pars=sim_pars,brts=sim_data[[i]],cond=1));test3

alphamax=30;testPc=rep(NA,alphamax)
for (alpha in 10:alphamax){
  testPc[alpha]=MBD:::calculate_conditional_probability(brts=brts,pars=pars,tips_interval=tips_interval,cond=cond,soc=soc,alpha=alpha,methode=methode,
                                                        abstol=abstol, reltol=reltol, minimum_multiple_births=minimum_multiple_births)$Pc
}
plot(testPc)

#experiment
# thresholds=rep(NA,2000)
thresholds[i]=testPc[alphamax];thresholds[!is.na(thresholds)]

########################################################################################################################
#data analysis
hist(  results[,1]-results[,2]  ,xlim = c(-1,1),breaks = 50)
abline(v=(sim_pars[1]-sim_pars[2]),col="red")
abline(v=median(results[,1]-results[,2]),col="blue")

fx = function (res,sim_pars=sim_pars,a,b,c,d,...){
 out=a*res[,1]+b*res[,2]+c*res[,3]*res[,4]^d
 hist(  out,...)
 abline(v=(a*sim_pars[1]+b*sim_pars[2]+c*sim_pars[3]*sim_pars[4]^d),col="red")
 abline(v=median(out),col="blue")
 ritorna = abs(median(out)-sim_pars)/sim_pars
 return(out)
}
fx(res=results,sim_pars = sim_pars,a=1,b=-1,c=1,d=2,xlim = c(-1,1),breaks = 50)

