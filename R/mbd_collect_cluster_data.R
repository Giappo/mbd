#HOW TO CREATE KEYS ON THE CLUSTER
# ssh-keygen -t rsa -C "glaudanno@gmail.com" -f "ir_rsa"
# mkdir ~/.ssh
# chmod 700 ~/.ssh
# cat id_rsa.pub >> ~/.ssh/authorized_keys

# @Giappo: add doc
#' Does something G
#' @export
get_data_from_cluster = function (sim_pars, s, account_name = "p274829" ){

  if(account_name=="cyrus" || account_name=="Cyrus"){account_name = "p257011"} #cyrus account_name = "257011"
  if(account_name=="giovanni" || account_name=="Giovanni"){account_name = "p274829"} #giovanni account_name = "274829"

  folder_name =  NULL
  for (i in 1:length(sim_pars)){
    folder_name = paste(folder_name,toString(sim_pars[i]),sep = '')
    if (i!=length(sim_pars)){folder_name = paste(folder_name,"-",sep = '')}
  }

  file_name = paste("mbd_MLE",s,".txt",sep = '')

  host_address = "peregrine.hpc.rug.nl"
  # remote_path = "/home/p274829/mbd_like2/sims/"
  remote_path = paste("/home/",account_name,"/mbd_like/sims/",sep = '')
  key_path = "C:/Users/P274829/Documents/keys/"
  remote_file = paste(remote_path,folder_name,"/",file_name,sep = '')

  out = tryCatch(
  RCurl::scp(host = host_address,path = remote_file,user = account_name,verbose = F,
                key = sprintf(c(
                  ssh.public.keyfile  = paste(key_path,"id_rsa.pub",sep = ''),
                  ssh.private.keyfile = paste(key_path,"id_rsa",sep = ''),
                  "rsa"
                ))
  ), error = function(e) return(NULL)
  )
  if (!is.null(out)){
    out = rawToChar(out) #does something!!!!
    out = strsplit(x = out, split = ",", fixed = T)
    out = as.numeric(unlist(out))
  }
  return(out)
}
# out=get_data_from_cluster(s=15);out

# @Giappo: add doc
#' Does something H
#' @export
get_all_data_from_cluster = function (sim_pars, max_sims = 1000, account_name = "p274829"){

  Npars = length(sim_pars)
  folder_name =  NULL
  for (i in 1:length(sim_pars)){
    folder_name = paste(folder_name,toString(sim_pars[i]),sep = '')
    if (i!=length(sim_pars)){folder_name = paste(folder_name,"-",sep = '')}
  }

  home_dir = substring(getwd(),1,21)
  data_folder = paste0(home_dir,"/Progress/RQ1-multiple-birth-death/results/new_results/")
  # data_folder = "F:/Dropbox/University/Progress/RQ1 - Multiple Births/Results + Reports/4_parameters/"
  folder_list = dir(path = data_folder )
  right_folder <- suppressWarnings( max( which( grepl(folder_name, folder_list) ) ) )

  if ( is.infinite(right_folder) ){
    new_folder_name = paste("[",Sys.Date(),"] ",folder_name,sep = '')
    dir.create(file.path(data_folder,new_folder_name ), showWarnings = FALSE)
    folder_list = dir(path = data_folder )
    right_folder <- suppressWarnings( max( which( grepl(folder_name, folder_list) ) ) )
  }

  local_path = paste0(data_folder, folder_list[right_folder])

  res_files = list.files(pattern=paste('[.]txt',sep = ''),path=local_path, full.names=TRUE)
  suppressWarnings(rm(targetTable))
  local_results = NULL
  if (length(res_files)>0){ #loading data from local
    for(i in 1:length(res_files)){
      fileData<-read.table(file=res_files[i],header=FALSE,sep=",")
      ifelse(exists("targetTable"),targetTable<-rbind(targetTable,fileData),targetTable<-fileData)
    }
    local_results = targetTable;
    N_outputs = dim(local_results)[2]
    ids = local_results[,N_outputs]
  }
  suppressWarnings(rm(targetTable))
  ifelse ( exists("ids"),missing_ids <- (1:max_sims)[-ids],missing_ids <- (1:max_sims) )

  flush.console()
  at_least_one_result_from_cluster = 0
  for (s in missing_ids){ #loading data from cluster
    print(s)
    out=get_data_from_cluster(s = s, sim_pars = sim_pars, account_name = account_name)
    if (!is.null(out)){
    at_least_one_result_from_cluster == 1
    write.table(matrix(out,ncol = length(out)),file = paste(local_path,"/mbd_MLE",s,".txt",sep = ''),append = F,row.names = F,col.names = F, sep = ",")
    ifelse(exists("targetTable"),targetTable<-rbind(targetTable,out),targetTable<-out)
    }
  }
  all_results <- NULL
  if(exists("targetTable")){all_results <- targetTable}

  suppressWarnings(rm(targetTable))
  all_results = rbind(all_results,local_results);
  if (Npars == 4){parnames <- c("lambda","mu","nu","q");}else if (Npars == 3){parnames <- c("nu","mu","q");}
  if (length(dim(all_results))>1){
  colnames(all_results) = (c(parnames,"LL","multiple_born","number_of_tips","percentage_multiple_species","tree_id"))
  }else{names(all_results) = (c(parnames,"LL","multiple_born","number_of_tips","percentage_multiple_species","tree_id"))}
  rownames(all_results) = NULL

  tidy_results=all_results[order(all_results[,dim(all_results)[2]]),]
  dimnames(tidy_results)[[2]]<-(c(parnames,"LL","multiple_born","number_of_tips","percentage_multiple_species","tree_id"))
  bad_results=tidy_results[rowSums(tidy_results[,1:(Npars+1)]==rep(-1,(Npars+1)))==(Npars+1),];#print(Nbad<-dim(bad_results)[1])
  results=tidy_results[rowSums(tidy_results[,1:(Npars+1)]==rep(-1,(Npars+1)))!=(Npars+1),];#print(N<-dim(results)[1])
  print(paste("There are ",Nbad<-dim(bad_results)[1]," bad results.",sep = ''))
  print(paste("There are ",N<<-dim(results)[1]," good results.",sep = ''))

  working_account = "noone"
  if (at_least_one_result_from_cluster){working_account=account_name}

  # path <<- local_path
  return(list(results=results,path=local_path,working_account=working_account))
}
