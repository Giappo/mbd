#' Imports data
#' No globals please
#' @inheritParams default_params_doc
# @param parnames sonething
# mbd_import_data <- function(
#   parnames = c("lambda","mu","nu","q")
# ) {
#   parnames <<- parnames
#   Npars <<- length(parnames);
#   base_path=dirname(dirname(getwd()))
#   if( !exists("path")  ){ #&& interactive()
#     path<<-choose.dir(paste(base_path,"/Progress/RQ1 - Multiple Births/Results + Reports/4_parameters/",sep = ''), "Choose a suitable folder")
#   }
#   res_files = list.files(pattern=paste('[.]txt',sep = ''),path=path, full.names=TRUE)
#   for(i in 1:length(res_files)){
#     fileData <- utils::read.table(file=res_files[i],header=FALSE,sep=",")
#     ifelse(exists("targetTable"),targetTable<-rbind(targetTable,fileData),targetTable<-fileData)
#   }
#   all_results = targetTable
#   
#   tidy_results=all_results[order(all_results[,dim(all_results)[2]]),]
#   dimnames(tidy_results)[[2]]<-(c(parnames,"LL","multiple_born","number_of_tips","percentage_multiple_species","tree_id"))
#   bad_results=tidy_results[rowSums(tidy_results[,1:(Npars+1)]==rep(-1,(Npars+1)))==(Npars+1),];#print(Nbad<-dim(bad_results)[1])
#   results=tidy_results[rowSums(tidy_results[,1:(Npars+1)]==rep(-1,(Npars+1)))!=(Npars+1),];#print(N<-dim(results)[1])
#   print(paste("There are ",Nbad<-dim(bad_results)[1]," bad results.",sep = ''))
#   print(paste("There are ",N<<-dim(results)[1]," good results.",sep = ''))
#   
#   load( as.character(paste(path,"/data/general_settings",sep = '')), envir = e <- globalenv() )
#   load( as.character(paste(path,"/data/sim_data",sep = '')), envir = e <- globalenv()  )
#   suppressWarnings( rm(targetTable,all_results,fileData,i) )
#   return(results)
# }

#' analysis function
#' @inheritParams default_params_doc
#' @param path something
#' @param titolo something
#' @param pdfname something
#' @param percentage_hidden_outliers something
#' @param openit something 
#' @param mother_folder something
correlation_analysis <- function(
  results,
  path,
  titolo = NULL,
  pdfname,
  sim_pars = sim_pars,
  percentage_hidden_outliers = 0.04, 
  openit = 0, 
  idparsopt, 
  mother_folder
) {
  Npars=length(sim_pars);#pars_interval=list(c(0,quantile(results[,1],.95)),c(0,quantile(results[,2],.95)),c(0,quantile(results[,3],.95)))
  if (missing(idparsopt)){estimated_pars=1:Npars}else{estimated_pars=idparsopt}
  par_names = colnames(results)[1:Npars]

  truevalues_color="red"; points_color = "azure4"; medians_color = "blue3"#"chartreuse3";
  medians_color_name = "blue"; truevalues_color_name = "red";

  medians=rep(0,Npars);for (idpar in 1:Npars){medians[idpar]=stats::median(results[,idpar])}
  medians_string=paste0( "MLE Medians (",medians_color_name,") = (",paste(signif(medians,2),sep = "''",collapse = ", "),")")
  truevalues_string=paste0( "True Values (",truevalues_color_name,") = (",paste(signif(sim_pars,2),sep = "''",collapse = ", "),")")
  axislimits=rep(NA,Npars)
  for (i in 1:Npars){
    axislimits[i] <- stats::quantile(
      results[,i],
      probs = 1 - percentage_hidden_outliers
    )
  }

  #pdf
  grDevices::pdf(file = paste(path,"/",pdfname,".pdf",sep=''));#plot.new();
  graphics::par(mfrow=c(length(estimated_pars),length(estimated_pars)))
  graphics::par(oma=c(0,0,2,0));
  for (i in estimated_pars){for (j in estimated_pars){
    good.lines = results[,i]<axislimits[i] & results[,j]<axislimits[j]
    ifelse(any(good.lines)>0,good.results <- results[good.lines,],good.results <- results)

    if (i==j){graphics::hist((good.results[,i]),main=NULL,xlab = paste(par_names[i]),breaks = 15); #,breaks = 15
      graphics::abline(v=sim_pars[i],col = truevalues_color)
      graphics::abline(v=medians[i],col = medians_color)
    }
    else{graphics::plot(good.results[,i]~good.results[,j],xlab=par_names[j],ylab=par_names[i],cex=0.3,col=points_color);
      graphics::points(x=sim_pars[j],y=sim_pars[i],col=truevalues_color,pch=10,cex=1.5)
      graphics::points(x=medians[j],y=medians[i],col=medians_color,pch=10,cex=1.5)
    }
  }}
  graphics::title(main=(titolo.pdf<-(paste("\n\n",titolo,"\n",medians_string,"\n",truevalues_string,sep = ''))),outer=T)
  grDevices::dev.off()

  if (openit==1){
    file.show(normalizePath(paste(path,"/",pdfname,".pdf",sep='')))
  }

  if (!missing(mother_folder))
    {
    grDevices::pdf(file = paste(mother_folder,"/",pdfname,".pdf",sep=''));#plot.new();
    graphics::par(mfrow=c(length(estimated_pars),length(estimated_pars)))
    graphics::par(oma=c(0,0,2,0));
    for (i in estimated_pars){for (j in estimated_pars){
      good.lines = results[,i]<axislimits[i] & results[,j]<axislimits[j]
      ifelse(any(good.lines)>0,good.results <- results[good.lines,],good.results <- results)

      if (i==j){graphics::hist((good.results[,i]),main=NULL,xlab = paste(par_names[i]),breaks = 15); #,breaks = 15
        graphics::abline(v=sim_pars[i],col = truevalues_color)
        graphics::abline(v=medians[i],col = medians_color)
      }
      else{graphics::plot(good.results[,i]~good.results[,j],xlab=par_names[j],ylab=par_names[i],cex=0.3,col=points_color);
        graphics::points(x=sim_pars[j],y=sim_pars[i],col=truevalues_color,pch=10,cex=1.5)
        graphics::points(x=medians[j],y=medians[i],col=medians_color,pch=10,cex=1.5)
      }
    }}
    graphics::title(main=titolo.pdf,outer = T)
    grDevices::dev.off()
    }

}

#' Percentiles function
#' @inheritParams default_params_doc
percentiles_function <- function(
  results, 
  sim_pars, 
  printit = 1, 
  quantiles_choice = c(.25, .50, .75)
) {
  quantiles_names <- format(round(quantiles_choice,2),nsmall = 2)
  Npars <- length(sim_pars);#pars_interval=list(c(0, stats::quantile(results[,1],.95)),c(0,quantile(results[,2],.95)),c(0,quantile(results[,3],.95)))
  parnames <-  colnames(results)[1:Npars]
  percentiles <- vector("list",3)
  for (idpar in 1:Npars){percentiles[[idpar]] <- stats::quantile(results[,idpar], quantiles_choice)}
  percentiles <- t(matrix(unlist(percentiles), nrow = 3, ncol = Npars));
  colnames(percentiles) <- quantiles_names; rownames(percentiles) <- parnames
  if (printit==1)
  {
    print(percentiles); print(sim_pars)
  }
  out <- percentiles
}

#' Converts a matrix in a dataframe that can be used with ggplot
#' @inheritParams default_params_doc
#' @param Matrix something
#' @param heatmap.name heatmap name
heatmap2dataframe = function(
  x,
  y,
  Matrix, 
  x.name="x", 
  y.name="y", 
  heatmap.name="HeatMap"
){
  Matrix2=matrix(Matrix,nrow = length(x)*length(y))
  df=data.frame(expand.grid(x=x,y=y),HeatMap=Matrix2)
  #names(df$x)=names(x);names(df$y)=names(y);names(df$HeatMap)=names(Matrix);
  colnames(df)[1] <- x.name;
  colnames(df)[2] <- y.name;
  colnames(df)[3] <- heatmap.name;
  return(df)
}

#' Converts branching times to "time intervals between branching times" and "birth at nodes" vectors
#' @inheritParams default_params_doc
brts2time_intervals_and_births <- function(brts
) {
  time_points=-unlist(unname(sort(abs(brts),decreasing = T)) )
  branching_times = -sort(abs(as.numeric(time_points)),decreasing = TRUE)
  births=c(0,unname(table(branching_times))[-1])
  unique_branching_times=as.numeric(names(table(branching_times)))
  time_intervals=c((diff(unique_branching_times)),(abs(utils::tail(unique_branching_times,1))) )
  births=births[-1]
  out=list(time_intervals=time_intervals,births=births)
  return(out)
}

#' Provides estimated number of species integrating the P-equation
#' @inheritParams default_params_doc
#' @param test_parameters something
#' @param output something
mbd_P_eq <- function(
  test_parameters,
  age=15,
  max_number_of_species = 2000, 
  precision = 50L,
  soc = 2,
  output = 0
){
  N0=soc
  
  # Use sim_pars for pars
  create_mbd_P_matrix <- function(
    pars,
    max_number_of_species
  )
  {
    lambda=pars[1];mu=pars[2];q=pars[3];nvec=0:max_number_of_species
    M=matrix(0,nrow=max_number_of_species+1,ncol=max_number_of_species+1)
    testit::assert(!"Do not call hyperA")
    # M = lambda * hyperA::hyperA(N = max_number_of_species,k = 0,q = q)
    # M[row(M) == col(M) - 1] = mu*nvec[2:(max_number_of_species+1)]
    # diag(M) = - mu*nvec - lambda*(1-(1-q)^nvec)
    # return(M)
  }

  nvec=0:max_number_of_species
  v0=rep(0,(max_number_of_species+1));v0[N0+1]=1
  transition_matrix = create_mbd_P_matrix(pars=test_parameters, max_number_of_species = max_number_of_species)
  vf=expoRkit::expv(v=v0,x=transition_matrix,t=age,m = precision)

  nmedio=sum(nvec*vf)
  std=sqrt( sum(nvec^2*vf)-nmedio^2 )
  if (output==1){
    graphics::plot(log(vf))
    print(paste("Sim pars are:",test_parameters[1],test_parameters[2],test_parameters[3],". Average n is",nmedio,"with std:",std))
  }
  return(list(avg_n=nmedio,std_n=std))
}

# @Giappo: add doc
#' Gets heatmap of my matrices.
#' @inheritParams default_params_doc
#' @param ... something
#' @export
myheatmap <- function(
  matrix,
  logs=1,
  colormap = grDevices::heat.colors(15),
  ...
){
  if (is.matrix(matrix)==F){matrix=as.matrix(matrix)}
  T=t(matrix[nrow(matrix):1,1:ncol(matrix)])
  if (logs ==1 ) {
  # image(log(T))
    graphics::image(log(T),col=rev(colormap),...)
  }else{graphics::image(T)}
}

# @Giappo: add doc
#' Does something J
#' @inheritParams default_params_doc
myheatmap2 <- function(x,y,z,x.name,y.name,z.name,x.splits,y.splits){

  if (missing(x.splits)){x.splits=round( (length(x))/10 )}
  if (missing(y.splits)){y.splits=round( (length(y))/10 )}
  if (missing(x.name)){x.name=NULL}
  if (missing(y.name)){y.name=NULL}
  if (missing(z.name)){z.name=NULL}

  lX=x; lY=y;
  pretty.X.at	<-	pretty(range(lX),n=x.splits)
  pretty.X.lab	<-	round(pretty.X.at,2)
  pretty.Y.at	<-	pretty(range(lY),n=y.splits)
  pretty.Y.lab	<-	round(pretty.Y.at,2)

  jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  graphics::filled.contour(t(z), color = jet.colors, nlevels=100,#asp = 1, #frame.plot = T, axes=F,
                 ylab = y.name, xlab = x.name,
                 plot.axes={
                   # axis(1,at=seq(0,1,length.out = x.splits+1),labels=pretty.X.lab)
                   # axis(2,at=seq(0,1,length.out = y.splits+1),labels=pretty.Y.lab)
                   graphics::axis(1,at=seq(0,1,length.out = length(pretty.X.lab) ),labels=pretty.X.lab)
                   graphics::axis(2,at=seq(0,1,length.out = length(pretty.Y.lab) ),labels=pretty.Y.lab)
                 },
                 key.title = graphics::title(main=z.name)
  )
}

#' Checks for NA, NaN or negative components in a vector (usually used for Qt)
#' @inheritParams default_params_doc
#' @param v sometjing
#' @param display_output something
negatives_correction <- function(
  v,
  pars,
  display_output = 0
){
  problems=0
  if ( any(is.na(v)) ){ problems=1;
  NA_components=which(is.na(v) & !is.nan(v) )
  NaN_components=which(is.nan(v))
  if (display_output==1){
    cat("There are non-numeric components for par values:",pars,"\n")
    if(length(NA_components)>0) {cat("NA component are:",NA_components)}
    if(length(NaN_components)>0) {cat("NaN component are:",NaN_components)}
  }
  }

  if ( any(v<0) & (problems)==0){
    negatives=v[v<0]
    maximum=max(v)
    v[v<0 & (abs(v)/abs(max(v)))<1e-10]=0
  }
  return(v)
}

#' Compares two functions that accept same arguments args=(...)
#' @inheritParams default_params_doc
#' @param function1 something 
#' @param function2 something 
#' @param output something 
#' @param ... something 
compare_functions <- function(
  function1,
  function2,
  iterations = 1,
  output = 1,
  ...
){
  Atime1=0;Atime2=0;
  for (i in 1:iterations){
    t0=proc.time()
    Atest1=( function1(...) )
    Atime1=Atime1+proc.time()-t0
    t1=proc.time()
    Atest2=( function2(...) )
    Atime2=Atime2+proc.time()-t1
  }
  if (output==1){print(list(Atest1,Atest2))}
  print(Atime1);print(Atime2);
  print( all.equal(Atest1,Atest2) )
  return(invisible(list(Atest1,Atest2)))
}

#' Checks function calls within a function
#' @inheritParams default_params_doc
called_functions <- function(
  function.name, 
  recursive = FALSE, 
  checked.functions = NULL
) {

    # Get the function's code:
    function.code <- deparse(get(function.name))

    # break code up into sections preceding left brackets:
    left.brackets <- c(unlist(strsplit(function.code, split="[[:space:]]*\\(")))

    called.functions <- unique(c(unlist(sapply(left.brackets, function(x) {
        # Split up according to anything that can't be in a function name.
        # split = not alphanumeric, not '_', and not '.'
        words <- c(unlist(strsplit(x, split="[^[:alnum:]_.]")))

        last.word <- utils::tail(words, 1)
        last.word.is.function <- tryCatch(is.function(get(last.word)),
                                      error=function(e) return(FALSE))
        return(last.word[last.word.is.function])
    }))))

    if (recursive){

        # checked.functions: We need to keep track of which functions
        # we've checked to avoid infinite loops.
        functs.to.check <- called.functions[!(called.functions %in% checked.functions)]

        testit::assert(!"Should never call listFunctions, whatever it is")
        # called.functions <- unique(c(called.functions,
        #     do.call(c, lapply(functs.to.check, function(x) {
        #         listFunctions(x, recursive = T, checked.functions = c(checked.functions, called.functions))
        #         }))))
    }
    return(called.functions)
}

#' Checks if some matrix entries are infinite, NaN, NA or (only in the lower triangle) negative
#' @inheritParams default_params_doc
#' @param Mlist something
#' @param sign_check something
matrix_check <- function(
  Mlist,
  sign_check = 0
) {
  max_k=length(Mlist)
  black_list=as.complex(rep(0,max_k))
  negative_list=rep(0,max_k)
  for (k in 1:max_k){
    matrix=as.matrix(Mlist[[k]])
    if( (max(is.infinite(matrix))>0) |
        (max(is.nan(matrix))>0)          |
        (max(is.na(matrix))>0)           )
    {black_list[k]=black_list[k]+1}
    if (max(diag(matrix))>0)
    {black_list[k]=black_list[k]+2i}
  }
  if (sum(Re(black_list))!=0)
  {print(paste("We have non-numeric values for k=",(1:length(black_list))[(Re(black_list)>0)]))}
  if (sum(Im(black_list))!=0)
  {print(paste("We have positive diagonal values for k=",(1:length(black_list))[(Im(black_list)>0)]))}
  if (sign_check==1){
    entries_sign=rep(1,max_k)
    for (k in 1:max_k){
      lower_triangle=(as.matrix(Mlist[[k]]))[lower.tri(Mlist[[k]])]
      if (min(lower_triangle)<0){entries_sign[k]=-1;negative_list[k]=-1}
    }
    if (min(entries_sign)<0){print(paste("There is a negative entry on k=",which((entries_sign)<0) ))}
  }
  return(list(black_list=black_list,negative_list=negative_list))
}

# RJCB: No, do not install missing packages here
#require(ape)
#require(abind)
#require(mvbutils)

#' Does something
#' @inheritParams default_params_doc
append_multiple <- function(
  x, 
  values, 
  after
) {
  if(!is.list(values) & length(values)!=length(after)) stop("In append_multiple: values must be a list unless it has the same length as after")

  x2 <- x
  for(i in 1:length(values)) {
    x2 <- append(x2,values[[i]],after[i]+length(x2)-length(x))
  }
  x2
}

#' Function returning branch lengths named by ending node
#' @inheritParams default_params_doc
branchLengths <- function(tr)
{
  if(length(tr$node)>0) {type='NAMED'} else{type='UNNAMED'}
  if(length(tr$root)>0) {root.length=tr$root} else{root.length=0}

  node.numbers <- sort(unique(tr$edge[,1]))
  BL1 <- as.numeric(
    stats::na.omit(
      c(
        root.length,
        tr$edge.length[match(node.numbers,tr$edge[,2])]
      )
    )
  )

  tip.numbers=1:length(tr$tip)
  BL2=tr$edge.length[match(tip.numbers,tr$edge[,2])]

  if(type=='NAMED')
  {
    res=data.frame(	Number=c(sort(node.numbers),tip.numbers),
                    Name=c(tr$node,tr$tip),
                    Lengths=c(BL1,BL2))
  }
  if(type=='UNNAMED')
  {
    res=data.frame(	Number=c(sort(node.numbers),tip.numbers),
                    Name=c(rep(NA,tr$Nnode),tr$tip),
                    Lengths=c(BL1,BL2))
  }

  return(res)
}

#' Does something
#' @inheritParams default_params_doc
summarize_beast_posterior <- function(
  INPUT.TREES,
  INPUT.XML = NULL,
  subsamp = NULL
){


  #If an XML with a starting tree is provided, get the node names from that tree
  node_names_xml=NULL
  if(!is.null(INPUT.XML)){
    subst.xml=scan(INPUT.XML,what="raw")
    where.starting.tree=which(subst.xml=="id=\"startingTree\">")+1
    start.tr <- ape::read.tree(text=subst.xml[where.starting.tree])
    node_names_xml=start.tr$node.label
  }

  branch_info_tot=tree_tot=list()
  subsamps=c()

  for(t in 1:length(INPUT.TREES)){

    #Read the trees
    subst.tr=scan(INPUT.TREES[t],what="raw")

    #Find the equal signs in tree file
    equals=(1:length(subst.tr))[subst.tr=="="]

    #Use position of equal signs to find position of trees and their names
    tree.positions=equals+2
    tree.names.pos=equals-2

    #How many trees?
    Ntrees=length(tree.positions)
    if(is.null(subsamp)) subsamp=Ntrees
    samp=sample(Ntrees,round(subsamp/length(INPUT.TREES)),replace=FALSE)

    #Get correspondance between numbers and taxa
    #How many taxa?
    split1=strsplit(subst.tr[5],"=")[[1]][2]
    ntax=as.numeric(strsplit(split1,";")[[1]])

    #Find the "Translate"
    Translate.pos=(1:length(subst.tr))[subst.tr=="Translate"]
    comma.pos=which(subst.tr==",")
    if(length(comma.pos)==0){
      numbers.pos=seq(Translate.pos+1,Translate.pos+2*ntax,2)
    } else{
      numbers.pos=seq(Translate.pos+1,Translate.pos+3*ntax,3)
    }
    taxa.pos=numbers.pos+1
    numbers=subst.tr[numbers.pos]
    taxa=array(,length(taxa.pos))
    for(i in 1:length(taxa.pos)) taxa[i]=strsplit(subst.tr[taxa.pos[i]],",")[[1]][1]


    for(i in 1:length(samp)){

      #Isolate the trees and tree names
      tree1=subst.tr[tree.positions[samp[i]]]

      #Break the character string
      tree_sub=substring(tree1,0:100000,0:100000)

      #Trim the end
      tree_sub2=tree_sub[tree_sub!=""]

      colons=which(tree_sub2==":")
      parenthesis_colons=(colons-1)[tree_sub2[colons-1]==")"]
      if(is.null(node_names_xml)){
        node_names=paste("Node",seq_along(c(NA,parenthesis_colons)), sep="")
      } else{node_names=node_names_xml}
      tree_sub2=append_multiple(tree_sub2,node_names,c(parenthesis_colons,length(tree_sub2)-1))

      #Find positions of opening and closing brackets
      croch1.pos=(1:length(tree_sub2))[tree_sub2=="["]
      croch2.pos=(1:length(tree_sub2))[tree_sub2=="]"]

      #Trim the brackets and the information in-between

      remove=c()
      for (k in 1:length(croch1.pos))
      {
        commas=which(tree_sub2[croch1.pos[k]:croch2.pos[k]]==",")
        commas=c(croch1.pos[k],commas+croch1.pos[k]-1,croch2.pos[k])
        node_name_k=tree_sub2[croch1.pos[k]-2]
        name_arg=value_arg=array(,length(commas)-1)
        for(j in 1:(length(commas)-1)){

          arg_j=tree_sub2[(commas[j]+1):(commas[j+1]-1)]
          equal=which(arg_j=="=")
          name_arg[j]=paste(arg_j[2:(equal-1)],collapse="")
          value_arg[j]=paste(arg_j[(equal+1):length(arg_j)],collapse="")

        }

        if(k==1){
          res=array(,c(length(croch1.pos),length(value_arg)+1))
          colnames(res)=c("name",name_arg)
        }

        res[k,]=c(node_name_k,value_arg)
        remove=c(remove,croch1.pos[k]:croch2.pos[k])
      }

      res[res[,"name"]%in%numbers,"name"]=taxa[match(res[res[,"name"]%in%numbers,"name"],numbers)]

      #Re-paste the new tree into a single character string
      tree.text=paste(tree_sub2[-remove],collapse="")

      tr <- ape::read.tree(text=tree.text)
      tr$tip.label=taxa[as.numeric(tr$tip.label)]

      #trees[[i]]=tr
      branch_info=branchLengths(tr)[-1,]
      res=res[match(branch_info[,"Name"],res[,"name"]),]
      branch_info=cbind(tr$edge[match(branch_info[,"Number"],tr$edge[,2]),],
                        #tr$edge.length[match(branch_info[,"Number"],tr$edge[,2])],
                        branch_info[,-1],res[,-1,drop=FALSE])
      colnames(branch_info)[1:2]=c("Ancestor","Descendant")

      if(i==1) {
        branch_infos=array(,c(dim(branch_info),length(samp)))
        dimnames(branch_infos)=list(NULL,dimnames(branch_info)[[2]],1:length(samp))
      }

      branch_infos[,,i]=as.matrix(branch_info)
      tree_tot[[length(tree_tot)+1]]=tr

    }

    branch_info_tot[[t]]=branch_infos
    subsamps[t]=length(samp)
  }


  branch_info_out <- abind::abind(branch_info_tot, along = 3)

  return(list(branch_info=branch_infos, trees=tree_tot))
}

#' Does something
#' @inheritParams default_params_doc
#' @param file_name file name
#' @param maxtree something
extract_posterior <- function(
  file_name = "simcophylo_1_RUN1.(time).trees.txt",
  maxtree = 10
) {

  names_i=file_name

  beast_posterior <- summarize_beast_posterior(
      INPUT.TREES = names_i, 
      subsamp=maxtree #Number of samples to keep
    )

  # beast_posterior$trees # a list of the trees in phylo format
  # beast_posterior$branch_info # a list of branch lengths and associated information
  brts=list()
  for (i in 1:maxtree){
    brts[[i]]=unname(ape::branching.times(beast_posterior$trees[[i]]) )
  }

  dist=vector("list")
  maxindex=max(unlist(lapply(brts,length)))

  for (i in 1:maxindex){
  dist[[i]]=unlist(lapply(brts,"[[",i))
  }

  return(list(brts=brts,dist=dist))
}

# # Optimizers taken from DDD #I commented them because I wanna know if i can use them from DDD directly since they are the same
# simplex = function(fun,trparsopt,optimpars,...)
# {
#   numpar = length(trparsopt)
#   reltolx = optimpars[1]
#   reltolf = optimpars[2]
#   abstolx = optimpars[3]
#   maxiter = optimpars[4]
#
#   ## Setting up initial simplex
#   v = t(matrix(rep(trparsopt,each = numpar + 1),nrow = numpar + 1))
#   for(i in 1:numpar)
#   {
#     parsoptff = 1.05 * trparsopt[i]/(1 - trparsopt[i])
#     trparsoptff = parsoptff/(1 + parsoptff)
#     fac = trparsoptff/trparsopt[i]
#     if(v[i,i + 1] == 0)
#     {
#       v[i,i + 1] = 0.00025
#     } else {
#       v[i,i + 1] = v[i,i + 1] * min(1.05,fac)
#     }
#   }
#
#   fv = rep(0,numpar + 1)
#   for(i in 1:(numpar + 1))
#   {
#     fv[i] = -fun(trparsopt = v[,i], ...)
#   }
#
#   how = "initial"
#   itercount = 1
#   string = itercount
#   for(i in 1:numpar)
#   {
#     string = paste(string, v[i,1]/(1 - v[i,1]), sep = " ")
#   }
#   string = paste(string, -fv[1], how, "\n", sep = " ")
#   cat(string)
#   flush.console()
#
#   tmp = order(fv)
#   if(numpar == 1)
#   {
#     v = matrix(v[tmp],nrow = 1,ncol = 2)
#   } else {
#     v = v[,tmp]
#   }
#   fv = fv[tmp]
#
#   ## Iterate until stopping criterion is reached
#   rh = 1
#   ch = 2
#   ps = 0.5
#   si = 0.5
#
#   v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))
#
#   while(itercount <= maxiter & ( ( is.nan(max(abs(fv - fv[1]))) | (max(abs(fv - fv[1])) - reltolf * abs(fv[1]) > 0) ) + ( (max(abs(v - v2) - reltolx * abs(v2)) > 0) | (max(abs(v - v2)) - abstolx > 0) ) ) )
#   {
#     ## Calculate reflection point
#
#     if(numpar == 1)
#     {
#       xbar = v[1]
#     } else {
#       xbar = rowSums(v[,1:numpar])/numpar
#     }
#     xr = (1 + rh) * xbar - rh * v[,numpar + 1]
#     fxr = -fun(trparsopt = xr, ...)
#
#     if(fxr < fv[1])
#     {
#       ## Calculate expansion point
#       xe = (1 + rh * ch) * xbar - rh * ch * v[,numpar + 1]
#       fxe = -fun(trparsopt = xe, ...)
#       if(fxe < fxr)
#       {
#         v[,numpar + 1] = xe
#         fv[numpar + 1] = fxe
#         how = "expand"
#       } else {
#         v[,numpar + 1] = xr
#         fv[numpar + 1] = fxr
#         how = "reflect"
#       }
#     } else {
#       if(fxr < fv[numpar])
#       {
#         v[,numpar + 1] = xr
#         fv[numpar + 1] = fxr
#         how = "reflect"
#       } else {
#         if(fxr < fv[numpar + 1])
#         {
#           ## Calculate outside contraction point
#           xco = (1 + ps * rh) * xbar - ps * rh * v[,numpar + 1]
#           fxco = -fun(trparsopt = xco, ...)
#           if(fxco <= fxr)
#           {
#             v[,numpar + 1] = xco
#             fv[numpar + 1] = fxco
#             how = "contract outside"
#           } else {
#             how = "shrink"
#           }
#         } else {
#           ## Calculate inside contraction point
#           xci = (1 - ps) * xbar + ps * v[,numpar + 1]
#           fxci = -fun(trparsopt = xci, ...)
#           if(fxci < fv[numpar + 1])
#           {
#             v[,numpar + 1] = xci
#             fv[numpar + 1] = fxci
#             how = "contract inside"
#           } else {
#             how = "shrink"
#           }
#         }
#         if(how == "shrink")
#         {
#           for(j in 2:(numpar + 1))
#           {
#
#             v[,j] = v[,1] + si * (v[,j] - v[,1])
#             fv[j] = -fun(trparsopt = v[,j], ...)
#           }
#         }
#       }
#     }
#     tmp = order(fv)
#     if(numpar == 1)
#     {
#       v = matrix(v[tmp],nrow = 1,ncol = 2)
#     } else {
#       v = v[,tmp]
#     }
#     fv = fv[tmp]
#     itercount = itercount + 1
#     string = itercount;
#     for(i in 1:numpar)
#     {
#       string = paste(string, v[i,1]/(1 - v[i,1]), sep = " ")
#     }
#     string = paste(string, -fv[1], how, "\n", sep = " ")
#     cat(string)
#     flush.console()
#     v2 = t(matrix(rep(v[,1],each = numpar + 1),nrow = numpar + 1))
#   }
#   if(itercount < maxiter)
#   {
#     cat("Optimization has terminated successfully.","\n")
#   } else {
#     cat("Maximum number of iterations has been exceeded.","\n")
#   }
#   out = list(par = v[,1], fvalues = -fv[1], conv = as.numeric(itercount > maxiter))
#   invisible(out)
# }
#
# optimizer = function(optimmethod = 'simplex',optimpars = c(1E-4,1E-4,1E-6,1000),fun,trparsopt, ...)
# {
#   if(optimmethod == 'simplex')
#   {
#     out = simplex(fun = fun,trparsopt = trparsopt,optimpars = optimpars,...)
#   }
#   if(optimmethod == 'subplex')
#   {
#     minfun = function(fun,trparsopt,...)
#     {
#       return(-fun(trparsopt = trparsopt,...))
#     }
#     out = subplex::subplex(par = trparsopt,fn = minfun,control = list(abstol = optimpars[3],reltol = optimpars[1],maxit = optimpars[4]),fun = fun,...)
#     out = list(par = out$par, fvalues = -out$value, conv = out$convergence)
#   }
#   return(out)
# }

###################### STUFF THAT I DON'T REALLY USE ######################
#HG shit
# hyperA_HannoX <- function(N,k,q) {# fast O(N)
#   j <- 0:k
#   A1 <- (1-q)^(2*k) * choose(k, j) * ((2-q)/(1-q))^j
#   # A1 <- (1-q)^(2*k) * choose(k, j) * (2)^j
#   A <- diag(A1[1], nrow=N+1, ncol=N+1)
#   A[1:(k+1),1] <- A1
#   for (dst in 2:N) {
#     src <- dst - 1
#     A[(dst+1):(N+1),dst] <- A[src:(N-1),src] + A[(src+1):N,src]
#   }
#   A * outer(0:N,0:N,function(m,n) { (q^(m-n))*(1-q)^(2*n-m) })
# }
#
# hyperA_HannoX2 <- function(N,k,q) {# fast O(N)
#   j <- 0:k
#   # A1 <- (1-q)^(2*k) * choose(k, j) * ((2-q)/(1-q))^j
#   A1 <- (1-q)^(2*k) * choose(k, j) * (2)^j
#   N <- N + 1
#   A <- diag(A1[1], nrow=N+2, ncol=N+2)
#   A[1:(k+1),1] <- A1
#   for (dst in 2:N) {
#     src <- dst - 1
#     s <- src:min(N,2*src+k-1)
#     A[s+2,dst] <- A[s,src] + A[s+1,src]
#     m <- s-1; n <- src-1;
#     A[s,src] <- A[s,src] * q^(m-n)*(1-q)^(2*n-m)
#   }
#   A[N,N] = A[N,N] * (1-q)^(N-1);
#   A[1:N,1:N]
# }
#
# hyperA_HannoX3 <- function(N,k,q) {# fast O(N)
#   j <- 0:k
#   # A1 <- (1-q)^(2*k) * choose(k, j) * ((2-q)/(1-q))^j
#   A1 <- (1-q)^(k) * choose(k, j) * (2)^j
#   N <- N + 1
#   A <- diag(A1[1], nrow=N+2, ncol=N+2)
#   A[1:(k+1),1] <- A1
#   for (dst in 2:N) {
#     src <- dst - 1
#     s <- src:min(N,2*src+k-1)
#     A[s+2,dst] <- A[s,src] + A[s+1,src]
#     m <- s-1; n <- src-1;
#     A[s,src] <- A[s,src] * q^(m-n)*(1-q)^(2*n-m)
#   }
#   A[N,N] = A[N,N] * (1-q)^(N-1);
#   A[1:N,1:N]
# }
#
# hyperB_HannoX <- function(N,k,q,b) {# fast O(N)
#   k=k-b
#   j <- 0:k
#   A1 <- choose(k, j) * ((2-q)/(1-q))^j
#   A <- diag(A1[1], nrow=N+1, ncol=N+1)
#   A[1:(k+1),1] <- A1
#   for (dst in 2:N) {
#     src <- dst - 1
#     A[(dst+1):(N+1),dst] <- A[src:(N-1),src] + A[(src+1):N,src]
#   }
#   A * outer(0:N,0:N,function(m,n) { (q^(m-n))*(1-q)^(2*n-m) })
# }
#
# hyperA_Hanno3=function(N,k,q){#faster version, performs the same calculations
#   nvec=0:N
#   A=matrix(0,nrow=length(nvec),ncol=length(nvec))
#   z=(2-q)/(1-q)
#   for (m in nvec){
#     for (n in 0:m){
#       out=0
#       j = 0:k
#       out = sum( choose(k,j)*choose(n,(m-n-j))*z^j )
#       A[m+1,n+1]=out*(q^(m-n))*( (1-q)^(2*n-m) )
#     }
#   }
#   A=A*(1-q)^(2*k) #this has to be applied later
#   return(A)
# }
#
# hyperA_Hanno5=function(N,k,q){#faster version, revised boundaries using binomial properties
#   nvec=0:N
#   A=matrix(0,nrow=length(nvec),ncol=length(nvec))
#   z=(2-q)/(1-q)
#   for (n in nvec){
#     upper_limit=min((2*n+k),N)
#     for (m in n:upper_limit){
#       out=0
#       j = 0:k
#       out = sum( choose(k,j)*choose(n,(m-n-j))*z^j )
#       A[m+1,n+1]=out*(q^(m-n))*(1-q)^(2*n-m)
#     }
#   }
#   A=A*(1-q)^(2*k)
#   return(A)
# }
#
# hyperA_Hanno6=function(N,k,q){#slower but more stable
#   nvec=0:N
#   A=matrix(0,nrow=length(nvec),ncol=length(nvec))
#   z=(2-q)/(1-q)
#   l0=600
#   for (n in nvec){
#     upper_limit=min((2*n+k),N)
#     for (m in n:upper_limit){
#       out=0
#       j = 0:k
#       out = l0 + log(sum(exp(lchoose(k,j)+lchoose(n,(m-n-j))+j*log(z)-l0)))
#       out = out + (m-n)*log(q) + (2*n-m)*log(1-q)
#       A[m+1,n+1]=exp(out)
#     }
#   }
#   A=A*(1-q)^(2*k)
#   return(A)
# }
#
# hyperA_Hanno7=function(N,k,q){#slower but more stable
#   nvec=0:N
#   A=matrix(-Inf,nrow=length(nvec),ncol=length(nvec))
#   z=(2-q)/(1-q)
#   l0=600
#   for (n in nvec){
#     upper_limit=min((2*n+k),N)
#     for (m in n:upper_limit){
#       out=0
#       j = 0:k
#       out = l0 + log(sum(exp(lchoose(k,j)+lchoose(n,(m-n-j))+j*log(z)-l0)))
#       out = out + (m-n)*log(q) + (2*n-m)*log(1-q)
#       A[m+1,n+1]=out
#     }
#   }
#   A=exp(A)
#   A=A*(1-q)^(2*k)
#   return(A)
# }
#
# hyperB_Hanno7=function(N,k,q,b){#slower but more stable
#   nvec=0:N
#   k=k-b
#   B=matrix(-Inf,nrow=length(nvec),ncol=length(nvec))
#   z=(2-q)/(1-q)
#   l0=600
#   for (n in nvec){
#     upper_limit=min((2*n+k),N)
#     for (m in n:upper_limit){
#       out=0
#       j = 0:k
#       out = l0 + log(sum(exp(lchoose(k,j)+lchoose(n,(m-n-j))+j*log(z)-l0)))
#       out = out + (m-n)*log(q) + (2*n-m)*log(1-q)
#       B[m+1,n+1]=out
#     }
#   }
#   B=exp(B)
#   # B=A*(1-q)^(2*k)
#   return(B)
# }
#
# hyperB_Hanno3=function(N,k,q,b){#faster version, performs the same calculations
#   nvec=0:N
#   k=k-b
#   B=matrix(0,nrow=length(nvec),ncol=length(nvec))
#   z=(2-q)/(1-q)
#   for (m in nvec){
#     for (n in 0:m){
#       out=0
#       j = 0:k
#       out = sum( choose(k,j)*choose(n,(m-n-j))*z^j )
#       B[m+1,n+1]=out*(q^(m-n))*(1-q)^(2*n-m)
#     }
#   }
#   # A=A*(1-q)^(2*k) #this has to be applied later
#   return(B)
# }

# hyperA_Hanno4=function(N,k,q){#faster version, performs the same calculations
#   nvec=0:N
#   A=matrix(0,nrow=length(nvec),ncol=length(nvec))
#   z=(2-q)/(1-q)
#   for (m in nvec){
#     for (n in 0:m){
#       out=0
#       j = 0:k
#       out = sum( choose(k,j)*choose(n,(m-n-j))*z^j )
#       A[m+1,n+1]=out*(q^(m-n))*(1-q)^(2*n-m)
#     }
#   }
#   # A=A*(1-q)^(2*k) #this has to be applied later
#   return(A)
# }

# hyperA <- function(m,n,q,k) {#works almost really well #, tol = .Machine$double.eps
#   z=(2-q)/(1-q)
#   j <- (m-2*n)*((m-2*n)>0)
#   A=(2*n-m+j); B=(m-n-j); C=(j); D=(k-j)
#   term_j0=j*log(z)+lchoose(n,(m-n-j))+lchoose(k,j)+(m-n)*log(q)+(2*k+2*n-m)*log(1-q)
#   term_j0=exp(term_j0)
#   term <- rep(0,length(m))
#   term[ ] <- term_j0
#   term <- term_j0
#   partial.sum <- term
#
#   while( term > 0 && (j <- j+1) <= (k) ){ #max(abs(term)) > tol
#     A <- A+1
#     C <- C+1
#     term <- term * z * (B * D)/(A * C)
#     B <- B-1
#     D <- D-1
#     partial.sum <- partial.sum + term
#   }
#   out=partial.sum
#   return(out)
# }
#
# # nvec=0:2000;k=90;q=0.2
# hyperA_Hanno=function(N,k,q){
#   nvec=0:N
#   A=matrix(0,nrow=length(nvec),ncol=length(nvec))
#   z=(2-q)/(1-q)
#   for (m in nvec){
#     for (n in 0:m){
#       out=0
#       for (j in 0:min((m-n),k)){
#         out = out + choose(k,j)*choose(n,(m-n-j))*z^j
#       }
#       A[m+1,n+1]=out*(q^(m-n))*(1-q)^(2*n-m)
#     }
#   }
#   A=A*(1-q)^(2*k)
#   return(A)
# }
#
# hyperA_Hanno2=function(N,k,q){
#   nvec=0:N
#   A=matrix(0,nrow=length(nvec),ncol=length(nvec))
#   z=(2-q)/(1-q)
#   for (m in nvec){
#     for (n in 0:m){
#       out=0
#       j = 0:min((m-n),k)
#       out = sum( choose(k,j)*choose(n,(m-n-j))*z^j )
#       A[m+1,n+1]=out*(q^(m-n))*(1-q)^(2*n-m)
#     }
#   }
#   A=A*(1-q)^(2*k)
#   return(A)
# }

# hyperB <- function(m,n,q,k,b) {#works really well. you need to add binomial(k,b) at the end
#   z=(2-q)/(1-q)
#   j <- (m-2*n)*((m-2*n)>0)
#   A=(2*n-m+j); B=(m-n-j); C=(j); D=(k-b-j)
#   term_j0=j*log(z)+lchoose(n,(m-n-j))+lchoose((k-b),j)+(m-n)*log(q)+(2*n-m)*log(1-q)
#   # term_j0=term_j0+(2*k-2*b)*log(1-q) #this is include at nodes in the main program
#   term_j0=exp(term_j0)
#   term <- rep(0,length(m))
#   term[ ] <- term_j0
#   partial.sum <- term
#
#   while( term > 0 && (j <- j+1) <= (k-b) ){
#     A <- A+1
#     C <- C+1
#     term <- term * z * (B * D)/(A * C)
#     B <- B-1
#     D <- D-1
#     partial.sum <- partial.sum + term
#   }
#   out=partial.sum#*choose(k,b)*q^b*(1-q)^(k-b) #this is already in the loglik script because it is faster
#   return(out)
# }
# hyperAredux <- function(m,n,q,k,l0) {#works almost really well #, tol = .Machine$double.eps
#   # l0=700
#   z=(2-q)/(1-q)
#   j <- (m-2*n)*((m-2*n)>0)
#   A=(2*n-m+j); B=(m-n-j); C=(j); D=(k-j)
#   term_j0=j*log(z)+lchoose(n,(m-n-j))+lchoose(k,j)+(m-n)*log(q)+(2*k+2*n-m)*log(1-q) + l0
#   term_j0=exp(term_j0)
#   term <- rep(0,length(m))
#   term[ ] <- term_j0
#   partial.sum <- term
#
#   while( term > 0 && (j <- j+1) <= (k) ){ #max(abs(term)) > tol
#     A <- A+1
#     C <- C+1
#     term <- term * z * (B * D)/(A * C)
#     B <- B-1
#     D <- D-1
#     partial.sum <- partial.sum + term
#   }
#   out=partial.sum
#   return(out)
# }
#
# hyperA2 <- function(m,n,q,k,l0=-710) {#works "almost" really well #, tol = .Machine$double.eps
#   z=(2-q)/(1-q)
#   j2 <- (m-2*n)*((m-2*n)>0)
#   A2=(2*n-m+j2); B2=(m-n-j2); C2=(j2); D2=(k-j2) #B and D must be non negative!
#   term_j0=j2*log(z)+lchoose(n,(m-n-j2))+lchoose(k,j2)+(m-n)*log(q)+(2*k+2*n-m)*log(1-q)
#   lj2 <- rep(0,length(m))
#   lj2[ ] <- term_j0
#   # l0 = 0
#   partial.sum2 <- exp( lj2 - l0 )
#
#   cond = ( term_j0==-Inf )
#   partial.sum = partial.sum2[!cond]
#   lj = lj2[!cond] - l0
#   A = A2[!cond]
#   B = B2[!cond]
#   C = C2[!cond]
#   D = D2[!cond]
#   j = j2[!cond]
#
#   while( lj > -Inf && (j <- j+1) <= (k) && length(j)>0 ){ #max(abs(term)) > tol
#     A <- A+1
#     C <- C+1
#     lj <- lj + (log(z) + log(B) + log(D) - log(A) - log(C) )
#     B <- B-1*(B!=0)
#     D <- D-1*(D!=0)
#     partial.sum <- partial.sum + exp( lj  )
#     # print(c(A,B,C,D))
#   }
#   partial.sum2[cond]=0
#   partial.sum2[!cond]=partial.sum
#   out=exp(l0) * partial.sum2;  return(out)}
#
# hyperAbrob <- function(m,n,q,k) {#works "almost" really well #, tol = .Machine$double.eps
#
#   # m=as.brob(m);n=as.brob(n);
#   z=(2-q)/(1-q)
#   j2 <- (m-2*n)*((m-2*n)>0)
#   A2=(2*n-m+j2); B2=(m-n-j2); C2=(j2); D2=(k-j2) #B and D must be non negative!
#   term_j0=j2*log(z)+lchoose(n,(m-n-j2))+lchoose(k,j2)+(m-n)*log(q)+(2*k+2*n-m)*log(1-q)
#   lj2 <- rep(0,length(m))
#   lj2[ ] <- term_j0
#
#   partial.sum2 <- brob( lj2 )
#
#   cond = ( term_j0==-Inf )
#   partial.sum = partial.sum2[!cond]
#   lj = lj2[!cond]
#   A = A2[!cond]
#   B = B2[!cond]
#   C = C2[!cond]
#   D = D2[!cond]
#   j = j2[!cond]
#
#   while( lj > -Inf && (j <- j+1) <= (k) && length(j)>0 ){ #max(abs(term)) > tol
#     A <- A+1
#     C <- C+1
#     lj <- lj + (log(z) + log(B) + log(D) - log(A) - log(C) )
#     B <- B-1*(B!=0)
#     D <- D-1*(D!=0)
#     partial.sum <- partial.sum + brob( lj  )
#     # print(c(A,B,C,D))
#   }
#   partial.sum2[cond]=0
#   partial.sum2[!cond]=partial.sum
#   out= (partial.sum2);  return(out)}

# hyperA3 <- function(m,n,q,k) {#works almost really well #, tol = .Machine$double.eps
#   z=(2-q)/(1-q)
#   j <- (m-2*n)*((m-2*n)>0)
#   A=(2*n-m+j); B=(m-n-j); C=(j); D=(k-j)
#   term_j0=j*log(z)+lchoose(n,(m-n-j))+lchoose(k,j)+(m-n)*log(q)+(2*k+2*n-m)*log(1-q)
#   term_j0=brob(term_j0)
#   # term <- rep(0,length(m))
#   term <- term_j0
#   partial.sum <- term_j0
#
#   # term=as.brob(term)
#   # partial.sum=as.brob(partial.sum)
#
#   while( term > 0 && (j <- j+1) <= (k) ){ #max(abs(term)) > tol
#     A <- A+1
#     C <- C+1
#     term <- term * z * (B * D)/(A * C)
#     B <- B-1
#     D <- D-1
#     partial.sum <- partial.sum + term
#   }
#   out=partial.sum
#   return(out)
# }
# pochhammer <- function(a, n, log)
# {
#   a.len <- length(a)
#   n.len <- length(n)
#
#   a.mat <- matrix(a, a.len, n.len)
#   n.mat <- matrix(n, a.len, n.len, byrow = TRUE)
#   out.log <- lgamma(a.mat+n.mat) - lgamma(a.mat)
#
#   if(log)
#   {
#     out <- out.log
#   }
#   else
#   {
#     out <- exp(out.log)
#   }
#
#   return(out)
# }
# ghypergeo <- function(a, b, z, k = 10)
# {
#   ## needs pochhammer function that you can find on top of this
#   ## evaluates pFq Hypergeometric Function
#   ## a "n-by-p matrix"
#   ## b "n-by-q matrix"
#   ## z "positive vector of length one or n"
#   ## k "positive integer" Indicates how long of the series should be used.
#
#   if(!is.matrix(a)|| !is.matrix(b) ||
#      dim(a)[1] != dim(b)[1])
#   {
#     stop("The arguments `a' and `b' must be matrices \n and should have the same row dimension.")
#   }
#
#   a.dim <- dim(a)
#   nObs <- a.dim[1]
#   a.nCol <-a.dim[2]
#   b.dim <- dim(b)
#   b.nCol <- b.dim[2]
#
#   n <- k+1
#   n.series <- 0:k
#
#   zpower.series <- matrix(n.series, nObs, n, byrow = TRUE)*
#     matrix(log(z), nObs, n)
#
#   nfact.series <- lgamma(n.series+1) # log(n!)
#   an.series <- pochhammer(t(a), n.series, log = TRUE)#
#   bn.series <- pochhammer(t(b), n.series, log = TRUE)#
#
#   an.ary <- array(t(an.series), c(n, a.nCol, nObs))
#   bn.ary <- array(t(bn.series), c(n, b.nCol, nObs))
#
#   an <- apply(an.ary, c(1, 3), sum)
#   bn <- apply(bn.ary, c(1, 3), sum)
#
#   nlog <- an - bn + t(zpower.series) - nfact.series
#   out <- matrix(apply(exp(nlog), 2, sum))
#
#   return(out)
# }
#
# regghypergeo <- function(a, b, z, k = 10)
# {
#   ## the regularized hypergeometric function pFq(a, b, z)/Prod[Gamma(b)]
#   gammaprod <- apply(gamma(b), 1, prod)
#   out <- ghypergeo(a = a, b = b, z = z, k = k)/gammaprod
#   return(out)
# }
#
#
# hyper <- function(A, B, C, z, tol = sqrt(.Machine$double.eps),maxit= 10000) {
#   term <- z
#   term[ ] <- 1
#   partial.sum <- term
#   n <- 0
#   while(max(abs(term)) > tol && (n <- n+1) < maxit) {
#     term <- term*A*B/C * z/n
#     partial.sum <- partial.sum + term
#     A <- A+1
#     B <- B+1
#     C <- C+1
#   }
#   if(n == maxit) warning("convergence may not have been achieved")
#   partial.sum
# } #hyper(2-4,-10,1+4+4,1.9/0.9);hyper(2:3-4:5,-10,1+4:5+2*(2:3),1.9/0.9)
#
# hyper2 <- function(A, B, C, z, tol = sqrt(.Machine$double.eps),maxit= 10000) {
#   # if (C<=0){return(0)}
#   # ans=sign(C)
#   # out[ans<0]=0
#   # if (min(ans<0)==1){return(out)}
#   term <- z
#   term[ ] <- 1
#   partial.sum <- term
#   n <- 0
#   A1=A;B1=B;C1=C
#   while(max(abs(term)) > tol && (n <- n+1) <= (maxit)) {
#     term <- term*A1*B1/C1 * z/n
#     partial.sum <- partial.sum + term
#     A1 <- A1+1
#     B1 <- B1+1
#     C1 <- C1+1
#   }
#   if(n == maxit) warning("convergence may not have been achieved")
#   out=partial.sum
#   # out[ans>=0]=partial.sum
#   return(out)
#   # partial.sum
# }
#CHOOSE FUNCTION FOR BIG NUMBERS (prototype)
# bchoose=function(m,n){
#   #remember to get as.numeric(out)
#   out=brob(lchoose(m,n))
#   return(out)
# }

# CORRELATION ANALYSIS EXTRA FUNCTIONS
# #data with no outliers
# fun_to_min = function(x,confidence_interval=0.95){
#   if (x<=0 || x>=1){return(Inf)}# x in [0,1]
#   cresults=results[,1:3]
#   cpars_interval=list(c(0,quantile(cresults[,1],x)),c(0,quantile(cresults[,2],x)),c(0,quantile(cresults[,3],x)))
#   goodrows=(cresults[,1]<=cpars_interval[[1]][2])&(cresults[,2]<=cpars_interval[[2]][2])&(cresults[,3]<=cpars_interval[[3]][2])
#   clean_data=results[goodrows,1:3]
#   out=abs(dim(clean_data)[1]/dim(results)[1]-confidence_interval)
# }
# exclude_outliers=function(results,confidence_interval=0.95){
#   out_quantile=subplex::subplex(par = 0.95,fn = fun_to_min,confidence_interval=confidence_interval)$par
#   pars_interval=list(c(0,quantile(results[,1],out_quantile)),c(0,quantile(results[,2],out_quantile)),c(0,quantile(results[,3],out_quantile)))
#   goodrows=(results[,1]<=pars_interval[[1]][2])&(results[,2]<=pars_interval[[2]][2])&(results[,3]<=pars_interval[[3]][2])
#   out_clean_data=results[goodrows,]
#   return( list(out_quantile=out_quantile,out_clean_data=out_clean_data) )
# }
# confidence_interval = 0.9
# no_outliers_results = exclude_outliers(results,confidence_interval = confidence_interval)$out_clean_data
# no_outliers_data_title = paste("\n CORRELATION ANALYSIS (First ",confidence_interval*100,"% of the distribution)\n", N ," trees (",max_sims-N," failed)",sep = '')
# no_outliers_data_pdfname = paste("correlation_analysis_first_",confidence_interval*100,"_percent_of_the_distribution",sep = '')
# percentiles_function(results = no_outliers_results,sim_pars = sim_pars)
# correlation_analysis(results = no_outliers_results,sim_pars = sim_pars,titolo=no_outliers_data_title,pdfname=no_outliers_data_pdfname,path=path)

#poster

# #poster function
# poster_analysis=function(results,path,titolo,pdfname,sim_pars=sim_pars,ranges){
#
#   Npars=length(sim_pars);#pars_interval=list(c(0,quantile(results[,1],.95)),c(0,quantile(results[,2],.95)),c(0,quantile(results[,3],.95)))
#   if(Npars==4){par_names = c(expression(lambda),expression(mu),expression(nu),expression(q))}
#   if(Npars==3){par_names = c(expression(lambda),expression(mu),expression(q))}
#
#   truevalues_color="coral"; medians_color = "darkviolet"; points_color = "azure4";
#   truevalues_color_name = "violet"; medians_color_name = "orange";
#
#   medians=rep(0,Npars);for (idpar in 1:Npars){medians[idpar]=median(results[,idpar])}
#   medians_string=paste0( "MLE Medians (",medians_color_name,") = (",paste(signif(medians,2),sep = "''",collapse = ", "),")")
#   truevalues_string=paste0( "True Values (",truevalues_color_name,") = (",paste(signif(sim_pars,2),sep = "''",collapse = ", "),")")
#
#   printa.la.figura=function(results,par_names,ranges,titolo,medians_string,truevalues_string,medians_color,truevalues_color){
#     par(mfrow=c(2,2));par(oma=c(0,0,2,0));
#     line_thickness=5
#
#     # for (i in 1:Npars)
#     # { graphics::hist((results[,i]),main=NULL,xlab = par_names[i], breaks = 15,xlim = c(0,ranges[i])); #,breaks = 15
#     #   abline(v=sim_pars[i],col = truevalues_color,lwd=line_thickness)
#     #   abline(v=medians[i],col = medians_color,lwd=line_thickness)
#     # }
#     # title(main=paste("\n",titolo,"\n     ",medians_string,"\n",truevalues_string,sep = ''),outer=T)
#     for (i in 1:Npars)
#     { graphics::hist((results[,i]),main=NULL,xlab = NULL, breaks = 15,xlim = c(0,ranges[i])); #,breaks = 15
#       graphics::abline(v=sim_pars[i],col = truevalues_color,lwd=line_thickness)
#       graphics::abline(v=medians[i],col = medians_color,lwd=line_thickness)
#     }
#   }
#
#   #pdf
#   pdf(file = paste(path,"/poster_",pdfname,".pdf",sep=''));#plot.new();
#   printa.la.figura(results=results,par_names=par_names,ranges=ranges,titolo=titolo,
#                    medians_string=medians_string,truevalues_string=truevalues_string,
#                    medians_color=medians_color,truevalues_color=truevalues_color)
#   grDevices::dev.off()
#
#   #high_resolution_image for ppt
#   tiff(file = paste(path,"/poster_",pdfname,".tiff",sep=''), height = 17, width = 17, units = 'cm',
#        compression = "lzw", res = 300)
#   printa.la.figura(results=results,par_names=par_names,ranges=ranges,titolo=titolo,
#                    medians_string=medians_string,truevalues_string=truevalues_string,
#                    medians_color=medians_color,truevalues_color=truevalues_color)
#   grDevices::dev.off()
#
#   #png
#   png(file = paste(path,"/poster_",pdfname,".png",sep=''));#plot.new();
#   printa.la.figura(results=results,par_names=par_names,ranges=ranges,titolo=titolo,
#                    medians_string=medians_string,truevalues_string=truevalues_string,
#                    medians_color=medians_color,truevalues_color=truevalues_color)
#   grDevices::dev.off()
# }

# #poster output
# poster_path = "F:/Dropbox/University/Progress/RQ1 - Multiple Births/Poster/figures"
# toobig = which ( ( mbd_res[,colnames(mbd_res)=="nu"] )>5 )
# mbd_res_cleaned = mbd_res[-toobig,]
# percentiles_function(results = mbd_res_cleaned,sim_pars = sim_pars)
# correlation_analysis(results = mbd_res_cleaned,path=path,titolo = "At least one multiple birth - cleaned dataset",pdfname = "correlation_analysis_at_least_one_multiple_birth_cleaned_dataset",sim_pars = sim_pars)
# poster_analysis(results = mbd_res_cleaned,path=poster_path,titolo = "\n 200 simulated trees having at least one multiple birth event",pdfname = "parameter_estimates",sim_pars = sim_pars,ranges=c(1,1,5,1))
# # poster_analysis(results = mbd_res_cleaned,path=poster_path,titolo = "",pdfname = "parameter_estimates",sim_pars = sim_pars,ranges=c(1,1,5,1))
# # dim(mbd_res_cleaned)
# #
# # i=1; par(mfrow=c(1,1));ranges=c(1,1,5,1);
# # graphics::hist((mbd_res_cleaned[,i]),main=NULL,xlab = pippo[i],xlim = c(0,ranges[i]));i=i+1
# # graphics::hist((mbd_res_cleaned[,i]),main=NULL,xlab = pippo[i],xlim = c(0,ranges[i]));i=i+1
# # graphics::hist((mbd_res_cleaned[,i]),main=NULL,xlab = pippo[i],xlim = c(0,ranges[i]));i=i+1
# # graphics::hist((mbd_res_cleaned[,i]),main=NULL,xlab = pippo[i],xlim = c(0,ranges[i]));i=i+1
# #
# # pippo=c(expression(lambda),expression(mu),expression(nu),expression(q))
