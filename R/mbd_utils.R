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
  n_pars <- length(sim_pars);#pars_interval=list(c(0, quantile(results[, 1],.95)), c(0, quantile(results[, 2],.95)), c(0, quantile(results[, 3],.95)))
  if (missing(idparsopt)){estimated_pars=1:n_pars}else{estimated_pars=idparsopt}
  par_names <- colnames(results)[1:n_pars]

  truevalues_color <- "red"; points_color = "azure4"; medians_color = "blue3"#"chartreuse3";
  medians_color_name <- "blue"; truevalues_color_name = "red";

  medians <- rep(0, n_pars)
  for (idpar in 1:n_pars) { 
    medians[idpar] <- stats::median(results[, idpar])
  }
  medians_string <- paste0( "MLE Medians (", medians_color_name, ") = (", paste(signif(medians, 2), sep = "''", collapse = ", "), ")")
  truevalues_string <- paste0( "True Values (", truevalues_color_name, ") = (", paste(signif(sim_pars, 2), sep = "''", collapse = ", "), ")")
  axislimits <- rep(NA, n_pars)
  for (i in 1:n_pars){
    axislimits[i] <- stats::quantile(
      results[, i],
      probs <- 1 - percentage_hidden_outliers
    )
  }

  #pdf
  grDevices::pdf(file = paste(path, "/", pdfname, ".pdf", sep=''));#plot.new();
  graphics::par(mfrow = c(length(estimated_pars), length(estimated_pars)))
  graphics::par(oma = c(0, 0, 2, 0));
  for (i in estimated_pars) { for (j in estimated_pars) {
    good.lines = results[, i] < axislimits[i] & results[, j]<axislimits[j]
    ifelse(any(good.lines)>0, good.results <- results[good.lines, ], good.results <- results)

    if ( i== j){graphics::hist((good.results[, i]), main=NULL, xlab = paste(par_names[i]), breaks = 15)
      graphics::abline(v = sim_pars[i], col = truevalues_color)
      graphics::abline(v = medians[i], col = medians_color)
    }
    else{graphics::plot(good.results[, i]~good.results[, j], xlab = par_names[j], ylab = par_names[i], cex=0.3, col=points_color);
      graphics::points(x = sim_pars[j], y = sim_pars[i], col = truevalues_color, pch = 10, cex = 1.5)
      graphics::points(x = medians[j], y = medians[i], col = medians_color, pch = 10, cex = 1.5)
    }
  }}
  graphics::title(main=(titolo.pdf<-(paste("\n\n", titolo, "\n", medians_string, "\n", truevalues_string, sep = ''))), outer=T)
  grDevices::dev.off()

  if (openit==1){
    file.show(normalizePath(paste(path, "/", pdfname, ".pdf", sep='')))
  }

  if (!missing(mother_folder))
    {
    grDevices::pdf(file = paste(mother_folder, "/", pdfname, ".pdf", sep=''));#plot.new();
    graphics::par(mfrow=c(length(estimated_pars), length(estimated_pars)))
    graphics::par(oma=c(0, 0, 2, 0));
    for (i in estimated_pars){for (j in estimated_pars){
      good.lines = results[, i]<axislimits[i] & results[, j]<axislimits[j]
      ifelse(any(good.lines)>0, good.results <- results[good.lines, ], good.results <- results)

      if (i==j){graphics::hist((good.results[, i]), main=NULL, xlab = paste(par_names[i]), breaks = 15)
        graphics::abline(v=sim_pars[i], col = truevalues_color)
        graphics::abline(v=medians[i], col = medians_color)
      }
      else{graphics::plot(good.results[, i]~good.results[, j], xlab=par_names[j], ylab=par_names[i], cex=0.3, col=points_color);
        graphics::points(x=sim_pars[j], y=sim_pars[i], col=truevalues_color, pch=10, cex=1.5)
        graphics::points(x=medians[j], y=medians[i], col=medians_color, pch=10, cex=1.5)
      }
    }}
    graphics::title(main=titolo.pdf, outer = T)
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
  quantiles_names <- format(round(quantiles_choice, 2), nsmall = 2)
  n_pars <- length(sim_pars);#pars_interval=list(c(0, stats::quantile(results[, 1],.95)), c(0, quantile(results[, 2],.95)), c(0, quantile(results[, 3],.95)))
  parnames <-  colnames(results)[1:n_pars]
  percentiles <- vector("list", 3)
  for (idpar in 1:n_pars){percentiles[[idpar]] <- stats::quantile(results[, idpar], quantiles_choice)}
  percentiles <- t(matrix(unlist(percentiles), nrow = 3, ncol = n_pars));
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
  x_name="x",
  y_name="y",
  heatmap.name="HeatMap"
){
  Matrix2=matrix(Matrix, nrow = length(x)*length(y))
  df=data.frame(expand.grid(x=x, y=y),HeatMap=Matrix2)
  #names(df$x)=names(x);names(df$y)=names(y);names(df$HeatMap)=names(Matrix);
  colnames(df)[1] <- x_name;
  colnames(df)[2] <- y_name;
  colnames(df)[3] <- heatmap.name;
  return(df)
}

#' Converts branching times to "time intervals between branching times" and "birth at nodes" vectors
#' @inheritParams default_params_doc
#' @noRd
brts2time_intervals_and_births <- function(
  brts
) {
  time_points=-unlist(unname(sort(abs(brts), decreasing = T)))
  branching_times = -sort(abs(as.numeric(time_points)), decreasing = TRUE)
  births=c(0, unname(table(branching_times))[-1])
  unique_branching_times=as.numeric(names(table(branching_times)))
  time_intervals=c((diff(unique_branching_times)),(abs(utils::tail(unique_branching_times, 1))))
  births=births[-1]
  out=list(time_intervals=time_intervals, births=births)
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
  testit::assert(max_number_of_species < Inf)
  init_n_lineages <- soc

  # Use sim_pars for pars
  create_mbd_P_matrix <- function(
    pars,
    max_number_of_species
  )
  {
    lambda <- pars[1]
    mu <- pars[2]
    q <- pars[3]
    nvec <- 0:max_number_of_species
    my_matrix <- matrix(
      0,
      nrow = max_number_of_species + 1,
      ncol = max_number_of_species + 1
    )
    testit::assert(!"Do not call hyperA")
    # my_matrix <- lambda * hyperA::hyperA(N = max_number_of_species, k = 0, q = q)
    # my_matrix[row(my_matrix) == col(my_matrix) - 1] = mu *nvec[2:(max_number_of_species+1)]
    # diag(my_matrix) = - mu * nvec - lambda * (1 - (1 - q) ^ nvec)
    # return(my_matrix)
  }

  nvec <- 0:max_number_of_species
  v0 <- rep(0,(max_number_of_species + 1))
  v0[init_n_lineages + 1] <- 1
  transition_matrix <- create_mbd_P_matrix(
    pars = test_parameters,
    max_number_of_species = max_number_of_species
  )
  vf <- expoRkit::expv(
    v = v0,
    x = transition_matrix,
    t = age, m = precision
  )

  nmedio <- sum(nvec *vf)
  std <- sqrt(sum(nvec^2 * vf) - nmedio^2)
  if (output == 1) {
    graphics::plot(log(vf))
    print(paste("Sim pars are:", test_parameters[1], test_parameters[2], test_parameters[3], ". Average n is", nmedio, "with std:", std))
  }
  list(avg_n = nmedio, std_n = std)
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
  if (is.matrix(matrix) == FALSE) {
    matrix <- as.matrix(matrix)
  }
  matrix_t=t(matrix[nrow(matrix):1, 1:ncol(matrix)])
  if (logs ==1 ) {
  # image(log(matrix_t))
    graphics::image(log(matrix_t), col=rev(colormap),...)
  }else{graphics::image(matrix_t)}
}

# @Giappo: add doc
#' Does something J
#' @inheritParams default_params_doc
myheatmap2 <- function(x, y, z, x_name, y_name, z_name, x_splits, y_splits){

  if (missing(x_splits)) {
    x_splits <- round( (length(x)) / 10)
  }
  if (missing(y_splits)) {
    y_splits <- round( (length(y)) / 10)
  }
  if (missing(x_name)) {
    x_name <- NULL
  }
  if (missing(y_name)) {
    y_name <- NULL
  }
  if (missing(z_name)) {
    z_name <- NULL
  }

  lX <- x
  lY <- y
  pretty_X.at	<-	pretty(range(lX), n = x_splits)
  pretty_X.lab	<-	round(pretty_X.at, 2)
  pretty_Y.at	<-	pretty(range(lY), n = y_splits)
  pretty_Y.lab	<-	round(pretty_Y.at, 2)

  jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  graphics::filled.contour(
    t(z),
    color = jet.colors,
    nlevels = 100,
    ylab = y_name,
    xlab = x_name,
    plot.axes = {
      graphics::axis(
        1, at = seq(0, 1, length.out = length(pretty_X.lab)),
        labels = pretty_X.lab
      )
      graphics::axis(
        2, at = seq(0, 1, length.out = length(pretty_Y.lab)),
        labels = pretty_Y.lab
      )
    },
    key_title = graphics::title(main = z_name)
  )
}

#' Checks for NA, NaN or negative components in a vector (usually used for Qt)
#' @inheritParams default_params_doc
#' @param v sometjing
#' @param display_output something
#' @noRd
negatives_correction <- function(
  v,
  pars,
  display_output = 0
){
  problems <- 0
  if (any(is.na(v))) {
    problems <- 1
    na_components <- which(is.na(v) & !is.nan(v))
    nan_components <- which(is.nan(v))
    if (display_output == 1) {
      cat("There are non-numeric components for par values:", pars, "\n")
      if (length(na_components) > 0) { cat("NA component are:", na_components)}
      if (length(nan_components) > 0) { cat("NaN component are:", nan_components)}
    }
  }

  if (any(v < 0) & (problems) == 0) {
    #negatives <- v[v < 0]
    #maximum <- max(v)
    v[v < 0 & (abs(v) / abs(max(v))) < 1e-10] <- 0
  }
  v
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
    Atest1=( function1(...))
    Atime1=Atime1+proc.time()-t0
    t1=proc.time()
    Atest2=( function2(...))
    Atime2=Atime2+proc.time()-t1
  }
  if (output==1){print(list(Atest1,Atest2))}
  print(Atime1);print(Atime2);
  print( all.equal(Atest1,Atest2))
  return(invisible(list(Atest1,Atest2)))
}

#' Checks function calls within a function
#' @inheritParams default_params_doc
called_functions <- function(
  function_name,
  recursive = FALSE,
  checked_functions = NULL
) {

    # Get the function's code:
    function.code <- deparse(get(function_name))

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

        # checked_functions: We need to keep track of which functions
        # we've checked to avoid infinite loops.
        functs.to.check <- called.functions[!(called.functions %in% checked_functions)]

        testit::assert(!"Should never call listFunctions, whatever it is")
        # called.functions <- unique(c(called.functions,
        #     do.call(c, lapply(functs.to.check, function(x) {
        #         listFunctions(x, recursive = T, checked_functions = c(checked_functions, called.functions))
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
  black_list=as.complex(rep(0, max_k))
  negative_list=rep(0, max_k)
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
    entries_sign=rep(1, max_k)
    for (k in 1:max_k){
      lower_triangle=(as.matrix(Mlist[[k]]))[lower.tri(Mlist[[k]])]
      if (min(lower_triangle)<0){entries_sign[k]=-1;negative_list[k]=-1}
    }
    if (min(entries_sign)<0){print(paste("There is a negative entry on k=", which((entries_sign)<0)))}
  }
  return(list(black_list=black_list, negative_list=negative_list))
}

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
    x2 <- append(x2, values[[i]], after[i]+length(x2)-length(x))
  }
  x2
}

#' Function returning branch lengths named by ending node
#' @inheritParams default_params_doc
branchLengths <- function(tr)
{
  if (length(tr$node) > 0) {
    type <- 'NAMED'
  } else {
    type <- 'UNNAMED'
  }
  if (length(tr$root) > 0) {
    root.length <- tr$root
  } else{
    root.length <- 0
  }

  node.numbers <- sort(unique(tr$edge[, 1]))
  b_l_1 <- as.numeric(
    stats::na.omit(
      c(
        root.length,
        tr$edge.length[match(node.numbers, tr$edge[, 2])]
      )
    )
  )

  tip.numbers <- 1:length(tr$tip)
  b_l_2 <- tr$edge.length[match(tip.numbers, tr$edge[, 2])]

  if (type=='NAMED') {
    res <- data.frame(
      Number = c(sort(node.numbers), tip.numbers),
      Name = c(tr$node, tr$tip),
      Lengths = c(b_l_1, b_l_2)
    )
  }
  if (type=='UNNAMED') {
    res <- data.frame(
      Number = c(sort(node.numbers), tip.numbers),
      Name = c(rep(NA, tr$Nnode), tr$tip),
      Lengths = c(b_l_1, b_l_2)
    )
  }

  return(res)
}

#' Does something
#' @inheritParams default_params_doc
summarize_beast_posterior <- function(
  input_trees_path,
  input_xml = NULL,
  subsamp = NULL
){


  #If an XML with a starting tree is provided, get the node names from that tree
  node_names_xml=NULL
  if(!is.null(input_xml)){
    subst.xml=scan(input_xml, what="raw")
    where.starting.tree=which(subst.xml=="id=\"startingTree\">")+1
    start.tr <- ape::read.tree(text=subst.xml[where.starting.tree])
    node_names_xml=start.tr$node.label
  }

  branch_info_tot=tree_tot=list()
  subsamps=c()

  for(t in 1:length(input_trees_path)){

    #Read the trees
    subst.tr=scan(input_trees_path[t], what="raw")

    #Find the equal signs in tree file
    equals=(1:length(subst.tr))[subst.tr=="="]

    #Use position of equal signs to find position of trees and their names
    tree.positions=equals+2
    tree.names.pos=equals-2

    #How many trees?
    n_trees=length(tree.positions)
    if(is.null(subsamp)) subsamp=n_trees
    samp=sample(n_trees, round(subsamp / length(input_trees_path)), replace=FALSE)

    #Get correspondance between numbers and taxa
    #How many taxa?
    split1=strsplit(subst.tr[5], "=")[[1]][2]
    ntax=as.numeric(strsplit(split1, ";")[[1]])

    #Find the "Translate"
    Translate.pos=(1:length(subst.tr))[subst.tr=="Translate"]
    comma.pos=which(subst.tr==", ")
    if(length(comma.pos)==0){
      numbers.pos=seq(Translate.pos+1,Translate.pos+2 *ntax, 2)
    } else{
      numbers.pos=seq(Translate.pos+1,Translate.pos+3 *ntax, 3)
    }
    taxa.pos=numbers.pos+1
    numbers=subst.tr[numbers.pos]
    taxa=array(, length(taxa.pos))
    for(i in 1:length(taxa.pos)) taxa[i]=strsplit(subst.tr[taxa.pos[i]], ", ")[[1]][1]


    for(i in 1:length(samp)){

      #Isolate the trees and tree names
      tree1=subst.tr[tree.positions[samp[i]]]

      #Break the character string
      tree_sub=substring(tree1, 0:100000, 0:100000)

      #Trim the end
      tree_sub2=tree_sub[tree_sub!=""]

      colons=which(tree_sub2==":")
      parenthesis_colons=(colons-1)[tree_sub2[colons-1]==")"]
      if(is.null(node_names_xml)){
        node_names=paste("Node", seq_along(c(NA, parenthesis_colons)), sep="")
      } else {
        node_names <- node_names_xml
      }
      tree_sub2=append_multiple(tree_sub2, node_names, c(parenthesis_colons, length(tree_sub2)-1))

      #Find positions of opening and closing brackets
      croch1.pos=(1:length(tree_sub2))[tree_sub2=="["]
      croch2.pos=(1:length(tree_sub2))[tree_sub2=="]"]

      #Trim the brackets and the information in-between

      remove=c()
      for (k in 1:length(croch1.pos))
      {
        commas=which(tree_sub2[croch1.pos[k]:croch2.pos[k]]==", ")
        commas=c(croch1.pos[k], commas+croch1.pos[k]-1, croch2.pos[k])
        node_name_k=tree_sub2[croch1.pos[k]-2]
        name_arg=value_arg=array(, length(commas)-1)
        for(j in 1:(length(commas)-1)){

          arg_j=tree_sub2[(commas[j]+1):(commas[j+1]-1)]
          equal=which(arg_j=="=")
          name_arg[j]=paste(arg_j[2:(equal-1)], collapse="")
          value_arg[j]=paste(arg_j[(equal+1):length(arg_j)], collapse="")

        }

        if(k==1){
          res=array(, c(length(croch1.pos), length(value_arg)+1))
          colnames(res)=c("name", name_arg)
        }

        res[k, ]=c(node_name_k, value_arg)
        remove=c(remove, croch1.pos[k]:croch2.pos[k])
      }

      res[res[, "name"]%in%numbers, "name"]=taxa[match(res[res[, "name"]%in%numbers, "name"], numbers)]

      #Re-paste the new tree into a single character string
      tree.text=paste(tree_sub2[-remove], collapse="")

      tr <- ape::read.tree(text=tree.text)
      tr$tip.label=taxa[as.numeric(tr$tip.label)]

      #trees[[i]]=tr
      branch_info=branchLengths(tr)[-1, ]
      res=res[match(branch_info[, "Name"], res[, "name"]), ]
      branch_info=cbind(tr$edge[match(branch_info[, "Number"], tr$edge[, 2]), ],
                        #tr$edge.length[match(branch_info[, "Number"], tr$edge[, 2])],
                        branch_info[, -1], res[, -1, drop=FALSE])
      colnames(branch_info)[1:2]=c("Ancestor", "Descendant")

      if(i==1) {
        branch_infos=array(, c(dim(branch_info), length(samp)))
        dimnames(branch_infos)=list(NULL, dimnames(branch_info)[[2]], 1:length(samp))
      }

      branch_infos[,, i]=as.matrix(branch_info)
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

  names_i <- file_name

  beast_posterior <- summarize_beast_posterior(
      input_trees_path = names_i,
      subsamp=maxtree #Number of samples to keep
    )

  # beast_posterior$trees # a list of the trees in phylo format
  # beast_posterior$branch_info # a list of branch lengths and associated information
  brts <- list()
  for (i in 1:maxtree) {
    brts[[i]] <- unname(ape::branching.times(beast_posterior$trees[[i]]))
  }

  dist <- vector("list")
  maxindex <- max(unlist(lapply(brts, length)))

  for (i in 1:maxindex){
    dist[[i]] <- unlist(lapply(brts, "[[", i))
  }

  return(list(brts=brts, dist=dist))
}

#' Imports data
#' No globals please
#' @inheritParams default_params_doc
# @param parnames sonething
# mbd_import_data <- function(
#   parnames = get_mbd_param_names()
# ) {
#   parnames <<- parnames
#   n_pars <<- length(parnames);
#   base_path=dirname(dirname(getwd()))
#   if( !exists("path")  ){ #&& interactive()
#     path<<-choose.dir(paste(base_path, "/Progress /RQ1 - Multiple Births /Results + Reports / 4_parameters /", sep = ''), "Choose a suitable folder")
#   }
#   res_files = list.files(pattern=paste('[.]txt', sep = ''), path=path, full.names=TRUE)
#   for(i in 1:length(res_files)){
#     fileData <- utils::read.table(file=res_files[i], header=FALSE, sep=", ")
#     ifelse(exists("targetTable"), targetTable<-rbind(targetTable, fileData), targetTable<-fileData)
#   }
#   all_results = targetTable
#
#   tidy_results=all_results[order(all_results[, dim(all_results)[2]]), ]
#   dimnames(tidy_results)[[2]]<-(c(parnames, "LL", "multiple_born", "number_of_tips", "percentage_multiple_species", "tree_id"))
#   bad_results=tidy_results[rowSums(tidy_results[, 1:(n_pars+1)]==rep(-1,(n_pars+1)))==(n_pars+1), ];#print(Nbad<-dim(bad_results)[1])
#   results=tidy_results[rowSums(tidy_results[, 1:(n_pars+1)]==rep(-1,(n_pars+1)))!=(n_pars+1), ];#print(N<-dim(results)[1])
#   print(paste("There are ",Nbad<-dim(bad_results)[1], " bad results.", sep = ''))
#   print(paste("There are ",N<<-dim(results)[1], " good results.", sep = ''))
#
#   load( as.character(paste(path, "/ data / general_settings", sep = '')), envir = e <- globalenv())
#   load( as.character(paste(path, "/ data / sim_data", sep = '')), envir = e <- globalenv()  )
#   suppressWarnings( rm(targetTable, all_results, fileData, i))
#   return(results)
# }
