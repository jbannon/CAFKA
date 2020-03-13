#####
#
#
#
# CAFKA: Causality-Informed Fisher kernels
#  @jamesBannon
#
#

#
####

get_parent_permutations <- function(nparents, max_parents = 4){
  # iterate over over parent list
  if(nparents >max_parents){
    stop("too many parents passed",nparents)
  }
  perm_mat <- permutations(2,nparents,v=c(0,1),repeats.allowed =TRUE)
  return(perm_mat)
}

strip_list<-function(tronco_model,issue_list,mode="indistinguishable"){
  
  if(!tolower(mode) %in% c("indistinguishable","ones","zeros","zeroes")){
    stop("Invalid Mode Passed: ",mode)
  }
  
  if(length(issue_list)==0){
    return(list(model=tronco_model,removed=list()))
  }
  
  
  if(tolower(mode)=="indistinguishable"){
    removed_events = c() #empty vector
    for(i in seq(1,length(issue_list))){
      for(j in seq(2,nrow(issue_list[[i]]))){
        temp_t <- issue_list[[i]][j,1]
        temp_g <- issue_list[[i]][j,2]
        tronco_model <- delete.event(tronco_model,gene=temp_g,type=temp_t)
        event_string <- paste0(temp_g,"-",temp_t)
        removed_events<-append(removed_events,event_string)
      }
      removed_events <- append(removed_events,"-1") # sentinel marking end of equivalence a class
    }
    
    results <- list(model = tronco_model, removed = removed_events)
    return(results)
    
  } else{
    removed_events = c() 
    for(i in seq(1,length(issue_list))){
      temp_t <- issue_list[[i]][1]
      temp_g <- issue_list[[i]][2]
      tronco_model <- delete.event(tronco_model,gene=temp_g,type=temp_t)
      event_string <- paste0(temp_g,"-",temp_t)
      removed_events<-append(removed_events,event_string)
    }
    results <- list(model = tronco_model, removed = removed_events)
    return(results)
  }
}

handle_consolidated_data <-function(model,pathologies,only=NULL){
  removed_events = list()
  cat("\nProcessing Indistinguishable Events")
  results <- strip_list(model,pathologies$indistinguishable,mode="indistinguishable")
  removed_events$indistinguishable <- results$removed
  model<-results$model  
  
  cat("\nRemoving Probability 1 Events")
  results <-strip_list(model,pathologies$ones,mode="ones")
  removed_events$ones <-results$removed
  model <- results$model
  
  cat("\nRemoving Probability 0 Events")
  
  results <-strip_list(model,pathologies$zeroes,mode="zeroes")
  removed_events$zeros <-results$removed
  model <- results$model
  
  saveRDS(removed_events,"removed_events_list.rds") #more record keeping
  
  final_results <- list(tronco_model = model,removed=removed_events)
  
  return(final_results)
  
}

compute_parameters<- function(data_set, node_names,parent_list,
                              smooth=TRUE, 
                              smoothing.constant = 1,
                              max_parents=4,
                              epsilon=0.01){
  
  params <- matrix(list(0),nrow(parent_list),ncol(parent_list))
  d<-2 # number of possible outcomes
  npar_toggle <-FALSE
  for(i in 1:nrow(parent_list)){
    parameter_vector <- numeric(0)
    parents <- parent_list[i,1][[1]]
    parent_names <- node_names[parents]
    nparents <- length(parents)
    gene_name <- rownames(parent_list)[i]
    all_names <- c(gene_name,parent_names)
    
    if(length(nparents)==1 && parents[1]==-1){ #rootlike node
      parameter<- (sum(data_set[,gene_name])+smoothing.constant)/(nrow(data_set)+smoothing.constant*d)
      parameter_vector <- c(parameter_vector,parameter)
    }else{
      
      possible_configurations <- get_parent_permutations(nparents,max_parents)
      
      relevant_data <- data.frame(data_set[,all_names])
      
      for(j in 1:nrow(possible_configurations)){
        joint_instances <- nrow(relevant_data[which(apply(relevant_data, 1, function(x) identical(unname(x), c(1,possible_configurations[j,])))),]) + smoothing.constant
        joint_instances <- ifelse(joint_instances==0,epsilon,joint_instances)
        marginal_likelihood <- nrow(relevant_data[which(apply(as.matrix(relevant_data[,parent_names]),1,function(x) identical(unname(x),possible_configurations[j,]))),])+smoothing.constant*d
        parameter <- ifelse(marginal_likelihood==0,0.5,joint_instances/marginal_likelihood)
        parameter_vector <- c(parameter_vector, parameter)
      }
    }
    
    params[i,1] <- list(parameter_vector)
  }
  
  return(params)
}


## gives raw fisher kernel in induced feature space with just gradient w.r.t params of Loglikelihood
## these samples can then be used later to construct 
compute_all_sample_gradients<- function(data_set, node_names,parent_list, 
                                        model_parameters,
                                        max_parents=4,
                                        normalize=FALSE){
  ## need to make sure this has same cannonical ordering
  ## takes in a data set which has patients as rows and events as columns
  ## we iterate over the samples and compute the fk gradient at each sample
  ##  (returned by as a list) which then we stack (rbind) into a n_patients x L matrix
  
  embedded_design_matrix <- matrix(0, nrow(data_set),length(unlist(model_parameters)))
  
  rownames(embedded_design_matrix)<-rownames(data_set)
  for(i in 1:nrow(data_set)){
    ## embedding of this sample
    sample <- data_set[i,]
    gradient_vector<- numeric(0)
    for(j in 1:nrow(parent_list)){
      parents <- parent_list[j,1][[1]]
      parent_names <- node_names[parents]
      nparents <- length(parents)
      
      gene_name <- rownames(parent_list)[j]
      
      all_names <- c(gene_name,parent_names)
      
      if(length(nparents)==1 && parents[1]==-1){ #rootlike node
        obtained_theta <- model_parameters[j,][[1]]
        gradient_subvector <- ((1/obtained_theta)^sample[gene_name])*((-1/(1-obtained_theta))^(1-sample[gene_name]))
      }else{
        possible_configurations <- get_parent_permutations(nparents,max_parents)
        gradient_subvector <- integer(nrow(possible_configurations))
        thetas <- model_parameters[j,][[1]]
        dex <- which(apply(possible_configurations, 1, function(x) return(all(x == sample[parent_names]))))
        obtained_theta <- thetas[dex]
        mutation_observed <- as.numeric(sample[gene_name])
        gradient_subvector[dex]<-((1/obtained_theta)^sample[gene_name])*((-1/(1-obtained_theta))^(1-sample[gene_name]))
      }
      gradient_vector <-unname(c(gradient_vector,gradient_subvector))
    }
      embedded_design_matrix[i,] <- gradient_vector
  }
  return(embedded_design_matrix)
}

rename_parent_lists <- function(event_matrix, parent_matrix){
  for(i in 1:nrow(parent_matrix)){
    gene_id <- event_matrix[rownames(parent_matrix)[i],"event"]
    rownames(parent_matrix)[i] <- gene_id
  }
  return(parent_matrix)
}

get_sbcn_info <- function(X_train,fit="bic",p.val=0.01){
  model <- import.genotypes(X_train)
  model_issues <- consolidate.data(model)
  model_issues_resolution <- handle_consolidated_data(model,model_issues)
  model_resolved <- model_issues_resolution$tronco_model
  
  SBCN <- tronco.capri(model_resolved,pvalue=p.val)
  
  if(fit=="bic"){
    parent_lists <- SBCN$model$capri_bic$parents.pos
  } else{
    parent_lists <- SBCN$model$capri_aic$parents.pos
  }
  
  parent_lists <- rename_parent_lists(as.events(SBCN),parent_lists)
  sbcn_params <- compute_parameters(as.matrix(select(X_train, as.genes(SBCN))),
                                    node_names = rownames(parent_lists),
                                    parent_list=parent_lists)
  
  network.info<-list(node.names=as.genes(SBCN),
                     parent.list=parent_lists)
  cat("\n")
  return(list(model=SBCN,
              network.info=network.info,
              network.parameters=sbcn_params))
}

genotypes_to_CFK <- function(genotypes, 
                             thresh=0.0, 
                             k=100, 
                             filter.mode="top",
                             fisher.mode="practical"){
  if(filter.mode=="top"){
    genotypes <- select(genotypes,names(sort(colSums(genotypes),
                                          decreasing = TRUE))[1:min(k,ncol(genotypes))])
    SBCN <- get_sbcn_info(genotypes)
    cfk<-compute_fisher_kernel(as.matrix(genotypes),
                              node_names=SBCN$network.info$node.names,
                              parent_list=SBCN$network.info$parent.list,
                              model_parameters = SBCN$network.parameters,mode = fisher.mode)
    return(cfk)
  }
}
