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
    
    if(normalize){
      norm <- twonorm(gradient_vector)
      if(norm==0 || is.na(norm)){
        embedded_design_matrix[i,] <- gradient_vector
      } else{
        embedded_design_matrix[i,] <- gradient_vector/twonorm(gradient_vector)
      }
    } else{
      embedded_design_matrix[i,] <- gradient_vector
    }
  }
  return(embedded_design_matrix)
}