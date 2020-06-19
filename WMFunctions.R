create.supraadjacency.matrix <- function(WholeNet, type, N, L, delta_zeta, is.weighted.graph){

  Graphs <- WholeNet[which(lapply(WholeNet, `[[`, "type") == type)]
  Graphs <- lapply(Graphs, `[[`, "graph")
  
  Idem_Matrix <- Diagonal(N, x = 1)
  
  Col_Node_Names <- character()
  Row_Node_Names <- character()
  
  t1 <- Sys.time()  
  if (getDoParWorkers()>1) registerDoSEQ()
  cl <- makeCluster(L, type = "FORK")
  registerDoParallel(cl, cores = L)
  SupraAdjacencyResult <- foreach (i=1:L, .packages=c('igraph', 'Matrix')) %dopar% {
    
    SupraAdjacencyMatrixPart <- Matrix(0,ncol=N*L,nrow=N,sparse = TRUE)
    
    Adjacency_Layer <-  as_adjacency_matrix(Graphs[[i]], attr="weight", sparse = TRUE)
    
    if (!is.weighted.graph){
      Adjacency_Layer <-  as_adjacency_matrix(Graphs[[i]], sparse = TRUE)
    }
    
    ## We order the matrix by the node name. This way all the matrix will have the same. Additionally we include a label with the layer number for each node name.
    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),order(colnames(Adjacency_Layer))]
    Layer_Row_Col_Names <- paste(colnames(Adjacency_Layer),i,sep="_")
    
    ## We fill the diagonal blocks with the adjacencies matrix of each layer.
    Position_ini_row <- 1
    Position_end_row <- N
    Position_ini_col <- 1 + (i-1)*N
    Position_end_col <- N + (i-1)*N
    SupraAdjacencyMatrixPart[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] <- (1-delta_zeta)*(Adjacency_Layer)
    
    # avoid division by zero for monoplex network
    L_mdfd <- L-1
    if (L == 1) L_mdfd <- 1
    
    ## We fill the off-diagonal blocks with the transition probability among layers.
    for (j in 1:L){
      Position_ini_col <- 1 + (j-1)*N
      Position_end_col <- N + (j-1)*N
      if (j != i){
        SupraAdjacencyMatrixPart[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] <- (delta_zeta/(L_mdfd))*Idem_Matrix
      }
    }
    return(list(SupraAdjacencyMatrixPart, Layer_Row_Col_Names))
  }
  
  stopCluster(cl)
  #cat(type, " SupraAdjacency Time for parallel part: ", format(Sys.time()-t1), "\n")
  t2 <- Sys.time()  
  #SupraAdjacencyMatrix <- Matrix(0,ncol=N*L,nrow=N*L,sparse = TRUE)
  SupraAdjacencyResult <- unlist(SupraAdjacencyResult, recursive = FALSE)
  
  
  # Row-Col names are even indexed
  Col_Names <- do.call('c',SupraAdjacencyResult[seq(2,2*L,by=2)])
  
  # Parallele parts of the SupraAdjacencyMatrix are odd indexed
  SupraAdjacencyMatrix <- do.call('rbind', SupraAdjacencyResult[seq(1,2*L,by=2)])
  
  #SupraAdjacencyMatrix <- rbind.fill.matrix(SupraAdjacencyResult[seq(1,2*L,by=2)])
  
  #cat("SupraAdjacency Time for rest: ", format(Sys.time()-t2), "\n")
  rownames(SupraAdjacencyMatrix) <- Col_Names
  colnames(SupraAdjacencyMatrix) <- Col_Names  
  
  return(SupraAdjacencyMatrix)
}
create.bipartite.matrix <- function(WholeNet, N, M, gene_pool_nodes_sorted, phenotype_pool_nodes_sorted, numCores, isWeighted){
  #WholeNet <- FullNet
  Gene_phenotype_Network <- WholeNet[which(lapply(WholeNet, `[[`, "type") == "bipartite")]
  Gene_phenotype_Network <- lapply(Gene_phenotype_Network, `[[`, "DF")[[1]]
  
  # Get the Subset of Gene-phenotype relations which have common genes in whole network
  Gene_phenotype_Network <- Gene_phenotype_Network[which(Gene_phenotype_Network$from %in% gene_pool_nodes_sorted), ]
  Gene_phenotype_Network <- Gene_phenotype_Network[which(Gene_phenotype_Network$to %in% phenotype_pool_nodes_sorted), ]
  
  Gene_phenotype_Network <- graph.data.frame(Gene_phenotype_Network, directed = FALSE) 
  
  el <- as_edgelist(Gene_phenotype_Network)
  value <- edge_attr(Gene_phenotype_Network, name = "weight")
  if (!isWeighted){
    value <- rep(1, nrow(el))
  }
  
  Bipartite_matrix <- Matrix(data=0, nrow=N, ncol=M)
  rownames(Bipartite_matrix) <- gene_pool_nodes_sorted
  colnames(Bipartite_matrix) <- phenotype_pool_nodes_sorted
  rindx <- unlist(mclapply(el[,1], function(x) which(rownames(Bipartite_matrix) %in% x), mc.cores=20))
  cindx <- unlist(mclapply(el[,2], function(x) which(colnames(Bipartite_matrix) %in% x), mc.cores=20))
  
  lenRindx <- length(rindx)
  partLen <- floor(lenRindx/numCores)
  
  
  if (partLen == 0){
    cat("numCores:", numCores, " - partLen: ", partLen, "\n")
    stop("Assigned numCores is greater than data length! Assign less!")
  }
  
  rindx_parts <- list()
  cindx_parts <- list()
  value_parts <- list()
  for(i in 1:numCores){
    stInd <- (i-1)*partLen + 1
    endInd <- i*partLen
    if (i==numCores){
      endInd <- lenRindx
    }
    rindx_parts[[i]] <- rindx[stInd:endInd]
    cindx_parts[[i]] <- cindx[stInd:endInd]
    value_parts[[i]] <- value[stInd:endInd]
  }
  
  cl <- makeCluster(numCores, type="FORK")
  registerDoParallel(cl, cores = numCores)
  Bipartite_matrix_result <- foreach (i=1:numCores, .packages='Matrix') %dopar% {
    for(j in 1:length(rindx_parts[[i]])){
      Bipartite_matrix[rindx_parts[[i]][j],cindx_parts[[i]][j]] <- value_parts[[i]][j]
    }
    return(Bipartite_matrix)
  }
  stopCluster(cl)
  Bipartite_matrix <- Reduce('+', Bipartite_matrix_result)
  
  return(Bipartite_matrix)
}
create.suprabibartite.matrix <- function(Bipartite_matrix, N, M, LG, LP){
  SupraBipartiteMatrix <- Matrix(0,nrow=N*LG,ncol=M*LP,sparse = TRUE)
  
  Row_Node_Names <- sprintf(paste0(rep(rownames(Bipartite_matrix), LG), "_%d"), 
                            rep(seq_len(LG), each=N))
  SupraBipartiteMatrix <- do.call(rbind, replicate(LG, Bipartite_matrix, simplify=FALSE))
  
  rownames(SupraBipartiteMatrix) <- Row_Node_Names
  
  
  Col_Node_Names <- sprintf(paste0(rep(colnames(Bipartite_matrix), LP), "_%d"), 
                            rep(seq_len(LP), each=M))
  
  SupraBipartiteMatrix <- do.call(cbind, replicate(LP, SupraBipartiteMatrix, simplify=FALSE))
  colnames(SupraBipartiteMatrix) <- Col_Node_Names  
  return(SupraBipartiteMatrix)
}
create.transition.matrix.gene_phenotype <- function(SupraBipartiteMatrix, N, M, LG, LP, lambda, no.cores){
  #### Transition Matrix for the inter-subnetworks links
  Transition_Gene_phenotype <- Matrix(0,nrow=N*LG,ncol=M*LP,sparse = TRUE)
  
  Col_Sum_Bipartite <- colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)
  
  if (getDoParWorkers()>1) registerDoSEQ()
  cl <- makeCluster(no.cores, type = "FORK")
  registerDoParallel(cl, cores = no.cores)
  
  # columnwise normalization for propability  
  Transition_Gene_phenotype2 <-foreach (j = 1:(M*LP), .packages=c('Matrix')) %dopar% {
    if (Col_Sum_Bipartite[j] != 0){
      Transition_Gene_phenotype[,j] <- (lambda*SupraBipartiteMatrix[,j]) /Col_Sum_Bipartite[j]
    }else{
      Transition_Gene_phenotype[,j] <- Transition_Gene_phenotype[,j]
    }
  }  
  stopCluster(cl)
  
  Transition_Gene_phenotype2 <- Matrix(unlist(Transition_Gene_phenotype2, recursive = FALSE),
                                    nrow=N*LG,ncol=M*LP,sparse = TRUE)  

  
  rownames(Transition_Gene_phenotype2) <- rownames(SupraBipartiteMatrix)
  colnames(Transition_Gene_phenotype2) <- colnames(SupraBipartiteMatrix)
  Transition_Gene_phenotype <- Transition_Gene_phenotype2  
  
  return(Transition_Gene_phenotype)
}
create.transition.matrix.phenotype_gene <- function(SupraBipartiteMatrix, N, M, LG, LP, lambda, no.cores){
  Transition_phenotype_Gene <- Matrix(0,nrow=M*LP,ncol=N*LG,sparse = TRUE)
  
  
  Row_Sum_Bipartite <- rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)
  
  
  if (getDoParWorkers()>1) registerDoSEQ()
  cl <- makeCluster(no.cores, type = "FORK")
  registerDoParallel(cl, cores = no.cores)
  
  Transition_phenotype_Gene2 <- foreach (i=1:(N*LG), .packages=c('Matrix')) %dopar% {
    if (Row_Sum_Bipartite[i] != 0){
      Transition_phenotype_Gene[,i] <- (lambda*SupraBipartiteMatrix[i,])/Row_Sum_Bipartite[i]
    }else{
      Transition_phenotype_Gene[,i] <- Transition_phenotype_Gene[,i]
    }
  }
  stopCluster(cl)
  #print(Transition_phenotype_Gene2)
  
  Transition_phenotype_Gene2 <- Matrix(unlist(Transition_phenotype_Gene2, recursive = FALSE),
                                    nrow=M*LP,ncol=N*LG,sparse = TRUE)
  
  #print(Transition_phenotype_Gene2)
  
  #print(dim(Transition_phenotype_Gene2))
  #print(dim(SupraBipartiteMatrix))
  rownames(Transition_phenotype_Gene2) <- colnames(SupraBipartiteMatrix)
  colnames(Transition_phenotype_Gene2) <- rownames(SupraBipartiteMatrix)
  Transition_phenotype_Gene <- Transition_phenotype_Gene2
  #printSpMatrix2(Transition_phenotype_Gene, col.names = TRUE)  
  
  return(Transition_phenotype_Gene)
}
create.gene.transition.multiplex.network <- function(SupraAdjacencyMatrix, SupraBipartiteMatrix, N, LG, lambda, numCores){
  #### Transition Matrix for the intra-subnetworks links
  Transition_Multiplex_Network <- Matrix(0,nrow=N*LG,ncol=N*LG,sparse = TRUE)
  
  rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
  colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)
  
  Col_Sum_Multiplex <- colSums(SupraAdjacencyMatrix,na.rm=FALSE,dims=1, sparseResult=FALSE)
  Row_Sum_Bipartite <- rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE) 
  
  
  partLen <- floor(N*LG/numCores)
  # below can only happen with toy samples
  if (partLen==0){
    stop("Assigned numCores is greater than data length! Assign less!")
  }
  
  parts <- vector('list',numCores)
  for(i in 1:numCores){
    stInd <- (i-1)*partLen + 1
    endInd <- i*partLen
    if (i==numCores){
      endInd <- N*LG
    }
    parts[[i]][['start']] <- stInd
    parts[[i]][['end']] <- endInd
  }
  
  cl <- makeCluster(numCores, type="FORK")
  
  registerDoParallel(cl, cores = numCores)
  
  Transition_Multiplex_Network_Result <- foreach (i=1:numCores, .packages='Matrix') %dopar% {  
    for (j in parts[[i]]["start"]:parts[[i]]["end"]){
      if(Row_Sum_Bipartite[j] != 0){
        Transition_Multiplex_Network[,j] <- ((1-lambda)*SupraAdjacencyMatrix[,j]) /Col_Sum_Multiplex[j]
      } else {
        Transition_Multiplex_Network[,j] <- SupraAdjacencyMatrix[,j] /Col_Sum_Multiplex[j]
      }
    }
    return(Transition_Multiplex_Network)
  }
  
  Transition_Multiplex_Network <- Reduce('+', Transition_Multiplex_Network_Result)
  stopCluster(cl)
  
  return(Transition_Multiplex_Network)
}
create.phenotype.transition.multiplex.network <- function(SupraAdjacencyMatrixphenotype, SupraBipartiteMatrix, M, LP, lambda, numCores){
  Transition_Multiplex_Network_phenotype <- Matrix(0,nrow=M*LP,ncol=M*LP,sparse = TRUE)
  
  rownames(Transition_Multiplex_Network_phenotype) <- rownames(SupraAdjacencyMatrixphenotype)
  colnames(Transition_Multiplex_Network_phenotype) <- colnames(SupraAdjacencyMatrixphenotype)
  
  Col_Sum_Multiplex <- colSums(SupraAdjacencyMatrixphenotype,na.rm=FALSE,dims=1, sparseResult=FALSE)
  Col_Sum_Bipartite <- colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE) 
  
  partLen <- floor(M*LP/numCores)
  # below can only happen with toy samples
  if (partLen==0){
    stop("Assigned numCores is greater than data length! Assign less!")
  }
  
  parts <- vector('list',numCores)
  for(i in 1:numCores){
    stInd <- (i-1)*partLen + 1
    endInd <- i*partLen
    if (i==numCores){
      endInd <- M*LP
    }
    parts[[i]][['start']] <- stInd
    parts[[i]][['end']] <- endInd
  }
  cl <- makeCluster(numCores, type="FORK")
  
  registerDoParallel(cl, cores=numCores)
  Transition_Multiplex_Network_phenotype_Result <- foreach (i=1:numCores, .packages='Matrix') %dopar% {  
    for (j in parts[[i]]["start"]:parts[[i]]["end"]){
      if(Col_Sum_Bipartite[j] != 0){
        Transition_Multiplex_Network_phenotype[,j] <- ((1-lambda)*SupraAdjacencyMatrixphenotype[,j]) /Col_Sum_Multiplex[j]
      } else {
        Transition_Multiplex_Network_phenotype[,j] <- SupraAdjacencyMatrixphenotype[,j] /Col_Sum_Multiplex[j]
      }
    }
    return(Transition_Multiplex_Network_phenotype)
  }
  
  Transition_Multiplex_Network_phenotype <- Reduce('+', Transition_Multiplex_Network_phenotype_Result)
  stopCluster(cl) 
  return(Transition_Multiplex_Network_phenotype)
}
create.WalkMatrix <- function(inputFileName, numCores){
  # inputFileName <- "input_files"
  # numCores <- 10
  network_range <<- c(0.001, 1) 
  t0 <- Sys.time()
  t1 <- Sys.time()
  Settings <- read.settings(inputFileName)
  Parameters <- Settings$Params
  FilesDF <- Settings$FilesDF
  LG <- Settings$LG
  LP <- Settings$LP

  cat("Creating Walk Matrix For: ", paste0(FilesDF$file_name[FilesDF$type%in%c("gene", "phenotype", "bipartite")], collapse = ","), "\n")
  
  FullNet <- read.network.layers(FilesDF, Parameters)
  gene_pool_nodes_sorted <- FullNet$gene_pool_nodes_sorted
  phenotype_pool_nodes_sorted <- FullNet$phenotype_pool_nodes_sorted
  FullNet <- FullNet$NetworkLayers
  
  N=length(gene_pool_nodes_sorted)
  M <- length(phenotype_pool_nodes_sorted)

  SupraAdjacencyMatrix <- create.supraadjacency.matrix(FullNet, "gene", N, LG, Parameters$delta, TRUE)

  
  SupraAdjacencyMatrixphenotype <- create.supraadjacency.matrix(FullNet, "phenotype", M, LP, Parameters$zeta, TRUE)
 
  BipartiteMatrix <- create.bipartite.matrix(FullNet, N, M, gene_pool_nodes_sorted, phenotype_pool_nodes_sorted, numCores, Parameters$weighted)

  SupraBipartiteMatrix <- create.suprabibartite.matrix(BipartiteMatrix, N, M, LG, LP)

  Transition_Gene_phenotype <- create.transition.matrix.gene_phenotype(SupraBipartiteMatrix, N, M, LG, LP, Parameters$lambda, numCores)

  Transition_phenotype_Gene <- create.transition.matrix.phenotype_gene(SupraBipartiteMatrix, N, M, LG, LP, Parameters$lambda, numCores)

  Gene_Transition_Multiplex_Network <- create.gene.transition.multiplex.network(SupraAdjacencyMatrix, 
                                                                                SupraBipartiteMatrix, 
                                                                                N, LG, Parameters$lambda, numCores)
  phenotype_Transition_Multiplex_Network <- create.phenotype.transition.multiplex.network(SupraAdjacencyMatrixphenotype, 
                                                                                    SupraBipartiteMatrix, 
                                                                                    M, LP, Parameters$lambda, numCores)

  Multiplex_Heterogeneous_Matrix <- rbind(cbind(Gene_Transition_Multiplex_Network, Transition_Gene_phenotype),
                                          cbind(Transition_phenotype_Gene, phenotype_Transition_Multiplex_Network))
  
  cat("Time to generate Walk Matrix : ", format(Sys.time()-t0), "\n")  
  # if there are still some paralllel workers stop and force sequential
  registerDoSEQ() 
  return(list(WM = Multiplex_Heterogeneous_Matrix, 
       genes = gene_pool_nodes_sorted,
       phenotypes = phenotype_pool_nodes_sorted,
       LG = LG,
       LP = LP,
       N = N,
       M = M))
}