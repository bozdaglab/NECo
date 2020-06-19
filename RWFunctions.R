read.settings <- function(input.file){
  
  files <-  read.table(paste0(input.file, ".txt"), header=TRUE, sep="\t", stringsAsFactors = FALSE)
  LG <- length(files$type[files$type=="gene"])
  LP <- length(files$type[files$type=="phenotype"])
  
  Parameters <- get.default.parameters(LG, LP)  
  
  output=list(
    FilesDF=files, 
    Params=Parameters, 
    LG=LG,
    LP=LP)   
}
read.network.layers <- function(filesDF, Parameters_File, Layers){
  # read the gene layer names
  filesDF <- filesDF[filesDF$type%in%c("gene", "phenotype", "bipartite"),]
  NetworkLayers <- vector("list",nrow(filesDF))
  j <- 1
  for(f in filesDF$file_name){
    NetworkLayers[[j]][["DF"]] <-  readRDS(paste0(f, ".rds"))
    NetworkLayers[[j]][["type"]] <- filesDF$type[j]
    
    NetworkLayers[[j]][["graph"]]  <- graph.data.frame(NetworkLayers[[j]][["DF"]], directed = FALSE)
    
    NetworkLayers[[j]][["name"]] <- f
    j <- j+1
  }
  
  gene_pool_nodes_sorted <- generate.pool.nodes(NetworkLayers, type = "gene")
  phenotype_pool_nodes_sorted <- generate.pool.nodes(NetworkLayers, type = "phenotype")
  
  
  idx <- which(lapply(NetworkLayers, `[[`, "type") == "gene")
  for(i in idx){
    NetworkLayers[[i]][["graph"]] <- add.missing.nodes.to.graph(NetworkLayers[[i]][["graph"]], 
                                                                "gene", gene_pool_nodes_sorted)  
  }
  
  idx <- which(lapply(NetworkLayers, `[[`, "type") == "phenotype")
  for(i in idx){
    NetworkLayers[[i]][["graph"]] <- add.missing.nodes.to.graph(NetworkLayers[[i]][["graph"]], 
                                                                "phenotype", phenotype_pool_nodes_sorted)  
  }  
  
  output=list(
    NetworkLayers=NetworkLayers, 
    gene_pool_nodes_sorted=gene_pool_nodes_sorted,
    phenotype_pool_nodes_sorted=phenotype_pool_nodes_sorted)     
  return(output)
}
add.missing.nodes.to.graph <- function(g, type, pool_nodes_sorted){
  ## We add to each layer the missing nodes of the total set of nodes, of the pool of nodes.
  Node_Names_Layer <- V(g)$name
  Missing_Nodes <- pool_nodes_sorted[which(!pool_nodes_sorted %in% Node_Names_Layer)]
  g <- add_vertices(g ,length(Missing_Nodes), name=Missing_Nodes)
  #print(vcount(g))
  return(g)
}
generate.pool.nodes <- function(FullNet, type){
  idx <- which(lapply(FullNet, `[[`, "type") == type)
  DFs <- lapply(FullNet[idx], `[[`, "DF")
  Node_Names_all <- unique(c(unlist(lapply(DFs, '[[', 'from')), unlist(lapply(DFs, '[[', 'to'))))  
  
  ## We remove duplicates and sort
  pool_nodes_sorted <- sort(Node_Names_all)  
  
  return(pool_nodes_sorted)
}
get.default.parameters <- function(LG, LP){
  r <- 0.7
  delta <- 0.5
  zeta <- 0.5
  tau <- rep(1, LG)
  phi <- rep(1, LP)
  lambda <- 0.5
  eta <- 0.5
  weighted <- 1
  parameters <- list(r, delta, zeta, tau, phi, lambda, eta, weighted)
  names(parameters) <- c("r", "delta", "zeta", "tau", "phi", "lambda", "eta", "weighted")
  return(parameters)
}
get.seed.scores <- function(GeneSeeds, phenotypeSeeds, eta, LG, LP, tau, phi) {
  
  n <- length(GeneSeeds)
  m <- length(phenotypeSeeds)
  
  if ((n != 0 && m!= 0)){
    
    Seed_Genes_Layer_Labeled <- paste0(rep(GeneSeeds,LG), sep="_",rep(seq(LG), length.out = n*LG,each=n))
    Seeds_Genes_Scores <- rep(((1-eta) * tau)/n,n)
    
    Seed_phenotypes_Layer_Labeled <- paste0(rep(phenotypeSeeds,LP), sep="_",rep(seq(LP), length.out = m*LP,each=m))
    Seeds_phenotypes_Scores <- rep((eta * phi)/m,m)    
    
  } else {
    eta <- 1
    if (n == 0){
      Seed_Genes_Layer_Labeled <- character()
      Seeds_Genes_Scores <- numeric()
      
      Seed_phenotypes_Layer_Labeled <- paste0(rep(phenotypeSeeds,LP), sep="_",rep(seq(LP), length.out = m*LP,each=m))
      Seeds_phenotypes_Scores <- rep((eta * phi)/m, m)      
    } else {
      Seed_Genes_Layer_Labeled <- paste0(rep(GeneSeeds,LG), sep="_",rep(seq(LG), length.out = n*LG,each=n))
      Seeds_Genes_Scores <- rep(tau/n, n)            
      
      Seed_phenotypes_Layer_Labeled <- character()
      Seeds_phenotypes_Scores <- numeric()        
    }
  }
  
  ### We prepare a data frame with the seeds.
  Seeds_Score <- data.frame(Seeds_ID = c(Seed_Genes_Layer_Labeled, Seed_phenotypes_Layer_Labeled),
                            Score = c(Seeds_Genes_Scores, Seeds_phenotypes_Scores)  ,stringsAsFactors = FALSE)
  return(Seeds_Score)  
}  
rank_genes <- function(Number_Genes, Number_Layers,Results,Seeds){
  ## We sort the score to obtain the ranking of Genes and Diseases.
  genes_rank <- data.frame(Node = character(length = Number_Genes), Score = 0)
  genes_rank$Node <- gsub("_1", "", row.names(Results)[1:Number_Genes])
  
  ## We calculate the Geometric Mean among the genes in the different layers.
  genes_rank$Score <- geometric.mean(as.vector(Results[,1]),Number_Layers,Number_Genes)
  
  genes_rank_sort <- genes_rank[with(genes_rank, order(-Score, Node)), ]
  
  ### We remove the seed genes from the Ranking
  genes_rank_sort_NoSeeds <- genes_rank_sort[which(!genes_rank_sort$Node %in% Seeds),]
  
  genes_rank_sort_NoSeeds$Rank <- seq(1, nrow(genes_rank_sort_NoSeeds))
  genes_rank_sort_NoSeeds <-genes_rank_sort_NoSeeds[,c("Rank","Node","Score")]   
  
  return(genes_rank_sort_NoSeeds)
}
rank_phenotypes <- function(Number_Genes,Num_Gene_Layers,Number_Diseases,Num_phenotype_Layers, Results,Seeds){
  
  ## rank_diseases
  diseases_rank <- data.frame(Node = character(length = Number_Diseases), Score = 0)
  diseases_rank$Node <- gsub("_1", "", row.names(Results)[(Number_Genes*Num_Gene_Layers+1):(Number_Genes*Num_Gene_Layers+Number_Diseases)])
  
  diseases_rank$Score <- geometric.mean(as.vector(Results[,1])[(Number_Genes*Num_Gene_Layers+1):nrow(Results)],Num_phenotype_Layers,Number_Diseases)
  
  diseases_rank_sort <- diseases_rank[with(diseases_rank, order(-Score, Node)), ]
  diseases_rank_sort_NoSeeds <- diseases_rank_sort[which(!diseases_rank_sort$Node %in% Seeds),]
  
  diseases_rank_sort_NoSeeds$Rank <- seq(1, nrow(diseases_rank_sort_NoSeeds))
  diseases_rank_sort_NoSeeds <-diseases_rank_sort_NoSeeds[,c("Rank","Node","Score")]   
  
  return(diseases_rank_sort_NoSeeds)
}
geometric.mean <- function(Scores, L, N) {
  
  FinalScore <- numeric(length = N)
  
  for (i in seq_len(N)){
    FinalScore[i] <- prod(Scores[seq(from = i, to = N*L, by=N)])^(1/L)
  }
  
  return(FinalScore)
}
Random_Walk_Restarts <- function(Walk_Matrix, r, Seeds_Score){
  ### We define the threshold and the number maximum of iterations for the randon walker.
  Threshold <- 1e-10
  if(r==0) Threshold <- 1e-6
  NetworkSize <- ncol(Walk_Matrix)
  
  ### We initialize the variables to control the flux in the RW algo.
  residue <- 1
  iter <- 1
  
  #### We define the prox_vector(The vector we will move after the first RW iteration. We start from The seed. We have to take in account
  #### that the walker with restart in some of the Seed genes, depending on the score we gave in that file).
  prox_vector <- Matrix(0,nrow = NetworkSize,ncol=1, sparse=TRUE)
  
  prox_vector[which(colnames(Walk_Matrix) %in% Seeds_Score[,1])] <- (Seeds_Score[,2])
  
  prox_vector  <- prox_vector/sum(prox_vector)
  restart_vector <-  prox_vector
  while(residue >= Threshold){
    old_prox_vector <- prox_vector
    prox_vector <- (1-r)*(Walk_Matrix %*% prox_vector) + r*restart_vector
    
    residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
    
    iter <- iter + 1; 
    #if(iter %% 10==0) print(residue)
    
  }
  
   # print("RWR-MH number of iteration: ")
   # print(iter-1)  
  return(prox_vector)
}
run.RWR.AllNodes <- function(WM, num.cores){
    genes <- WM$genes
    phenotypes <- WM$phenotypes 
    LG <- WM$LG
    LP <- WM$LP
    N <- WM$N
    M <- WM$M
    
    WM <- WM$WM
    
    t0 <- Sys.time()
 
    nodes <- c(genes, phenotypes)
    params <- get.default.parameters(LG, LP)

    cat("Number of nodes to be processed: ", length(nodes), "\n")
    cl <- makeCluster(num.cores)
    registerDoParallel(cl, cores = num.cores)    
    
    rankDF <- foreach(seed = nodes, .export=as.vector(lsf.str(.GlobalEnv)),
                  .packages=c('Matrix', 'data.table')) %dopar% {
        
        if (startsWith(seed, 'g')){ # if seed is gene
          GeneSeeds <- c(seed)
          phenotypeSeeds <- vector()
          
        }else{ # if seed is phenotype
          GeneSeeds <- vector()
          phenotypeSeeds <- c(seed)
        }
 
        Seeds_Score <- get.seed.scores(GeneSeeds, phenotypeSeeds, params$eta, LG, LP, 
                                       params$tau/LG, params$phi/LP)
                      
        if (sum(Seeds_Score$Score)!=1) stop("ERROR: Seeds Scores don't add up to 1!")
        
        # Run RWR
        Random_Walk_Results <- Random_Walk_Restarts(WM, params$r, Seeds_Score)
    
        RWGeneRankDF <- rank_genes(N, LG, Random_Walk_Results, GeneSeeds)
        RWphenotypeRankDF <- rank_phenotypes(N, LG, M, LP, Random_Walk_Results, phenotypeSeeds)
        dfRank <- rbindlist(list(RWGeneRankDF, RWphenotypeRankDF))
        
        dfRank <- dfRank[, c(2,3)]
        colnames(dfRank) <- paste0(c(seed, seed), c(".n", ".s"))
        
        return(as.data.frame(t(dfRank), stringsAsFactors=FALSE))
    }
    stopCluster(cl)
    rw_names <- sapply(rankDF, function(x){rownames(x)})
    
    rankDF <- rbindlist(rankDF, use.names = F)
    rownames(rankDF) <- rw_names  
    cat("Time to generate neighborhood file : ", format(Sys.time()-t0), "\n") 
    
    return(rankDF)
}
