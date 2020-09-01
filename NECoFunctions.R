suppressMessages(library(igraph))
library(Matrix)
library(data.table)
suppressMessages(library(foreach))
suppressMessages(library(doParallel))

source("RWFunctions.R")
source("WMFunctions.R")
generateNbhd <- function(inputFileName, numCores){
    WM <- create.WalkMatrix(inputFileName, numCores)
    nbhd <- run.RWR.AllNodes(WM, numCores)
    return(nbhd)
}
writeDFWithNaToFile <- function(df, f_name){
    lst <- lapply(split(df, row.names(df)), function(x){ v <- unname(unlist(x)); v[!is.na(v)] })
    
    df <- data.frame(a=I(unlist(lapply(lst, paste,collapse="\t"))))
    write.table(df, file=f_name, quote=FALSE, row.names=FALSE, col.names = FALSE)
}
trimTopNbhd <- function(dfNbhd, file_name, nbhd_types = c("Gg", "Gp","Pg", "Pp"), 
                        top_Ggs=c(150, 250), top_Gps=c(150, 250), 
                        top_Pgs=c(150, 250), top_Pps=c(150, 250)){
    
    t0 <- Sys.time()

    # to preserve the row names of node ranks
    # also remove the ".n" part from row name to keep gene names only
    rw_names <- sapply(strsplit(rownames(dfNbhd)[seq(1, nrow(dfNbhd), 2)], "\\."), `[[`, 1)
    
    no_of_genes <- sum(startsWith(rw_names, "g"))
    
    no_of_phenotypes <- ncol(dfNbhd) - no_of_genes + 1
    cat("Number of genes in the nbhd file:", no_of_genes, "\n")
    cat("Number of phenotypes in the nbhd file:", no_of_phenotypes, "\n")
    
    dfNbhd <- as.data.frame(dfNbhd)
    
    dfNbhdScr <- dfNbhd[seq(2, nrow(dfNbhd), 2), ]
    dfNbhd <- dfNbhd[seq(1, nrow(dfNbhd), 2), ]
    
    # convert score df to numeric
    dfNbhdScr <- as.data.frame(lapply(dfNbhdScr, function(x) as.numeric(x)))
    
    # some rows (for unconnected nodes in the graph)  have 0 values after 2nd column
    # those ranks will be removed, because they create an unnecessary noise
    rws_zero <- unlist(lapply(1:nrow(dfNbhdScr), function(i){
    if(sum(dfNbhdScr[i,1:ifelse(ncol(dfNbhdScr)>100, 100, ncol(dfNbhdScr))]==0)>1){
      return(i)
    }
    }))
    
    dfNbhdScr[dfNbhdScr==0] <- NA
    
    # make the corresponding NA values in nbhd as NA
    dfNbhd[is.na(dfNbhdScr)] <- NA
    
    rws_zero_gene <- rws_zero[rws_zero <= no_of_genes]
    rws_zero_phenotypes <- rws_zero[rws_zero > no_of_genes]
    
    if (length(rws_zero) > 0){
        # cat("NA/zero valued rows after 2nd column:", length(rws_zero), "\n")
       
        dfNbhdScr <- dfNbhdScr[-rws_zero,]
        dfNbhd <- dfNbhd[-rws_zero,]
        rownames(dfNbhd) <- rw_names[-rws_zero]
        # add the seed nodes as first node for removed rws_zero case
        dfNbhd <- cbind(rw_names[-rws_zero], dfNbhd, stringsAsFactors =F)
        # duplicate the first column, since we have seeds in the first column of dfNbhd
        dfNbhdScr <- cbind(dfNbhdScr[,1], dfNbhdScr)
        
    }else{
        dfNbhd <- cbind(rw_names, dfNbhd, stringsAsFactors =F) 
        dfNbhdScr <- cbind(dfNbhdScr[,1], dfNbhdScr)
    }
    
    g_rws <- c(1:(no_of_genes - length(rws_zero_gene)))
    s_rws <- c((no_of_genes - length(rws_zero_gene) + 1):nrow(dfNbhd))
    
    cl <- makeCluster(length(nbhd_types))
    registerDoParallel(cl, cores = length(nbhd_types))   
    
    foreach(i = nbhd_types, .export=as.vector(lsf.str(.GlobalEnv))) %dopar% {    
        if(i=="Gg"){
            top_Ns <- top_Ggs
        }else if(i=="Gp"){
            top_Ns <- top_Gps
        }else if(i=="Pg"){
            top_Ns <- top_Pgs
        }else if(i=="Pp"){
            top_Ns <- top_Pps
        }
        
        for (no_top in top_Ns){
            list_nbhd <- list()
            list_topN <- list()
            list_topN[[i]] <- no_top
            
            # subset top_Gg
            if (i == "Gg"){
                gg_cols <- c(1:list_topN[[i]])
                list_nbhd[[i]] <- dfNbhd[g_rws, gg_cols]
            }else if (i == "Gp"){
                if (list_topN[[i]] > no_of_phenotypes) list_topN[[i]] <- no_of_phenotypes
                gp_cols <- c(1, (no_of_genes + 1):(no_of_genes + list_topN[[i]] - 1))
                list_nbhd[[i]]  <- dfNbhd[g_rws, gp_cols]
                gp_cols <- gp_cols[-1] # won't use the first column of gene for score
                gp_cols <- c(gp_cols[1], gp_cols) # instead duplicate the first score column of strain      
            }else if (i == "Pg"){
                sg_cols <- c(2:(list_topN[[i]] + 1)) # exlude first column which includes phenotypes
                list_nbhd[[i]] <- dfNbhd[s_rws, sg_cols]
            }else if (i == "Pp"){
                if (list_topN[[i]] > no_of_phenotypes) list_topN[[i]] <- no_of_phenotypes
                    
                pp_cols <- c(1, (no_of_genes + 1):(no_of_genes + list_topN[[i]]) - 1)
                list_nbhd[[i]] <- dfNbhd[s_rws, pp_cols]
            }
          
            nbhd_file_name <- paste0(file_name, "_", i, "_Top", list_topN[[i]], ".txt")
            writeDFWithNaToFile(list_nbhd[[i]], nbhd_file_name)
            cat(nbhd_file_name, " created!\n")
          
            rm(list_nbhd, list_topN) # remove the variables after done
        }
    }
    stopCluster(cl)
    cat("Time to generate top N neighborhoods: ", format(Sys.time()-t0), "\n") 
}
