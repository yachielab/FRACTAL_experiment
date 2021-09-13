library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggtree)
library(tidyr)
library(igraph)
library(ape)
library(treeio)
library(hash)
library(aplot)
library(ggthemes)
library(ggforce)

get_root_node <- function(tree_tbl){
  tmp <- subset(tree_tbl,parent==node)
  root_node=tmp$node
  return(root_node)
}

get_upstream_tree2 <- function(tree_tbl,lev,parent_vector,out_df=data.frame(),cnt=0){
  if(cnt>=lev){
    return(out_df)
  }else{
    cnt <- cnt+1
    out_df <- rbind(out_df,subset(tree_tbl,parent %in% parent_vector)) %>% unique
    return(Recall(tree_tbl,lev,out_df$node,out_df=out_df,cnt))
  }
}

gen_parent_to_child <- function(tree_tbl){
  tree_tbl <- filter(tree_tbl,parent!=node)
  tree_tbl$parent <- as.character(tree_tbl$parent)
  tree_tbl$node <- as.character(tree_tbl$node)
  parent_to_child_hash <- split(tree_tbl$node,tree_tbl$parent)
  parent_to_child_hash <- hash(parent_to_child_hash)
  return(parent_to_child_hash)
}

get_upstream_tree3 <- function(tree_tbl,node_now,parent_to_child,node_to_n_tips,thresh_n_tips,out_tbl,todo_vector,visited=c()){
  if(length(todo_vector)==0){
    return(out_tbl)
  }
  
  visited <- c(node_now,visited)
  if(node_to_n_tips[[node_now]]<=thresh_n_tips){
    next_node <- todo_vector[1]
    todo_vector <- todo_vector[-1]
    out_tbl <- rbind(out_tbl,filter(tree_tbl,node==next_node))
    return(Recall(tree_tbl,next_node,parent_to_child,node_to_n_tips,thresh_n_tips,out_tbl,todo_vector=todo_vector,visited=visited))
  }else{
    childs <- parent_to_child[[node_now]]
    todo_vector <- c(childs,todo_vector)
    for(child_node in todo_vector){
      todo_vector <- todo_vector[-1]
      out_tbl <- rbind(out_tbl,filter(tree_tbl,node==child_node))
      return(Recall(tree_tbl,child_node,parent_to_child,node_to_n_tips,thresh_n_tips,out_tbl,todo_vector=todo_vector,visited=visited))
    }
  }
  
}

get_upstream_tree_dfs <- function(tree_tbl,root_node,parent_to_child,node_to_n_tips,thresh_n_tips){
  node_list <- c()
  stack_list <- c(as.character(root_node))
  
  print("Starting Depth First Search...")
  starting.time <- Sys.time()
  while(TRUE){
    if(length(stack_list)==0){break}
    
    node_now <- stack_list[1]
    child_nodes <- parent_to_child[[node_now]]
    stack_list <- stack_list[-1]
    node_list <- c(node_now,node_list)
    
    if(node_to_n_tips[[node_now]]>thresh_n_tips){
      stack_list <- c(child_nodes,stack_list)
    }
  }
  print(Sys.time()-starting.time)
  
  return(node_list)
}

fnc_dfs_Ntips <- function(tree_data,root_node,vec_ntips,node_completed=c(),reverse=F){
  root_node <- as.character(root_node)
  
  if(length(node_completed)==(Ntip(tree_data)+Nnode(tree_data))){
    return(vec_ntips)
  }
  
  if(length(child(tree_data,root_node))==0){
    vec_ntips[root_node]=1
    node_completed <- c(node_completed,root_node)
    return(Recall(tree_data,root_node=parent(tree_data,root_node),vec_ntips,node_completed=node_completed,reverse=T))
    
  }else{
    child_nodes <- child(tree_data,root_node) %>% as.character
    
    if(all(!is.na(vec_ntips[child_nodes])) && !is.null(vec_ntips) && reverse==T){
      vec_ntips[root_node] <- sum(vec_ntips[child_nodes])
      node_completed <- c(node_completed,root_node)
      return(Recall(tree_data,root_node=parent(tree_data,root_node),vec_ntips,node_completed=node_completed,reverse=T))
      
    }else{
      for(node in child_nodes){
        if(!node %in% node_completed){return(Recall(tree_data,root_node=node,vec_ntips,node_completed=node_completed,reverse=F))}
      }
    }
  }
}

pop_vec <- function(elem,vec){
  vec <- vec[!is.element(vec,elem)]
  return(vec)
}

fnc_dfs_Ntips_stack <- function(tree_data,root_node){
  print("Building edge hash...")
  edge_df <- tree_data$edge %>% as_tibble()
  edge_df$V1 <- as.character(edge_df$V1)
  edge_df$V2 <- as.character(edge_df$V2)
  parent_to_child_hash <- split(edge_df$V2,edge_df$V1)
  parent_to_child_hash <- hash(parent_to_child_hash)
  print("done.")
  stack_list <- c(as.character(root_node))
  cnt <- 0
  
  print("Starting Depth First Search...")
  starting.time <- Sys.time()
  dfs_list <- character(Nnode(tree_data))
  while(TRUE){
    if(length(stack_list)==0){break}
    
    node_now <- stack_list[1]
    child_nodes <- parent_to_child_hash[[node_now]]
    stack_list <- stack_list[-1]
    
    if(length(child_nodes)!=0){
      cnt <- cnt+1
      dfs_list[cnt] <- node_now
      stack_list <- c(child_nodes,stack_list)
    }
  }
  dfs_list <- rev(dfs_list)
  print(Sys.time()-starting.time)
  
  print("Starting node to Ntips hash generation...")
  tips <- edge_df$V2[! edge_df$V2 %in% edge_df$V1]
  num_tips <- hash()
  num_tips[tips] <- 1
  starting.time <- Sys.time()
  for(node_now in dfs_list){
    child_nodes <- parent_to_child_hash[[node_now]]
    if(length(child_nodes)==2){
      num_tips[node_now] <- num_tips[[child_nodes[1]]] + num_tips[[child_nodes[2]]]
    }else{
      num_tips[node_now] <- 0
      for(i in child_nodes){
        num_tips[node_now] <- num_tips[[i]] + num_tips[[node_now]]
      }
    }
  }
  print(Sys.time()-starting.time)
  return(num_tips)
}

fill_label <- function(tree_tbl){
  tree_tbl$label[!tree_tbl$node %in% tree_tbl$parent] <- paste0("Clade",1:length(tree_tbl$label[!tree_tbl$node %in% tree_tbl$parent]))
  tree_tbl$label[tree_tbl$node %in% tree_tbl$parent] <- NA
  tree_tbl$branch.length[tree_tbl$node==tree_tbl$parent] <- NA
  tree_tbl <- tree_tbl[order(tree_tbl$label),]
  return(tree_tbl)
}

reassign_node <- function(tree_tbl){
  conv_vec <- 1:nrow(tree_tbl)
  names(conv_vec) <- tree_tbl$node
  
  tree_tbl$node <- 1:nrow(tree_tbl)
  tree_tbl$parent <- conv_vec[as.character(tree_tbl$parent)]
  
  ret_vec <- as.numeric(names(conv_vec))
  names(ret_vec) <- conv_vec
  
  L <- list()
  L$tree <- tree_tbl
  L$conversion_vector <- ret_vec
  
  return(L)
}

get_original_node_to_clade <- function(tree_list){
  df_use <- na.omit(tree_list$tree)
  clade_vec <- df_use$label
  names(clade_vec) <- tree_list$conversion_vector[as.character(df_use$node)]
  return(clade_vec)
}

gen_vector_root_to_child <- function(anc_df,root_node,tree_tbl,node_now){
  final_vec <- c(root_node)
  n_vec <- anc_df$node; names(n_vec) <- anc_df$parent
  
  p <- root_node
  if(nrow(anc_df)!=0){
    while(TRUE){
      chld <- n_vec[as.character(p)] %>% unname
      final_vec <- c(final_vec,chld)
      if(is.na(n_vec[as.character(chld)])){
        break
      }
      p <- chld
    }
    
    if(chld!=parent(tree_tbl,node_now)$node){
      stop()
    }
  }
  
  return(final_vec)
}

gen_ancestor_hash <- function(node_now,tree_tbl,root_node){
  # tree_tbl <- filter(tree_tbl,parent!=node)
  # print("hoge")
  anc_df <- ancestor(tree_tbl,node_now)
  anc_df <- anc_df[anc_df$parent!=anc_df$node,]
  v <- gen_vector_root_to_child(anc_df,root_node,tree_tbl,node_now)
  return(rev(v))
}

path_to_root <- function(tree_data,root_node){
  print("Building edge hash...")
  edge_df <- tree_data$edge %>% as_tibble()
  edge_df$V1 <- as.character(edge_df$V1)
  edge_df$V2 <- as.character(edge_df$V2)
  child_to_parent_hash <- hash()
  child_to_parent_hash[edge_df$V2] <- edge_df$V1
  print("done.")
  
  tips <- edge_df$V2[! edge_df$V2 %in% edge_df$V1]
  
  out_hash <- hash()
  cnt <- 0
  for(t in tips){
    cnt <- cnt+1
    out_v <- c()
    node <- t
    while(T){
      out_v <- c(out_v,node)
      if(length(child_to_parent_hash[[node]])==0){
        break
      }else{
        node <- child_to_parent_hash[[node]]
      }
    }
    out_hash[paste0("tip",t)] <- out_v
    if(cnt%%10000==0){
      print(paste(cnt,"nodes were processed."))
    }
  }
  return(out_hash)
}

get_mrca <- function(tip1,tip2,hash){
  tip1 <- paste0("tip",tip1)
  tip2 <- paste0("tip",tip2)
  
  shared_branch <- intersect(hash[[tip1]],hash[[tip2]])
  mrca <- shared_branch[1]
  return(mrca)
}

dist_norm_mrca <- function(tip,mrca,hash){
  tip <- paste0("tip",tip)
  
  tip_to_root <- hash[[tip]]
  # dist_mrca <- 0
  # for(i in tip_to_root){
  #   if(i==mrca){
  #     break
  #   }
  #   dist_mrca <- dist_mrca+1
  # }
  dist_mrca <- which(tip_to_root==mrca)-1
  dist_norm <- dist_mrca/length(tip_to_root)
  # dist_norm <- dist_mrca
  return(dist_norm)
}

calc_dist <- function(tip_merged,hash){
  tip1 <- strsplit(tip_merged,":")[[1]][1]
  tip2 <- strsplit(tip_merged,":")[[1]][2]
  
  mrca <- get_mrca(tip1,tip2,hash)
  dist_norm1 <- dist_norm_mrca(tip1,mrca,hash)
  dist_norm2 <- dist_norm_mrca(tip2,mrca,hash)
  
  return(dist_norm1+dist_norm2)
}

calc_dist2 <- function(t,hash){
  tip1 <- strsplit(tip_merged,":")[[1]][1]
  tip2 <- strsplit(tip_merged,":")[[1]][2]
  
  mrca <- get_mrca(tip1,tip2,hash)
  dist_norm1 <- dist_norm_mrca(tip1,mrca,hash)
  dist_norm2 <- dist_norm_mrca(tip2,mrca,hash)
  
  return(dist_norm1+dist_norm2)
}

calc_dist_wrapper <- function(v,hash){
  v <- sapply(v,calc_dist,hash=hash)
  # return(mean(v))
  return(v)
}

calc_dist_ntips <- function(tip_merged,tip_to_root,node_to_ntips){
  tip1 <- strsplit(tip_merged,":")[[1]][1]
  tip2 <- strsplit(tip_merged,":")[[1]][2]
  
  mrca <- get_mrca(tip1,tip2,tip_to_root)
  dist_ntips <- node_to_ntips[[mrca]]
  return(dist_ntips)
}


eq_return <- function(x){return(x)}

fnc_entropy <- function(x){
  E <- 0
  for(i in x){if(i!=0){E <- E+i*log2(i)}}
  return(-1*E)
}