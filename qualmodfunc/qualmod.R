library(tidyverse)
library(reshape2)
library(igraph)
library(DiagrammeR)

initialise_foodweb <- function(positive_edges_list, negative_edges_list){
    # Accepts list of positive edges and a list of negative edges. 
    # It returns a DiagrammeR graph object with positive edges 
    # having a 'color' green and 'sign' of +1, and negative edges
    # having a 'color' red and 'sign' -1.
  
    # convert the positive and negative edge list into dataframes
    neg_edges_df <- 
        negative_edges_list %>%
        melt() %>%
        setNames(., c("labelfrom", "labelto")) %>%
        mutate_if(is.factor, as.character)
    
    pos_edges_df <- 
        positive_edges_list %>%
        melt() %>%
        setNames(., c("labelfrom", "labelto")) %>%
        mutate_if(is.factor, as.character)
    
    # create an igraph object from above edges and get the names of vertices (spp_list)
    spp_list <-
        bind_rows(neg_edges_df, pos_edges_df) %>% 
        graph_from_data_frame(., directed = TRUE) %>%
        get.vertex.attribute(.) %>% 
        get("name", .) %>% 
        sort()
    
    # construct NDF from spp_list and set attributes for the nodes
    NDF <-
          create_node_df(
            n = length(spp_list),
            label = spp_list,
            color = "black",
            fillcolor = "none",
            shape = "ellipse"
          )
    
    dict <- setNames(seq(length(spp_list)), spp_list)
    
    # construct EDF from positive and negative edge dataframes and set atributes of the edges 
    neg_EDF <-
        create_edge_df(
          from = unname(dict[neg_edges_df$labelfrom]),
          to = unname(dict[neg_edges_df$labelto]),
          rel = "-1",
          color = "#990033",
          arrowhead = "dot"
        )
    
    pos_EDF <- 
        create_edge_df(
          from = unname(dict[pos_edges_df$labelfrom]),
          to = unname(dict[pos_edges_df$labelto]),
          rel = "1",
          color = "#40e0d0",
          arrowhead = "vee"
        )
    
    EDF <- combine_edfs(neg_EDF, pos_EDF)
    
    # create the DiagrammeR graph object
    G <- create_graph(
        nodes_df = NDF,
        edges_df = EDF,
        directed = TRUE,
        attr_theme = NULL)
    
    return(G)
}

qualitative_community_matrix <- function(web){
    # Takes a DiagrammeR graph of the food web, returns 
    # the qualitative community matrix, with positive and 
    # negative interactions marked by positive and negative 1. 
    # It also returns two named vectors mapping species labels 
    # to indices of matrix and vice versa.
  
    NDF <- get_node_df(web)
    EDF <- get_edge_df(web)
    
    length <- nrow(NDF)
    Mq <- matrix(0, nrow = length, ncol = length)
    
    to_idx <- EDF[,"to"]
    from_idx <- EDF[, "from"]
    Mq[cbind(to_idx, from_idx)] <- as.numeric(EDF[,"rel"])
    
    labelToIndex <- setNames(NDF[, "id"], NDF[, "label"])
    indexToLabel <- setNames(NDF[, "label"], as.character(NDF[, "id"]))
    
    result <- list(Mq=Mq, labelToIndex=labelToIndex, indexToLabel=indexToLabel)
    return(result)
}

get_condition_df <- function(labelToIndex, validation) { 
  
    # Takes validation criteria which is a list of species responses, and returns
    # a dataframe with three columns: species' names (label), their corresponding conditions 
    # and their corresponding indices in the qualitative community matrix.
  
    if (hasArg(validation)) {
      
      conditions_idx_df <- 
          validation %>%
          melt() %>%
          setNames(., c("conditions", "speciesNames")) %>%
          mutate_if(is.factor, as.character) %>% 
          mutate(idx = unname(labelToIndex[speciesNames]))
    }
    
    else {
      
      conditions_idx_df <- data.frame(
          conditions = numeric(),
          speciesNames = character(),
          idx = numeric()
      )
    }
    
    # returns a dataframe with each species' name (label), index and condition of the validation criteria
    return(conditions_idx_df)
}

foodweb_neat_plot <- function(G, title = "foodweb_neat") {
    
    # Returns a neat plot of the food web, which is a DiagrammeR graph object.
  
    edf <- get_edge_df(G) %>%
        filter(from != to) %>%
        mutate(smaller = pmin(from, to), larger = pmax(to, from))
    
    # subset unidirectional interactions
    df_uni <-
        edf[!(duplicated(edf[, c("smaller", "larger")]) | 
              duplicated(edf[, c("smaller", "larger")], fromLast = TRUE)), ]
    
    edge_df_uni <- create_edge_df(
        from = df_uni[ ,"from"],
        to = df_uni[ ,"to"],
        color = df_uni[ ,"color"],
        arrowhead = df_uni[ ,"arrowhead"],
        arrowtail = "none",
        dir = "both")
    
    # subset bi-directional interactions
    df_bi <- edf %>%
        select(smaller, larger) %>%
        duplicated(fromLast = TRUE) %>%
        subset(edf, .) %>%
        select(from, to, rel, arrowhead, smaller, larger)
    
    df_bi_dup <-
        edf %>%
        select(smaller, larger) %>%
        duplicated() %>%
        subset(edf, .) %>%
        select(from, to, rel, arrowhead, smaller, larger)
    
    df_both <-
        merge(df_bi, df_bi_dup, by = c("smaller", "larger")) %>%
        select(-smaller, -larger) %>%
        mutate(color = case_when(rel.x == "-1" & rel.y == "-1" ~ "red",
                                 rel.x == "1" & rel.y == "1" ~ "#006600",
                                 rel.x != rel.y ~ "#000099"))  %>%
        mutate(from = pmax(from.x, from.y)) %>%
        mutate(to = pmin(to.x, to.y)) %>%
        mutate(arrowhead = case_when(from == from.x ~ arrowhead.x,
                                     from == from.y ~ arrowhead.y)) %>%
        mutate(arrowtail = case_when(to == to.x ~ arrowhead.y,
                                     to == to.y ~ arrowhead.x))
    
    edge_df_bi <- create_edge_df(
        from = df_both[ ,"from"],
        to = df_both[ ,"to"],
        color = df_both[ ,"color"],
        arrowhead = df_both[ ,"arrowhead"],
        arrowtail = df_both[ ,"arrowtail"],
        dir = "both")
    
    edges_df <- combine_edfs(edge_df_bi, edge_df_uni)
    
    edges_df <- edges_df[order(edges_df$from, edges_df$to),] %>%
        mutate(id = 1: length(id))
    
    G.neat <- create_graph(
        nodes_df = get_node_df(G),
        edges_df = edges_df,
        directed = TRUE,
        attr_theme = NULL) 
    
    render_graph(G.neat, title = title)
}
