library(tidyverse)
library(R.utils)
library(DiagrammeR)
library(igraph)

getUnobservedInts <- function(df, desiredResponsesMask, boolLen, str4true) {
    # This function first turns the dataframe of parameter-sweep results 
    # (i.e. observed responses), written in 'pos' and 'neg' format, into 
    # a vector of binary values (each row is represented by a binary value), 
    # and then convert the binary values into a list of integers 
    # which represents the observed integers. After discarding observed integers 
    # from all possible integers, the funciton returns a list of unobserved integers.
  
    observedInts <- 
      df %>%
      select(., desiredResponsesMask) %>%   # select desired responses
      apply(., 2, function(x) ifelse(x == str4true, "1", "0")) %>%  # pos as 1 & neg as 0
      as.data.frame() %>%
      apply(., 1, paste0, collapse="") %>%  # binary strings representing responses
      strtoi(., base = 2L)  # convert the binary string into integers
    
    unobservedInts <-
      c(0: (2**boolLen-1)) %>%  # all possible responses as integers 
      discard(., . %in% observedInts)   # discard observed
    
    return(unobservedInts)
}

getUnobservedInts2 <- function(df, desiredResponsesMask, boolLen) {
  
    # This function is similar to getUnobservedInts, but takes the dataframe
    # of parameter-sweep results written in '0' and '1' format.
    
    observedInts <- 
      df %>%
      select(., desiredResponsesMask) %>%   # select desired responses
      as.data.frame() %>%
      apply(., 1, paste0, collapse="") %>%  # binary strings representing responses
      strtoi(., base = 2L)  # convert the binary string into integers
    
    unobservedInts <-
      c(0: (2**boolLen-1)) %>%  # all possible responses as integers 
      discard(., . %in% observedInts)   # discard observed
    
    return(unobservedInts)
}

getUnobservedBooldf <- function(unobservedInts, desiredResponses){
    # The purpose of this function is to turn a list of integers into a
    # dataframe in Boolean formats (species responses written as 0s and 1s)
    # with an additional output column that meets the requirements of logicopt()
    # for Boolean minimization.
  
    unobservedBooldf <-
        append((2**length(desiredResponses))-1, unobservedInts) %>%
        intToBin(.) %>%   # convert the integer back into binary strings
        str_split(., "") %>%  # chop strings to 0s and 1s
        do.call(rbind, .) %>%
        as.data.frame() %>%
        setNames(., desiredResponses) %>%
        slice(., -1) %>%
        mutate(., unob = "1") # add an output column to meet the requirements of logicopt()
  
    return(unobservedBooldf)
}

getPCUList <- function(optEqn, str4true, str4flase, desiredResponses){
  
    # This function converts the optimized equations (i.e. results of Boolean minimization) 
    # into a list of PCUs.
    
    # split the string of equation and get PCUs as a list.
    unobservedList <- 
        str_split(optEqn, " [=] ", simplify = TRUE)[2] %>%
        str_split(., " [+] ") %>%
        unlist(.) %>%
        str_split(., "[*]")
      
    # convert the binary variabe into descriptive strings, e.g.
    # RABBITS_ALBATROSSES->"posrabbits_albatrosses" and
    # rabbits_albatrosses->"negrabbits_albatrosses"
    
    dictionary <- setNames(
        append(paste0(str4true, desiredResponses), paste0(str4flase, desiredResponses)),
        append(str_to_upper(desiredResponses), str_to_lower(desiredResponses))
    ) # create a named vector as dictionary
    
    PCUList_unordered <- 
        unobservedList %>%
        lapply(., function(i) dictionary[i]) %>% # match and replace use dictionary
        lapply(., function(i) unname(i)) # unname each list
    
    PCUList <- PCUList_unordered[order(sapply(PCUList_unordered,length))] # sort PCUs
    
    return(PCUList) # show the PCUList
}


get_edgelist_singleAnte <- function(PCUList) { 
    
    # Get the edgelist for the simplest kind of implication network, 
    # where every logical implication statement in its single-antecedent form is included.

    andCounter = 0
    orCounter = 0
    edgesList = list()
    
    for (PCU in PCUList[]){
      
        for (p in PCU[]){
          
            if (length(PCU) == 1) {
              
                ante = 'True'
                
                notpp = ifelse(grepl('pos', p), str_replace(p, 'pos', 'neg'),
                               str_replace(p, 'neg', 'pos'))
                cons = notpp
              
            } else {
              
                ante = p
                cons = NULL
                
                qs = PCU[!PCU %in% p]
                
                for (q in qs) {
                  
                    notqq = ifelse(grepl('pos', q), str_replace(q, 'pos', 'neg'),
                                   str_replace(q, 'neg', 'pos'))
                    cons = append(cons, notqq)
                }
            }
            
            # Create the edges list
            
            if (length(cons) == 1) {
              
                edgesList[[length(edgesList)+1]] = append(ante, cons)
              
            } else {
              
                orNode = paste('or', as.character(orCounter), sep = "")
                
                edgesList[[length(edgesList)+1]] = append(ante, orNode)
                
                orCounter = orCounter + 1
              
                for (c in cons) {
                    edgesList[[length(edgesList)+1]] = append(orNode, c)
                }
            }
        }
    }
    return(edgesList)
}

get_edgelist_certainAnte <- function(PCUList, alwaysAnteList) {
    
    # Get the edgelist for the implication network where certain responses are 
    # always treated as the antecedent (alwaysAnteList)

    andCounter = 0
    orCounter = 0
    edgesList = list()
    
    for (PCU in PCUList[]){
      
        # split into a list of antecedents and consequents
        anteList = list()
        consList_raw = list()
        for(p in PCU[]) {
            ifelse(p %in% alwaysAnteList, anteList[[length(anteList)+1]] <- p, 
                  consList_raw[[length(consList_raw)+1]] <- p)
        }
        
        # The consequents need to be negated
        consList = list()
        for (q in consList_raw) {
            notqq = ifelse(grepl('pos', q), str_replace(q, 'pos', 'neg'),
                           str_replace(q, 'neg', 'pos'))
            consList = append(consList, notqq)
        }
        
        # if there are no antecedents, the consequents will be attached to the True node
        if (length(anteList) == 0) {anteList = append(anteList, 'True')}
        # if there are no consequents, the antecedents will be attached to the False node
        if (length(consList) == 0) {consList = append(consList, 'False')}
        
        # identify the upper and lower node, which depends on how many antecedents and consequents
        # do upper and antecedents
        if (length(anteList) < 2) {
          
            upperNode <- anteList[[1]]
        
        } else {
            
            # the upper node is an And node
            andCounter = andCounter + 1
            upperNode <- paste('and', as.character(andCounter), sep = "")
            
            # and every antecedent has an edge with this And node
            for (p in anteList) {
              edgesList[[length(edgesList)+1]] = append(p, upperNode)
            }
        }
        
        # do lower and consequents
        if (length(consList) < 2) {
          
            lowerNode = consList[[1]]
          
        } else {
          
            # the lower node is an Or node
            orCounter = orCounter + 1
            lowerNode <- paste('or', as.character(orCounter), sep = "")
            
            # and every consequent has an edge with this Or node
            for (q in consList) {
              edgesList[[length(edgesList)+1]] = append(lowerNode, q)
            }
        }
        # link the upper and lower node
        edgesList[[length(edgesList)+1]] = append(upperNode, lowerNode)
    }
    return(edgesList)
}

draw_implication_network <- function(edgesList, niceNames){
  
    # Draw the implication nework from the edgelist get from the funtion 
    # get_edgelist_singleAnte() or get_edgelist_certainAnte().
      
    controlSymbol = '&darr; ';
  
    # Use edgelist to create an igraph object and get the nodelist from the igraph object
    edges_list <- do.call(rbind, edgesList)
    nodes_list <- 
        edges_list %>% graph_from_edgelist(., directed = TRUE) %>%
        get.vertex.attribute(.) %>% get("name", .)
    
    # set attributes for response nodes
    nodes_df_resp <- 
        # get response nodes
        data.frame(names = nodes_list) %>%
        mutate_if(is.factor, as.character) %>%
        filter(grepl('pos|neg', names)) %>% 
      
        # set sign, fillcolor and node shape 
        mutate(respSign = case_when(grepl('pos', names) ~ "+", TRUE ~ "&#8210;"),  
               fillcolor = case_when(grepl('pos', names) ~ "white", TRUE ~ "gray"), 
               shape = "box") %>% 
        # separate the strings into two parts: control and response species
        mutate(dropSign = substring(names, 4)) %>%
        separate(dropSign, c("contSpp", "respSpp"), "_") 
    
    # check if niceNames are defined
    if (hasArg(niceNames)){      
        contSppNiceName = unname(niceNames[nodes_df_resp$contSpp])
        respSppNiceName = unname(niceNames[nodes_df_resp$respSpp])
    } else {
        contSppNiceName = nodes_df_resp$contSpp
        respSppNiceName = nodes_df_resp$respSpp
    }
    
    # set labels for response nodes
    nodes_df_resp <- nodes_df_resp %>%
        mutate(contSppNiceName = contSppNiceName,
               respSppNiceName = respSppNiceName) %>%
        mutate(label = paste0('< <font point-size="10">', controlSymbol, contSppNiceName, '</font>',
                              '<br align="left"/> &nbsp; &nbsp; ', respSppNiceName, ' ', respSign, ' >')) %>%
        
        # select node attribute columns 
        select(names, shape, fillcolor, label) 
    
    # set attributes for Boolean nodes and
    # combine response node dataframe and Boolean node dataframe into one dataframe
    nodes_df <- 
        # get Boolean nodes
        data.frame(names = nodes_list) %>%
        mutate_if(is.factor, as.character) %>%
        filter(!grepl('pos|neg', names)) %>% 
        
         # set attributes for Boolean nodes 
        mutate(shape = "circle",
               fillcolor = "white",
               label = case_when(grepl('or', names) ~ "or", 
                                 grepl('and', names) ~ "\"&\"", #NOTE!!!!!!!!!
                                 grepl('False', names) ~ "False", 
                                 TRUE ~ "True")) %>%
         
        # combine the two dataframe
        bind_rows(., nodes_df_resp)  
    
    ## construct NDF and EDF for DiagrammeR   
    ndf <- create_node_df(
        n = nrow(nodes_df),
        names = nodes_df[, "names"],
        label = nodes_df[ , "label"],
        shape = nodes_df[, "shape"],
        fillcolor = nodes_df[, "fillcolor"])
    
    
    edges_df <- edges_list %>%
        as.data.frame() %>% setNames(., c("labelfrom", "labelto")) %>%
        mutate_if(is.factor, as.character) %>%
        left_join(., ndf[, c("id", "names")], by = c("labelfrom"="names")) %>% rename(from = id) %>%
        left_join(., ndf[, c("id", "names")], by = c("labelto"="names")) %>% rename(to = id)
    
    edf <- create_edge_df(
        from = edges_df[ , "from"],
        to = edges_df[ , "to"])
    
    G <- create_graph(
        nodes_df = ndf,
        edges_df = edf,
        directed = TRUE,
        attr_theme = NULL) %>% 
        add_global_graph_attrs(
          attr = "style",
          value = '\"rounded, filled\"',
          attr_type = "node") %>%
        add_global_graph_attrs(
          attr = "width",
          value = 0,
          attr_type = "node") %>%
        add_global_graph_attrs(
          attr = "margin",
          value = 0,
          attr_type = "node")
      
      dot <- gsub("\\'","",generate_dot(G))
      
      DiagrammeR(diagram = dot, type = "grViz")
}

