#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



library(BiocManager)
#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.15/bioc"))
options(repos = BiocManager::repositories())


#library(shiny)
library(DT)
library(DiagrammeR)
library(HydeNet)
library(igraph)


#first define local read and write filenames and functions  
#outputDir <- "DN_Model"
outputDir <- paste0("DN_Model", Sys.Date())
#outputDir <-  file.path(getwd(), outputDir)
if (!dir.exists(outputDir)) {dir.create(outputDir)}


#make a function to create freq tables
make_freq_tables <- function(BN_glob = BN_glob, outputDir){
  #now loop through xpar
  for(x in BN_glob$nodes){     #iterate over all  nodes
    #print(paste('node is: ', x))
    xpar = eval(parse(text = paste("BN_glob$parents$",x,sep = ''))) #get parents
    
    if (is.null(xpar) ) { #treat the root nodes
      xfact = eval(parse(text = paste("BN_glob$factorLevels$",x,sep = '')))
      print(paste('factors root node,', x, ' are: ', xfact))
      tbl <- eval(parse(text = paste('as.data.frame(xtabs(~',x,',data = BN_glob$data))',sep = '')))
      tbl$Freq <- tbl$Freq/sum(tbl$Freq) # make freq into a prob
      eval(parse(text = paste('levels(tbl$', x, ') <- BN_glob$factorLevels$', x ,sep = '')))
      print(tbl)
      fileName = paste0(x,".csv")
      fileDirName = file.path(outputDir, fileName)
      write.csv(tbl, fileDirName, row.names = FALSE)
      #also record table in global env
      assign(x, tbl, envir = .GlobalEnv) 
    }
    
    if (!is.null(xpar) ) { #treat the cpt nodes
      print(paste('parent are: ', xpar))
      cpt_tbl = eval(parse(text = paste('as.data.frame(ftable(cpt(BN_glob$nodeFormula$',x,', BN_glob$data)))')))
      eval(parse(text = paste('levels(cpt_tbl$', x, ') <- BN_glob$factorLevels$', x ,sep = '')))
      for(par in xpar){ 
        xfact = eval(parse(text = paste("BN_glob$factorLevels$",par,sep = '')))
        print(paste('factors of,', par, ' are: ', xfact))
        eval(parse(text = paste('levels(cpt_tbl$', par, ') <- BN_glob$factorLevels$', par ,sep = '')))
        
      }
      print(cpt_tbl)
      fileName = paste0(x,".csv")
      fileDirName = file.path(outputDir, fileName)
      write.csv(cpt_tbl, fileDirName, row.names = FALSE)
      #also record table in global env
      assign(x, cpt_tbl, envir = .GlobalEnv)
    }
  } #end loop
} #end function


#function to update model-list after updating CPTs
updateBN <- function(BN_glob, outputDir){ 
  new_df <- BN_glob$data #initialize new data with old data
  for(cn in colnames(new_df)){ #put in the level names for humans
    xfact = eval(parse(text = paste("BN_glob$factorLevels$",cn,sep = '')))
    #print(paste('factors of,', cn, ' are: ', list(xfact)))
    eval(parse(text = paste('levels(new_df$', cn, ') <- BN_glob$factorLevels$', cn ,sep = '')))
  }
  
  # now loop though and input make new data based on new probs and cndprobs
  
  un_finished_nodes <- BN_glob$nodes #list of undone nodes to cross off as they are processed
  while(length(un_finished_nodes) > 0){ #keep looping through unfinished nodes
    for(x in un_finished_nodes){     #iterate over all unfinished nodes
      #print(paste('beginning for loop child node is: ', x))
      xpar = eval(parse(text = paste("BN_glob$parents$",x,sep = ''))) #get parents
      #cat(paste('beginning for loop child node is: ', x, '\n its parents are: ', toString(xpar), '\n'))
      n=nrow(new_df) #no of data points -- keep the same number made in the skeleton
      
      if (is.null(xpar) ){ #treat the  root nodes unconditionally
        cat('made into root part \n')
        #read the new file ... need to generalize the storage folder
        # fileName = paste0(x,".csv")
        # fileDirName = file.path(outputDir, fileName)
        # eval(quote(dftest <- read.csv(fileDirName)))
        
        #get cpt and root pt directly from the global envir
        dftest <- eval(parse(text = paste("dftest <- ",x,sep = '')), envir = .GlobalEnv)
        print(dftest)
        N <- nrow(dftest) #number of factors
        #print(paste('there are ', N, 'levels'))
        
        txtFL <- eval(parse(text = paste("BN_glob$factorLevels$",x,sep = '')))
        #print(paste('its levels are ', toString(txtFL)))
        
        
        #now input the read probabilities into Probs to construct the DF
        Probs  <- dftest$Freq 
        print(paste('adding column: ',x))
        #txt = paste0("df <- cbind(df, ",x," <- as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
        txt = paste0("new_df$",x," <- as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs))")
        eval(parse(text = txt))#, envir = globalenv())
        eval(parse(text = paste('levels(new_df$',x, ') <- BN_glob$factorLevels$',x,sep = '')))
        
        # cross this node off the list
        un_finished_nodes = un_finished_nodes[!(un_finished_nodes %in% x)]   
        
        
        
      } # end of root treatments
      
      else{ #treat nodes that are not a root ( are conditioned or have  parents)
        
        if(any(xpar %in% un_finished_nodes)) {next} #only proceed if all parents are already in new_df
        cat('Made into conditional part. \n')
        
        # fileName = paste0(x,".csv")
        # fileDirName = file.path(outputDir, fileName)
        # eval(quote(dftest <- read.csv(fileDirName)))
        
        #get cpt and root pt directly from the global envir
        dftest <- eval(parse(text = paste("dftest <- ",x,sep = '')), envir = .GlobalEnv)
        cat("dftest is: \n")
        print(dftest)
        
        
        #***** calc cond prob fucntion
        #cat(paste('parents are: ', list(xpar), '\n'))
        # for each par combo of a node, e.g. payoff has (decision, child) for parents ...
        # i need to subset the new_df for each combinations of values of the parents and find
        # the correspond probs of the child in dftest and call multinomial for these probs
        # to choose values for the child for each combo of the parents values.
        
        #1) get subsetting argument for the probs
        arguPR = ""
        for (i in xpar) {  arguPR = paste0(arguPR,paste0("dftest$",i," ,"))}
        
        #2) get subsetting argument for the data table
        arguDT = ""
        for (i in xpar) {  arguDT = paste0(arguDT,paste0("new_df$",i," ,"))}
        
        #3) get the combo of parent levels
        txt = paste0("combos <- interaction(", arguPR, " sep = ',')")
        combos = levels(eval(parse(text = txt)))
        
        #4) find parents levels in the new_df 
        txtFL <- eval(parse(text = paste("BN_glob$factorLevels$",x,sep = '')))
        #print(paste('the child levels are ', toString(txtFL)))
        N <- length(unlist(txtFL))
        #print(paste('there are ', N, 'child levels'))
        
        
        arguPR = unlist(strsplit(arguPR, ','))
        arguDT = unlist(strsplit(arguDT, ','))
        
        for(j in strsplit(as.character(combos), ' ')){ #loop over combos
          yy = unlist(strsplit(j[1], ',')) #split to a list
          #5) get probs from input cpt 
          arguPR2 = ""
          for(k in range(1,length(arguPR)))
          {arguPR2 = paste0(arguPR2
                            ,paste0(eval(parse(text = paste0('arguPR[',k,']'))),"== '"
                                    , eval(parse(text = paste0('yy[',k,']')))
                                    ,"' &"))
          }
          arguPR2 = substring(arguPR2, 1, nchar(arguPR2)-1) #get rid of last &
          txt2 <- paste('Probs <- dftest[',arguPR2, ' ,]$Freq')
          Probs = eval(parse(text = txt2))
          print(paste('Probs are: ', list(Probs)))
          
          #6) use Probs from input CPT to remake child values for each subset of combo levels
          arguDT2 = ""
          for(k in range(1,length(arguDT)))
          {arguDT2 = paste0(arguDT2
                            ,paste0(eval(parse(text = paste0('arguDT[',k,']'))),"== '"
                                    , eval(parse(text = paste0('yy[',k,']')))
                                    ,"' &"))
          }
          arguDT2 = substring(arguDT2, 1, nchar(arguDT2)-1) #get rid of last &
          txt2 <- paste('Subset_DT <- new_df[',arguDT2, ' ,]')
          Subset_DT <- eval(parse(text = txt2))
          
          #7) generate new child values into subset table
          n = nrow(Subset_DT) #no. of values to generate
          txt = paste0("chld <- as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs))")
          chld <- eval(parse(text = txt))
          #put in factor leveLs
          eval(parse(text = paste('levels(chld) <- BN_glob$factorLevels$', x ,sep = '')))
          
          
          
          #7) put new child values into subset table and then full table
          txt = paste0("Subset_DT$",x," <- chld")
          eval(parse(text = txt))
          #cat("Subset_DT is: \n")
          #print(Subset_DT)
          #update the original dataframe... will be completed on last pass
          new_df[rownames(Subset_DT),] <- Subset_DT
          
        } #end of for loop over subset combos to update cpt
        
        # cross this node off the list
        un_finished_nodes = un_finished_nodes[!(un_finished_nodes %in% x)]
        
      }   #end of conditioned node treatment
      
      
      
      
      
    } #end of for loop over un_finished_nodes
    cat('End of for loop. \n')
    cat(paste('unfinished nodes are: ', list(un_finished_nodes),'\n \n'))
    
  } # end of while loop when un_finished_nodes > 0
  
  #put back the original factor levels
  for(cn in colnames(new_df)){ 
    eval(parse(text = paste('levels(new_df$', cn, ') <- levels(BN_glob$data$', cn, ')' ,sep = '')))
  }
  
  return(new_df)
  
} #end of function


#make a function to create utility tables from previously created freq table i.e, combos of input already created
make_util_tables <- function(BN_glob = BN_glob, outputDir){
  #now loop through xpar
  for(x in BN_glob$nodes){     #iterate over all  nodes
    print(paste('node is: ', x))
    xpar = eval(parse(text = paste("BN_glob$parents$",x,sep = ''))) #get parents
    Util <- eval(parse(text = paste("BN_glob$nodeUtility$",x,sep = ''))) #test if node is a utility node
    
    if (!is.null(xpar) & Util) { #treat the utility nodes
      print(paste('parent are: ', xpar))
      #first read table written previously as a CPT 
      cat('Made into utility part. \n')
      # fileName = paste0(x,".csv")
      # fileDirName = file.path(outputDir, fileName)
      # eval(quote(util_tbl <- read.csv(fileDirName)))
      # 
      
      #get cpt and rooot pt directly from the global envir
      util_tbl <- eval(parse(text = paste("util_tbl <- ",x,sep = '')), envir = .GlobalEnv)
      
      cat("util_tbl is: \n")
      print(util_tbl)
      
      #AND NOW MODIFY THIS TABLE
      if("Freq" %in% colnames(util_tbl)){ #only update if it needs updating
        
        util_tbl = subset(util_tbl, select = -Freq) #delete the Freq column
        
        
        #the column name with the payoff
        P = noquote(colnames(util_tbl)[length(colnames(util_tbl))])
        # and cut the table to one level of the payoff
        util_tbl <- eval(parse(text = paste0("subset(util_tbl, ",P," == unique(subset(util_tbl,  select=",P,"))[1,1])")))
        #and now put a dummy value into the payoff
        eval(parse(text = paste0("util_tbl[,'",P,"'] = 1.0")))
        

        print(util_tbl)       
        # fileName = paste0(x,".csv")
        # fileDirName = file.path(outputDir, fileName)
        # write.csv(util_tbl, fileDirName, row.names = FALSE)
        #write payoff to global envir
        assign(x, util_tbl, envir = .GlobalEnv)
        print('yes')
      }
      else{
        print('you already updated the table')
      }
    }    
    
  } #end loop
} #end function

#FUNCTION TO MAKE UTILITY TABLE
updateDN <- function(BN_glob = BN_glob, outputDir){
  new_df <- BN_glob$data #initialize new data with old data
  for(cn in colnames(new_df)){ #put in the level names for humans
    xfact = eval(parse(text = paste("BN_glob$factorLevels$",cn,sep = '')))
    #print(paste('factors of,', cn, ' are: ', list(xfact)))
    eval(parse(text = paste('levels(new_df$', cn, ') <- BN_glob$factorLevels$', cn ,sep = '')))
  }
  
  #GET TABLE NAMES
  ute = names(BN_glob$nodeUtility[BN_glob$nodeUtility == TRUE])
  NON_ute = names(BN_glob$nodeUtility[BN_glob$nodeUtility == FALSE])
  
  # fileName = paste0(ute,".csv")
  # fileDirName = file.path(outputDir, fileName)
  # eval(quote(util_tbl <- read.csv(fileDirName)))
  
  #get cpt and root pt directly from the global envir
  util_tbl <- eval(parse(text = paste("util_tbl <- ",ute, sep = '')), envir = .GlobalEnv)
  
  
  #MAKE NEW MODEL LIST
  BN_glob2 <- HydeNetwork(BN_glob$network_formula
                          , data = new_df[ ,NON_ute])
  
  #MAKE THE UTILITY TABLE WITH NEW VALUES
  yy ='' 
  
  for (i in 1:length(util_tbl[[1]])) {
    for (j in 2:length(util_tbl)){
      if( j == 2){
        xx=paste0(" ifelse(",colnames(util_tbl)[j]," == '",util_tbl[i,noquote(paste0(colnames(util_tbl)[j]))],"'")
        #print(paste('j=2: ', xx))
      }
      else if(i != length(util_tbl[[1]]) && j == length(util_tbl)){
        xx=paste0(", ",util_tbl[i,noquote(paste0(colnames(util_tbl)[j]))],", ")
        #print(paste('last j non end: ', xx))
      }
      else if(i == length(util_tbl[[1]]) && j == length(util_tbl)){
        xx=paste0(", ",util_tbl[i,noquote(paste0(colnames(util_tbl)[j]))],", 0 ", paste(rep(")", length(util_tbl[[1]]) + 1), collapse = ""))
        #print(paste('last j end: ', xx))
      }
      
      else{
        xx=paste0(" && ",colnames(util_tbl)[j]," == '",util_tbl[i,noquote(paste0(colnames(util_tbl)[j]))],"'")
        #print(paste('middle j: ', xx))
      }
      yy <- paste(yy,xx)
      
    }
  }
  print(yy)
  
  pyoff =colnames(util_tbl)[length(util_tbl)]
  txt <- paste("BN_glob2 <- setNode(BN_glob2, ", pyoff, ", 'determ', define=fromFormula(), nodeFormula = ", pyoff, " ~ ",yy )
  print(txt)
  eval(parse(text = txt))
  assign("BN_glob2", BN_glob2, envir = .GlobalEnv)
  
} #END FUNCTION




#FUNCTION TO OUTPUT NEW RECONSTRUCED BN AND FAKE DATA DF
out <- function(ndata, skel, levls){ #INPUT NO. OF DATA ROWS AND RAW SKELETON and wheher the first iteration
  first = TRUE # set first time for 
  n=ndata
  make_df = TRUE #flag for making the first df
  fmula <- ''  # reconstruct a formula without errors
  fmula1 <- '' #the first formula without any '+'
  fmula_par <- ''
  fmula_ch <- ''
  for(x in skel$nodes){     #iterate over all nodes
    xpar = eval(parse(text = paste("skel$parents$",x,sep = ''))) #get parents
    if( !is.null(xpar)){xpar <- toString(xpar) #put '*' between parent names
    while (grepl(', ', xpar)){xpar = gsub(', ', '*', xpar)} }
    if(first == TRUE){ #first time through: set dataframe and formula 
      if (is.null(xpar) ){ #treat the first root nodes
        fmula_par = c(fmula_par, x) 
        
        #print(paste('FE=T parent node  ',x))
        #print(paste('the parent formula is: ', toString(fmula_par)))
        txtFL <- eval(parse(text = paste("levls$",x,sep = '')))
        #print(paste('its levels are ', txtFL))
        
        N <- length(unlist(txtFL))
        if(N==0){N=3} #3 levels if no level input
        #make uniform dist.
        #print(paste('there are ', N, 'levels'))
        if(make_df ==TRUE){ # make the initial df
          #print(paste('making df for column',x))
          #txt = paste0("df <- data.frame(",x,"= as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
          txt = paste0("df <- data.frame(",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
          
          df <- eval(parse(text = txt))#, envir = globalenv())
          make_df = FALSE 
        }
        else{ # cbind to df
          #txt = paste0("df <- cbind(df, ",x,"= as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
          txt = paste0("df <- cbind(df, ",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
          #print(paste('adding column',x))
          df <- eval(parse(text = txt))#, envir = globalenv())
        }
        
      }
      else{ #treat the first node that is not a root (or has  parents)
        fmula_ch = c(fmula_ch, x, ' | ', xpar)
        #first = FALSE
        #print(paste('FE=T child node',x))
        #print(paste('the child formula is: ', toString(fmula_ch)))
        txtFL <- eval(parse(text = paste("levls$",x,sep = '')))
        #print(paste('its levels are ', toString(txtFL)))
        #now replace the factorlevels
        #txt2 = paste0('lst2 <- unlist(strsplit(txtFL," "))')
        #eval(parse(text = txt2))
        #txt = paste0('skel$factorlevls$',Ch_name(), '<<- lst2')  #input$vec2')
        #txt = paste0('skel$factorLevels$',x, ' <<- lst2') #put into 
        #eval(parse(txt))
        
        
        N <- length(unlist(txtFL))
        print(paste('there are ', N, 'levels'))
        if(N==0){N=3} #3 levels if no level input
        Probs  <- c(rep(1/N,N)) #make uniform dist.
        if(make_df ==TRUE){ # make the initial df
          
          #txt = paste0("df <- cbind(df,",x,"= as.factor(1:",N," %*% rmultinom(n,1,prob=Probs)))")
          
          #print(paste('making df for column',x))
          txt = paste0("df <- data.frame(",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
          df <- eval(parse(text = txt))#, envir = globalenv())
          make_df = FALSE 
        }
        else {# cbind to df
          #print(paste('adding column',x))
          #txt = paste0("df <- cbind(df, ",x,"= as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
          txt = paste0("df <- cbind(df, ",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
          df <- eval(parse(text = txt))#, envir = globalenv())
        }
        
        
      }
      fmula1 <- c(fmula_par, fmula_ch) # record first pass without '+'
      fmula_par <- ''
      fmula_ch <- '' #re-initialize
      
      first = FALSE #next only pass through the else below
      
    } #end of first edge
    else { #treat the non first nodes/edge or add more edges
      #this part treats a new edge input with 3 cases: 
      # 1) new parent with existing child in BN blob,
      # 2) new child with existing parent in BN blob and 
      # 3) adding a new edge in the existing BN blob between Parent and child.
      # this implies we can only add new nodes as parent or child of existing nodes in the BN.
      if (is.null(xpar) ) { #treat case 1
        fmula_par = c(fmula_par, ' + ', x) 
        #print(paste('FE=F parent node  ',x))
        #print(paste('the parent formula is: ', toString(fmula_par)))
        txtFL <- eval(parse(text = paste("levls$",x,sep = '')))
        #print(paste('its levels are ', txtFL))
        
        N <- length(unlist(txtFL))
        if(N==0){N=3} #3 levels if no level input
        Probs  <- c(rep(1/N,N)) #make uniform dist.
        #print(paste('there are ', N, 'levels'))
        if(make_df ==TRUE){ # make the initial df
          #print(paste('making df for column',x))
          #txt = paste0("df <- data.frame(",x,"= as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
          txt = paste0("df <- data.frame(",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
          
          df <- eval(parse(text = txt))#, envir = globalenv())
          make_df = FALSE 
        }
        else{ # cbind to df
          #txt = paste0("df <- cbind(df, ",x,"= as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
          txt = paste0("df <- cbind(df, ",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
          #print(paste('adding column',x))
          df <- eval(parse(text = txt))#, envir = globalenv())
        } 
      } #end of add more edges root nodes section                      
      else{ #treat case 2 or not a root or with parents
        fmula_ch = c(fmula_ch, ' + ', x, ' | ', xpar)
        #print(paste('FE=T child node',x))
        #print(paste('the child formula is: ', toString(fmula_ch)))
        txtFL <- eval(parse(text = paste("levls$",x,sep = '')))
        #print(paste('its levels are ', toString(txtFL)))
        N <- length(unlist(txtFL))
        #print(paste('there are ', N, 'levels'))
        if(N==0){N=3} #3 levels if no level input
        Probs  <- c(rep(1/N,N)) #make uniform dist.
        if(make_df ==TRUE){ # make the initial df
          #print(paste('making df for column',x))
          txt = paste0("df <- data.frame(",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
          df <- eval(parse(text = txt))#, envir = globalenv())
          make_df = FALSE 
        }
        else {# cbind to df
          #print(paste('adding column',x))
          #txt = paste0("df <- cbind(df, ",x,"= as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
          txt = paste0("df <- cbind(df, ",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
          df <- eval(parse(text = txt))#, envir = globalenv())
        }
      } #end of add more edges with parents
    } #end of add more edges section
    
  } #end of for loop
  #if(first == TRUE){ 
  fmula = c('~', fmula1, fmula_par, fmula_ch)
  print(paste('the final out function formula is: ', toString(fmula)))
  #}
  # and now remake the network with fake data and a proper formula
  fmula1 = toString(fmula)
  while (grepl(',', fmula1)){fmula1 = gsub(',', '', fmula1)} #get rid of commas
  txt = paste0("HydeNetwork(",fmula1,", data = df)")
  Bnet <- eval(parse(text = txt))
  
  Bnet$nodeDecision <-skel$nodeDecision #and don' forget to update the nodes
  Bnet$nodeUtility <- skel$nodeUtility
  
  Bnet$factorLevels <- levls
  
  list(Bnet,df) # need them both for the good bays net and the table 
} 



# Define UI for decision app
ui <- fluidPage(

#    FIRST BUILD A SKELETON AND GET VARIABLES
    tabsetPanel(
        tabPanel("MAKE BAYES NET SKELETON"
                  ,h2("enter a new  parent node name")
                  ,column(width = 6
                          ,textInput("P_node_name3", "Node", "Parent")
                          ,radioButtons("rb1", "Choose node type:",
                                       choiceNames = list(
                                         "random",
                                         "decision",
                                         "utility"
                                       ),
                                       choiceValues = list(
                                         "random", "decision", "utility")
                          )
                          , textOutput("txt1")
                          ,    selectizeInput(
                            "vec1"
                            , "Enter the names of the possible events, decisions, or outcomes (e.g. 'heads tails' or 'success failure' "
                            , choices = list('heads', 'tails', 'success', 'failure', 'up','down', 'left', 'right','forward', 'backward', 'good', 'bad', 'red', 'green', 'blue', 'hot', 'cold', 'medium', 'middle')
                            , multiple = TRUE
                            , options = list(create = TRUE)
                          )
                          , verbatimTextOutput("oid1")
                          
                  )
                  
                  ,h2("enter a new child node name")
                  ,column(width = 6
                          ,textInput("Ch_node_name3", "Node", "Child")
                          ,radioButtons("rb2", "Choose node type:",
                                        choiceNames = list(
                                          "random",
                                          "decision",
                                          "utility"
                                        ),
                                        choiceValues = list(
                                          "random", "decision", "utility")
                          )
                          , textOutput("txt2")
                          ,    selectizeInput(
                            "vec2"
                            , "Enter the names of the possible events, decisions, or outcomes (e.g. 'heads tails' or 'success failure' "
                            , choices = list('heads', 'tails', 'success', 'failure', 'up','down', 'left', 'right','forward', 'backward', 'good', 'bad', 'red', 'green', 'blue', 'hot', 'cold', 'medium', 'middle')
                            , multiple = TRUE
                            , options = list(create = TRUE)
                          )
                          , verbatimTextOutput("oid2")
                  )
              #     ,column(width = 12
              #             ,actionButton("update_gr1", "make initial edge")
              #             ,hr()  
              #     )
              # , h3("see graph")
              # , grVizOutput("Skel3")
                 ,column(width = 12
                          ,actionButton("update_gr", "Update skeleton")
                          ,hr()
                  )
                 
            , h3("see new node")
            , grVizOutput("Skel3")
            ,h3(" Shows the formula. ")
            ,verbatimTextOutput("Formula")
        
            )
        
#NEXT USE BN SKELETON TO GENERATE SIMULATED DATA AND RECONSTRUCT BN
        ,tabPanel("MAKE BAYES NET"
                  ,column(width = 12
                          ,actionButton("update_BN", "Update BaysNet")
                          ,hr()
                  )
                  
                  ,fluidRow(column(width = 12, h3("BNs Bayes net (BN)")
                                   , grVizOutput("BN")
                  )
                  )
                  ,fluidRow(column(width = 12, h3("BNs Bayes data.")
                                   , h4("The data used to construct the cpt's of BNs BN.")
                                   ,dataTableOutput("BNTable")
                  )
                  ,helpText("This is the raw data that is used to calulate the conditional probability tables used in the BNS Bayes Net.")
                  )
                  ,fluidRow(column(width = 12, h3('View of the BNs model object.')
                                   ,h4(" Shows the details of the BNs model. ")
                                   ,verbatimTextOutput("BNObj"))
                  )
                  ,column(width = 12
                          ,actionButton("Write_BN", "Write the Bayes Net to disk")
                          ,hr()
                  )
            )

#NEXT READ BN AND UPDATE BAYES NET CPT CSV'S AND ROOT PT CSV'S
        ,tabPanel("UPDATE BAYES NET"
                  # ,fluidRow(column(width = 12
                  #                  , h3("First get a CPT")
                  #                  , fileInput("fileBN", "Choose CSV File", accept = ".csv") 
                  #                  ,checkboxInput("header", "Header", TRUE)
                  # ))
                  
                  
                  
                  ,fluidRow(column(width = 12
                                   , h3("First get a CPT")
                                   , uiOutput("selectUIBN")
                  ))
                  ,fluidRow(column(width = 12
                                   , h3("Uploaded file")
                                   , tableOutput("contentsBN")
                  ))
                  
                  ,fluidRow(column(width = 12
                                   , h3("Update and Save Probability Table")
                                   , DT::dataTableOutput("probs.df_data_BN")
                                   
                  ))
                ,fluidRow(column(width = 12
                                   , h3("Update the model")
                                   ,actionButton("update_model", "Update the Bayes Net")
                                   
                  ))

              )

#NEXT READ UDATED BN AND MAKE DECISION NET
#MODIFY UTILITY NODE CSV'S, MAKE SURE DECISION NODE ARE LABELED, MAKE THE DN, AND 
#CALCULATE MEU
      ,tabPanel("UPDATE DECISION NET"  
                ,fluidRow(column(width = 12
                                 , h3("Update the payoff node (initially its built as a CPT)")
                                 ,actionButton("update_UT_node", "Update the utility node")

                ))
                # ,fluidRow(column(width = 12
                #                  , h3("Now update the payoff node")
                #                  , fileInput("fileDN", "Choose  File", accept = ".csv")
                #                  ,checkboxInput("header", "Header", TRUE)
                # ))
                ,fluidRow(column(width = 12
                                 , h3("First get the payoff")
                                 , uiOutput("selectUIDN")
                ))
                

                ,fluidRow(column(width = 12
                                 , h3("Uploaded file")
                                 , tableOutput("contentsDN")
                ))

                ,fluidRow(column(width = 12
                                 , h3("Update and Save Payoff Table")
                                 , DT::dataTableOutput("probs.df_data_DN")

                ))
                ,fluidRow(column(width = 12
                                 , h3("Update the DN and calculate the Maximum Expected Utility")
                                 ,actionButton("updateDN", "Update the Decision Net")

                ))
                , h3("see the final DN")
                , grVizOutput("BN_glob2")
                ,h3(" Shows the MEU. ")
                ,verbatimTextOutput("MEU")
              )  
            )
        )

# Define server oto make a decision net
server <- function(input, output, session) {
  
    #all the following are updated when update_gr button is clicked
    dia_func2 <- eventReactive(input$update_gr,{a <- #initial formula
        paste('~', input$P_node_name3, 
              '+', input$Ch_node_name3 ,'|', input$P_node_name3, sep = ' ') }) 
  
    dia_func3 <- eventReactive(input$update_gr,{a <-  #updated formula for new parent, child, or both
        paste('.~.', '+', input$Ch_node_name3 ,'|', input$P_node_name3, sep = ' ') })
    
    Par_name <- eventReactive(input$update_gr,{a <- input$P_node_name3 })
    Ch_name <- eventReactive(input$update_gr,{a <- input$Ch_node_name3 })
    
    Par_fac_lev <- eventReactive(input$update_gr,{a <- input$vec1 })
    Ch_fac_lev <- eventReactive(input$update_gr,{a <- input$vec2 })
    
    RB_1 <- eventReactive(input$update_gr,{a <- input$rb1 })
    RB_2 <- eventReactive(input$update_gr,{a <- input$rb2 })
    
    value <- reactiveValues(first = TRUE) #1st skeleton render flag
    #need to pass the factor levels and decision and utility types past the skeleton
    levls <- list()
    dec_nodes<- reactiveValues()
    util_nodes<- reactiveValues()
    
    skel <- HydeNetwork(~ x + y|x) #make an initial network object
    
    
    
    output$Skel3 <- renderDiagrammeR({ #build skeleton and display
      
        if(value$first == FALSE){

            skel <<- update(skel,  as.formula(dia_func3()))
            #set decision and utility nodes in plot
            if(RB_1() == 'decision'){ 
              txt = paste0("setDecisionNodes(skel, ",Par_name(),")")
              skel <<- eval(parse(text = txt))
              #now put this value into dec_nodes for later processing (update forgets these values)
              txt1 = paste0('lst1 <- unlist(strsplit(Par_name()," "))')
              eval(parse(text = txt1))
              txt = paste0('dec_nodes$',Par_name(), ' <<- lst1') 
              eval(parse(text = txt))
              
              
            }  
            if(RB_1() == 'utility'){ 
              txt = paste0("setUtilityNodes(skel, ",Par_name(),")")
              skel <<- eval(parse(text = txt))
              #now put this value into util_nodes for later processing (update forgets these values)
              txt1 = paste0('lst1 <- unlist(strsplit(Par_name()," "))')
              eval(parse(text = txt1))
              txt = paste0('util_nodes$',Par_name(), ' <<- lst1') 
              eval(parse(text = txt))
              
            }  
            #set decision and utility nodes in plot for child nodes
            if(RB_2() == 'decision'){ 
              txt = paste0("setDecisionNodes(skel, ",Ch_name(),")")
              skel <<- eval(parse(text = txt))
              #now put this value into dec_nodes for later processing (update forgets these values)
              txt1 = paste0('lst1 <- unlist(strsplit(Ch_name()," "))')
              eval(parse(text = txt1))
              txt = paste0('dec_nodes$',Ch_name(), ' <<- lst1') 
              eval(parse(text = txt))
              
            }  
            if(RB_2() == 'utility'){ 
              txt = paste0("setUtilityNodes(skel, ",Ch_name(),")")
              skel <<- eval(parse(text = txt))
              #now put this value into util_nodes for later processing (update forgets these values)
              txt1 = paste0('lst1 <- unlist(strsplit(Ch_name()," "))')
              eval(parse(text = txt1))
              txt = paste0('util_nodes$',Ch_name(), ' <<- lst1') 
              eval(parse(text = txt))
              
            }
            
            # now set the factor levels 1st the parent
            txt1 = paste0('lst1 <- unlist(strsplit(Par_fac_lev()," "))')
            eval(parse(text = txt1))
            #print(paste('VECT1 IS', lst1))
            txt = paste0('skel$factorLevels$',Par_name(), ' <<- lst1') #put into skel for diag
            eval(parse(text = txt))
            txt2 = paste0('levls$',Par_name(), ' <<- lst1') #put levls for later processin
            eval(parse(text = txt2))
            
            # now set the factor levels 2nd the child           
            txt2 = paste0('lst2 <- unlist(strsplit(Ch_fac_lev()," "))')
            eval(parse(text = txt2))
            #print(paste('VECT2 IS', lst2))
            txt = paste0('skel$factorLevels$',Ch_name(), '<<- lst2')  #input$vec2')
            eval(parse(text = txt))
            txt2 = paste0('levls$',Ch_name(), ' <<- lst2') #put levels for later processin
            eval(parse(text = txt2))
            

        }
        else{ # first time through
          skel <<- HydeNetwork(as.formula(dia_func2())) #make a new skel object 1st time
          value$first<-FALSE

        }
      plot(skel)       
    })
    
        output$Formula <- renderPrint({
            paste('text: ', dia_func3()
                  , 'formula: ',toString(skel$network_formula), sep=' ')
        })  
        
       
# now reconstruct the bnet with fake data

  # use reactiveValues to store out() output
        
        BN_reactive <- reactiveValues()
          # BNDF <- out(ndata = 50, skel=skel, first = TRUE) #call once
          # BN_reactive$BN = BNDF[[1]] #this will initialize for ~ x+x|y
          # BN_reactive$DF = BNDF[[2]]
          #skel <- BN_reactive$BN
          
        
  # update when the update button is clicked presumably after the BN skeleton is updated

        observeEvent(input$update_BN,{
          BNDF <- out(ndata = 50, skel=skel, levls=levls) #call to update data
          BN_reactive$BN = BNDF[[1]] #update new edges
          BN_reactive$DF = BNDF[[2]]
          assign("BN_glob", BN_reactive$BN, envir = .GlobalEnv)
          skel <<- BN_reactive$BN

  # and the diagrapmR
        output$BN <- renderDiagrammeR({
            plot(BN_reactive$BN)
        })
        
  # can render the dataframe with factors/char      
        output$BNTable <- renderDT({  
          Btable <- as.data.frame(BN_reactive$DF)
        })

  #show the model    
        output$BNObj <- renderPrint({ 
            writeNetworkModel(BN_reactive$BN, pretty=TRUE)

        })
        
        }) # end of update_BN
        
  # write the decision net to disk
        observeEvent(input$Write_BN,{
          make_freq_tables(BN_glob, outputDir) #BN_glob is in the global env
          
          
        })
        

  # this is the experimental radio output
        output$txt1 <- renderText({
          paste("You chose", input$rb1)
        })
        output$txt2 <- renderText({
          paste("You chose", input$rb2)
    
        })
        # this is the syntax to set the nodes to decision or utility
        # net <- setDecisionNodes(net, hit1, hit2, hit3)
        # net <- setUtilityNodes(net, payoff)
        # 
        # c(net$nodeDecision$hit2, net$nodeUtility$payoff)
        # 
        # ## [1] TRUE TRUE
  
        
        output$oid1 <- renderPrint({
          
          req(input$vec1)
          
          cat("As string:\n")
          cat(input$vec1)
          #skel$factorLevels < input$vec1
          
        })
        
        output$oid2 <- renderPrint({
          
          req(input$vec2)
          
          cat("As string:\n")
          cat(input$vec2)
          #skel$factorLevels 
          
        })
        
# now treat the csv file I/O and editing for the BN        
        # dataBN <- reactive({
        #   
        #   file <- input$fileBN #
        #   
        #   ext <- tools::file_ext(file$datapath)
        #   
        #   req(file)
        #   validate(need(ext == "csv", "Please upload a BN csv file"))
        #   
        #   d1 <- read.csv(file$datapath, header = input$header)
        # 
        #   assign("d1", d1, envir = .GlobalEnv) 
        #   
        #   d1
        # })
        
 #select the global env prob table df to edit       
        output$selectUIBN<-renderUI({
          req(BN_glob$data)
          selectInput(inputId='selectprobBN', label='select probability table', choices = names(BN_glob$data))
        })
        
        dataBN <- reactive({
          
          txt = paste0('d1 <- ',input$selectprobBN)
          eval(parse(text = txt), envir = .GlobalEnv)
          
          
          #assign("d1", d1, envir = .GlobalEnv) 
          
          d1


        })
        
        output$contentsBN <- renderTable({
          dataBN()
        })
        
        output$probs.df_data_BN <- renderDataTable({
          
          d1 <-dataBN()
          
          df <- datatable(
            d1,
            selection = 'none', editable = TRUE, 
            rownames = TRUE,
            extensions="Buttons",
            
            options = list(
              paging = TRUE,
              searching = TRUE,
              fixedColumns = TRUE,
              autoWidth = TRUE,
              ordering = TRUE,
              dom = 'Bfrtip',
              buttons = c('csv', 'excel')
            ),
            
            class = "display"
          )
          
          #txt = paste0(input$selectprob, '<- df ')
          #eval(parse(text = txt))
          #assign(input$selectprob, d1, envir = .GlobalEnv)
          
          return(df)
        })
        
        # Every rendered DT will create a input_cell_edit object that contains the row and column index of the edit.
        observeEvent(input$probs.df_data_BN_cell_edit, {
          d1[input$probs.df_data_BN_cell_edit$row,input$probs.df_data_BN_cell_edit$col] <<- input$probs.df_data_BN_cell_edit$value
          
          assign(input$selectprobBN, d1, envir = .GlobalEnv) #put the updated df in the global envir
          
        })
        
        # update the bayes net when new cpts are input
        observeEvent(input$update_model,{
          
          new_df <- updateBN(BN_glob, outputDir)
          BN_glob$data <- new_df
          
          assign("BN_glob", BN_glob, envir = .GlobalEnv)#BN_glob is in the global env
          
        })
        
        observeEvent(input$update_UT_node,{
          
          make_util_tables(BN_glob, outputDir)

        })
        
# now treat the csv file I/O and editing for the utility node of the DN        
        # dataDN <- reactive({
        # 
        #   file <- input$fileDN #
        # 
        #   ext <- tools::file_ext(file$datapath)
        # 
        #   req(file)
        #   validate(need(ext == "csv", "Please upload a DN utility csv file"))
        # 
        #   d1 <- read.csv(file$datapath, header = input$header)
        # 
        #   assign("d1", d1, envir = .GlobalEnv)
        # 
        #   d1
        # })
        # 
        # #select the global env prob table df of the payoff to edit
        output$selectUIDN<-renderUI({
          req(BN_glob$data)
          selectInput(inputId='selectprobDN', label='select payoff table', choices = names(BN_glob$data))
        })

        dataDN <- reactive({

          txt = paste0('d1 <- ',input$selectprobDN)
          eval(parse(text = txt), envir = .GlobalEnv)


          #assign("d1", d1, envir = .GlobalEnv)

          d1


        })
        output$contentsDN <- renderTable({
          dataDN()
        })
        
        output$probs.df_data_DN <- renderDataTable({
          
          d1 <-dataDN()
          
          df <- datatable(
            d1,
            selection = 'none', editable = TRUE, 
            rownames = TRUE,
            extensions="Buttons",
            
            options = list(
              paging = TRUE,
              searching = TRUE,
              fixedColumns = TRUE,
              autoWidth = TRUE,
              ordering = TRUE,
              dom = 'Bfrtip',
              buttons = c('csv', 'excel')
            ),
            
            class = "display"
          )
          
          return(df)
        })
        
        # Every rendered DT will create a input_cell_edit object that contains the row and column index of the edit.
        observeEvent(input$probs.df_data_DN_cell_edit, {
          d1[input$probs.df_data_DN_cell_edit$row,input$probs.df_data_DN_cell_edit$col] <<- input$probs.df_data_DN_cell_edit$value
          
          assign(input$selectprobDN, d1, envir = .GlobalEnv) #put the updated df in the global envir
          
        })
          
        # update the DECISION NET when AFTER THE UTILITY NODE IS UPDATED
        observeEvent(input$updateDN,{
            
            BN_glob2 <- updateDN(BN_glob, outputDir)
            
            
            #and now iterate over the nodes to set util and decision and get names   
              
            for(x in BN_glob$nodes){     #iterate over all  nodes
              print(paste('node is: ', x))
              xpar = eval(parse(text = paste("BN_glob$parents$",x,sep = ''))) #get parents
              Util <- eval(parse(text = paste("BN_glob$nodeUtility$",x,sep = ''))) #test if node is a utility node
              dec <- eval(parse(text = paste("BN_glob$nodeDecision$",x,sep = ''))) #test if node is a utility node
                  
                  
              if ( dec) { #treat the decision nodes
                dec_name <- x
                print(paste('decision node is: ', x))
                txt = paste0("BN_glob2 <- setDecisionNodes(BN_glob2,",x,")")
                eval(parse(text = txt))
               }
                  
              if (!is.null(xpar) & Util) { #treat the utility nodes
                util_name <- x
                print(paste('utility node is: ', x))
                txt = paste0("BN_glob2 <- setUtilityNodes(BN_glob2,",x,")")
                eval(parse(text = txt))
                }    
                  
                print(paste('parent are: ', xpar))
                  
            } #end loop
    
            assign("BN_glob2", BN_glob2, envir = .GlobalEnv)  #BN_glob is in the global env
                
            #4) calculate the MEU
                
            writeNetworkModel(BN_glob2, pretty=TRUE)
            compiledNet <- compileJagsModel(BN_glob2, #data = evidence,
                                             n.chains = 3,
                                             n.adapt = 5000)
      
            post <- HydeSim(compiledNet,
                            variable.names = BN_glob$nodes,
                            n.iter=10000)
                
            dplyr::sample_n(post, 20)
                
                
            compiledNets <- compileDecisionModel(BN_glob2) 
                
            samples <- lapply(compiledNets,
                              HydeSim,
                              variable.names = BN_glob$nodes,
                              n.iter=10000)
                
            lapply(samples, head)
            txt = paste0("lapply(samples, function(l) mean(l$",util_name,"))")
            eval(parse(text = txt))
    
          output$BN_glob2 <- renderDiagrammeR({ #build skeleton and display
          plot(BN_glob2)       
          })
              
          output$MEU <- renderPrint({
            
            txt <- paste('dec_values <- BN_glob$factorLevels$', dec_name)
            eval(parse(text = txt))
            #dec_values <- BN_glob$factorLevels$Decision
            #paste('the decisions are: ', toString(dec_values))
            txt <- paste('meu <- toString(lapply(samples, function(l) mean(l$', util_name, ')))')
            eval(parse(text = txt))
            msg <- paste('the decisions are: ', toString(dec_values), '. \n',
                         'The expected value of each decision: ', meu, sep=' ')
            strsplit(msg, "\n")[[1]]
          })
          
          
      })  

}

#one more line
# Run the application add
shinyApp(ui = ui, server = server)
