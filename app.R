#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(igraph)
library(DiagrammeR)
library(HydeNet)
library(datasets)

# Define UI for application that draws a histogram
ui <- fluidPage(

#    FIRST BUILD A SKELETON AND GET VARIABLES
    tabsetPanel(
        tabPanel("Hydenet Skeleton"
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
                            , choices = list('d','f')
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
                            , choices = list('a','s')
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
        ,tabPanel("BN Net"
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
                  
                  
                  
        ) 
        
    )
)

# Define server oto make a decision net
server <- function(input, output) {
  
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
            print(paste('VECT1 IS', lst1))
            txt = paste0('skel$factorLevels$',Par_name(), ' <<- lst1') #put into skel for diag
            eval(parse(text = txt))
            txt2 = paste0('levls$',Par_name(), ' <<- lst1') #put levls for later processin
            eval(parse(text = txt2))
            
            # now set the factor levels 2nd the child           
            txt2 = paste0('lst2 <- unlist(strsplit(Ch_fac_lev()," "))')
            eval(parse(text = txt2))
            print(paste('VECT2 IS', lst2))
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
        
    #FUNCTION TO OUTPUT NEW RECONSTRUCED BN AND FAKE DATA DF
        out <- function(ndata, skel){ #INPUT NO. OF DATA ROWS AND RAW SKELETON and wheher the first iteration
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
                      
                        print(paste('FE=T parent node  ',x))
                        print(paste('the parent formula is: ', toString(fmula_par)))
                        txtFL <- eval(parse(text = paste("levls$",x,sep = '')))
                        print(paste('its levels are ', txtFL))
                        
                        N <- length(unlist(txtFL))
                        if(N==0){N=3} #3 levels if no level input
                        #make uniform dist.
                        print(paste('there are ', N, 'levels'))
                        if(make_df ==TRUE){ # make the initial df
                          print(paste('making df for column',x))
                          #txt = paste0("df <- data.frame(",x,"= as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
                          txt = paste0("df <- data.frame(",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
                          
                          df <- eval(parse(text = txt))#, envir = globalenv())
                          make_df = FALSE 
                        }
                        else{ # cbind to df
                          #txt = paste0("df <- cbind(df, ",x,"= as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
                          txt = paste0("df <- cbind(df, ",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
                          print(paste('adding column',x))
                          df <- eval(parse(text = txt))#, envir = globalenv())
                        }
                        
                    }
                    else{ #treat the first node that is not a root (or has  parents)
                        fmula_ch = c(fmula_ch, x, ' | ', xpar)
                        #first = FALSE
                        print(paste('FE=T child node',x))
                        print(paste('the child formula is: ', toString(fmula_ch)))
                        txtFL <- eval(parse(text = paste("levls$",x,sep = '')))
                        print(paste('its levels are ', toString(txtFL)))
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
                          
                          print(paste('making df for column',x))
                          txt = paste0("df <- data.frame(",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
                          df <- eval(parse(text = txt))#, envir = globalenv())
                          make_df = FALSE 
                        }
                        else {# cbind to df
                          print(paste('adding column',x))
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
                      print(paste('FE=F parent node  ',x))
                      print(paste('the parent formula is: ', toString(fmula_par)))
                      txtFL <- eval(parse(text = paste("levls$",x,sep = '')))
                      print(paste('its levels are ', txtFL))
                      
                      N <- length(unlist(txtFL))
                      if(N==0){N=3} #3 levels if no level input
                      Probs  <- c(rep(1/N,N)) #make uniform dist.
                      print(paste('there are ', N, 'levels'))
                      if(make_df ==TRUE){ # make the initial df
                        print(paste('making df for column',x))
                        #txt = paste0("df <- data.frame(",x,"= as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
                        txt = paste0("df <- data.frame(",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
                        
                        df <- eval(parse(text = txt))#, envir = globalenv())
                        make_df = FALSE 
                        }
                      else{ # cbind to df
                        #txt = paste0("df <- cbind(df, ",x,"= as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
                        txt = paste0("df <- cbind(df, ",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
                        print(paste('adding column',x))
                        df <- eval(parse(text = txt))#, envir = globalenv())
                        } 
                    } #end of add more edges root nodes section                      
                    else{ #treat case 2 or not a root or with parents
                      fmula_ch = c(fmula_ch, ' + ', x, ' | ', xpar)
                      print(paste('FE=T child node',x))
                      print(paste('the child formula is: ', toString(fmula_ch)))
                      txtFL <- eval(parse(text = paste("levls$",x,sep = '')))
                      print(paste('its levels are ', toString(txtFL)))
                      N <- length(unlist(txtFL))
                      print(paste('there are ', N, 'levels'))
                      if(N==0){N=3} #3 levels if no level input
                      Probs  <- c(rep(1/N,N)) #make uniform dist.
                      if(make_df ==TRUE){ # make the initial df
                        print(paste('making df for column',x))
                        txt = paste0("df <- data.frame(",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
                        df <- eval(parse(text = txt))#, envir = globalenv())
                        make_df = FALSE 
                      }
                      else {# cbind to df
                        print(paste('adding column',x))
                        #txt = paste0("df <- cbind(df, ",x,"= as.factor(1:3 %*% rmultinom(",n,",1,prob=c(.5,.3,.2))))")
                        txt = paste0("df <- cbind(df, ",x,"= as.factor(1:",N," %*% rmultinom(",n,",1,prob=Probs)))")
                        df <- eval(parse(text = txt))#, envir = globalenv())
                      }
                    } #end of add more edges with parents
                } #end of add more edges section
                
            } #end of for loop
            #if(first == TRUE){ 
              fmula = c('~', fmula1, fmula_par, fmula_ch)
              print(paste('the final formula is: ', toString(fmula)))
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
# now reconstruct the bnet with fake data

  # use reactiveValues to store out() output
        
        BN_reactive <- reactiveValues()
          # BNDF <- out(ndata = 50, skel=skel, first = TRUE) #call once
          # BN_reactive$BN = BNDF[[1]] #this will initialize for ~ x+x|y
          # BN_reactive$DF = BNDF[[2]]
          #skel <- BN_reactive$BN
          
        
  # update when the update button is clicked presumably after the BN skeleton is updated

        observeEvent(input$update_BN,{
          BNDF <- out(ndata = 50, skel=skel) #call to update data
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
        
}

# Run the application add
shinyApp(ui = ui, server = server)
