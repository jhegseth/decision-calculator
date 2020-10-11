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
                          ,textInput("P_node_name3", "Node", "w")
                  )
                  ,h2("enter a new child node name")
                  ,column(width = 6
                          ,textInput("Ch_node_name3", "Node", "z")
                  )
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

# Define server logic required to draw a histogram
server <- function(input, output) {



    dia_func2 <- eventReactive(input$update_gr,{a <- 
        paste('~', input$P_node_name3, 
              '+', input$Ch_node_name3 ,'|', input$P_node_name3, sep = ' ') }) 
    
    dia_func3 <- eventReactive(input$update_gr,{a <- 
        paste('.~.', #'+',input$P_node_name3, 
              '+', input$Ch_node_name3 ,'|', input$P_node_name3, sep = ' ') })
    
    value <- reactiveValues(first = TRUE)
    skel <- HydeNetwork(~ x + y|x)
    
    
    output$Skel3 <- renderDiagrammeR({
        #skel<- HydeNetwork(as.formula(dia_func2()))
        if(value$first == TRUE){
            #skel <<- update(skel, .~. + p|q)
            #skel <<- update(skel, .~. - y|x)
            skel <<- HydeNetwork(as.formula(dia_func2()))#~ p + q|p)
            value$first<-FALSE
            }
        else{
            #get the current formula
            #now fix the formula
            while (grepl('\\+ \\+', sf)){sf = gsub('\\+ \\+', '\\+', sf)} # make sure there are no spaces beteen +s
            while (grepl('\\+\\+', sf)){sf = gsub('\\+\\+', '\\+', sf)} #get rid of repeated +s
            skel <<- update(skel,  as.formula(dia_func3()))
            }
        plot(skel)
    })
    
        output$Formula <- renderPrint({
            paste('text: ', dia_func3()
                  , 'formula: ',toString(skel$network_formula), sep=' ')
        })  
        
    #FUNCTION TO OUTPUT NEW RECONSTRUCED BN AND FAKE DATA DF
        out <- function(ndata, skel){ #INPUT NO. OF DATA ROWS AND RAW SKELETON
            n=ndata
            first_node = TRUE
            fmula <- '~'  # reconstruct a formula without errors
            for(x in skel$nodes){     #iterate over all nodes
                xpar = eval(parse(text = paste("skel$parents$",x,sep = ''))) #get parents
                if( !is.null(xpar)){xpar <- toString(xpar) #put '*' between parent names
                    while (grepl(', ', xpar)){xpar = gsub(', ', '*', xpar)} }
                if(first_node == TRUE){ #first time through: set dataframe and formula 
                    if (is.null(xpar)){ #treat the first root nodes
                        fmula = c(fmula, x) 
                        first_node = FALSE 
                        print(paste('first  node is a root node  ',x))
                        txt = paste0("data.frame(",x,"= as.factor(1:3 %*% rmultinom(n,1,prob=c(.5,.3,.2))))")
                        df <- eval(parse(text = txt))
                    }
                    else{ #treat the first node that is not a root (or has  parents)
                        fmula = c(fmula,  x, ' | ', xpar)
                        first_node = FALSE
                        print(paste('first node is with children ',x))
                        txt = paste0("data.frame(",x,"= as.factor(1:3 %*% rmultinom(n,1,prob=c(.5,.3,.2))))")
                        df <- eval(parse(text = txt))
                        #.GlobalEnv$data <- df
                    }
                }
                else { #treat the non first nodes
                    if (is.null(xpar)){ #treat secondary root nodes
                        fmula = c(fmula, ' + ', x)
                        print(paste('secondary root node ',x))
                        txt = paste0("cbind(df,",x,"= as.factor(1:3 %*% rmultinom(n,1,prob=c(.5,.3,.2))))")
                        df <- eval(parse(text = txt))
                    }
                    else{ #treat the non first node that is not a root or with parents
                        fmula = c(fmula, ' + ',  x, ' | ', xpar)
                        print(paste('secondary node with children ',x))
                        txt = paste0("cbind(df,",x,"= as.factor(1:3 %*% rmultinom(n,1,prob=c(.5,.3,.2))))")
                        df <- eval(parse(text = txt))
                        # for (prnt in xpar) {print(paste0('parent(s) of ',x, ' are ',prnt))
                        #   fmula = c(fmula, ' + ', x, ' | ', xpar)} 
                    }
                }
                
            }
            # and now remake the network with fake data and a proper formula
            fmula1 = toString(fmula)
            fmula1 = gsub(',','',fmula1)
            txt = paste0("HydeNetwork(",fmula1,", data = df)")
            Bnet <- eval(parse(text = txt))
            
            list(Bnet,df) # need them both for the good bays net and the table 
        }        
# now reconstruct the bnet with fake data

  # use reactiveValues to store out() output
        BN_reactive <- reactiveValues(
          BN = out(ndata = 50, skel=skel)[[1]] #this will initialize for ~ x+x|y
          ,DF = out(ndata = 50, skel=skel)[[2]]
        )
  # update when the update button is clicked presumably after the BN skeliton is updated
        observeEvent(input$update_BN,{ 
          BN_reactive$BN <- out(ndata = 50, skel=skel)[[1]]
          BN_reactive$DF <- out(ndata = 50, skel=skel)[[2]]
        })

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

  
}

# Run the application 
shinyApp(ui = ui, server = server)
