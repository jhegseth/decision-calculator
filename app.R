#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(igraph)
library(DiagrammeR)
library(HydeNet)


#and now the two graphs
#first the cars graph
#make the diagram

ui <- shinyUI(fluidPage(titlePanel("Explore several exmples of BN's and ID's.")
    ,navbarPage("Choose view", inverse = TRUE
        ,tabPanel("Two examples"
            ,fluidRow(column(width = 6, h4("Cars Bayes net")
                              , h3("The BN")
                              , grVizOutput("Car")
                             ),
                        column(width = 6, h4("Blackjack Influence diagram")
                               , h3("The Influence diagram")
                               , grVizOutput("Bjack")
                             ),
                            helpText("Both of these diagrams show models of the Bayes net type.")
                        )
                    )
        ,navbarMenu("Make skeleton"
            ,tabPanel("Make a node"
                ,fluidPage(h2('enter a new node name')
                    , fluidRow(column(width = 6
                                ,textInput("node_name", "Node", "Node1")
                                ,checkboxInput("Show_gr", "Show graph", FALSE),
                                    )
                            ,column(width = 6
                                , h3("see new node")
                                , grVizOutput("Skel")
                                    )
                                )
                    ,helpText("This is the raw data that is used to calulate the conditional probability tables used in the CARS Bayes Net.")
                            )
                ,fluidPage(h2('cars model object')
                    , fluidRow(column(width = 6
                            , textOutput("skelObj")
                                    )
                            ,column(width = 6
                                , h3("Summary")
                                , verbatimTextOutput("This is the CARS Bayes net model object")
                                    )
                                )
                            )
                    )
            ,tabPanel("make a dependency"
                ,h2('View of the blackjack model object')
                ,helpText(" Show many aspects of the model object. ")
                ,fluidPage(h2('model object'))
                ,fluidRow(column(10
                    ,textOutput("BjackObj"))
                         )
                       )
               )

               # navbarMenu("Show data",
               #            tabPanel("CARS",
               #                     fluidPage(h2('cars table'),
               #                               fluidRow(
               #                                   column(width = 12,
               #                                          dataTableOutput("CarsTable")
               #                                   ),
               #                                   column(width = 5, h3("Summary"),
               #                                          verbatimTextOutput("This is the CARS data table.")
               #                                   )
               # 
               # 
               #                               ),
               #                               helpText("This is the raw data that is used to calulate the conditional probability tables used in the CARS Bayes Net.")
               #                     )
               #                     ,
               #                     fluidPage(h2('cars model object'),
               #                               fluidRow(
               #                                   column(width = 12,
               #                                          textOutput("CarsObj")
               #                                   ),
               #                                   column(width = 5, h3("Summary"),
               #                                          verbatimTextOutput("This is the CARS Bayes net model object")
               #                                   )
               # 
               # 
               #                               )
               #                     )
               # 
               # 
               #            ),
               #            tabPanel("Blackjack",
               #                     h2('View of the blackjack model object'),
               #                     helpText(" Show many aspects of the model object. "),
               # 
               #                     fluidPage(h2('model object')),
               # 
               #                     fluidRow(
               #                         column(12
               #                                ,textOutput("BjackObj"))
               #                     )
               # 
               #            )
               # ),

        ,tabPanel("About"
            ,fluidRow(h3("On the ferriswheel in Paris"))
            ,fluidRow(column(6, imageOutput("image1"))
                ,column(6, includeMarkdown("jh_paragraph.Rmd"))
                      )
            ,fluidRow(h3("Check out his resume and contact him at:"
                         ,a(href="https://www.linkedin.com/pub/john-hegseth/6/864/b20", "jhegseth")
                        )
                     )
                 )
    )   
))


#now make the server side app

server <- shinyServer(function(input, output) {
    # first assemble examples to be displayed
    mtcars2 <- transform(mtcars,
                         cyl = factor(cyl),
                         gear=factor(gear),
                         am = factor(am))
    
    carNet <- HydeNetwork(~ cyl
                          + disp | cyl
                          + hp | disp
                          + wt
                          + gear
                          + mpg | disp*hp*wt*gear,
                          data=mtcars2)
    data(BlackJack, package="HydeNet")    
    
    
    # and the diagrapmR
    output$Car <- renderDiagrammeR({
        plot(carNet)   
    })
    output$CarsTable = renderDataTable({
        mtcars
    })
    
    output$CarsObj <- renderText({
        writeNetworkModel(carNet, pretty=FALSE)
        
    })
    #BUILD SKELETON
    
    observeEvent(input$show_gr,{
        a <- paste('~', input$node_name, '+ disp' ,'|', input$node_name, sep = ' ')
        firstN <-paste('~',as.character(input$node_name),'+', 'disp', sep = ' ') #, '+', 'disp', '|',input$node_name, sep = ' ')
        
            print(input$a)
    })
    
    output$Skel <- renderDiagrammeR({
        
        #this works at the console
        #a <- "~ cyl+ disp | cyl" this is a charactr string
        #> carNet <- HydeNetwork(as.formula(a)) convert to a formula
        #> plot(carNet)
        # skelNet <- HydeNetwork(~ input$node_name
        #                        + disp | cyl
        #                        + hp | disp
        #                        + wt
        #                        + gear
        #                        + mpg | disp*hp*wt*gear,
        #                        #data=mtcars2
        # )
        # rule 1st node gets ~, all others get +
        # , first dep of a node gets | all others get *
        a <- paste('~', input$node_name, '+ disp' ,'|', input$node_name, sep = ' ')
        firstN <-paste('~',as.character(input$node_name),'+', 'disp', sep = ' ') #, '+', 'disp', '|',input$node_name, sep = ' ')
        input$Show_gr
        skelNet <- HydeNetwork(as.formula(a))
        plot(skelNet)   
    })
    # and see the object made so far
    output$skelObj <- renderText({
        a <- paste('~', input$node_name, '+ disp' ,'|', input$node_name, sep = ' ')
        firstN <-paste('~',as.character(input$node_name),'+', 'disp', sep = ' ') #, '+', 'disp', '|',input$node_name, sep = ' ')
        skelNet <- HydeNetwork(as.formula(a))
        writeNetworkModel(skelNet, pretty=FALSE)
    })
    
    #and the second graph
    #the blackjack example
    #data(BlackJack, package="HydeNet")
    
    output$Bjack <- renderDiagrammeR({
        #data(BlackJack, package="HydeNet")
        plot(BlackJack)
    })
    
    output$BjackObj <- renderText({
        #data(BlackJack, package="HydeNet")
        writeNetworkModel(BlackJack, pretty=FALSE)
    })
    # image1 sends pre-rendered images stored on server
    output$image1 <- renderImage({
        list(
            src = "jh_paris_ferriswheel.jpg",
            width = "600px",
            height = "400px",
            alt = "this is me"
        )}, deleteFile = FALSE)
    
    
    
})


shinyApp(ui = ui, server = server)