#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
source("utils.R")
library(shiny)
library(shinythemes)
library(bslib)
library("latex2exp")
library(ggplot2)
require(gridExtra)

# Define UI for application that draws a histogram
ui <- 
    navbarPage("Targeting Sample Size", selected = "Home", collapsible = TRUE, inverse = TRUE, theme = shinytheme("readable"),
               tags$head(
                   # Note the wrapping of the string in HTML()
                   tags$style(HTML("
      @import url('https://fonts.googleapis.com/css?family=Google+Sans|Noto+Sans|Castoro');
     
      body {
        font-family: 'Google Sans', sans-serif!important;
      }"
                                   ))
               ),
               tabPanel("Home", fluidPage(
                   
                   # Application title
                   titlePanel(column(12, h1("A Sample Size Calculation for Training and Certifying Targeting Policies"), align="center")),
                   titlePanel(column(12, h3( tags$a(href="https://www.kellogg.northwestern.edu/faculty/directory/timoshenko_artem/", "Artem Timoshenko", target="_new")), align="center")),
                   titlePanel(column(12, h5("Kellogg School of Management, Northwestern University", align="center"))),
                              
                   br(),
                   
                   column( 12,
                       br(),br(),
                       p("Lorem ipsum dolor sit amet, consectetuer adipiscing elit. Donec odio. Quisque volutpat mattis eros. Nullam malesuada erat ut turpis. Suspendisse urna nibh, viverra non, semper suscipit, posuere a, pede."),
                       p("Donec nec justo eget felis facilisis fermentum. Aliquam porttitor mauris sit amet orci. Aenean dignissim pellentesque felis.
Morbi in sem quis dui placerat ornare. Pellentesque odio nisi, euismod in, pharetra a, ultricies in, diam. Sed arcu. Cras consequat."),
                   br() ,align="center"
                   
                   ),
                   column(6,shiny::a(icon("file"), "Paper",href="https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4228297",
                                                 class = "btn btn-primary btn-lg", target="_new"
                                            
                   ),align="right"),
                   column(6,shiny::a(icon("github"),"Code",href="https://github.com/JaloteEshi/Training_Certifying_Targeting_Policies/tree/master/Code",
                                                 class = "btn btn-default btn-lg",  target="_new"
                                    
                   ),align="left"),
                  
               )
                        
                        
                        
                        
                        
                        ),
               tabPanel("Analyze your data",
                        fluidPage(
                            
                            # Application title
                            titlePanel("Upload your file and the script will run"),
                            
                            # Sidebar with a slider input for number of bins 
                            sidebarLayout(
                                sidebarPanel(
                                    h5("The graph will show. Please wait for atleast a minute."),
                                    # Input: Select a file ----
                                    fileInput("file1", "Choose CSV File",
                                              accept = ".csv", placeholder = "No file selected"),
                                    
                                    # Horizontal line ----
                                    tags$hr(),
                                    
                                    # Input: Checkbox if file has header ----
                                    checkboxInput("header", "Header", TRUE),
                                    
                                    # Input: Select separator ----
                                    radioButtons("sep", "Separator",
                                                 choices = c(Comma = ",",
                                                             Semicolon = ";",
                                                             Tab = "\t"),
                                                 selected = ","),
                                    
                                    # Input: Select quotes ----
                                    radioButtons("quote", "Quote",
                                                 choices = c(None = "",
                                                             "Double Quote" = '"',
                                                             "Single Quote" = "'"),
                                                 selected = '"'),
                                    
                                    radioButtons("type", "Graph Type",
                                                 choices = c(
                                                             "First Graph" = '1',
                                                             "Second Graph" = '2',
                                                             "Third Graph" = '3',
                                                             "Fourth Graph" = '4'
                                                             
                                                             ),
                                                 selected = '1'),
                                    actionButton("graph", "Submit"),
                                ),
                                
                                
                                # Show a plot of the generated distribution
                                mainPanel(
                                    span(textOutput("errorMessage"), style="color:red"),
                                    plotOutput("distPlot")
                                )
                            )
                        )), 
               tabPanel("Sample Python Notebook", fluidPage(
                   
                   # Application title
                   titlePanel("Sample Python Notebook"),
                   
                   column( 12,
                           br(),br(),
                           p("Lorem ipsum dolor sit amet, consectetuer adipiscing elit. Donec odio. Quisque volutpat mattis eros. Nullam malesuada erat ut turpis. Suspendisse urna nibh, viverra non, semper suscipit, posuere a, pede."),
                         
                           br() ,align="center"
                           
                   ),
                   column( 12,  
                       # Show a plot of the generated distribution
                   # https://jaloteeshi.github.io/ArtemPythonCode/lab/index.html
                           tags$iframe(src = "https://jaloteeshi.github.io/ArtemTargetSample/",
                                       height = '5000px',
                                       width = '100%', frameborder="0"),
                   ),
                       
                   )
               )
    )

    
    

# Define server logic required to draw a histogram
server <- function(input, output) {

    # # Define any Python packages needed for the app here:
    # PYTHON_DEPENDENCIES = c('matplotlib', 'seaborn','pandas', 'numpy', 'scipy','tqdm')
    # 
    # # ------------------ App virtualenv setup (Do not edit) ------------------- #
    # 
    # virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
    # python_path = Sys.getenv('PYTHON_PATH')
    # 
    # # Create virtual env and install dependencies
    # reticulate::virtualenv_create(envname = virtualenv_dir, python = python_path)
    # reticulate::virtualenv_install(virtualenv_dir, packages = PYTHON_DEPENDENCIES, ignore_installed=TRUE)
    # reticulate::use_virtualenv(virtualenv_dir, required = T)
    
    
    # Import python functions to R
    #reticulate::source_python('test.py')
  observeEvent(input$graph, {
    output$distPlot <- renderPlot({
      req(input$file1)
      
      df <- data.frame(read.csv(input$file1$datapath,
                                header = input$header,
                                sep = input$sep,
                                quote = input$quote))
      
      para_python <-input_priors_csv(df)
      
      s_c <- para_python$s_c
      w_c <- para_python$w_c
      mu_c <- para_python$mu_c
      sigma_c <- para_python$sigma_c
      
      
      errorMessage <- sanity_check(para_python)
      
      
      output$errorMessage <- renderText({ 
        paste(errorMessage)})
      
      
      
      DELTA <- 1.5 # Expected Performance Requirement (see Section 3)
      GAMMA <- 0.3 # Probability Requirement (see Section 3.5)
      
      ALPHA <- 0.05 # (1-ALPHA) indicates confidence in a statistical test (see Section 4)
      BETA <- 0.3 # (1-BETA) indicates power in a statistical test (see Section 4)
      
      N_MAX <- 100000 # Maximum experiment size (Algorithms 2-5)
      B <- 10000 # Simulation iterations (Algorithms 2 and 4)
      N_PARAL <- 1000 # Values of Ntr for parallel evaluation (memory constraint: B x N_PARAL x C)
      
      plotSelected <-{    
        
      if(input$type=='1'){
        plot1(s_c ,
              w_c ,
              mu_c,
              sigma_c,
              DELTA ,
              GAMMA,
              
              ALPHA ,
              BETA ,
              
              N_MAX,
              B ,
              N_PARAL)
        
      }else if(input$type=='2'){
        
        plot2(s_c ,
              w_c ,
              mu_c,
              sigma_c,
              DELTA ,
              GAMMA,
              N_MAX,
              B ,
              N_PARAL
              
        )
      }else if(input$type=='3'){
        
        plot3(s_c ,
              w_c ,
              mu_c,
              sigma_c,
              DELTA ,
              GAMMA,
              
              ALPHA ,
              BETA ,
              
              N_MAX,
              B ,
              N_PARAL)
        
      }else if(input$type=='4'){
        
        
        plot4(s_c ,
              w_c ,
              mu_c,
              sigma_c,
              DELTA ,
              GAMMA,
              
              ALPHA ,
              BETA ,
              
              N_MAX,
              B ,
              N_PARAL)
      }
      
        
        }
      
      output$distPlot <- renderPlot({
        print(grid.arrange(
          plotSelected,
          ncol=1))
        
      })
      

      
     
      #plot_function(df)
      # list(src = "plot.png")
    })
  })
  
    
    
}

plot1 <- function(s_c ,
                w_c ,
                mu_c,
                sigma_c,
                DELTA ,
                GAMMA,
                
                ALPHA ,
                BETA ,
                
                N_MAX,
                B ,
                N_PARAL ,
                
                ntr_set
                
                
                ){
  
  ntr_set <- seq(0,N_MAX,(N_MAX)/100)
  output <- stats_g_prob_d_expected_improvement(0,0.5,mu_c,sigma_c,s_c,w_c,ntr_set,0)
  
  perf_mean <- output$quant
  
  ################################
  ## equation 13 viz
 plot1<- ggplot() +
    geom_line(aes(x=ntr_set, y=perf_mean))+
    ylab("Expected Profit Improvement")+
    xlab("Size of Training Data")+
    ggtitle(TeX(r"(Equation 13 in Mean($V_{\pi_t - \pi_0}^{post}$))", bold=TRUE))+
    theme(plot.title = element_text(hjust = 0.5))
  gc() 
  
  
  return (plot1)
}

plot2 <- function(s_c ,
                w_c ,
                mu_c,
                sigma_c,
                DELTA ,
                GAMMA,
                N_MAX,
                B ,
                N_PARAL
                
){
  
  
  # Algorithm 1
  perf_no_training <- get_expected_perf(mu_c, sigma_c, s_c, w_c, 0)
  perf_no_training <- perf_no_training$mu_delta
  perf_limit <- get_expected_perf_limit(mu_c, sigma_c, s_c, w_c)
  
  print(paste("Expected policy performance with no training data:  ", round(perf_no_training,2)))
  print(paste("Expected policy performance limit (Equation 15):  ", round(perf_limit,2)))
  
  sample_d_exp <- solve_d_expected_improvement(DELTA,mu_c,sigma_c,s_c,w_c) 
  
  print(paste("Sample size for $",round(DELTA, 2), "-expected improvement"))
  print(paste("   Analytical Solution (exact):  ", sample_d_exp))
  
  ################################
  # Equation 24 Viz2
  # Quality of the Analytical Approximation (Figure 2 in the paper)
  ntr_set1 <- seq(1,N_MAX,(N_MAX)/100)
  
  prob_aprx <- stats_g_prob_d_expected_improvement(DELTA,GAMMA,mu_c,sigma_c,s_c,w_c,ntr_set1,0)
  prob_aprx <- prob_aprx$prob
  
  prob_sim <- stats_g_prob_d_expected_improvement(DELTA,GAMMA,mu_c,sigma_c,s_c,w_c,ntr_set1,B)
  prob_sim <- prob_sim$prob
  
  
  plot2<-ggplot()+
    geom_line(aes(x=ntr_set1, y=prob_sim, color='Simulation'),  linetype=1)+
    geom_line(aes(x=ntr_set1, y=prob_aprx, color='Approximation'),  linetype=2)+
    scale_color_manual(name = "Method", values = c("Simulation" = "red", "Approximation" = "darkblue"))+
    xlab("Size of Training Data")+
    ylab(TeX(sprintf(r"(Prob( $V_{\pi_t - \pi_0}^{post}\geq$$\DELTA = %.2f$))", DELTA)))+
    ggtitle("Equation 24: Analytical Approximation")+
    theme(plot.title = element_text(hjust = 0.5))
  gc() 
  
  
  return (plot2)
}

plot3<-function(s_c ,
                w_c ,
                mu_c,
                sigma_c,
                DELTA ,
                GAMMA,
                
                ALPHA ,
                BETA ,
                
                N_MAX,
                B ,
                N_PARAL){
  
  
 # sample_g_prob_d_epx_sim <- solve_g_prob_d_expected_improvement(DELTA,GAMMA,mu_c,sigma_c,s_c,w_c,N_MAX,N_PARAL,B)
  
  print(paste("Sample size for ", round(DELTA,2), "-expected ", round((1-GAMMA)*100),"%-probable improvement"))
  #print(paste("  Simulation-based solution:  ", sample_g_prob_d_epx_sim))
  
  
  
  # Algorithm 5
  #sample_g_prob_d_epx_approx <- solve_g_prob_d_expected_improvement(DELTA,GAMMA,mu_c,sigma_c,s_c,w_c,N_MAX) 
  
  #print(paste("  Analytical approximation-based solution: ", sample_g_prob_d_epx_approx))  
  
  ################################
  # Equation 37 Viz
  # Choose a reasonable BETA (for a given prior)
  ntr_set2 <- seq(1,N_MAX,(N_MAX)/100)
  power_set_aprx <- list()
  power_set_sim <- list()
  for (ntr in ntr_set2){
    power_set_aprx <- append(power_set_aprx, get_prob_ab_cert_aprx(ntr, N_MAX-ntr, ALPHA, mu_c, sigma_c, s_c, w_c))
    #     power_set_aprx.append(get_prob_ab_cert_aprx(ntr, N_MAX-ntr, ALPHA, mu_c, sigma_c, s_c, w_c,taylor=1))
    power_set_sim <- append(power_set_sim, get_prob_ab_cert_sim(ntr, N_MAX-ntr, ALPHA, mu_c, sigma_c, s_c, w_c, B))
  }
  
  plot3 <-ggplot()+
    geom_line(aes(x=ntr_set2, y=as.numeric(power_set_sim), color='Simulation'),  linetype=1)+
    geom_line(aes(x=ntr_set2, y=as.numeric(power_set_aprx), color='Approximation'),  linetype=2)+
    scale_color_manual(name = "Method", values = c("Simulation" = "red", "Approximation" = "darkblue"))+
    ggtitle("Equation 24: Analytical Approximation")+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Equation 37: Analytical Approximation")+
    ylab(TeX(sprintf(r"(Prob( $T_{ce}\geq z_{$ %.2f$}$ ))", 1-ALPHA)))
  gc() 
  
 return (plot3)
  
  }

plot4<-function(s_c ,
                w_c ,
                mu_c,
                sigma_c,
                DELTA ,
                GAMMA,
                
                ALPHA ,
                BETA ,
                
                N_MAX,
                B ,
                N_PARAL){
  
  
  
  ################################
  # Approximation for alpha, beta Policy Certification
  # Quality of the Analytical Approximation (Figure 5 in the paper)
  ntr_set3 <- seq(1, N_MAX, (N_MAX)/100)
  nce_set_aprx <- solve_nce_aprx(ntr_set3,ALPHA,BETA,mu_c,sigma_c,s_c,w_c)
  nce_set_sim <- rep(length(ntr_set3))
  for (n in seq(1, length(ntr_set3))){
    nce_set_sim[n] = solve_nce_sim(ntr_set3[n],N_MAX,ALPHA,BETA,mu_c,sigma_c,s_c,w_c,B)
  }
  nce_set_aprx[nce_set_aprx<0] <- N_MAX*2
  nce_set_sim[nce_set_sim<0] <- N_MAX*2
  
  
  plot4<-ggplot()+
    geom_line(aes(x=ntr_set3, y=as.numeric(nce_set_sim), color='Simulation'),  linetype=1)+
    geom_line(aes(x=ntr_set3, y=as.numeric(nce_set_aprx), color='Approximation'),  linetype=2)+
    scale_color_manual(name = "Method", values = c("Simulation" = "red", "Approximation" = "darkblue"))+
    ggtitle(TeX(r"(Approximation for $(\alpha,\beta)$ Policy Certification)"))+
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("Size of the Training Sample")+
    ylab("Size of the Certification Sample")+
    coord_cartesian(xlim = c(0, N_MAX), ylim = c(0, N_MAX))
  gc() 
  # 
  # 
  # 
  # ################################
  # # Algorithm 3
  # 
  result <- solve_ab_certification(ALPHA,BETA,mu_c,sigma_c,s_c,w_c,N_MAX)
  ntr_aprx <- result$ntr_opt
  nce_aprx <- result$nce_opt
  
  # result <- solve_ab_certification(ALPHA,BETA,mu_c,sigma_c,s_c,w_c,N_MAX,B=0,taylor=1)
  # ntr_taylor <- result$ntr_opt
  # nce_taylor <- result$nce_opt
  # 
  # result <- solve_ab_certification(ALPHA,BETA,mu_c,sigma_c,s_c,w_c,N_MAX,B)
  # ntr_sim <- result$ntr_opt
  # nce_sim <- result$nce_opt
  # 
   print(paste("Sample size for (a,b)-certification, confidence=", (1-ALPHA)*100, "%, power=", (1-BETA)*100, "%"))
  # 
  # print(paste("  Simulation-based solution:  N_tr=", ntr_sim, ", N_ce=", nce_sim, ", Total=", ntr_sim+nce_sim))
  # 
   print(paste("  Analytical approximation (main):  N_tr=", ntr_aprx, ", N_ce=", nce_aprx, ", Total=", ntr_aprx+nce_aprx))
  # 
return(plot4)
  }
# Run the application 
shinyApp(ui = ui, server = server)
