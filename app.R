## Author: Italo Duran
## duran01@bu.edu
## BU BF591
## Assignment 7

# Welcome to R Shiny. All that glitters is not gold.
library(tidyverse)
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(dplyr)
library(ggrepel)

# Define UI for application that draws a histogram

ui <- fluidPage(theme = shinythemes::shinytheme("slate"),
                h1("BF591 Assignment 7"),
                h3("Italo Duran"),
                h4("Visualize differential expression results"),
                tags$h5(HTML("To use this application, download the CSV <b>deseq_res.csv</b> from the data folder in my github:")),
                tags$a(href="https://github.com/BF591-R/bf591-assignment-7-imd9/tree/main/data",target="_blank",rel="noopener noreferrer","Click here for 'My Github' CSV sample file!"),
                HTML("<br><br>"),
                sidebarLayout(sidebarPanel(
                                fileInput("file1","Load differential expression results", accept = ".csv", placeholder = "deseq_res.csv"),
                                HTML(paste(rep("<p>A volcano plot can be generated with <b>'log<sub>2</sub> fold-change'</b> on the x-axis and <b>'p-adjusted'</b> on the y-axis.</p>"), collapse = "")),
                                radioButtons("button1", "Choose the column for the x-axis",choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"),selected = "log2FoldChange"),
                                radioButtons("button2", "Choose the column for the y-axis",choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"),selected = "padj"),
                                colourInput("col1", label = "Base point color" ,"#138086"),
                                colourInput("col2", label = "Highlight point color" ,"#EEB462"),
                                sliderInput(inputId = "slider_p", min = -300, max = 0,label = "Select the magnitude of the p adjusted coloring:", value = -150, step = 1),
                                submitButton("Plot", icon("r-project",class="fab fa-r-project fa-1x"), width = '100%')),
                              mainPanel(tabsetPanel(tabPanel("Plot", plotOutput("volcano")),tabPanel("Table", tableOutput("table"))))))

##################### Define server logic required to draw a histogram########################

server <- function(input, output, session) {
    #####' load_Data ##########
    #' @details Okay this one is a little weird but bear with me here. This is 
    #' still a "function", but it will take no arguments. The `reactive({})` bit 
    #' says "if any of my inputs (as in, input$...) are changed, run me again". 
    #' This is useful when a user clicks a new button or loads a new file. In 
    #' our case, look for the uploaded file's datapath argument and load it with 
    #' read.csv. Return this data frame in the normal return() style.
    
    load_data <- reactive({
        req(input$file1) #this function makes the program not run until a file is given 
        f=input$file1
        if (is.null(f)){return(NULL)} #if there's no file then return nothing
        else{datafile=read.csv(f$datapath, header= TRUE, sep=",")}%>% 
            rename(Gene = X)%>% #once the line above is executed, this renames the X values as Gene 
            return() }) #datafile
        
    ########' Volcano plot ################
    #' @param dataf The loaded data frame.
    #' @param x_name The column name to plot on the x-axis
    #' @param y_name The column name to plot on the y-axis
    #' @param slider A negative integer value representing the magnitude of
    #' p-adjusted values to color. Most of our data will be between -1 and -300.
    #' @param color1 One of the colors for the points.
    #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
    #' @return A ggplot object of a volcano plot
    #' @details I bet you're tired of these plots by now. Me too, don't worry.
    #' This is _just_ a normal function. No reactivity, no bells, no whistles. 
    #' Write a normal volcano plot using geom_point, and integrate all the above 
    #' values into it as shown in the example app. The testing script will treat 
    #' this as a normal function.
    ######' !!sym() may be required to access column names in ggplot aes().#########
    #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
    
    volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
        p<-ggplot(data = dataf,aes(x = !!sym(x_name), y=-log10(!!sym(y_name)))) + 
            geom_point(aes(color = padj< 1*10^(slider)))+
            labs( color = str_glue('{y_name} 1 x 10^ {slider}'))+ #for labels to follow the selected button,like f' in python.
            scale_color_manual(values = c(color1, color2)) + theme_classic() +
            theme(panel.background = element_rect(fill = "grey23"),
                  panel.border = element_blank(),
                  panel.grid.major = element_line(color = "grey35"),
                  panel.grid.minor = element_line(color = "grey25"),
                  plot.background = element_rect(fill = "grey75",colour="grey75"),
                  plot.margin = margin(0.5, 0.5, 0.3, 0.5, "cm"),
                  legend.key = element_rect(fill = "grey81", colour="grey81"),
                  legend.background = element_rect(fill = "grey75"),
                  legend.position = "bottom" ) # colour = "black",size = 25, size = 1,theme_classic(), plot.margin = margin(2, 2, 2, 2, "cm")
                                               # plot.background = element_rect(fill = "grey30")
        return(p)
    }
    #############' Draw and filter table ##############
    #' @param dataf Data frame loaded by load_data()
    #' @param slider Negative number, typically from the slider input.
    #' @return Data frame filtered to p-adjusted values that are less than 
    #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits displayed.
    #' @details Same as above, this function is a standard R function. Tests will 
    #' evaluate it normally. Not only does this function filter the data frame to 
    #' rows that are above the slider magnitude, it should also change the format 
    #' of the p-value columns to display more digits. This is so that it looks 
    #' better when displayed on the web page. I would suggest the function 
    #' `formatC()`
    #' @examples draw_table(deseq_df, -210)
    #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
    #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
    #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
    draw_table <- function(dataf, slider) {
        dataf %>%
        arrange(pvalue) %>% #took dataframe and arranged the pvalues
        filter(padj < 10^slider) %>% #filtered for the p-adjusted with the slider
        mutate(pvalue = formatC(.$pvalue, digits = 3, format = "g" ), # Used mutate to format both pvals & p_adjusted. digits(for numbr of decimal spaces)
               padj = formatC(.$padj, digits = 3, format = "g")) %>% # format can also use: d(integers), f,e,E,G,fg(real numbers), s(strings) 
        return() 
    }
    #' These outputs aren't really functions, so they don't get a full skeleton, 
    #' but use the renderPlot() and renderTabel() functions to return() a plot 
    #' or table object, and those will be displayed in your application.
    #' here we load the fields from the server and UI to connect it between them so it gives us an ouput on the app.
    output$volcano <- renderPlot(volcano_plot(dataf = load_data(), slider = input$slider_p, x_name=input$button1, y_name=input$button2, color1=input$col1, color2=input$col2)) # replace this NULL
    # Same here, just return the table as you want to see it in the web page
    output$table <- renderTable(draw_table(dataf = load_data(), slider = input$slider_p)) # replace this NULL
}
# Run the application
shinyApp(ui = ui, server = server)
