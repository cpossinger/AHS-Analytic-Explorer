library(skimr)
library(shinyFiles)
library(plyr)
library(dplyr)
library(magrittr)
library(broom)
library(shiny)
library(stringr)
library(survival)
library(survminer)
#library(plotly)
library(purrr)
library(ggplot2)
#library(skimr)
#library(survivalAnalysis)
library(shinythemes)
#library(tm)
library(rms)

make_regression_result_table <- function(results, x_values = NULL, retain_cols=c("base_name"), hr_w_ci = FALSE){
  
  var_betas = diag(results$var) # var is a covariates matrix of the beta coefficients
  
  z_betas = results$coefficients / sqrt(var_betas)
  
  pvalue = ( 1 - pnorm(abs(z_betas)) ) * 2 # two-tailed
  
  hazard_ratios         = exp( results$coefficients)
  #inverse_hazard_ratios = exp(-results$est_betas)
  
  alpha = 0.05
  
  conf_betas <- cbind(
    lcl = (results$coefficients-qnorm(1-alpha/2)*sqrt(var_betas) ),
    ucl = (results$coefficients+qnorm(1-alpha/2)*sqrt(var_betas) )
  )
  
  if( !is.null(x_values) && (length(x_values) == nrow(conf_betas)) )
    hazard_ratio_cis = exp(conf_betas*x_values)
  else
    hazard_ratio_cis = exp(conf_betas)
  
  sig =
    ifelse(pvalue < 0.001, 0.001,
           ifelse(pvalue < 0.01 , 0.01,
                  ifelse(pvalue < 0.05 , 0.05,
                         ifelse(pvalue < 0.1  , 0.1,
                                1))))
  
  output_table <- data.frame(
    term = names(results$coefficients),
    beta = round(results$coefficients, digits = 4),
    se   = sqrt(var_betas) %>% round(4), 
    beta.lcl = round(conf_betas[,"lcl"], digits = 3),
    beta.ucl = round(conf_betas[,"ucl"], digits = 3),
    pvalue = round(pvalue,              digits = 3),
    sig = sig,
    hr = round(hazard_ratios,        digits = 3),
    hr.lcl = round(hazard_ratio_cis[,"lcl"], digits = 3),
    hr.ucl = round(hazard_ratio_cis[,"ucl"], digits = 3),
    stringsAsFactors = F) %>% mutate(
      hr.sig  = ifelse( (hr.lcl < 1 & hr.ucl < 1) | (hr.lcl >= 1 & hr.ucl >= 1), "Yes", "No")
    ) 
  output_table %<>% mutate(hr_w_ci = format_ci(hr, hr.lcl, hr.ucl))
  
  for(colname in retain_cols){
    output_table[[colname]] <- results$model[[colname]]
  }
  
  rownames(output_table) <- output_table$term
  
  output_table %<>% left_join(
    lapply(names(results$assign), function(varname){
      data.frame(variable=varname, term = names(results$coefficients)[results$assign[[varname]]])
    }) %>% {do.call(rbind, .)}
  ) %>% select(variable, everything())
  
  return(output_table)
  
}

format_ci <- function(value, lcl, ucl){
  ci <- sprintf("%0.2f (%0.2f-%0.2f)", value, lcl, ucl)
  return(ci)
}


make_tidy_regression_result <- function(result, retain_cols = c("base_name"), custom_cols = NULL){
  tidy_result <- result %>% 
    make_regression_result_table %>%
    mutate_all(as.character) %>%  
    tidyr::pivot_longer(cols = names(.) %>% setdiff(c("term", "variable")), 
                        names_to = "result_type", 
                        values_to = "value") #%>%
  #mutate(name = paste(result$model$name_prefix, term, result_type, sep="."))
  for(colname in retain_cols){
    tidy_result[[colname]] <- result$model[[colname]]
  }
  if(!is.null(custom_cols)){
    tidy_result %<>% { cbind(., sapply(custom_cols, rep, nrow(.))) }
  }
  tidy_result$valuename <- tidy_result %>% select(one_of(c(retain_cols, names(custom_cols))), term, result_type) %>% apply(1, paste0, collapse=".") 
  
  return(tidy_result)
  
}

# Expects a data frame with variable, term, hr, hr.lcl, hr.ucl, hr_w_ci
make_forest_plot <- function(fp.data, breaks = c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)){
  fp <- 
    ggplot(data=fp.data, aes(x=term, y=hr, ymin=hr.lcl, ymax=hr.ucl)) +
    geom_errorbar(aes(col=variable)) + 
    geom_point(aes(col=variable)) + 
    geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1
    xlab("Covariates") + 
    ylab("HR (95% CI)") +
    theme_minimal() +
    theme(
      plot.margin = unit(c(1,10,1,1), "lines"),
      #legend.position = c(1.4,0.5),
      legend.position = "none",
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    
    annotate("segment", y = min(breaks),  yend = max(breaks), x = 0.5, xend=0.5) + 
    scale_y_continuous(breaks = breaks, expand = ) +  
    #geom_label(aes(y = max(breaks)*1.05, label=hr_w_ci, hjust=0), # use a white background
    #           size=3.33, label.padding = unit(0.5,"lines"), color="white", fill="white") + 
    geom_text(aes(y = max(breaks)*1.05, label=hr_w_ci, hjust=0), size=3.33) + 
    coord_flip(clip = "off") 
  
  return(fp)
}



# Define UI #### 
ui <- fluidPage(theme = shinytheme("flatly"), 
                
                # Application title
                titlePanel("Adventist Health Study Analytic Explorer"),
                
                tabsetPanel(
                  tabPanel("Cox Analysis",
                           # Sidebar with a slider input for number of bins 
                           sidebarLayout(
                             sidebarPanel(
                               # . . . Load Projects Button ####
                               fileInput("project_file","Load Project",accept = ".survproj"),
                               textInput("project_name", "Project Name:"),
                               actionButton("save_project", "Save Project"),
                               tags$hr(),
                               tabsetPanel(
                                 tabPanel("Load Data",
                                          # . . . Load Data Button ####
                                          textInput("dataset_path", label = "Datafile Path:"),
                                          actionButton("load_data", "Load Data File"),
                                 ),
                                 tabPanel("Select File",
                                          shinyFilesButton("select_data_file","Select Data File",title = "",multiple = FALSE,style = "margin-top: 15px"))
                               ),
                               # . . . File Path Output ####
                               fluidRow(
                                 column(6,
                                        tags$br(), tags$strong("Loaded Datafile: "), textOutput("loaded_dataset_path")),
                                 # . . . Number of Rows In Dataset ####
                                 column(6,
                                        tags$br(),tags$strong(textOutput("nrows"))),
                               ), 
                               tags$hr(),
                               # . . . Saved Models SI ####
                               selectInput("saved_models","Saved Models",choices = c()),
                               # . . . Load Model Button ####
                               actionButton("load_model","Load Model"),
                               # . . . Delete Model Button ####
                               actionButton("delete_model","Delete Model"),
                               fluidRow(
                                 column(6,
                                        # . . . Selected Model Formula####
                                        tags$br(), tags$strong("Selected Formula: "),textOutput("selected_formula")),
                                 column(6,
                                        # . . . Selected Filtering ####
                                        tags$br(), tags$strong("Selected Filtering: "),textOutput("selected_filtering"))                    
                               )),
                             mainPanel(
                               fluidRow(
                                 wellPanel(
                                   flowLayout(
                                     column(12,
                                            # . . . Dataset Variables####
                                            selectInput("cox_dataset_vars", label = "Dataset Variables", choices = c(), multiple = TRUE,selectize = FALSE,size = 10)),
                                     column(10,
                                            verticalLayout(
                                              # . . . Outcome Button ####
                                              actionButton("cox_outcome_button",label = "Outcome",style = "margin-top: 45px"),
                                              
                                              # . . . Start Time Button ####
                                              actionButton("cox_start_time_button",label = "Start Time",style = "padding:12px"),
                                              
                                              # . . . End Time Button ####
                                              actionButton("cox_end_time_button",label = "End Time"),
                                              
                                              # . . . Predictor Button ####
                                              actionButton("cox_predictor_button",label = "Predictor"))),
                                     
                                     column(12,
                                            # . . . Model Variables ####
                                            selectInput("cox_model_vars", label = "Model Variables", choices = c(), multiple = TRUE,selectize = FALSE,size = 10)),
                                     column(10,
                                            verticalLayout(
                                              
                                              # . . . Interact Button ####
                                              actionButton("cox_interact_button",label = "Interact",style = "margin-top: 45px"),
                                              
                                              # . . . Strata Button ####
                                              actionButton("cox_strata_button",label = "Strata",style = "padding:9px 21px"),
                                              
                                              # . . . Spline Button ####
                                              actionButton("cox_spline_button",label = "Spline",style = "padding:9px 21px"),
                                              
                                              # . . . Remove Button ####
                                              actionButton("cox_remove_button",label = "Remove")))
                                   ))),
                               
                               tags$hr(), 
                               fluidRow(
                                 splitLayout(
                                   column(6,
                                          # . . . Save Filtering ####
                                          textInput("filter",label = tags$strong("Filter Dataset"),value = ""),
                                          actionButton("filter_button","Filter"),
                                          actionButton("unfilter_button","Remove Filtering")),      
                                   column(6,
                                          # . . . Model Name Input ####
                                          textInput("cox_model_name",label = "Model Name"),
                                          # . . . Run Model ####
                                          actionButton("run_cox_model","Run Model"),
                                          # . . . Save Model ####
                                          actionButton("save_cox_model","Save Model")),
                                   
                                   
                                   column(6,
                                          tags$strong("Model:"),
                                          # . . . Model String####
                                          textOutput("cox_model_str"),
                                          tags$style(type="text/css", "#cox_model_str {white-space: pre-wrap;}"),
                                          textOutput("cox_model_error"),
                                          tags$style(type="text/css", "#cox_model_error {white-space: pre-wrap;}")
                                          
                                   ))
                               ),
                               tags$hr(),
                               tabsetPanel(
                                 tabPanel("Model Summary",
                                          tags$strong("Model Summary"),
                                          # . . . Model Summary Table ####
                                          tableOutput("cox_model_summary")
                                 ), 
                                 tabPanel("Forest Plot",
                                          # . . . Forest Plot ####
                                          plotOutput("cox_forest_plot")
                                 ), 
                                 tabPanel("Proportional Hazard Assumption",
                                          # . . . PHA Test Table ####
                                          tableOutput("cox_fit_test"),
                                          
                                          # . . . PHA Test Plot ####
                                          plotOutput("cox_fit_plot")
                                 ),
                                 tabPanel("Martingale Residuals",
                                          # . . . Residual Plot #### 
                                          plotOutput("cox_residual_plot"))), 
                             )
                           ) 
                  ),
                  tabPanel("Data Prep",
                           textAreaInput("data_prep",label = "Enter Data Prep Function",value = "function(d){\nreturn(d)\n}",width = "600px",height = "600px"),
                           actionButton("apply_data_prep",label = "Apply Data Prep Function"),
                           actionButton("remove_data_prep",label = "Remove Data Prep Function")),
                  # . . . Variables Tab ####
                  tabPanel(title = "Variables",
                           selectInput("varname", label = "Variable", choices = c(), multiple = FALSE),
                           textOutput("var_desc"),
                           plotOutput("var_hist"),
                           DT::dataTableOutput("var_stats")
                  )
                )
)


# Define server logic required to draw a histogram

server <- function(input, output, session) {
  # Cox Proportional Survival Analysis Tab ####
  # : Data Prep Tab ####
  # . . . Text Area Input Function ####
  
  # . . . Apply Data Prep Button ####
  observeEvent(input$apply_data_prep,{
    proj$data_prep <- input$data_prep
  })
  
  # . . . Remove Data Prep Button ####
  observeEvent(input$remove_data_prep,{
    updateTextAreaInput(session,"data_prep",value = "function(d){\nreturn(d)\n}")
    proj$data_prep <- "function(d){\nreturn(d)\n}"
  }) 
  # : Sidebar ####
  
  # : Projects ####
  proj <- reactiveValues()
  
  # . . . Load project file #### 
  observeEvent(input$project_file, {
    message("Loading project.")
    new_proj <- readRDS(input$project_file$datapath)
    proj$loaded_dataset_path <- ""
    updateTextInput(session = session, inputId = "dataset_path", value = new_proj$dataset_path)
    updateTextInput(session = session, inputId = "project_name", value = new_proj$name)
    updateTextAreaInput(session,"data_prep",value = new_proj$data_prep)
    proj$data_prep <- new_proj$data_prep
    message("  Loaded dataset_path:", proj$loaded_dataset_path)
    proj$saved_models <- (new_proj$saved_models)
    message("  Loaded models:", names(proj$saved_models))
    message("Project attributes",new_proj)
  })
  
  # . . . Save Project Button #### 
  observeEvent(input$save_project, {
    message("Save Project Button function")
    
    if(length(input$project_name) > 0){
      dest_path <- file.path(getwd(), paste0(input$project_name, ".survproj"))
      message("Saving project to:", dest_path)
      saveRDS( list(
        name = input$project_name,
        dataset_path = proj$loaded_dataset_path,
        saved_models = proj$saved_models,
        data_prep    = proj$data_prep
      ), dest_path)
    } else {
      message("Project file name required.")
    }
  })
  
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  
  shinyFileChoose(input,"select_data_file",roots = volumes,session = session)
  
  data_file <- reactive({
    if(!is.null(input$select_data_file))
      parseFilePaths(volumes,input$select_data_file)
  })
  
  observe({
    if(data_file()[,"datapath"] %>% as.character != "character(0)"){
      updateTextInput(session,"dataset_path",value = data_file()[,"datapath"] %>% as.character)
    }
  })
  
  # . . . Dataset Object #### 
  dataset <- reactiveVal()
  
  # . . . Load Dataset Button ####
  observeEvent({input$load_data},{
    message("Load Data Button Pressed or Changed Data_prep event")
    #selected_file_path <- data_file()[,"datapath"] %>% as.character  
    selected_file_path <- input$dataset_path
    if(file.exists(selected_file_path)){
      message("New data path exists: saving.")
      output$loaded_dataset_path <- renderText("Loading file ...")
      proj$loaded_dataset_path <- selected_file_path 
    } else {
      message("New data path not found:", selected_file_path, ". Doing nothing.")
      output$loaded_dataset_path <- renderText("File not found.")
    }
    if(!is.null(proj$loaded_dataset_path) && proj$loaded_dataset_path != ""){
      message("Loading data file:", proj$loaded_dataset_path)
      # TODO/FIXME wrap this in try()
      result <- try(readr::read_csv(proj$loaded_dataset_path))
      if(class(result)[1] == "try-error"){
        message("THERE WAS A TRY_ERROR")
        message(result)
        output$loaded_dataset_path <- renderText("Error loading file.")
      } else {
        message("Finished loading data file.")
        output$loaded_dataset_path <- renderText(input$dataset_path)
        message("Parsing and applying data prep function.")
        if(length(proj$data_prep) == 0){
          proj$data_prep <- "function(d){\nreturn(d)\n}"
        }
        message("dataprep text:", proj$data_prep)
        data_prep_func <- eval(parse(text = proj$data_prep))
        message("evaled data_prep_func:")
        print(data_prep_func)
        prepped_result <- data_prep_func(result)
        dataset(prepped_result)
        print(dataset() %>% head)
      }
    } else {
      
      output$dataset_path <- renderText("")
      return(NULL)
    }
  })
  
  
  
  # . . . File # of Rows ####
  observe({
    output$nrows <- renderText({
      paste0("Number of Rows: ",dataset() %>% nrow)
    })
  })
  
  # . . . Filter Button: set filtered_dataset ####
  
  observeEvent(input$filter_button,{
    #selected_model <- proj$saved_models[[input$saved_models]]
    #selected_model$filtering <- input$filter
    model_list$filtering <- input$filter
    
  })
  observeEvent(input$unfilter_button,{
    #selected_model <- proj$saved_models[[input$saved_models]]
    #selected_model$filtering <- "" 
    model_list$filtering <- NULL 
    updateTextInput(session,"filter",value = "")
  })
  
  # : Saved Models ####
  # . . . Model Names ####
  proj$saved_models <- list()
  
  # : Model List Variables ####
  # . . . Update Cox Dataset Vars ####
  observeEvent(dataset(),{
    dataset_vars <- names(dataset()) 
    if(is.null(dataset_vars)) dataset_vars <- character()
    message("dataset variables: ", paste0(dataset_vars, collapse=", "))
    updateSelectInput(session,"cox_dataset_vars",choices=dataset_vars) 
  }, ignoreNULL = FALSE)
  
  # : Model Attributes ####
  model_list <- reactiveValues(outcome = NULL,start_time = NULL,end_time = NULL,predictor = list(),data = NULL, filtering = list())
  
  # . . . Update Model Vars ####
  observeEvent(reactiveValuesToList(model_list), {
    message("Update Model Vars function")
    
    model_display <- c(
      paste0("Outcome: ", model_list$outcome),
      paste0("Start Time: ", model_list$start_time),
      paste0("End Time: ", model_list$end_time),
      paste0("Predictor: ", model_list$predictor)
    )
    
    message("Model Display:", str(model_display))
    
    updateSelectInput(session,"cox_model_vars",
                      choices=model_display)
  })
  
  # . . . Outcome Button
  observeEvent(input$cox_outcome_button,{
    
    message("Outcome button pressed.")
    
    model_list$outcome <- input$cox_dataset_vars[1]
    
    message("Outcome variable:", model_list$outcome)
    
  })
  
  # . . . Start Time Button ####
  observeEvent(input$cox_start_time_button,{
    model_list$start_time <- input$cox_dataset_vars[1]
  })
  
  # . . . End Time Button ####
  observeEvent(input$cox_end_time_button,{
    model_list$end_time <- input$cox_dataset_vars[1]
  })
  
  # . . . Predictor Button ####
  observeEvent(input$cox_predictor_button,{
    model_list$predictor %<>% c(input$cox_dataset_vars) %>% unique
  })
  
  # . . . Interact Button ####
  observeEvent(input$cox_interact_button,{ 
    
    model_list$predictor %<>% 
      c( paste0( input$cox_model_vars %>% str_remove(".*: ") %>% setdiff(""), collapse=":") ) %>% 
      setdiff("") %>% unique
    
  })
  
  # . . . Strata Button ####
  observeEvent(input$cox_strata_button,{
    message("Strata button pressed.")
    
    # Replace selected rows with strata()-fied values
    # e.g. c("cat","dog") with c("strata(cat)", "strata(dog)")
    # using str_replace("^(cat|dog/)$", "strata(\\1)")
    var_pattern <- paste0("^(",paste0(input$cox_model_vars %>% str_remove(".*: "), collapse="|"),")$")
    
    message("str_replace pattern:", var_pattern)
    
    model_list$predictor %<>%  str_replace(var_pattern, "strata(\\1)")
    
  })
  
  # . . . Spline Button ####
  observeEvent(input$cox_spline_button,{
    message("Spline button pressed.")
    
    # Replace selected rows with rcs()-fied values
    # e.g. c("cat","dog") with c("rcs(cat)", "rcs(dog)")
    # using str_replace("^(cat|dog)$", "rcs(\\1)")
    var_pattern <- paste0("^(",paste0(input$cox_model_vars %>% str_remove(".*: "), collapse="|"),")$")
    
    message("str_replace pattern:", var_pattern)
    
    model_list$predictor %<>%  str_replace(var_pattern, "rcs(\\1)")
    
  })
  
  
  # . . . Remove Button #### 
  observeEvent(input$cox_remove_button,{
    message("Remove button pressed.")
    
    for(selection in input$cox_model_vars){
      message("Removing:", str(selection))
      if(selection %>% str_detect("^Outcome: ")){
        model_list$outcome %<>% setdiff(selection %>% str_remove("^.+: "))
      } else if (selection %>% str_detect("^Predictor: ")){
        model_list$predictor %<>% setdiff(selection %>% str_remove("^.+: "))
      } else if (selection %>% str_detect("^Start Time: ")){
        model_list$start_time %<>% setdiff(selection %>% str_remove("^.+: "))
      } else if (selection %>% str_detect("^End Time: ")){
        model_list$end_time   %<>% setdiff(selection %>% str_remove("^.+: "))
      }
    }
    
  })
  
  # : Make Model Str and Object ####
  
  # . . . Model Formula ####
  
  make_coxph_formula <- function(model){
    paste0("Surv(",paste0(c(model$start_time, model$end_time,model$outcome), collapse=", "), ")",
           " ~ ",paste0(model$predictor,collapse = " + "))
    
  }
  # . . . Model Call ####
  make_coxph_call <- function(model){
    
    formula_str <- make_coxph_formula(model)
    
    # When we code the ability to load a dataset, this should just become data_str <- "dataset"
    data_str    <- "dataset()"
    
    # paste0 stuff to make coxph() function call, including filtering.
    
    if(length(model$filtering) == 0){
      filter_str <- ""
    } else {
      filter_str <- paste0("%>%","filter(",model$filtering,")")
    }
    message("Filtering Str: ",filter_str %>% length) 
    call_str <- paste0("coxph(", formula_str, ",", "data=", data_str, filter_str ,")")
    
    message("coxph() call string: ", call_str) 
    return(call_str)
    
  }
  # . . . Set Model Str to Model Object #### 
  observeEvent(reactiveValuesToList(model_list),{
    message("model list changed")
    
    if(length(model_list$outcome) == 0){ proj$model_str <- "Error: Outcome required"}
    else if(length(model_list$end_time) == 0){ proj$model_str <- "Error: End time required."}
    else if(length(model_list$predictor) == 0){ proj$model_str <- "Error: Predictor required."}
    else{
      model_formula <- make_coxph_formula(model_list) 
      model_call <-make_coxph_call(model_list)  
      
      message("Cox model string: ", model_call)
      
      proj$model_str <- model_call 
    }
  })
  
  
  # exmaple: coxph(Surv(time,event) ~ x1  + x2 + x3, data = hobbs)
  # . . . Render Model String Text ####
  observeEvent(reactiveValuesToList(model_list),{
    output$cox_model_str <- renderText({
      if(is.null(proj$model_str)) return(NULL)
      proj$model_str
    })
    message("cox model string render: ",proj$model_str)
  })
  
  
  # . . . Update Model Selection List ####
  
  observeEvent(proj$saved_models,{
    message("Update Model Selection List function")
    updateSelectInput(session,"saved_models",choices=names(proj$saved_models))
  })
  
  
  
  # . . . Save Button #### 
  observeEvent(input$save_cox_model,{
    
    if(input$cox_model_name == ""){
      message("No Name")
      showModal(
        modalDialog(
          title = "Missing Model Name", 
          "Model name required.",
          footer = tagList(
            modalButton("Ok")
          )
        )
      )
      
    }
    
    if(input$cox_model_name %in% names(proj$saved_models)){
      message("Common Name Found")
      showModal(
        modalDialog(
          title = "Overwrite existing model?", 
          'Do you wish to overwrite the model named "', input$cox_model_name, '"?',
          footer = tagList(
            actionButton("overwrite_model", "Yes"),
            modalButton("Cancel")
          )
        )
      )
    }
    
    cox_save_model()
    
    
  })
  
  observeEvent(input$overwrite_model,{
    removeModal()
    cox_save_model()
  })
  
  
  
  # . . . Save Model function ####
  cox_save_model <- function(){
    message("cox_save_model function")
    # Find position of model name, or return position for new model
    pos <- which(names(proj$saved_models) == input$cox_model_name)
    if(length(pos) == 0) {pos <- length(proj$saved_models) + 1}
    
    message("pos:", pos)
    message("proj$saved_models:", str(proj$saved_models))
    
    proj$saved_models[[pos]] <- list(
      cox_model_name = input$cox_model_name,
      outcome = model_list$outcome,
      start_time = model_list$start_time,
      end_time = model_list$end_time,
      predictor = model_list$predictor,
      filtering = model_list$filtering
    )
    make_coxph_call(proj$saved_models[[pos]]) 
    # Name the list element
    names(proj$saved_models)[pos] <- input$cox_model_name
  }
  
  # . . . Set Selected Model Strs ####    
  observeEvent(input$saved_models,{
    if(length(proj$saved_models) == 0 ){
      formula_str <- "No selected model."
      filter_str  <- "No selected model."
    } else {
      selected_model <- proj$saved_models[[input$saved_models]]
      formula_str <- make_coxph_formula(proj$saved_models[[input$saved_models]]) 
      filter_str  <- selected_model$filtering
    }
    output$selected_formula   <- renderText({ paste0(formula_str) })
    output$selected_filtering <- renderText({ paste0(filter_str) })
    
  }) 
  # . . . Delete Button #### 
  observeEvent(input$delete_model,{
    if(input$available_projects == ""){
      proj$saved_models[[input$available_models]] <- NULL
    }else{ 
      cox_projects$available_projects[[input$available_projects]][[input$available_models]] <- NULL
      print(cox_projects$available_projects[[input$available_projects]])
    }
    formula_str <- "No selected model."
    filter_str  <- "No selected model."
    output$selected_formula   <- renderText({ paste0(formula_str) })
    output$selected_filtering <- renderText({ paste0(filter_str) })
  })
  
  cox_model_obj <- reactiveVal()
  
  # . . . Run Cox Model Button ####
  observeEvent(input$run_cox_model,{
    model.obj <- try(eval(parse(text = proj$model_str)))
    if(class(model.obj) == "try-error"){
      output$cox_model_error <- renderText(model.obj)
    } else {
      output$cox_model_error <- renderText("")
      cox_model_obj(model.obj) 
    }
  })
  
  # . . . Load Model Button #### 
  observeEvent(input$load_model,{
    selected_model <- proj$saved_models[[input$saved_models]]
    
    message("Loading model:",selected_model$cox_model_name)
    
    updateTextInput(session = session, inputId = "cox_model_name", value = selected_model$cox_model_name)
    proj$model_str <- make_coxph_call(selected_model)
    model_list$outcome <- selected_model$outcome
    model_list$start_time <- selected_model$start_time
    model_list$end_time <- selected_model$end_time
    model_list$predictor <- selected_model$predictor
    model_list$filtering <- selected_model$filtering
    updateTextInput(session,"filter",value = selected_model$filtering)
  })
  
  
  
  # :   Summary Tabs ####    
  
  # . . . Model Summary ####
  output$cox_model_summary <- renderTable({
    message("Render model summary function")
    if(is.null(cox_model_obj())) return(NULL)
    cox_model_obj() %>% 
      tidy(conf.int=TRUE) %>% 
      mutate(across(c("conf.low", "conf.high"), exp), hr=exp(estimate)) %>% 
      select(term, beta=estimate, SE=std.error, hr, hr.lcl=conf.low, hr.ucl=conf.high, p.value)
  })
  
  # . . . Test Model Proportional Hazards Assumption ####
  cox_fit_test <- reactive({
    if(is.null(cox_model_obj())) return(NULL)
    cox.zph(cox_model_obj())
  })
  # . . . Render Test Model Proportional Hazards Assumption ####
  output$cox_fit_test <- renderTable({
    if(is.null(cox_model_obj())) return(NULL)
    cox_fit_test()$table %>% as_tibble
  })
  # . . . Plot Test Model Proportional Hazards Assumption ####
  output$cox_fit_plot <- renderPlot({
    if(is.null(cox_model_obj())) return(NULL)
    ggcoxzph(cox_fit_test() ) 
  })
  # . . . Predictions vs Martingale Residuals ####
  output$cox_residual_plot <- renderPlot({
    if(is.null(cox_model_obj())) return(NULL)
    ggcoxdiagnostics(cox_model_obj())
  }) 
  # . . . Model Hazard Ratios Plot ####
  output$cox_forest_plot <- renderPlot({
    if(is.null(cox_model_obj())) return(NULL)
    result_table <- cox_model_obj() %>% 
      make_regression_result_table %>% 
      filter(!str_detect(term, "^rcs\\("))
    if(nrow(result_table) > 0)
      result_table %>% make_forest_plot 
  })
  
  # : Variables Tab #### 
  # . . Variable Selection ####
  observe({
    updateSelectInput(session,"varname",choices=names(dataset()))
  })
  output$var_hist <- renderPlot({
    if(!is.null(dataset())){
      if(dataset()[[input$varname]] %>% is.numeric == TRUE){
        print(dataset()$input$varname %>% is.factor())
        print("if plot")
        var_plot <- ggplot(data = dataset(), 
                           aes(x = get(input$varname))) + 
          geom_histogram(aes(y = ..density..),fill = "gray82") +  
          xlab(input$varname)+
          geom_density(alpha = 0.2,fill = "#00FF78")
      } else{
        print("else plot")
        var_plot <- ggplot(data = dataset(), 
                           aes(x = get(input$varname))) + 
          geom_bar(fill = "#99ffcc")+ 
          xlab(input$varname)
        
      } 
      
      var_plot
      #ggplotly(var_plot)
    }
  })
  
  output$var_stats <- DT::renderDataTable(selection = "single",options = list(dom = "t"),{
    dataset() %>% skim %>% select(-complete_rate,-character.max,-character.min,-character.empty,-character.whitespace)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

# observeEvent(cox_projects$available_projects,{
#     message("Update Project Selection List function")
#     message("names avaiable projects: ", names(cox_projects$available_projects))
#     updateSelectInput(session,"available_projects",choices=names(cox_projects$available_projects))
#                      
#     
# })
# observeEvent(input$save_project,{
#     message("length of available projects",input$available_projects)
#     if(input$available_projects == ""){
#         message("No Name")
#         showModal(
#             modalDialog(
#                 title = "Save Project", 
#                 "Name Project",
#                 footer = tagList(
#                     textInput("name_project",label = ""),
#                     modalButton("Done")
#                 )
#             )
#         )
#         
#     }
#     
#     else{
#         message("Common Name Found")
#         showModal(
#             modalDialog(
#                 title = "Save Project", 
#                 'Do you wish to overwrite the project named "', input$available_projects,"or make a new project?", 
#                 footer = tagList(
#                     actionButton("overwrite_project", "Overwrite Project"),
#                     actionButton("new_project","New Project"),
#                     modalButton("Cancel")
#                 )
#             )
#         )
#     }
# }) 
# observeEvent(input$name_project,{
#    cox_projects$available_projects[[1]] <- cox_models$available_models
#    print(cox_models$available_models)
#    print("seperate")
#     
#    names(cox_projects$available_projects)[1] <- input$name_project
#    print(cox_projects$available_projects[[1]])   
# })
# observeEvent(input$overwrite_project,{
#     pos <- which(names(cox_projects$available_projects) == input$available_projects)
#     message("cox models: ", names(cox_projects$available_projects))
#     message("input projects: ", input$available_projects)
#     
#     message("pos: ",pos)
#     
#     cox_projects$available_projects[[pos]] <- cox_models$available_models
#     
# })
# observeEvent(input$new_project,{
#     showModal(
#         modalDialog(
#             title = "Save New Project",
#             "Name your new project",
#             footer = tagList(
#                 textInput("name_new_project",label = ""),
#                 modalButton("Done")
#             )
#         )
#     )
#     
# })
# observeEvent(input$name_new_project,{
#     
#     pos <- which(names(cox_projects$available_projects) == input$available_projects)
#     cox_projects$available_projects[[pos + 1]] <- cox_models$available_models
#     
#     names(cox_projects$available_projects)[pos+1] <- input$name_new_project 
# })
# observeEvent({input$load_project;input$delete_model},{
#     print(names(cox_projects$available_projects[[input$available_projects]]))
#     updateSelectInput(session,"available_models",choices=names(cox_projects$available_projects[[input$available_projects]]))
# })