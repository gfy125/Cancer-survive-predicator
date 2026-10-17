library(shiny)
library(shinydashboard)
library(DT)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(glmnet)
library(stringr)
library(tools)
library(biomaRt)
library(plotly)

# Load TCGA data and pretrained model
load("SKCM.rda")       # contains query2
load("cox_trials.RData")   # contains cv_fit, lambda_min, etc.
clinical_followup <- read.csv("clinical_followup.csv")

options(shiny.maxRequestSize = 200*1024^2)

# UI Section
ui <- dashboardPage(
  skin = "yellow",
  
  dashboardHeader(
    title = "Skin Cancer Survival Predictor",
    titleWidth = 350
  ),
  
  dashboardSidebar(
    width = 250,
    sidebarMenu(
      menuItem("About", tabName = "about", icon = icon("info-circle")),
      menuItem("Calculator", tabName = "calculator", icon = icon("calculator")),
      menuItem("Model Explanation", tabName = "model", icon = icon("chart-line")),
      menuItem("Dataset Explorer", tabName = "explore", icon = icon("database"))
    )
  ),
  
  dashboardBody(
    # Custom CSS styling
    tags$head(
      tags$style(HTML("
        /* Custom branding and UI elements */
        .skin-yellow .body {
        background-color: #fff6e8;
        }
        .skin-yellow .main-header .logo {
          background-color: #ffbd59;
          color: #000000;
          font-weight: bold;
          font-size: 18px;
        }
        .skin-yellow .main-header .navbar {
          background-color: #ffbd59;
        }
        .skin-yellow .main-sidebar {
          background-color: #ffbd59 !important;
        }
        .skin-yellow .main-sidebar .sidebar-menu > li > a {
          color: #000000;
          font-weight: 500;
        }
        .skin-yellow .main-sidebar .sidebar-menu > li.active > a {
          background-color: #fdce87;
          color: black;
        }
        .skin-yellow .main-sidebar .sidebar-menu > li:hover > a {
          background-color: #fdce87;
          color: black;
        }
        
        /* Content area styling */
        .content-wrapper {
          background-color: #fff6e8;
        }
        
        /* Box styling */
        .box.box-solid.box-primary > .box-header {
          color: #000000;
          background-color: #ffbd59;
        }
        .box.box-solid.box-primary {
          border-color: #ffbd59;
          background-color: #ffeeb6;
        }
        
        /* Button styling */
        .btn-primary {
          background-color: #ffbd59 !important;
          border-color: #ffbd59 !important;
          color: #000000 !important;
          font-weight: bold;
          transition: all 0.3s ease;
        }
        .btn-primary:hover {
          background-color: #fdce87 !important;
          border-color: #fdce87 !important;
        }
        
        /* Text elements */
        h1, h2, h3, h4 {
          font-weight: bold;
          color: #333;
        }
        
        /* Results section styling */
        .results-highlight {
          background-color: #ffeeb6;
          padding: 15px;
          border-left: 5px solid #ffbd59;
          margin-bottom: 15px;
        }
        
        /* Risk level indicators */
        .risk-low {
          color: #4CAF50;
          font-weight: bold;
        }
        .risk-medium {
          color: #FF9800;
          font-weight: bold;
        }
        .risk-high {
          color: #F44336;
          font-weight: bold;
        }
        .progress-bar {
          background-color: #fdce87 !important;  /* Lighter blue */
          color: black;
        }
        .black-bordered-tables table, 
        .black-bordered-tables th, 
        .black-bordered-tables td {
          border: 1px solid black !important;
          border-collapse: collapse !important;
        }
      "))
    ),
    
    tabItems(
      # About tab
      tabItem(tabName = "about",
              fluidRow(
                column(
                  width = 12,
                  box(
                    width = NULL,
                    status = "primary",
                    solidHeader = TRUE,
                    title = "Welcome to the Skin Cancer Survival Predictor",
                    p("Skin cancer is the abnormal growth of skin cells, often caused by excessive exposure to ultraviolet (UV) radiation from the 
                    sun or tanning beds. It is the most common type of cancer and typically develops on areas of the skin that are frequently 
                    exposed to sunlight, such as the face, neck, and arms. There are three main types: basal cell carcinoma, squamous cell 
                    carcinoma, and melanoma, the most dangerous form, which can spread to other parts of the body if not detected early. 
                    Preventative measures include using sunscreen, wearing protective clothing, and avoiding sun exposure during peak hours. 
                    Early detection and treatment significantly improve outcomes."),
                    p("This calculator utilises melanoma patient demographics and gene expression data to predict survival probability and provides 
                    healthcare professionals with personalised insights that may help improve patient outcomes and overall well-being."),
                    
                    h4("How to Use This Tool", style = "margin-top: 30px; border-bottom: 2px solid #ffbd59; padding-bottom: 10px;"),
                    tags$ol(
                      tags$li("Navigate to the ", tags$b("Calculator"), " tab."),
                      tags$li("Input patient details, including:",
                              tags$ul(
                                tags$li("Age"),
                                tags$li("Gender"),
                                tags$li("Cancer Stage")
                              )
                      ),
                      tags$li("Upload the patient's gene expression data (.csv or .RData file) or use default values."),
                      tags$li("Click on the ", tags$b("Predict"), " button."),
                      tags$li("View the results including survival curve, risk assessment, and recommended lifestyle changes.")
                    ),
                    
                    h4("About the Model", style = "margin-top: 30px; border-bottom: 2px solid #ffbd59; padding-bottom: 10px;"),
                    p("This application uses a Cox Proportional Hazards model with LASSO regularisation trained on The Cancer Genome Atlas (TCGA) 
                    melanoma dataset. The model integrates both clinical and genomic features to provide a comprehensive risk assessment."),
                    p("For more detailed information about the model, including feature importance and performance metrics, 
                    please visit the ", tags$b("Model"), " tab.")
                  )
                )
              )
      ),
      
      # Calculator tab
      tabItem(tabName = "calculator",
              h2("Calculate Survival Probability", style = "margin-bottom: 25px;"),
              fluidRow(
                # Left column for demographics
                column(width = 6,
                       box(
                         title = "Patient Demographics",
                         solidHeader = TRUE,
                         width = NULL,
                         status = "primary",
                         height = "auto",
                         collapsible = TRUE,
                         
                         # Demographics input fields
                         div(
                           class = "form-group",
                           selectInput("gender", "Gender", 
                                       choices = c("Male", "Female"),
                                       selected = "Male",
                                       width = "100%"
                           )
                         ),
                         div(
                           class = "form-group",
                           numericInput("age", "Age (Years)", 
                                        value = 60, 
                                        min = 0, 
                                        max = 120,
                                        step = 1,
                                        width = "100%"
                           )
                         ),
                         div(
                           class = "form-group",
                           selectInput("stage", "Cancer Stage", 
                                       choices = c("Stage 0", "Stage I", "Stage IA", 
                                                   "Stage IB", "Stage II", "Stage IIA",
                                                   "Stage IIB", "Stage IIC", "Stage III", 
                                                   "Stage IIIA", "Stage IIIB", "Stage IIIC", 
                                                   "Stage IV", "Unknown"),
                                       selected = "Stage II",
                                       width = "100%"
                           )
                         )
                       ),
                       
                       # Key gene expression inputs
                       uiOutput("geneExpressionInputs")
                ),
                
                # Right column for file upload and predict button
                column(width = 6,
                       box(
                         title = "Gene Expression Upload",
                         solidHeader = TRUE,
                         status = "primary",
                         width = NULL,
                         height = "auto",
                         collapsible = TRUE,
                         
                         p("Upload a file containing gene expression values for additional genes. The file can be in CSV format or RData format."),
                         
                         # File upload input
                         fileInput("file", "Upload Gene Expression Data", 
                                   accept = c(".csv", ".txt", ".RData", ".rda", ".Rdata"),
                                   width = "100%"),
                         
                         tags$div(
                           class = "alert",
                           style = "background-color: #fcbb65; color: black; border-color: #fcbb65;",
                           HTML("<strong>Note:</strong> If no file is uploaded, default values will be used for genes.")
                         ),
                         
                         # Template download link
                         div(
                           style = "margin-top: -15px;",
                           downloadLink("download_template", "Download template file")
                         ),
                         
                         # Predict button
                         div(
                           style = "text-align: center; margin-top: 20px;",
                           actionButton("predict", "Predict Survival", 
                                        class = "btn-primary",
                                        icon = icon("calculator"),
                                        style = "width: 80%; font-size: 16px; padding: 10px;")
                         )
                       )
                )
              ),
              
              # Results section (conditional on predict button being clicked)
              conditionalPanel(
                condition = "input.predict > 0",
                h3("Survival Prediction Results", style = "margin-top: 30px; margin-bottom: 20px;"),
                
                fluidRow(
                  column(width = 8,
                         box(
                           title = "Survival Probability Curve",
                           solidHeader = TRUE,
                           status = "primary",
                           width = NULL,
                           plotlyOutput("survivalPlot", height = "400px")
                         )
                  ),
                  
                  column(width = 4,
                         box(
                           title = "Risk Assessment",
                           solidHeader = TRUE,
                           status = "primary",
                           width = NULL,
                           
                           div(class = "results-highlight",
                               h4("Risk Score:"),
                               textOutput("riskScore"),
                               p(HTML("<span style='color: #b34700; font-weight: bold;'>See Below for Risk Distribution</span>")),
                               
                               h4("Risk Group:"),
                               htmlOutput("riskGroup"),
                               
                               h4("5-Year Survival Probability:"),
                               textOutput("survivalProb")
                           )
                         )
                  )
                ),
                
                # New separate section for Risk Distribution
                fluidRow(
                  column(width = 12,
                         box(
                           title = "Risk Distribution",
                           solidHeader = TRUE,
                           status = "primary",
                           width = NULL,
                           
                           plotlyOutput("riskBoxPlot", height = "250px"),
                           p("This shows where the patient's risk score (red dot) falls within the overall distribution of risk scores from the training cohort. The dashed lines indicate thresholds for risk groups.")
                         )
                  )
                ),
                
                fluidRow(
                  column(width = 12,
                         box(
                           title = "Lifestyle Recommendations",
                           solidHeader = TRUE,
                           status = "primary",
                           width = NULL,
                           
                           uiOutput("recommendations")
                         )
                  )
                )
              )
      ),
      
      # Model tab
      tabItem(tabName = "model",
              h2("Model Explanation", style = "margin-bottom: 25px;"),
              
              fluidRow(
                column(width = 12,
                       box(
                         title = "About the Cox Proportional Hazards Model",
                         solidHeader = TRUE,
                         status = "primary",
                         width = NULL,
                         
                         p("This survival calculator uses gene expression and clinical data from melanoma patients (SKCM) to build a survival prediction model using the Cox proportional hazards model. 
                         The goal is to identify genes and clinical features associated with patient survival, and to predict a risk score for each patient."),
                         
                         p("The Cox model is a semi-parametric regression model commonly used in medical research for analyzing survival data. It estimates the effect of different variables on the time to event (in this case, time to death or last follow-up)."),
                         
                         p("LASSO regularisation was used to select the most predictive features while avoiding overfitting. This approach helps identify the key genes and clinical variables that most strongly associate with survival outcomes.")
                       )
                )
              ),
              
              fluidRow(
                column(width = 12,
                       box(
                         title = "Feature Importance Analysis",
                         solidHeader = TRUE,
                         status = "primary",
                         width = NULL,
                         
                         h4("Top Genes"),
                         plotlyOutput("geneImportancePlot", height = "500px", width = "auto"),
                         p("This bar plot visualises the most influential genes identified by the LASSO-regularised Cox proportional hazards model, 
                           based on their non-zero coefficients. Positive coefficients (red bars) indicate genes whose higher expression is 
                           associated with increased risk (worse prognosis), while negative coefficients (blue bars) indicate genes whose 
                           higher expression is associated with decreased risk (better prognosis)."),
                         
                         h4("Clinical Features"),
                         plotlyOutput("clinicalImportancePlot", height = "400px"),
                         p("This bar plot visualises the most influential clinical demographics identified by the LASSO-regularised Cox proportional 
                           hazards model. Positive coefficients indicate features associated with increased risk, while negative coefficients 
                           indicate features associated with decreased risk.")
                       )
                )
              ),
              
              fluidRow(
                column(width = 6,
                       box(
                         title = "Risk Score Distribution",
                         solidHeader = TRUE,
                         status = "primary",
                         width = NULL,
                         plotlyOutput("riskDistributionPlot", height = "350px"),
                         p("This histogram shows the distribution of risk scores across all patients in the training cohort. 
                           The risk score is a linear predictor from the Cox proportional hazards model, with higher values indicating greater risk.")
                       )
                ),
                
                column(width = 6,
                       box(
                         title = "Model Performance Metrics",
                         solidHeader = TRUE,
                         status = "primary",
                         width = NULL,
                         
                         div(
                           style = "padding: 15px;",
                           tableOutput("metricsTable")
                         ),
                         
                         div(
                           p("This metric indicates the performance of our Cox proportional hazards model:"),
                           tags$ul(
                             tags$li(strong("Concordance Index (C-index):"), " Measures the model's ability to correctly rank patients by their risk of an event. Values range from 0.5 (no discrimination) to 1.0 (perfect discrimination).")
                           )
                         )
                       )
                )
              ),
              
              # Simplified PCA Plot Section with just one plot
              fluidRow(
                column(width = 12,
                       box(
                         title = "PCA Analysis of Gene Expression Data",
                         solidHeader = TRUE,
                         status = "primary",
                         width = NULL,
                         
                         tabsetPanel(
                           id = "pca_tab",
                           tabPanel(HTML("<span style='color: #b34700; font-weight: bold;'>PC1 vs PC2</span>"), value = "pc1_2"),
                           tabPanel(HTML("<span style='color: #b34700; font-weight: bold;'>PC2 vs PC3</span>"), value = "pc2_3"),
                           tabPanel(HTML("<span style='color: #b34700; font-weight: bold;'>PC3 vs PC4</span>"), value = "pc3_4")
                         ),
                         
                         plotlyOutput("pcaPlot", height = "500px"),
                         
                         p("This Principal Component Analysis (PCA) plot visualizes how melanoma samples cluster based on their gene expression profiles.
                            Each point represents a patient, with the axes representing the principal components that capture the most variation in the data.
                            Clustering patterns can reveal potential subgroups within melanoma patients that might correlate with clinical outcomes.")
                       )
                )
              ),
              # 
              # fluidRow(
              #   column(width = 12,
              #          box(
              #            title = "Information",
              #            solidHeader = TRUE,
              #            status = "warning",
              #            width = NULL,
              #            collapsible = TRUE,
              #            collapsed = TRUE,
              #            verbatimTextOutput("debugOutput")
              #          )
              #   )
              # )
      ),
      
      # Dataset Explorer tab
      tabItem(tabName = "explore",
              h2("TCGA Melanoma Dataset Explorer", style = "margin-bottom: 25px;"),
              
              fluidRow(
                column(width = 12,
                       box(
                         title = "Dataset Overview",
                         solidHeader = TRUE,
                         status = "primary",
                         width = NULL,
                         
                         p("This application uses data from The Cancer Genome Atlas (TCGA) Skin Cutaneous Melanoma (SKCM) project, 
                           which provides comprehensive molecular profiling of melanoma patients."),
                         
                         h4("Summary Statistics", style = "margin-top: 20px;"),
                         
                         div(class = "black-bordered-tables",
                             
                             h5(HTML("<span style=' font-weight: bold;'>Gender Distribution</span>")),
                             tableOutput("genderDist"),
                             
                             h5(HTML("<span style=' font-weight: bold;'>Stage Distribution</span>")),
                             tableOutput("stageDist"),
                             
                             h5(HTML("<span style=' font-weight: bold;'>Age Statistics (Years)</span>")),
                             tableOutput("ageStats")
                         )
                       )
                       
                )
              ),
              
              fluidRow(
                tabBox(
                  title = "Browse Dataset",
                  width = 12,
                  id = "datasetTabs",
                  
                  tabPanel("Clinical Data",
                           div(style = "margin-top: 10px;"),
                           DTOutput("clinicalDataTable")
                  ),
                  
                  tabPanel("Gene Expression Preview",
                           div(style = "margin-top: 10px;"),
                           p("Sample of gene expression data (log2 transformed FPKM values):"),
                           DTOutput("expressionDataTable")
                  )
                )
              )
      )
    )
  )
)

# Server Section
server <- function(input, output, session) {
  # CRITICAL FIX: Identify model's feature structure immediately on startup
  model_features <- reactive({
    if (exists("cv_fit")) {
      # Try to get the exact features used in the model
      model_features <- try({
        if (!is.null(cv_fit$glmnet.fit$beta)) {
          rownames(cv_fit$glmnet.fit$beta)
        } else {
          NULL
        }
      }, silent = TRUE)
      
      # If that fails, try coefficients
      if (inherits(model_features, "try-error") || is.null(model_features)) {
        model_features <- rownames(coef(cv_fit, s = lambda_min))
      }
      
      model_features
    } else {
      NULL
    }
  })
  
  # Extract top genes for display and input
  top_genes <- reactive({
    if (exists("coef_df")) {
      # Extract non-zero coefficients for genes
      nonzero_genes <- coef_df %>%
        filter(coefficient != 0) %>%
        filter(str_detect(feature, "^ENSG"))
      
      # Sort by coefficient magnitude and take top 5
      top_genes <- nonzero_genes %>%
        arrange(desc(abs(coefficient))) %>%
        head(5)
      
      # Try to get gene symbols if not already available
      if (!"gene_symbol" %in% colnames(top_genes)) {
        tryCatch({
          # Remove version suffix from Ensembl IDs
          top_genes <- top_genes %>%
            mutate(feature_clean = str_replace(feature, "\\..*$", ""))
          
          # Connect to Ensembl database
          ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
          
          # Extract unique cleaned Ensembl IDs
          ensembl_ids_clean <- unique(top_genes$feature_clean)
          
          # Query biomaRt for gene symbols
          gene_map <- getBM(
            attributes = c("ensembl_gene_id", "hgnc_symbol"),
            filters = "ensembl_gene_id",
            values = ensembl_ids_clean,
            mart = ensembl
          ) %>%
            rename(ensembl_id = ensembl_gene_id, gene_symbol = hgnc_symbol)
          
          # Join gene symbols back
          top_genes <- top_genes %>%
            left_join(gene_map, by = c("feature_clean" = "ensembl_id"))
          
          # Add label
          top_genes <- top_genes %>%
            mutate(label = ifelse(is.na(gene_symbol) | gene_symbol == "", 
                                  feature, 
                                  paste0(gene_symbol, " (", feature, ")")))
        }, error = function(e) {
          # If biomaRt fails, just use the feature names
          top_genes$label <- top_genes$feature
        })
      } else if (!"label" %in% colnames(top_genes)) {
        # If gene_symbol exists but no label, create label
        top_genes <- top_genes %>%
          mutate(label = ifelse(is.na(gene_symbol) | gene_symbol == "", 
                                feature, 
                                paste0(gene_symbol, " (", feature, ")")))
      }
      
      return(top_genes)
    } else {
      # Return empty dataframe if model not loaded
      data.frame(
        feature = character(),
        coefficient = numeric(),
        label = character()
      )
    }
  })
  
  # Download template file
  output$download_template <- downloadHandler(
    filename = function() {
      "gene_expression_template.csv"
    },
    content = function(file) {
      if (exists("selected_genes") && length(selected_genes) > 0) {
        # Create template with selected genes
        gene_template <- data.frame(matrix(ncol = length(selected_genes) + 1, nrow = 1))
        colnames(gene_template) <- c("sample_id", selected_genes)
        gene_template[1, ] <- c("PATIENT_1", rep(2.5, length(selected_genes)))
        
        write.csv(gene_template, file, row.names = FALSE)
      } else {
        # Default template if model isn't loaded
        genes <- top_genes()
        if (nrow(genes) > 0) {
          gene_template <- data.frame(matrix(ncol = nrow(genes) + 1, nrow = 1))
          colnames(gene_template) <- c("sample_id", genes$feature)
          gene_template[1, ] <- c("PATIENT_1", rep(2.5, nrow(genes)))
        } else {
          gene_template <- data.frame(sample_id = "PATIENT_1")
        }
        write.csv(gene_template, file, row.names = FALSE)
      }
    }
  )
  
  # Function to load data from various file formats
  loadExpressionData <- function(file_path) {
    tryCatch({
      file_ext <- tolower(tools::file_ext(file_path))
      
      if (file_ext %in% c("rdata", "rda")) {
        # Load RData/RDA file
        temp_env <- new.env()
        load(file_path, envir = temp_env)
        
        # Get loaded objects
        obj_names <- ls(temp_env)
        
        gene_data <- NULL
        
        # First priority: look for data frames or matrices
        for (obj in obj_names) {
          obj_data <- get(obj, envir = temp_env)
          if (is.data.frame(obj_data) || is.matrix(obj_data)) {
            if (exists("selected_genes")) {
              matching_genes <- intersect(colnames(obj_data), selected_genes)
              if (length(matching_genes) > 0) {
                gene_data <- obj_data
                break
              }
            } else {
              gene_data <- obj_data
              break
            }
          }
        }
        
        return(gene_data)
        
      } else if (file_ext == "csv") {
        # Handle CSVs only
        return(read.csv(file_path, stringsAsFactors = FALSE))
      } else {
        stop("Unsupported file type. Please upload a .csv, .rda, or .RData file.")
      }
    }, error = function(e) {
      warning(paste("Error loading expression data from file:", e$message))
      return(NULL)
    })
  }
  
  # Prediction logic
  observeEvent(input$predict, {
    # Show progress indicator
    withProgress(message = 'Generating prediction...', value = 0, {
      
      # CRITICAL FIX: Use the exact structure needed for the Cox model
      # Create a complete patient data frame
      # Stage 1: Initialise with zeros
      features <- model_features()
      
      if (is.null(features) || length(features) == 0) {
        showNotification("Error: Could not determine model features", type = "error")
        return(NULL)
      }
      
      # Create patient data with correct structure
      patient_data <- data.frame(matrix(0, nrow = 1, ncol = length(features)))
      colnames(patient_data) <- features
      rownames(patient_data) <- "new_patient"
      
      incProgress(0.2, detail = "Processing demographics...")
      
      # Step 2: Fill in demographics
      # Age
      if ("age_diagnosis" %in% features) {
        patient_data$age_diagnosis <- input$age
      }
      
      # Gender (one-hot encoded)
      if ("gendermale" %in% features) {
        patient_data$gendermale <- ifelse(tolower(input$gender) == "male", 1, 0)
      }
      if ("genderfemale" %in% features) {
        patient_data$genderfemale <- ifelse(tolower(input$gender) == "female", 1, 0)
      }
      
      # Cancer stage (one-hot encoded)
      for (feature in features) {
        if (grepl("^ajcc_pathologic_stage", feature)) {
          stage_name <- sub("^ajcc_pathologic_stage", "", feature)
          patient_data[, feature] <- ifelse(input$stage == stage_name, 1, 0)
        }
      }
      
      incProgress(0.4, detail = "Processing gene expression data...")
      
      # Step 3: Fill in gene expression values from inputs
      genes <- top_genes()
      for (i in 1:nrow(genes)) {
        gene_name <- genes$feature[i]
        input_name <- paste0("gene_", gene_name)
        
        if (gene_name %in% features && !is.null(input[[input_name]])) {
          patient_data[, gene_name] <- input[[input_name]]
        }
      }
      
      # Step 4: Add gene expression values from uploaded file
      if (!is.null(input$file)) {
        gene_data <- loadExpressionData(input$file$datapath)
        
        if (!is.null(gene_data) && nrow(gene_data) > 0) {
          # Add genes from the uploaded file to patient data
          for (gene in features) {
            # Skip non-gene features and those already handled
            if (!grepl("^ENSG", gene) || gene %in% genes$feature) {
              next
            }
            
            if (gene %in% colnames(gene_data)) {
              patient_data[, gene] <- gene_data[[gene]][1]
            }
          }
        }
      }
      
      incProgress(0.6, detail = "Calculating risk score...")
      
      # Convert to matrix for prediction
      patient_matrix <- as.matrix(patient_data)
      
      # Predict risk score
      tryCatch({
        risk_score <- as.numeric(predict(cv_fit, s = lambda_min, newx = patient_matrix, type = "link"))
        
        # Make sure it's a valid number
        if (is.na(risk_score) || !is.finite(risk_score)) {
          risk_score <- 0  # Default to neutral risk
          showNotification("Warning: Could not calculate a valid risk score. Using default.", type = "warning")
        }
        
        incProgress(0.8, detail = "Generating survival curve...")
        
        # Use the baseline survival to calculate survival curve
        baseline_surv <- baseline_surv  # From loaded data
        surv_times <- baseline_surv$time
        surv_probs <- as.vector(baseline_surv$surv) ^ exp(risk_score)
        
        # Create risk groups based on quantiles from training data
        risk_cutoffs <- quantile(cox_data$linear_risk, probs = c(0, 0.33, 0.66, 1), na.rm = TRUE)
        
        if (risk_score <= risk_cutoffs[2]) {
          risk_group <- "Low"
          risk_color <- "risk-low"
        } else if (risk_score <= risk_cutoffs[3]) {
          risk_group <- "Medium"
          risk_color <- "risk-medium"
        } else {
          risk_group <- "High"
          risk_color <- "risk-high"
        }
        
        # Create survival curve data frame for plotting
        surv_df <- data.frame(
          Time = surv_times,
          Probability = surv_probs
        )
        
        incProgress(1.0, detail = "Rendering results...")
        
        # Display results
        output$survivalPlot <- renderPlotly({
          p <- ggplot(surv_df, aes(x = Time, 
                              y = Probability,
                              text = paste0("Time: ", round(Time, 2), "<br>",
                                            "Probability: ", round(Probability, 3)))) +
            geom_line(aes(group = 1), size = 1.2, color = "#fdce87") +
            geom_ribbon(aes(group = 1, ymin = 0, ymax = Probability), alpha = 0.2, fill = "#fdce87") +
            geom_hline(yintercept = 0.5, linetype = "dashed", color = "#F44336") +
            labs(
              subtitle = paste("Patient age:", input$age, "| Gender:", input$gender, "| Stage:", input$stage),
              x = "Time (Years)",
              y = "Survival Probability"
            ) +
            scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
            scale_x_continuous(breaks = seq(0, max(surv_df$Time), by = 2)) +
            theme_minimal() +
            theme(
              plot.title = element_text(face = "bold", size = 16),
              plot.subtitle = element_text(size = 12),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10),
              panel.grid.minor = element_blank()
            )
          
          ggplotly(p, tooltip = "text")
        })
        
        output$riskScore <- renderText({
          paste(round(risk_score, 3))
        })
        
        output$riskGroup <- renderUI({
          HTML(paste0("<span class='", risk_color, "'>", risk_group, "</span>"))
        })
        
        output$survivalProb <- renderText({
          # Calculate 5-year survival probability
          if (max(surv_df$Time) >= 5) {
            five_year_prob <- approx(surv_df$Time, surv_df$Probability, xout = 5)$y
          } else {
            # If we don't have data out to 5 years, use the last point
            five_year_prob <- tail(surv_df$Probability, 1)
          }
          
          # Ensure value is between 0 and 100%
          five_year_prob <- pmax(0.01, pmin(0.99, five_year_prob))
          
          paste(round(five_year_prob * 100, 1), "%")
        })
        
        # NEW: Risk Score Box Plot
        output$riskBoxPlot <- renderPlotly({
          p <- ggplot() +
            geom_boxplot(data = cox_data, aes(x = factor(0), y = linear_risk,
                                              text = paste0("Linear Risk: ", round(linear_risk, 3))),
                         fill = "#ffbd59", alpha = 0.7, width = 0.5) +
            geom_jitter(data = cox_data, aes(x = factor(0), y = linear_risk,
                                             text = paste0("Linear Risk: ", round(linear_risk, 3))),
                        alpha = 0.2, width = 0.2, color = "#fdce87") +
            geom_point(aes(x = factor(0), y = risk_score,
                           text = paste0("Current Patient\nRisk Score: ", round(risk_score, 3))),
                       color = "red", size = 5) +
            geom_hline(yintercept = risk_cutoffs[2], linetype = "dashed", color = "#4CAF50") +
            geom_hline(yintercept = risk_cutoffs[3], linetype = "dashed", color = "#F44336") +
            geom_text(aes(x = factor(0.15), y = risk_cutoffs[2] + 0.15, 
                          label = "Low Risk Threshold"), hjust = 0, vjust = -0.5, 
                      color = "#4CAF50", size = 3.5) +
            geom_text(aes(x = factor(0.15), y = risk_cutoffs[3] + 0.15, 
                          label = "High Risk Threshold"), hjust = 0, vjust = 1.5, 
                      color = "#F44336", size = 3.5) +
            labs(y = "Risk Score") +
            theme_minimal() +
            theme(
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.minor = element_blank()
            )
          
           ggplotly(p, tooltip = "text")
        })
        
        # Generate lifestyle recommendations based on risk group
        output$recommendations <- renderUI({
          if (risk_group == "Low") {
            HTML('<div class="panel panel-success">
                  <div class="panel-heading">
                    <h3 class="panel-title"><i class="fa fa-check-circle"></i> Low Risk Recommendations</h3>
                  </div>
                  <div class="panel-body">
                    <ul>
                      <li><strong>Monitoring:</strong> Continue with regular skin examinations every 12 months.</li>
                      <li><strong>Sun Protection:</strong> Maintain sun protection practices including wearing SPF 30+ sunscreen daily.</li>
                      <li><strong>Self-Checks:</strong> Monitor any changes in existing moles using the ABCDE method.</li>
                      <li><strong>Nutrition:</strong> Maintain a healthy diet rich in antioxidants and omega-3 fatty acids.</li>
                      <li><strong>Exercise:</strong> Regular physical activity for at least 150 minutes per week.</li>
                    </ul>
                  </div>
                </div>')
          } else if (risk_group == "Medium") {
            HTML('<div class="panel panel-warning">
                  <div class="panel-heading">
                    <h3 class="panel-title"><i class="fa fa-exclamation-triangle"></i> Medium Risk Recommendations</h3>
                  </div>
                  <div class="panel-body">
                    <ul>
                      <li><strong>Monitoring:</strong> Increase frequency of skin examinations to every 6 months.</li>
                      <li><strong>Imaging:</strong> Consider more comprehensive imaging such as dermoscopy during checkups.</li>
                      <li><strong>Sun Protection:</strong> Strengthen sun protection with SPF 50+ sunscreen, reapplied every 2 hours when outdoors.</li>
                      <li><strong>Timing:</strong> Avoid sun exposure during peak hours (10am-4pm).</li>
                      <li><strong>Nutrition:</strong> Consider increasing intake of foods rich in selenium and vitamins D and E.</li>
                      <li><strong>Stress Management:</strong> Implement stress reduction techniques such as mindfulness or yoga.</li>
                    </ul>
                  </div>
                </div>')
          } else if (risk_group == "High") {
            HTML('<div class="panel panel-danger">
                  <div class="panel-heading">
                    <h3 class="panel-title"><i class="fa fa-exclamation-circle"></i> High Risk Recommendations</h3>
                  </div>
                  <div class="panel-body">
                    <ul>
                      <li><strong>Consultation:</strong> Schedule immediate consultation with oncologist to discuss treatment options.</li>
                      <li><strong>Monitoring:</strong> Arrange for frequent monitoring every 3 months with comprehensive skin examinations.</li>
                      <li><strong>Genetic Counseling:</strong> Consider genetic counseling to assess hereditary risk factors.</li>
                      <li><strong>Sun Protection:</strong> Ensure strict sun protection protocols with physical barriers and high SPF sunscreen.</li>
                      <li><strong>Treatment Options:</strong> Discuss specific immunotherapy and targeted therapy options based on molecular profile.</li>
                      <li><strong>Nutrition Consultation:</strong> Consider nutritional consultation to support immune function and overall health.</li>
                      <li><strong>Support Network:</strong> Engage with cancer support groups and mental health resources.</li>
                    </ul>
                  </div>
                </div>')
          } else {
            HTML('<div class="panel panel-info">
                  <div class="panel-heading">
                    <h3 class="panel-title"><i class="fa fa-info-circle"></i> General Recommendations</h3>
                  </div>
                  <div class="panel-body">
                    <p>No specific risk group could be determined. Here are general recommendations for melanoma patients:</p>
                    <ul>
                      <li>Consult with your healthcare provider for personalised advice.</li>
                      <li>Practice regular sun protection and skin self-examinations.</li>
                      <li>Maintain a healthy lifestyle with balanced nutrition and regular exercise.</li>
                    </ul>
                  </div>
                </div>')
          }
        })
      }, error = function(e) {
        # Handle errors in prediction
        showNotification(paste("Error in prediction:", e$message), type = "error")
        
        # Output debug info
        output$debugOutput <- renderPrint({
          cat("Error in prediction:", e$message, "\n\n")
          cat("Patient matrix dimensions:", dim(patient_matrix), "\n")
          cat("First few feature names:", head(colnames(patient_matrix)), "\n\n")
          cat("Model features count:", length(features), "\n")
          cat("First few model features:", head(features), "\n")
          cat("Missing features:", setdiff(features, colnames(patient_matrix)), "\n")
          cat("Extra features:", setdiff(colnames(patient_matrix), features), "\n")
        })
        
        # Set default values so the UI doesn't break
        output$riskScore <- renderText({ "Error" })
        output$riskGroup <- renderUI({ HTML("<span>Unable to determine</span>") })
        output$survivalProb <- renderText({ "Unable to calculate" })
        
        # Create empty plot
        output$survivalPlot <- renderPlotly({
          ggplot() + 
            annotate("text", x = 0.5, y = 0.5, label = paste("Error: Unable to generate survival curve\n", e$message)) +
            theme_void() +
            xlim(0, 1) + ylim(0, 1)
        })
        
        # Generic recommendations
        output$recommendations <- renderUI({
          HTML('<div class="panel panel-info">
                <div class="panel-heading">
                  <h3 class="panel-title"><i class="fa fa-info-circle"></i> Error in Prediction</h3>
                </div>
                <div class="panel-body">
                  <p>We could not generate specific recommendations due to an error: ', e$message, '</p>
                  <p>Here are general recommendations for melanoma patients:</p>
                  <ul>
                    <li>Consult with your healthcare provider for personalised advice.</li>
                    <li>Practice regular sun protection and skin self-examinations.</li>
                    <li>Maintain a healthy lifestyle with balanced nutrition and regular exercise.</li>
                  </ul>
                </div>
              </div>')
        })
      })
    })
  })
  
  # Model tab visualisations
  observe({
    # 1. Gene Importance Plot
    output$geneImportancePlot <- renderPlotly({
      if (!exists("coef_df") || is.null(coef_df)) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, label = "No coefficient data available") +
                 theme_void() +
                 xlim(0, 1) + ylim(0, 1))
      }
      
      # Filter gene features (Ensembl IDs)
      nonzero_genes <- coef_df %>%
        filter(coefficient != 0) %>%
        filter(str_detect(feature, "^ENSG"))
      
      if (nrow(nonzero_genes) == 0) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, label = "No genes with non-zero coefficients found") +
                 theme_void() +
                 xlim(0, 1) + ylim(0, 1))
      }
      
      # Get gene symbols if we can
      if (!"plot_label" %in% colnames(nonzero_genes)) {
        tryCatch({
          # Remove version suffix from Ensembl IDs
          nonzero_genes <- nonzero_genes %>%
            mutate(feature_clean = str_replace(feature, "\\..*$", ""))
          
          # Connect to Ensembl database
          ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
          
          # Extract unique cleaned Ensembl IDs
          ensembl_ids_clean <- unique(nonzero_genes$feature_clean)
          
          # Query biomaRt for gene symbols
          gene_map <- getBM(
            attributes = c("ensembl_gene_id", "hgnc_symbol"),
            filters = "ensembl_gene_id",
            values = ensembl_ids_clean,
            mart = ensembl
          ) %>%
            rename(ensembl_id = ensembl_gene_id, gene_symbol = hgnc_symbol)
          
          # Join gene symbols back
          nonzero_genes <- nonzero_genes %>%
            left_join(gene_map, by = c("feature_clean" = "ensembl_id"))
          
          # Create plot label
          nonzero_genes <- nonzero_genes %>%
            mutate(plot_label = ifelse(is.na(gene_symbol) | gene_symbol == "", 
                                       feature, 
                                       gene_symbol))
        }, error = function(e) {
          # If biomaRt fails, just use the feature names
          nonzero_genes$plot_label <- nonzero_genes$feature
        })
      }
      
      # Plot top genes by absolute coefficient
      top_genes <- nonzero_genes %>%
        arrange(desc(abs(coefficient))) %>%
        head(20)
      
      p <- ggplot(top_genes, aes(x = reorder(plot_label, coefficient), 
                            y = coefficient, 
                            fill = coefficient > 0,
                            text = paste0("Coefficient: ", round(coefficient, 3))
                            )) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = c("#2196F3", "#F44336"),
                          labels = c("Protective", "Risk"),
                          name = "Effect") +
        coord_flip() +
        labs(subtitle = "Based on Cox model coefficients",
             x = "Gene",
             y = "Coefficient Value (log hazard ratio)") +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 14),
          legend.position = "top",
          panel.grid.minor = element_blank()
        )
      
      ggplotly(p, tooltip = "text")
    })
    
    # 2. Clinical Feature Importance Plot
    output$clinicalImportancePlot <- renderPlotly({
      if (!exists("coef_df") || is.null(coef_df)) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, label = "No coefficient data available") +
                 theme_void() +
                 xlim(0, 1) + ylim(0, 1))
      }
      
      # Filter clinical features (not ENSG IDs)
      all_clinical <- coef_df %>%
        filter(coefficient != 0) %>%
        filter(!str_detect(feature, "^ENSG"))
      
      if (nrow(all_clinical) == 0) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, label = "No clinical features with non-zero coefficients found") +
                 theme_void() +
                 xlim(0, 1) + ylim(0, 1))
      }
      
      # Clean up feature names
      all_clinical <- all_clinical %>%
        mutate(feature = case_when(
          str_starts(feature, "ajcc_pathologic_stage") ~ str_remove(feature, "ajcc_pathologic_stage"),
          str_starts(feature, "gender") ~ str_remove(feature, "gender"),
          feature == "age_diagnosis" ~ "Age at Diagnosis",
          TRUE ~ feature
        ),
        feature = str_trim(feature))
      
      # Plot clinical features
      p <- ggplot(all_clinical, aes(
            x = reorder(feature, coefficient), 
            y = coefficient, 
            fill = coefficient >= 0,
            text = paste0("Coefficient: ", round(coefficient, 3))
          )) +
        geom_bar(stat = "identity") + 
        scale_fill_manual(values = c("#2196F3", "#F44336"),
                          labels = c("Protective", "Risk"),
                          name = "Effect") +
        coord_flip() +
        labs(subtitle = "Based on Cox model coefficients",
             x = "Feature",
             y = "Coefficient Value (log hazard ratio)") +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 14),
          legend.position = "top",
          panel.grid.minor = element_blank()
        )
      
      ggplotly(p, tooltip = "text")
    })
    
    # 3. Risk Distribution Plot
    output$riskDistributionPlot <- renderPlotly({
      if (!exists("cox_data") || is.null(cox_data) || !"linear_risk" %in% colnames(cox_data)) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, label = "No risk score data available") +
                 theme_void() +
                 xlim(0, 1) + ylim(0, 1))
      }
      
      p <- ggplot(cox_data, aes(x = linear_risk)) +
        geom_histogram(
          bins = 30, fill = "#ffbd59", color = "white",
          aes(text = paste0("Count: ", ..count.., "<br>",
                            "Linear Risk: ", round(..x.., 2)))
        ) +
        labs(subtitle = "From LASSO-regularised Cox model",
             x = "Linear Risk Score",
             y = "Count") +
        theme_minimal() +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          axis.title = element_text(size = 12),
          panel.grid.minor = element_blank()
        )
      
      ggplotly(p, tooltip = "text")
    })
    
    # 4. Metrics Table
    output$metricsTable <- renderTable({
      if (!exists("c_index") || is.null(c_index)) {
        return(data.frame(
          Metric = "Concordance Index",
          Value = "Not available"
        ))
      }
      
      data.frame(
        Metric = "Concordance Index",
        Value = round(c_index$c.index, 3)
      )
    }, align = 'lr', striped = TRUE, hover = TRUE, bordered = TRUE)
    
    output$pcaPlot <- renderPlotly({
      if (!exists("query2")) {
        return(ggplot() + 
                 annotate("text", x = 0.5, y = 0.5, label = "Gene expression data not available") +
                 theme_void() + xlim(0, 1) + ylim(0, 1))
      }
      
      tryCatch({
        fpkm <- assay(query2, "fpkm_uq_unstrand")
        log_expr <- log2(fpkm + 1)
        
        gene_var <- apply(log_expr, 1, var)
        top_var_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(1000, length(gene_var))]
        expr_for_pca <- log_expr[top_var_genes, ]
        
        pca_result <- prcomp(t(expr_for_pca), scale = TRUE, center = TRUE)
        
        clinical_pca <- colData(query2) %>% 
          as.data.frame() %>%
          dplyr::select(gender, ajcc_pathologic_stage)
        
        pca_data <- data.frame(
          Sample = colnames(expr_for_pca),
          Gender = clinical_pca$gender,
          Stage = ifelse(is.na(clinical_pca$ajcc_pathologic_stage), "Unknown", clinical_pca$ajcc_pathologic_stage),
          pca_result$x
        )
        
        var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
        
        gender_levels <- unique(pca_data$Gender)
        gender_colors <- "#999999"  # neutral color, symbol shows legend
        gender_shapes <- c("male" = "circle", "female" = "square", "unknown" = "triangle-up")
        
        # Data for legend traces - Stage
        stage_levels <- c("Stage 0", "Stage I", "Stage IA", "Stage IB",
                          "Stage II", "Stage IIA", "Stage IIB", "Stage IIC",
                          "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                          "Stage IV", "Unknown")
        stage_colors <- c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4",
                          "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", "#008080",
                          "#e6beff", "#9a6324", "#800000")
        names(stage_colors) <- stage_levels
        
        # Determine PC axes based on selected tab
        pc_selection <- input$pca_tab
        pc_x <- "PC1"; pc_y <- "PC2"  # Default
        if (pc_selection == "pc2_3") {
          pc_x <- "PC2"; pc_y <- "PC3"
        } else if (pc_selection == "pc3_4") {
          pc_x <- "PC3"; pc_y <- "PC4"
        }
        
        # Start plotly plot
        p <- plot_ly()
        
        # Add "Gender" group header dummy trace (no markers shown)
        p <- add_trace(p,
                       x = NA, y = NA,
                       type = 'scatter',
                       mode = 'markers',
                       marker = list(opacity = 0),
                       name = "Gender",
                       showlegend = TRUE,
                       hoverinfo = 'none')
        
        # Add one trace per gender level (with symbol and same color)
        for (g in gender_levels) {
          display_name <- tools::toTitleCase(g)
          p <- add_trace(p,
                         data = filter(pca_data, Gender == g),
                         name = display_name,
                         x = ~get(pc_x),
                         y = ~get(pc_y),
                         type = 'scatter',
                         mode = 'markers',
                         marker = list(symbol = gender_shapes[g], color = gender_colors, size = 8, opacity = 0.7),
                         name = g,
                         showlegend = TRUE,
                         text = ~paste("Sample:", Sample,
                                       "<br>Stage:", Stage,
                                       "<br>Gender:", Gender),
                         hoverinfo = 'text')
        }
        
        # Add "Stage" group header dummy trace (no markers shown)
        p <- add_trace(p,
                       x = NA, y = NA,
                       type = 'scatter',
                       mode = 'markers',
                       marker = list(opacity = 0),
                       name = "Stage",
                       showlegend = TRUE,
                       hoverinfo = 'none')
        
        # Add one trace per stage level (with color, fixed symbol)
        for (s in stage_levels) {
          if (s %in% pca_data$Stage) {
            p <- add_trace(p,
                           data = filter(pca_data, Stage == s),
                           x = ~get(pc_x),
                           y = ~get(pc_y),
                           type = 'scatter',
                           mode = 'markers',
                           marker = list(color = stage_colors[s], symbol = "circle", size = 8, opacity = 0.7),
                           name = s,
                           showlegend = TRUE,
                           text = ~paste("Sample:", Sample,
                                         "<br>Stage:", Stage,
                                         "<br>Gender:", Gender),
                           hoverinfo = 'text')
          }
        }
        
        # Final layout
        p <- layout(p,
                    title = list(text = "PCA of Gene Expression<br><sub>Top 1000 most variable genes</sub>"),
                    xaxis = list(title = paste0(pc_x, " (", var_explained[as.numeric(gsub("PC", "", pc_x))], "% variance)")),
                    yaxis = list(title = paste0(pc_y, " (", var_explained[as.numeric(gsub("PC", "", pc_y))], "% variance)")),
                    legend = list(title = list(text = "<b>Legend</b>")),
                    margin = list(t = 60))
        
        
      }, error = function(e) {
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = paste("Error calculating PCA:", e$message)) +
          theme_void() +
          xlim(0, 1) + ylim(0, 1)
      })
    })
    
  # Debug output
  output$debugOutput <- renderPrint({
    cat("Model Features:\n")
    print(head(model_features(), 20))
    
    cat("\nTop Genes:\n")
    print(top_genes())
    
    if (exists("cv_fit")) {
      cat("\nModel Information:\n")
      cat("Lambda min:", cv_fit$lambda.min, "\n")
      cat("Has lambda:", !is.null(cv_fit$lambda), "\n")
      cat("Has beta matrix:", !is.null(cv_fit$glmnet.fit$beta), "\n")
      
      if (!is.null(cv_fit$glmnet.fit$beta)) {
        cat("Beta matrix dimensions:", dim(cv_fit$glmnet.fit$beta), "\n")
      }
    } else {
      cat("\nModel not loaded properly.")
    }
  })
  
  # Dataset Explorer outputs
  output$genderDist <- renderTable({
    data.frame(
      Gender = c("Female", "Male"),
      Count = as.integer(c(180, 293))
    )
  })
  
  output$stageDist <- renderTable({
    data.frame(
      Stage = c("Stage 0", "Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB", "Stage IIC", 
                "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV"),
      Count = as.integer(c(1,4,3,4,6,2,10,30,14,2,20,17,10))
    )
  })
  
  output$ageStats <- renderTable({
    values <- c(20.35, 50.25, 61.38, 60.94, 73.29, 90.06, 132)
    
    # Format all values with 2 decimals except the last one as integer
    formatted_values <- c(
      sprintf("%.2f", values[-length(values)]),
      as.character(as.integer(values[length(values)]))
    )
    
    age_stats <- data.frame(
      Statistic = c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "NA's"),
      Value = formatted_values,
      stringsAsFactors = FALSE
    )
  })
  
  
  
  output$clinicalDataTable <- renderDT({
    if (!exists("query2")) {
      return(NULL)
    }
    
    clinical <- colData(query2) %>% 
      as.data.frame() %>%
      dplyr::select(gender, days_to_death, ajcc_pathologic_stage, age_at_diagnosis) %>%
      mutate(
        `Age (years)` = age_at_diagnosis / 365,
        Gender = gender,
        `Cancer Stage` = ajcc_pathologic_stage,
        `Days to Death` = days_to_death
      ) %>%
      dplyr::select(`Gender`, `Age (years)`, `Cancer Stage`, `Days to Death`)
    
    datatable(clinical, 
              options = list(pageLength = 10, 
                             scrollX = TRUE,
                             autoWidth = TRUE),
              filter = 'top',
              rownames = FALSE)
  })
  
  output$expressionDataTable <- renderDT({
    if (!exists("query2")) {
      return(NULL)
    }
    
    # Get expression data
    fpkm <- assay(query2, "fpkm_uq_unstrand")
    log_expr <- log2(fpkm + 1)
    
    # Get a sample of gene expression data
    sample_expr <- log_expr[1:min(20, nrow(log_expr)), 1:min(10, ncol(log_expr))]
    sample_expr_df <- as.data.frame(sample_expr)
    
    # Try to get gene symbols if possible
    gene_ids <- rownames(sample_expr_df)
    gene_ids_clean <- str_replace(gene_ids, "\\..*$", "")
    
    tryCatch({
      ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
      gene_map <- getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters = "ensembl_gene_id",
        values = gene_ids_clean,
        mart = ensembl
      )
      
      # Create a column for gene symbols
      gene_df <- data.frame(
        Ensembl_ID = gene_ids,
        Gene_Symbol = NA
      )
      
      for(i in 1:nrow(gene_df)) {
        clean_id <- str_replace(gene_df$Ensembl_ID[i], "\\..*$", "")
        matching_row <- which(gene_map$ensembl_gene_id == clean_id)
        
        if(length(matching_row) > 0 && !is.na(gene_map$hgnc_symbol[matching_row[1]]) && 
           gene_map$hgnc_symbol[matching_row[1]] != "") {
          gene_df$Gene_Symbol[i] <- gene_map$hgnc_symbol[matching_row[1]]
        } else {
          gene_df$Gene_Symbol[i] <- gene_df$Ensembl_ID[i]
        }
      }
      
      rownames(sample_expr_df) <- gene_df$Gene_Symbol
    }, error = function(e) {
      # Keep original rownames if biomaRt fails
      warning("Could not retrieve gene symbols:", e$message)
    })
    
    datatable(sample_expr_df, 
              options = list(pageLength = 10, 
                             scrollX = TRUE,
                             autoWidth = TRUE),
              caption = "Sample of log2(FPKM+1) gene expression values",
              rownames = TRUE)
    })
  })
}

# Run the application
shinyApp(ui = ui, server = server)
