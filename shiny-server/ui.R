library(shiny)
library(ggplot2)
# library(Cairo)   # For nicer ggplot2 output when deployed on Linux

shinyUI <- fluidPage(
  # Some custom CSS
  tags$head(
    tags$style(HTML("
      "))
  ),
  fluidRow(
    uiOutput("error_message"),
    div(
      uiOutput("ui_reference_samples"),
      uiOutput("ui_sample"),
      hr()
      # style = "font-size:60%"
    )
  ),
  fluidRow(
    div(
      column(3,
             textInput("gene_id", "", value = "TP53"),
             actionButton("gene_id_find", "",  icon = icon("search"))
      ),
      column(9, 
             tags$h5("Targeted Raw Counts Summary"),
               verbatimTextOutput("sample_summary"),             
               hr()
        )
    )
  ),
  fluidRow(
    div(
      column(2,
              sliderInput('cn_max_val', 'Copy Number Scale Max', 
                          min=4, max=500,
                          value=4, 
                          step=1, round=0)
      ),
      column(2,
             uiOutput("min_probe_coverage")
#              sliderInput('min_probe_coverage', 'Minimum Targeted Read Depth', 
#                          min=0, max=500,
#                          value=0, 
#                          step=10, round=0)
      ),
      column(2,
             uiOutput("display_wg_data"),
             uiOutput("display_targeted_data"),
             uiOutput("display_aux_target_data")
#              checkboxInput(inputId = "display_wg_data", label = "WG",value = FALSE),
#              checkboxInput(inputId = "display_targeted_data", label = "Targeted",value = TRUE),
#              checkboxInput(inputId = "display_aux_target_data", label = "Off Target",value = FALSE)
      ),
      column(2,
             uiOutput("use_reference_wg"),
             uiOutput("use_reference_targeted"),
             uiOutput("use_reference_aux_target")
#              checkboxInput(inputId = "use_reference_wg", label = "Use reference for WG",value = TRUE),
#              checkboxInput(inputId = "use_reference_targeted", label = "Use reference for Targeted",value = TRUE),
#              checkboxInput(inputId = "use_reference_aux_target", label = "Use reference for Off Target",value = TRUE)
      ),
      column(2,
             uiOutput("gc_correct_wg"),
             uiOutput("gc_correct_targeted"),
             uiOutput("gc_correct_aux_target"),
             uiOutput("display_error_bars"),
             uiOutput("display_baf"),
             uiOutput("display_fusion"),
             uiOutput("display_segmentation")
             #              checkboxInput(inputId = "gc_correct_wg", label = "GC correct WG",value = FALSE),
#              checkboxInput(inputId = "gc_correct_targeted", label = "GC correct Targeted",value = TRUE),
#              checkboxInput(inputId = "gc_correct_aux_target", label = "GC correct Off Target",value = FALSE)
      ),
      column(2,
#             actionButton(inputId = "recalculate_reference", label = "Recalculate Reference", icon = icon("cogs")),
             checkboxInput(inputId = "display_all_gene_names", label = "Show all annotations",value = FALSE),
             checkboxInput(inputId = "display_cytobands", label = "Show Cytobands",value = TRUE),
             uiOutput("recalculate_reference")
      )
    )
  ),
  fluidRow(
    plotOutput("plot1", 
      height = 256,
      click = "plot1_click",
      dblclick = "plot1_dblclick",
      hover = hoverOpts(
        id = "plot1_hover",
      ),
      brush = brushOpts(
        id = "plot1_brush",
        resetOnNew = TRUE
      )
    )
  ),
  fluidRow (
    # h4("Left plot controls right plot"),
    fluidRow(
      plotOutput("plot2", 
        height = 256,
        click = "plot2_click",
        dblclick = "plot2_dblclick",
        hover = hoverOpts(
          id = "plot2_hover",
        ),
        brush = brushOpts(
          id = "plot2_brush",
          resetOnNew = TRUE
        )
      )
    )
  ),
  fluidRow(
    plotOutput("plot3", 
      height = 256,
      # click = "plot3_click",
      dblclick = "plot3_dblclick",
      brush = brushOpts(
        id = "plot3_brush",
        resetOnNew = TRUE
      )
    )
  ),
  fluidRow(
    div(
      column(3),
      column(6, 
             tags$h5("Summary of N in selected region."),
             verbatimTextOutput("hilighted_data_summary"),             
             hr()
      ),
      column(3)
      # column(2, verbatimTextOutput("hilighted_data_summary"))
      # style = "font-size:60%"
    )
  ),
  fluidRow(
    div(
      column(2),
      column(8,dataTableOutput('hilighted_data')),
      column(2)
      # column(2, verbatimTextOutput("hilighted_data_summary"))
      # style = "font-size:60%"
    )
  )
)
