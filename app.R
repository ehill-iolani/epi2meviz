library(shiny)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(vegan)
library(shinycssloaders)

# Accessory Function(s) -------------------------------------------------------
bc_sample <- function(x, y) {
  x <- x[sample(nrow(x), y, replace = FALSE), ]
  x <- x[!duplicated(x$species), ]
  x <- length(unique(x$species))
  return(x)
}

# Barcode index 1 - 24 --------------------------------------------------------
bci <- c("barcode01", "barcode02", "barcode03", "barcode04", "barcode05",
         "barcode06", "barcode07", "barcode08", "barcode09", "barcode10",
         "barcode11", "barcode12", "barcode13", "barcode14", "barcode15",
         "barcode16", "barcode17", "barcode18", "barcode19", "barcode20",
         "barcode21", "barcode22", "barcode23", "barcode24")

# Define the UI ---------------------------------------------------------------
ui <- shinyUI(fluidPage(
  titlePanel("Post EPI2ME Analysis"),
  tabsetPanel(
    tabPanel("Upload EPI2ME .csv",
             titlePanel("Uploading Files and Setting Parameters"),
             sidebarLayout(
               sidebarPanel(
                 fileInput("file1", "Choose EPI2ME .csv File",
                           accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv"),
                           multiple = TRUE),
                 actionButton("button", "Submit"),
                 radioButtons("sep", "Separator",
                              c(Comma = ","),
                              ","),
                 textInput("filt", "Average EPI2ME Accuracy",
                           value = "80"),
                 checkboxGroupInput("barcodes", "Barcodes to Analyze",
                           choices = bci, selected = "all"),
               ),
               mainPanel(
                 withSpinner(tableOutput("contents"))
               )
             )
    ),
    tabPanel("Rarefaction Curve",
             fluidPage(
               headerPanel("Rarefaction Curves"),
               downloadButton("rare_butt", "Download PDF"),
               actionButton("rare_help", "Help"),
               mainPanel(
                 withSpinner(plotOutput("rarefaction"))
               )
             )
    ),
    tabPanel("Relative Abundance Plot",
             fluidPage(
               headerPanel("Relative Abundance"),
               downloadButton("relab_butt", "Download PDF"),
               actionButton("relab_help", "Help"),
               mainPanel(
                 withSpinner(plotOutput("relative_abundance"))
               )
             )
    ),
    tabPanel("Relative Abundance Table",
             fluidPage(
               headerPanel("Relative Abundance Table"),
               selectInput("relab_barcode", "Select Barcode",
                           choices = bci),
               downloadButton("relabf_butt", "Download CSV"),
               actionButton("relabf_help", "Help"),
               mainPanel(
                 withSpinner(tableOutput("relab"))
               )
             )
    ),
    tabPanel("Bray Curtis PCoA Plot",
             fluidPage(
               headerPanel("Bray Curtis PCoA"),
               downloadButton("braycurtis_butt", "Download PDF"),
               actionButton("braycurtis_help", "Help"),
               mainPanel(
                 withSpinner(plotOutput("bray_curtis"))
               )
             )
  )
)
)
)

# Define the server -----------------------------------------------------------
server <- function(input, output, session) {
  # Set large upload size limit (server side)
  options(shiny.maxRequestSize = 70 * 1024^2)

  data <- eventReactive(input$button, {
    # Filtering and cleaning the data -----------------------------------------

    # require that at least 1 input is available
    req(input$file1)

    # Stores the input file as a variable
    inFile <- input$file1

    # Parses the input file
    dat <- read.csv(inFile$datapath, header = TRUE, sep = input$sep)

    # Filters the data based on user input
    dat <- dat[dat$barcode %in% input$barcodes, ]

    # First pass filtering
    dat <- dat[dat$accuracy > input$filt, ]
    dat <- dat[dat$species != "unclassified", ]
    dat <- dat[dat$barcode != "unclassified", ]

    # Cleans up species and genus names
    dat$species <- word(dat$species, 1, 2, sep = " ")
    dat$species <- gsub("\\[|\\]", "", dat$species)
    dat <- dat[dat$species != "causative agent", ]
    dat$genus <- word(dat$species, 1, sep = " ")

    # Finds unqiue barcodes in data
    bc <- unique(dat$barcode)

    # Removes species occuring less than 5 times
    rdat <- dat %>% group_by(species) %>% tally()
    rdat <- rdat[rdat$n < 5, ]
    index <- dat$species %in% rdat$species
    rdat <- dat[!index, ]

    # Grabs dimensions of each barcode to set up rarefaction
    bcc <- data.frame()
    for (i in bc) {
      temp <- rdat[rdat$barcode == i, ]
      temp <- dim(temp)
      bcc <- rbind(bcc, temp)
    }
    rn <- mean(bcc[, 1])
    if (rn > 10000) {
      rnn <- 100
    } else {
      rnn <- 10
    }

    # Rareify all barcodes ----------------------------------------------------
    rare <- data.frame()

    for (i in bc) {
       temp <- rdat[rdat$barcode == i, ]
       reads <- seq_along(temp$species)
       urare <- data.frame()
       for (j in seq(1, nrow(temp), rnn)) {
              temp2 <- mean(replicate(2, bc_sample(temp, j)))
              temp2 <- as.data.frame(temp2)
              temp2$reads <- j
              urare <- rbind(urare, temp2)
       }
       urare$barcode <- i
       rare <- rbind(rare, urare)
    }

    # Renames columns of processed data
    names(rare) <- c("unique_species", "reads_sampled", "barcode")

  # Calculates Relative abundance ---------------------------------------------

    # Calculates relative abundance per barcode
    adat <- dat %>% group_by(genus, barcode) %>% tally()
    adat <- adat[adat$n > 5, ]
    relab <- data.frame()

    for (i in bc) {
      temp <- adat[adat$barcode == i, ]
      temp$rel_ab <- temp$n / sum(temp$n)
      relab <- rbind(relab, temp)
    }

    # Cleans table for presentation
    relabf <- data.frame(relab$genus, relab$barcode, relab$rel_ab)
    names(relabf) <- c("Genus", "Barcode", "rel_ab")
    relabf$rel_ab <- relabf$rel_ab * 100
    relabf <- relabf[with(relabf, order(Barcode, rel_ab, decreasing = TRUE)), ]
    names(relabf) <- c("Genus", "Barcode", "Relative Abundance (%)")

  # Conducts Bray Curtis PCoA -------------------------------------------------

    # Bray Curtis PCoA
    if (length(bc) < 3) {
      print("Not enough barcodes to perform Bray Curtis PCoA")
    } else {
      brayc <- dat[, c("species", "barcode")]
      brayc <- brayc %>% group_by(species, barcode) %>% tally()
      brayc <- brayc[brayc$n > 5, ]
      brayc <- brayc %>% pivot_wider(names_from = species, values_from = n)

      # Replaces NA with 0
      brayc[is.na(brayc)] <- 0
      bray_out <- vegdist(brayc[5:dim(brayc)[2]], method = "bray")

      # Performs PCoA and renames output
      bray_curtis_pcoa <- cmdscale(bray_out, k = 2, eig = TRUE, add = TRUE)
      bray_curtis_pcoa_dat <- as.data.frame(bray_curtis_pcoa$points)
      bray_curtis_pcoa_dat <- cbind(bray_curtis_pcoa_dat,
                                  brayc$barcode)
      names(bray_curtis_pcoa_dat) <- c("PCoA_1", "PCoA_2", "Barcode")
    }

  # Returns analysis outputs --------------------------------------------------

    # Returns processed data
    if (exists("bray_curtis_pcoa_dat") == FALSE) {
      combo <- list(rdat = rdat, rare = rare, relab = relab,
                    relabf = relabf)
      return(combo)
    } else {
      combo <- list(rdat = rdat, rare = rare, relab = relab,
                    relabf = relabf,
                    bray_curtis_pcoa_dat = bray_curtis_pcoa_dat)
      return(combo)
    }
  })

  # Visualizing/Downloading outputs -------------------------------------------

  # Display the table of processed data
  output$contents <- renderTable({
    head(data()$rdat)
  })

  # Display the rarefaction curve
  output$rarefaction <- renderPlot({
    ggplot(data = data()$rare, aes(x = reads_sampled,
                                     y = unique_species,
                                     color = barcode)) +
    geom_point(size = 2, aes(group = barcode)) +
    labs(x = "Reads Sampled", y = "Unique Species", fill = "Genus") +
    theme_bw()
  })

  # Download the rarefaction curve
  output$rare_butt <- downloadHandler(
    filename = function() {
      paste("rarefaction", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      print(ggplot(data = data()$rare, aes(x = reads_sampled,
                                           y = unique_species,
                                           color = barcode)) +
      geom_point(size = 2, aes(group = barcode)) +
      labs(x = "Reads Sampled", y = "Unique Species", fill = "Genus") +
      theme_bw())
      dev.off()
    }
  )

  # Display the rarefaction curve help
  observeEvent(input$rare_help, {
    showModal(modalDialog(
      title = "Rarefaction Curve Help",
      "Rarefaction curves are used to determine if enough reads have been
      sampled to capture the diversity of a sample. The x-axis represents the
      number of reads sampled and the y-axis represents the number of unique
      species found. If the curve plateaus, then enough reads have been sampled
      to capture the diversity of the sample. If the curve does not plateau,
      then more reads need to be sampled to capture the diversity of the
      sample.",
      easyClose = TRUE,
      footer = NULL
    ))
  })

  # Display the relative abundance plot
  output$relative_abundance <- renderPlot({
    ggplot(data = data()$relab, aes(x = barcode,
                                     y = rel_ab,
                                     fill = genus)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample ID", y = "Relative Abundance", fill = "Genus") +
    theme_bw()
  })

  # Display the relative abundance plot help
  observeEvent(input$relab_help, {
    showModal(modalDialog(
      title = "Relative Abundance Plot Help",
      "Relative abundance plots are used to visualize the relative abundance
      of each genus in each sample. The x-axis represents the sample ID and
      the y-axis represents the relative abundance of each genus. The legend
      represents the genus of each species.",
      easyClose = TRUE,
      footer = NULL
    ))
  })

  # Download the relative abundance plot
  output$relab_butt <- downloadHandler(
    filename = function() {
      paste("relative_abundance", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      print(ggplot(data = data()$relab, aes(x = barcode,
                                           y = rel_ab,
                                           fill = genus)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(x = "Sample ID", y = "Relative Abundance", fill = "Genus") +
      theme_bw())
      dev.off()
    }
  )

  # Display the relative abundance table help
  observeEvent(input$relabf_help, {
    showModal(modalDialog(
      title = "Relative Abundance Table Help",
      "The relative abundance table is used to visualize the relative abundance
      of each genus in each sample. The table is sorted by barcode and then
      relative abundance. The table can be downloaded as a .csv file.",
      easyClose = TRUE,
      footer = NULL
    ))
  })

  # Display the relative abundance table
  output$relab <- renderTable({
    data()$relabf[data()$relabf$Barcode == input$relab_barcode, ]
  })

  # Download the relative abundance table
  output$relabf_butt <- downloadHandler(
    filename = function() {
      paste("relative_abundance", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data()$relabf[data()$relabf$Barcode == input$relab_barcode, ],
                file, row.names = FALSE)
    }
  )

  # Display the Bray Curtis PCoA plot
    output$bray_curtis <- renderPlot({
      req(data()$bray_curtis_pcoa_dat)
      ggplot(data = data()$bray_curtis_pcoa_dat, aes(x = PCoA_1,
                                                     y = PCoA_2,
                                                     color = Barcode)) +
      geom_point(size = 4) +
      theme_bw()
    })

  # Download the Bray Curtis PCoA plot
  output$braycurtis_butt <- downloadHandler(
    filename = function() {
      paste("bray_curtis", ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      print(ggplot(data = data()$bray_curtis_pcoa_dat, aes(x = PCoA_1,
                                                           y = PCoA_2,
                                                           color = Barcode)) +
      geom_point(size = 4) +
      theme_bw())
      dev.off()
    }
  )

  # Display the Bray Curtis PCoA plot help
  observeEvent(input$braycurtis_help, {
    showModal(modalDialog(
      title = "Bray Curtis PCoA Plot Help",
      "Bray Curtis PCoA plots are used to visualize the similarity between
      samples based on species present. The x-axis represents the first 
      principal coordinate and they-axis represents the second principal 
      coordinate. The legend representsthe sample ID.",
      easyClose = TRUE,
      footer = NULL
    ))
  })
}

# Run the app -----------------------------------------------------------------
shinyApp(ui, server)
