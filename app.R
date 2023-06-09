library(shiny)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(vegan)
library(shinycssloaders)
library(rlang)
library(bslib)
library(plotly)

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
  theme = bs_theme(bootswatch = "sandstone"),
  navbarPage(
    title = "EPI2MEviz",
    tabPanel("Upload EPI2ME .csv",
      p(HTML("This app is designed to process the classifications 
      of reads from the Oxford Nanopore Technologies 16S EPI2ME analysis. 
      The app will analyze the data based on user input accuracy of 
      reads and the barcodes selected. The app will then generate a 
      rarefaction curve, a relative abundance plot, and a Bray Curtis
      PCoA plot. The app will also generate a table of the data for 
      each of the plots. The inputs for the app are the EPI2ME .csv file,
      the average accuracy of the reads, and the barcodes to analyze 
      recovered from the EPI2ME dashboard. You may also provide an 
      optional metadata file to add additional filtering options. The 
      metadata file must be a .csv file with the ''barcode'' in the first 
      column and ''id'' in the second column. The remaing metdata can be
      input into the following columns. 
      <br/> <br/> 
      If you do not remember which barcodes you used, you can leave the
      selection blank and the app will use all barcodes it finds in the 
      data set.
      <br/> <br/> 
      I have included 2 example datasets in the repository for you to test 
      the app. The example data is called ''example_data.csv'' and the 
      example metadata is called ''example_metadata.csv''. You can use 
      the example metadata file as a template for your own metadata file. 
      <br/> <br/> 
      If you have any questions or comments please
      contact me at ehill@iolani.org or create an issue on the
      github repository: 
      https://github.com/ehill-iolani/epi2meviz/issues")),
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
          textInput("filt", "Average EPI2ME Accuracy (1-99)",
            value = "80"),
          checkboxGroupInput("barcodes", HTML("Barcodes to Analyze"),
            choices = bci, inline = TRUE),
          actionButton("bar_help", "Help"),
          fileInput("metadata", "Choose Metadata File",
            accept = c("text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"),
            multiple = FALSE),
          actionButton("meta_help", "Help")
        ),
        mainPanel(
          withSpinner(tableOutput("contents"))
        )
      )
    ),
    navbarMenu("Plots",
      tabPanel("Rarefaction Curve",
        fluidPage(
          headerPanel("Rarefaction Curves"),
          uiOutput("rare_meta"),
          downloadButton("rare_butt", "Download PDF"),
          actionButton("rare_help", "Help"),
          mainPanel(
            withSpinner(plotlyOutput("rarefaction"))
          )
        )
      ),
      tabPanel("Relative Abundance Plot",
        fluidPage(
          headerPanel("Relative Abundance"),
          uiOutput("relab_meta"),
          downloadButton("relab_butt", "Download PDF"),
          actionButton("relab_help", "Help"),
          mainPanel(
            withSpinner(plotlyOutput("relative_abundance"))
          )
        )
      ),
      tabPanel("Bray Curtis PCoA Plot",
        fluidPage(
          headerPanel("Bray Curtis PCoA"),
          uiOutput("pcoa_meta"),
          downloadButton("braycurtis_butt", "Download PDF"),
          actionButton("braycurtis_help", "Help"),
          mainPanel(
            withSpinner(plotlyOutput("bray_curtis"))
          )
        )
      )
    ),
    navbarMenu("Tables",
      tabPanel("Rarefaction Table",
        fluidPage(
          headerPanel("Rarefaction Table"),
          selectInput("rare_barcodes", "Select Barcode",
            choices = ""),
          downloadButton("raref_butt", "Download CSV"),
          actionButton("raref_help", "Help"),
          mainPanel(
            withSpinner(tableOutput("rare"))
          )
        )
      ),
      tabPanel("Relative Abundance Table",
        fluidPage(
          headerPanel("Relative Abundance Table"),
          selectInput("relab_barcode", "Select Barcode",
            choices = ""),
          downloadButton("relabf_butt", "Download CSV"),
          actionButton("relabf_help", "Help"),
          mainPanel(
            withSpinner(tableOutput("relab"))
          )
        )
      ),
      tabPanel("Bray Curtis PCoA Table",
        fluidPage(
          headerPanel("Bray Curtis PCoA Table"),
          downloadButton("braycurtisf_butt", "Download CSV"),
          actionButton("braycurtisf_help", "Help"),
          mainPanel(
            withSpinner(tableOutput("bray_curtis_tab"))
          )
        )
      )
    )
  )
))

# Define the server -----------------------------------------------------------
server <- function(input, output, session) {
  # Set large upload size limit (server side)
  options(shiny.maxRequestSize = 70 * 1024^2)

  # Input metadata ------------------------------------------------------------
  observeEvent(input$meta_help, {
    showModal(modalDialog(
      title = "Metadata Input",
      HTML("Metadata can be added to the data set by uploading a .csv file with
      the barcodes and their corresponding information. The first column of the 
      .csv file must be the barcode ID the additional columns must be the 
      metadata value. The barcode ID must be the same as the barcode ID in the 
      data. The metadata will be added to the data set and can be used to
      filter the data after the analysis is complete."),
      easyClose = TRUE,
      footer = NULL
    ))
  })

  # Logic for processing EPI2ME data ------------------------------------------

  data <- eventReactive(input$button, {
    # Filtering and cleaning the data -----------------------------------------

    # require that at least 1 input is available
    req(input$file1)

    # Stores the input file as a variable
    inFile <- input$file1

    # Parses the input file(s)
    if (is.list(inFile) == TRUE) {
      dat <- do.call(rbind, lapply(inFile$datapath, function(x) {
                                   read.csv(x, header = TRUE,
                                   sep = input$sep)}))
    } else {
      dat <- read.csv(inFile$datapath, header = TRUE, sep = input$sep)
    }

    # Adds metadata to the data if it is uploaded
    if (length(input$metadata) > 0) {
      meta <- read.csv(input$metadata$datapath, header = TRUE)
      meta <- na.omit(meta)
      dat <- left_join(dat, meta, by = "barcode")
    }

    # Filters the data based on user input; if no barcodes are entered
    # the script will scrap the unique barcodes from the input data
    if (length(input$barcodes) == 0) {
      bcii <- unique(dat$barcode)
      dat <- dat[dat$barcode %in% bcii, ]
    } else {
    dat <- dat[dat$barcode %in% input$barcodes, ]
    }

    # First pass filtering; if the input is <1% or >99% accuracy
    # it will default to 80% accuracy
    if (input$filt < 1 | input$filt > 99) {
      ifilt <- 80
      dat <- dat[dat$accuracy > ifilt, ]
    } else {
      dat <- dat[dat$accuracy > input$filt, ]
      dat <- dat[dat$species != "unclassified", ]
      dat <- dat[dat$barcode != "unclassified", ]
    }

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
    # If metdata is not present
    if (dim(dat)[2] < 12) {
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
    names(rare) <- c("unique_species", "reads_sampled", "barcode")
    names(rare)[3] <- "Barcode"
    rare_meta <- names(rare)[3]
    } else {
    rare <- data.frame()
    for (i in bc) {
      temp <- rdat[rdat$barcode == i, ]
      reads <- seq_along(temp$species)
      urare <- data.frame()
      for (j in seq(1, nrow(temp), rnn)) {
        temp2 <- mean(replicate(2, bc_sample(temp, j)))
        temp2 <- as.data.frame(temp2)
        temp2$reads <- j
        temp2$id <- temp$id[1]
        urare <- rbind(urare, temp2)
      }
      urare$barcode <- i
      rare <- rbind(rare, urare)
    }
    names(rare) <- c("unique_species", "reads_sampled", "id", "barcode")
    names(rare)[4] <- "Barcode"
    rare_meta <- names(rare)[c(4, 3)]
    }

    # Update rarefaction table choices based on barcodes present in data ------
    updateSelectizeInput(session, "rare_barcodes",
                         choices = unique(rare$Barcode)
                                   [order(unique(rare$Barcode))])

  # Calculates Relative abundance ---------------------------------------------

    # Calculates relative abundance per barcode
    # If metdata is not present
    if (dim(dat)[2] < 12) {
    adat <- dat %>% group_by(genus, barcode) %>% tally()
    adat <- adat[adat$n > 1, ]
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
    names(relabf)[3] <- c("Relative Abundance (%)")
    relab_meta <- names(relabf)[2]
    } else {
    # If metadata is present
    adat <- dat[, c(10, 3, 12)]
    adat <- adat %>% group_by(adat[1], adat[2], adat[3]) %>% tally()
    adat <- adat[adat$n > 1, ]
    relab <- data.frame()
    for (i in bc) {
      temp <- adat[adat$barcode == i, ]
      temp$rel_ab <- temp$n / sum(temp$n)
      relab <- rbind(relab, temp)
    }
    # Cleans table for presentation
    relabf <- data.frame(relab$genus, relab$barcode, relab$id, relab$rel_ab)
    names(relabf)[1:4] <- c("Genus", "Barcode", "ID", "rel_ab")
    relabf$rel_ab <- relabf$rel_ab * 100
    relabf <- relabf[with(relabf, order(Barcode, rel_ab, decreasing = TRUE)), ]
    names(relabf)[4] <- c("Relative Abundance (%)")
    relab_meta <- names(relabf)[2:3]
    }

    # Update relative abundance table choices based on barcodes present in data
    updateSelectizeInput(session, "relab_barcode",
                         choices = unique(relabf$Barcode)
                                   [order(unique(relabf$Barcode))])

  # Conducts Bray Curtis PCoA -------------------------------------------------

    # Bray Curtis PCoA
    if (length(bc) < 3) {
      print("Not enough barcodes to perform Bray Curtis PCoA")
    # Logic if metadata is provided
    } else if (dim(dat)[2] > 11) {
      brayc <- dat[, c(7, 3, 12:ncol(dat))]
      brayci <- length(12:ncol(dat)) + 2
      brayc <- brayc %>%
               group_by(brayc[1], brayc[2], brayc[3:ncol(brayc)]) %>%
               tally()
      brayc <- brayc[brayc$n > 1, ]
      brayc <- brayc %>% pivot_wider(names_from = species, values_from = n)

      # Replaces NA with 0
      brayc[is.na(brayc)] <- 0
      bray_out <- vegdist(brayc[brayci:dim(brayc)[2]], method = "bray")

      # Performs PCoA and renames output
      bray_curtis_pcoa <- cmdscale(bray_out, k = 2, eig = TRUE, add = TRUE)
      bray_curtis_pcoa_dat <- as.data.frame(bray_curtis_pcoa$points)
      bray_curtis_pcoa_dat <- cbind(bray_curtis_pcoa_dat,
                                    brayc[1:brayci - 1])
      names(bray_curtis_pcoa_dat)[1:3] <- c("PCoA_1", "PCoA_2", "Barcode")
      pcoa_meta <- names(bray_curtis_pcoa_dat)[3:length(bray_curtis_pcoa_dat)]
    # Logic if no metadata is provided
    } else {
      brayc <- dat[, c("species", "barcode")]
      brayc <- brayc %>% group_by(species, barcode) %>% tally()
      brayc <- brayc[brayc$n > 1, ]
      brayc <- brayc %>% pivot_wider(names_from = species, values_from = n)

      # Replaces NA with 0
      brayc[is.na(brayc)] <- 0
      bray_out <- vegdist(brayc[2:dim(brayc)[2]], method = "bray")

      # Performs PCoA and renames output
      bray_curtis_pcoa <- cmdscale(bray_out, k = 2, eig = TRUE, add = TRUE)
      bray_curtis_pcoa_dat <- as.data.frame(bray_curtis_pcoa$points)
      bray_curtis_pcoa_dat <- cbind(bray_curtis_pcoa_dat,
                                    brayc$barcode)
      names(bray_curtis_pcoa_dat) <- c("PCoA_1", "PCoA_2", "Barcode")
      pcoa_meta <- names(bray_curtis_pcoa_dat)[3:length(bray_curtis_pcoa_dat)]
    }

  # Returns analysis outputs --------------------------------------------------

    # Returns processed data
    if (exists("bray_curtis_pcoa_dat") == FALSE) {
      combo <- list(rdat = rdat,
                    rare = rare,
                    rare_meta = rare_meta,
                    relab = relab,
                    relabf = relabf,
                    relab_meta = relab_meta)
      return(combo)
    } else {
      combo <- list(rdat = rdat,
                    rare = rare,
                    rare_meta = rare_meta,
                    relab = relab,
                    relabf = relabf,
                    relab_meta = relab_meta,
                    bray_curtis_pcoa_dat = bray_curtis_pcoa_dat,
                    pcoa_meta = pcoa_meta)
      return(combo)
    }
  })

  # Help/Visualizing/Downloading outputs --------------------------------------

  # Display the table of processed data
  output$contents <- renderTable({
    temp <- data()$rdat[order(data()$rdat$barcode), ]
    temp[!duplicated(temp$barcode), ]
  })

  # Display barcode help
  observeEvent(input$bar_help, {
    showModal(modalDialog(
      title = "Barcode Selection",
      HTML("Barcodes are the unique identifiers for each sample. They are
      typically 12-16 base pairs long and are used to identify the
      sample in the sequencing run. Barcodes are used to sort the
      sequences into their respective samples. <br/> <br/> 
      If you leave this field blank the app will process all barcodes found 
      in the dataset."),
      easyClose = TRUE,
      footer = NULL
    ))
  })

  # Display rarefaction plot options
  output$rare_meta <- renderUI({
    selectInput(
      "rare_meta",
      "Select Metadata",
      data()$rare_meta
    )
  })

  # Display the rarefaction curve
  output$rarefaction <- renderPlotly({
    ggplotly(ggplot(data = data()$rare, aes(x = reads_sampled,
                                     y = unique_species)) +
    geom_point(size = 2, aes(color = .data[[input$rare_meta]])) +
    labs(x = "Reads Sampled", y = "Unique Species") +
    theme_bw())
  })

  # Download the rarefaction curve
  output$rare_butt <- downloadHandler(
    filename = function() {
      paste("rarefaction_", input$rare_meta, ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      print(ggplot(data = data()$rare, aes(x = reads_sampled,
                                           y = unique_species)) +
      geom_point(size = 2, aes(color = .data[[input$rare_meta]])) +
      labs(x = "Reads Sampled", y = "Unique Species", fill = "Genus") +
      theme_bw())
      dev.off()
    }
  )

  # Display the rarefaction curve help
  observeEvent(input$rare_help, {
    showModal(modalDialog(
      title = "Rarefaction Curve Help",
      HTML("Rarefaction curves are used to determine if enough reads have been
      sampled to capture the diversity of a sample. The x-axis represents the
      number of reads sampled and the y-axis represents the number of unique
      species found. If the curve plateaus, then enough reads have been sampled
      to capture the diversity of the sample. If the curve does not plateau,
      then more reads need to be sampled to capture the diversity of the
      sample."),
      easyClose = TRUE,
      footer = NULL
    ))
  })

  # Display the rarefaction table
  output$rare <- renderTable({
    data()$rare[data()$rare$Barcode == input$rare_barcodes, ]
  })

  # Download the rarefaction table
  output$raref_butt <- downloadHandler(
    filename = function() {
      paste("rarefaction", "_", input$rare_barcodes, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data()$rare[data()$rare$Barcode == input$rare_barcodes, ],
                file, row.names = FALSE)
    }
  )

  # Display the rarefaction table help
  observeEvent(input$raref_help, {
    showModal(modalDialog(
      title = "Rarefaction Table Help",
      HTML("The rarefaction table is the data used to generate rarefaction
      curve seen on the preivous panel. It scales the subsamples of the data
      based the size of the data set. The table logs the subsampling intervale
      and the number of unique species recovered for each step in the 
      subsampling process. The table can be downloaded as a .csv file."),
      easyClose = TRUE,
      footer = NULL
    ))
  })

  # Display relative abundance plot options
  output$relab_meta <- renderUI({
    selectInput("relab_meta",
                "Select Metadata",
                choices = data()$relab_meta)
  })

  # Display the relative abundance plot
  output$relative_abundance <- renderPlotly({
    ggplotly(ggplot(data = data()$relabf, aes(x = .data[[input$relab_meta]],
                                     y = `Relative Abundance (%)`,
                                     fill = Genus)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Sample ID", y = "Relative Abundance (%)", fill = "Genus") +
    theme_bw())
  })

  # Display the relative abundance plot help
  observeEvent(input$relab_help, {
    showModal(modalDialog(
      title = "Relative Abundance Plot Help",
      HTML("Relative abundance plots are used to visualize the relative 
      abundance of each genus in each sample. The x-axis represents the 
      sample ID and the y-axis represents the relative abundance of each 
      genus. The legend represents the genus of each species."),
      easyClose = TRUE,
      footer = NULL
    ))
  })

  # Download the relative abundance plot
  output$relab_butt <- downloadHandler(
    filename = function() {
      paste("relative_abundance_", input$relab_meta, ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      print(ggplot(data = data()$relabf, aes(x = .data[[input$relab_meta]],
                                           y = `Relative Abundance (%)`,
                                           fill = Genus)) +
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
      HTML("The relative abundance table is used to visualize the relative 
      abundance of each genus in each sample. The table is sorted by barcode 
      and then relative abundance. The table can be downloaded as a 
      .csv file."),
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
      paste("relative_abundance", "_", input$relab_barcode, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data()$relabf[data()$relabf$Barcode == input$relab_barcode, ],
                file, row.names = FALSE)
    }
  )

  # Display the Bray Curtis PCoA plot options
  output$pcoa_meta <- renderUI({
    selectInput(
      "pcoa_meta",
      "Select Metadata",
      data()$pcoa_meta)
  })

  # Display the Bray Curtis PCoA plot w/ options
  output$bray_curtis <- renderPlotly({
    ggplotly(ggplot(data = data()$bray_curtis_pcoa_dat, aes(x = PCoA_1,
                                                     y = PCoA_2)) +
      geom_point(size = 4, aes(color = .data[[input$pcoa_meta]])) +
    theme_bw())
  })

  # Download the Bray Curtis PCoA plot
  output$braycurtis_butt <- downloadHandler(
    filename = function() {
      paste("bray_curtis_", input$pcoa_meta, ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file)
      print(ggplot(data = data()$bray_curtis_pcoa_dat, aes(x = PCoA_1,
                                                           y = PCoA_2)) +
      geom_point(size = 4, aes(color = .data[[input$pcoa_meta]])) +
      theme_bw())
      dev.off()
    }
  )

  # Display the Bray Curtis PCoA plot help
  observeEvent(input$braycurtis_help, {
    showModal(modalDialog(
      title = "Bray Curtis PCoA Plot Help",
      HTML("Bray Curtis PCoA plots are used to visualize the similarity between
      samples based on species present. The x-axis represents the first 
      principal coordinate and the y-axis represents the second principal 
      coordinate. The legend represents the sample ID. <br/> <br/>
      This will not run if there are less than 3 barcodes."),
      easyClose = TRUE,
      footer = NULL
    ))
  })

  # Display the Bray Curtis PCoA table
  output$bray_curtis_tab <- renderTable({
    data()$bray_curtis_pcoa_dat
  })

  # Download the Bray Curtis PCoA table
  output$braycurtisf_butt <- downloadHandler(
    filename = function() {
      paste("bray_curtis", "_", "pcoa", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data()$bray_curtis_pcoa_dat, file, row.names = FALSE)
    }
  )

  # Display the Bray Curtis PCoA table help
  observeEvent(input$braycurtisf_help, {
    showModal(modalDialog(
      title = "Bray Curtis PCoA Table Help",
      HTML("The Bray Curtis PCoA table is used to visualize the 
      similarity between samples based on species present. The table 
      can be downloaded as a .csv file.  <br/> <br/>
      This will not run if there are less than 3 barcodes."),
      easyClose = TRUE,
      footer = NULL
    ))
  })
}

# Run the app -----------------------------------------------------------------
shinyApp(ui, server)
