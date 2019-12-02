library(shiny)
library(ape)
library(ggplot2)

#Define UI

ui <-fluidPage(
  tags$head(tags$script(src = "message-handler.js")),
  
  #Application Title
  titlePanel("Sequence Overrepresentation (SOR) Index "),
  
  #Sidebar Layout
  sidebarLayout(  
    sidebarPanel(
      helpText("Determine probabilities of overrepresentation of identical
               sequence rakes in a dataset."),
    fileInput("file1", "Choose Fasta File", accept=c('.fa', '.fas', '.txt', '.fasta'), multiple = F),
    selectInput("model", "Substitution Model:",
                list("Hamming"="raw",
                     "Tamura-Nei (TN93)"="TN93",
                     "Jukes-Cantor (JC69)"="JC69",
                     "Kimura 2-param (K80)"="K80")),
    actionButton("submit", "Submit"),
    helpText("Choose a nucleotide sequence file using the \"Choose Fasta 
             File\" option above to input an aligned file. Make sure the
             file is devoid of hypermutants and has no consensus/outgroup
             as these sequences will skew the data. Perform the analysis 
             by uploading your file, choosing a distance model (default:
             absolute number of differences <Hamming>), and click submit.
             "),
    helpText("Output table will have each set of identical sequences, the
             names of the sequences in the set, and the probability that
             a set that size is found given the overall diversity of the
             supplied sequence file. For mathematical details of the
             statistical test, please refer to the citation supplied.")
    ),
    mainPanel(img(src='DRPlogo.png', align="right", height=80, width=290),
              plotOutput("plot"),
              hr(),
              hr(),
              hr(),
              tableOutput("table"),
              hr(),
              hr(),
              textOutput("text"),
              textOutput("text2")
              )
  )
)

server <- function(input, output, session) {
  
  data <- reactiveValues(sequences=NULL,
                         num_seqs=NULL,
                         dist_matrix=NULL,
                         avg_normdist=NULL,
                         dist_table=NULL,
                         avg_len=NULL,
                         prob_0dist=NULL,
                         plot_vals=NULL,
                         rakes=NULL,
                         print_rakes=NULL
                         )
  
  #Analysis Body
  
  observeEvent(input$submit, {
    #Get Sequence File
    
    inFile <- input$file1
    if(is.null(inFile))
      return(NULL)
    data$sequences <- ape::read.dna(inFile$datapath, format="fasta", as.character = T, as.matrix=F)
    if(length(data$sequences) == 0){
      session$sendCustomMessage(type='error', message = 'Please supply valid Fasta File')
      return(NULL)
    }
    
    #Get average length of sequences
    data$num_seqs <- length(data$sequences)
    tmp_len <- 0
    for(i in 1:data$num_seqs){
      tmp_len <- tmp_len + length(data$sequences[[i]]) - sum(data$sequences[[i]] == '-')
    }
    data$avg_len <- tmp_len/data$num_seqs
    
    #Calculate Pairwise Distance and Average Pairwise Distance
    data$dist_matrix <- dist.dna(as.DNAbin(data$sequences), model=input$model, pairwise.deletion = T, as.matrix = T)
    
    avg_dists <- 0
    for(i in 1:(data$num_seqs-1)){
      for(j in (i+1):data$num_seqs){
        avg_dists <- avg_dists + data$dist_matrix[i,j]
      }
    }
    data$avg_normdist <- avg_dists/(i*j/2)
    avg_dist <- data$avg_normdist * data$avg_len
    data$dist_matrix <- data$dist_matrix * data$avg_len
    
    #Isolate Sequence Rakes
    rake_list <- list()
    
    for(i in 2:data$num_seqs){
      in_rake = FALSE
      if(sum(data$dist_matrix[i-1,i:data$num_seqs] <0.5) == 0)
        next
      else{
        if(length(rake_list) > 0){
          for(j in 1:length(rake_list)){
            if(names(data$sequences)[i-1] %in% rake_list[[j]]) {
              in_rake = TRUE
              break
            }
          }
          if(in_rake) next
        } 
        rake_list[[paste0("rake",as.character(length(rake_list)+1))]] <- names(data$sequences)[data$dist_matrix[(i-1),1:data$num_seqs] == 0]
      }
    }
    data$rakes <- rake_list
    
    #Calculate Probabilities
    if(length(data$rakes) >0){  
      data$prob_0dist <- dpois(0, lambda=avg_dist)
      tmp_print <- data.frame(matrix(nrow=0, ncol=3))
      colnames(tmp_print) <- c("rake", "num", "ids")
      for(name in 1:length(names(data$rakes))){
        tmp_print[name,] <- c(names(data$rakes)[name], length(data$rakes[[name]]), paste0(data$rakes[[name]], sep='', collapse=', '))
      }
      data$print_rakes <- tmp_print[order(as.numeric(tmp_print$num), decreasing = T),]
      npairs <- data$num_seqs *(data$num_seqs-1)/2
      for(i in 1:length(data$print_rakes$num)){
        ipairs <- choose(as.numeric(data$print_rakes$num[i]), 2)
        data$print_rakes$prob[i] <- format(round(1 - pbinom(ipairs-1, npairs, data$prob_0dist, lower.tail = T), 16), digits=2)
        if(as.numeric(data$print_rakes$prob[i]) < 10^(-16)){ data$print_rakes$prob[i] <- "<10^-16" }
        else if(as.numeric(data$print_rakes$prob[i]) == 1) { data$print_rakes$prob[i] <- ">0.99" }
      }
      colnames(data$print_rakes) <- c("Rake Number", "Number of Sequences in Rake", "Sequence IDs", "Probability")
    }
    

    #Formatting output
    tmp_table<- data.frame(matrix(ncol = 2, nrow=0))
    colnames(tmp_table) <- c("HD", "count")
    for(i in 0:round(max(data$dist_matrix,0))){
      count <- sum(data$dist_matrix >= i-0.5 & data$dist_matrix < i+0.5)
      tmp_table[i+1,] <- c(i, count)
    }
    tmp_table[1,2] <- tmp_table[1,2]-data$num_seqs
    tmp_table[1:(round(max(data$dist_matrix),0)+1),2] <- tmp_table[1:(round(max(data$dist_matrix),0)+1),2]/2
    freq_vector<- as.vector(rep(tmp_table$HD, tmp_table$count))
    data$dist_table <- as.data.frame(freq_vector)
    colnames(data$dist_table) <- "dist"
    
    data$plot_vals <- TRUE
  })
  

  
  output$plot <- renderPlot({
    if(is.null(data$plot_vals)) return()
    else if(data$plot_vals)
      caption_gen=paste0(
                         "Figure: Graph of sequence pair distances using model: ", input$model, " as substitution model. ",
                         "Sequences pulled had an average p-distance of ", round(data$avg_normdist, 3),
                         " and average length of ", round(data$avg_len, 0), "."
      )
      ggplot(data$dist_table, aes(x=dist)) + geom_histogram()+
      labs(title="Inter-sequence distance distribution",
           caption=caption_gen,
           subtitle=paste0("Average p-Distance: ", as.character(round(data$avg_normdist, 3))),
           x="Genetic Distance (substitutions)",
           y="Number of Sequence Pairs")+
      theme_bw(base_size=18)+
      theme(plot.caption = element_text(hjust=0))
  })
  output$table <- renderTable({
    if(is.null(data$plot_vals)) return()
    data$print_rakes
  })
  output$text <- renderText({
    "Citation: Patro SC, Brandt LD, Bale MJ, Halvas EK, Joseph KW, Shao W,
     Wu X, Guo S, Murrell B, Wiegand A, Spindler J, Raley C, Hautman C,
     Sobolewski M, Fennessey CM, Hu WS, Luke B, Hasson JM, Niyongabo A,
     Capoferri AA, Keele BF, Milush J, Hoh R, Deeks SG, Maldarelli F,
     Hughes SH, Coffin JM, Rausch JW, Mellors JW, Kearney MF. 2019.
     Combined HIV-1 sequence and integration site analysis informs viral
     dynamics and allows reconstruction of replicating viral ancestors.
     Proc Natl Acad Sci U S A doi:10.1073/pnas.1910334116."
    
  })
  output$text2 <- renderText({
    "Please email michael.bale@nih.gov regarding any concerns/questions/bugs."
  })
  
  
  
  
  
}

shinyApp(ui, server)