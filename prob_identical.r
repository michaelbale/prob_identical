require(ape)
require(ggplot2)

SORindex <- function(rake_list, p_values, plot) {
  #S3 constructor for SORindex object
  
  #Attributes
  value <- list(rakes = rake_list, #Rake list
                p.values = p_values, #p-values of test on rakes
                dist.plot = plot #ggplot object
  )
  
  attr(value, "class") <- "SORindex"
  value
}

readSeqSet <- function(inPath, format="fasta") {
  if(!file.exists(inPath)) { stop("Supplied File does not exist") }
  
  sequences <- ape::read.dna(inPath, format=format, as.char = T, as.mat = F)
  
  sequences
}

getRakeList <- function(distMatrix) {
 

  
  #Determine identical rake sizes using getRakeList function
  rake_list <-  list()    
  numSeqs <- nrow(distMatrix)
  #Main Looping structure
  for(i in 2:numSeqs){
    
    #rake dummy variable
    in_rake = FALSE
    
    #If sequence is unique, go to next sequence
    if(sum(distMatrix[i-1,i:numSeqs] <0.5) == 0)
      next
    else{
      
      if(length(rake_list) > 0){
        #Create new Rake
        for(j in 1:length(rake_list)){
          if(dimnames(distMatrix)[[1]][i-1] %in% rake_list[[j]]) {
            in_rake = TRUE
            break
          }
        }
        if(in_rake) next
      }
      #add sequence to Rake ID
      rake_list[[paste0("rake",as.character(length(rake_list)+1))]] <- dimnames(distMatrix)[[1]][distMatrix[(i-1),1:numSeqs] == 0]
    }
  }
  if(length(rake_list) == 0) { stop("Sequences are all unique") }
  rake_list
}

calcSOR <- function(seqSet, model="raw") {
  
  #Assumed all sequences are of same length in the set 
  #Can be changed to get average length or pairwise lengths
  seqLen <- length(seqSet[[1]])
  numSeqs <- length(seqSet)
  #Get pairwise Hamming Distance Matrix
  distMatrix <- dist.dna(as.DNAbin(seqSet), model=model, pairwise.deletion = T, as.matrix = T)
  avgDist <- sum(distMatrix)/(choose(numSeqs,2)*2)
  distMatrix <- distMatrix * seqLen
  
  
  
  rake_list <- getRakeList(distMatrix)
  
  p_values <- calcP(rake_list, numSeqs, avgDist * seqLen)
  
  dist.plot <- makePlot(seqLen, avgDist, distMatrix)
  
  
  out <- SORindex(rake_list, p_values, dist.plot)
  
  out
}

calcP <- function(rake_list, numSeqs, avgDist) {
  
  prob_0dist <- dpois(0, lambda = avgDist)
  npairs <- choose(numSeqs, 2)
  getP <- function(y, npairs. = npairs, prob_0dist. = prob_0dist) {
    size <- length(y)
    ipairs <- choose(size, 2)
    p <- signif(1 - pbinom(ipairs - 1, npairs., prob_0dist., lower.tail = T), 3)
  }
  
  p.vec <- sapply(rake_list, getP)
  p.vec
}

makePlot <- function(seqLen, avgDist, distMatrix) {
  tmp.freq <- as.data.frame(table(round(distMatrix,0)))
  colnames(tmp.freq) <- c("Distance", "Frequency")
  tmp.freq[1,2] <- tmp.freq[1,2] - nrow(distMatrix)
  tmp.freq[, 2] <- tmp.freq[, 2] / 2
  freq.count <- unlist(apply(tmp.freq, 1, function(x) rep(x[1], x[2])))
  #freq.count <- unname(freq.count)
  freq.count <- as.data.frame(as.numeric(freq.count))
  colnames(freq.count) <- c("dist")
  caption <- paste0("Figure: Graph of sequence pair distances. ",
                    "Sequences pulled had an average p-distance of ", round(avgDist, 3),
                    "\nand length of ", round(seqLen, 0), "nt.")
  
  plotDists <- ggplot(data=freq.count, aes(x=dist))+
    geom_histogram()+
    labs(title="Inter-sequence distance distribution",
         caption=caption,
         subtitle=paste0("Average p-Distance: ", as.character(round(avgDist,3))),
         x="Genetic Distance (Substitutions)",
         y="Number of Sequence Pairs")+
    theme_bw(base_size=18)+
    theme(plot.caption = element_text(hjust=0))
  
  plotDists
}


getSOR <- function(inPath, model="raw", format="fasta") { 
  seqSet <- readSeqSet(inPath, format)
  out.SOR <- calcSOR(seqSet, model)
  
  }
