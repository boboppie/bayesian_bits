
# Rfunctions.R
# Written by Avril Coghlan, a.coghlan@ucc.ie

dnasmithwaterman <- function(S1, S2, gapopen, gapextend, mymatch, mymismatch)
{
           library(Biostrings)
           library(seqinr)
           S1 <- s2c(S1) # Convert the string to a vector of characters
           S2 <- s2c(S2) # Convert the string to a vector of characters
           lengthS1 <- length(S1)
           lengthS2 <- length(S2)
           # Define the DNA scoring matrix:
           s1 <- nucleotideSubstitutionMatrix(match = mymatch,
mismatch = mymismatch, baseOnly = TRUE)
           # match is usually +ve, mismatch -ve
           s <- s1[S2, S1]
           d <- gapopen # The gap open penalty is normally -ve
           d2 <- gapextend # The gap extension penalty is normally -ve

           # Make the table T with lengthS1+1 columns and lengthS2+1 rows:
           # eg. if S1 = VIVADAVIS and S2 = VIVALASVEGAS, then we will
           # have a matrix T with 10 columns and 13 rows
           T <- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)
           mycolnames <- c("NA",S1)
           myrownames <- c("NA",S2)
           colnames(T) <- mycolnames
           rownames(T) <- myrownames

           # Make a matrix for storing the traceback:
           T2 <- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)
           colnames(T2) <- mycolnames
           rownames(T2) <- myrownames

           # Initialise the first column and row of matrix T:
           rowlength <- length(T[1,])
           for (k in 1:rowlength)
           {
              if          (k == 1) { T[1,k] <- 0                    }
              else if (k == 2)     { T[1,k] <- 0                    }
              else                 { T[1,k] <- 0                    }
           }
           collength <- length(T[,1])
           for (k in 1:collength)
           {
              if          (k == 1) { T[k,1] <- 0                    }
              else if (k == 2)     { T[k,1] <- 0                    }
              else                 { T[k,1] <- 0                    }
           }

           # Initialise the first column and row of matrix T2:
           T2[1,] <- rep('-',lengthS1+1) # Fill in the first row
           T2[,1] <- rep('|',lengthS2+1) # Fill in the first column
           # print(T2)

           # Carry out the recurrence relation to fill in matrix T:
           for (j in 2:(nrow(T))) # Go through each row (j) at a time
           {
              for (i in 2:(ncol(T))) # Go through each column (i) at a time
              {
                 print(paste("Calculating for col i=",i," row j=",j))
                 # Set the value of T[i,j] ie. in column i and row j:
                 diag <- T[j-1,i-1]+s[j-1,i-1]
                 if (T2[j-1,i] == "|")      { up <- T[j-1,i] + d2   }
                 else                       { up <- T[j-1,i] + d + d2    }
                 # Still use extension penalty for first position
                 if (T2[j,i-1] == "-")      { left <- T[j,i-1] + d2  }
                 else                       { left <- T[j,i-1] + d  + d2 }
                 # Still use extension penalty for first position
                 T[j,i] <- max(c(diag,up,left,0))
                 if (diag < 0 && up < 0 && left < 0)
                 {
                    T2[j,i] <- "+"
                 }
                 else
                 {
                    if (diag > up && diag > left)
                    {
                       T2[j,i] <- ">"
                    }
                    else if (up > diag && up > left)
                    {
                       T2[j,i] <- "|"
                    }
                    else if (left > diag && left > up)
                    {
                        T2[j,i] <- "-"
                    }
                    else if (left == diag && up == diag)
                    {
                        T2[j,i] <- "*"
                    }
                    else if (up == left && up > diag && left > diag)
                    {
                        T2[j,i] <- "L"
                    }
                    else if (up == diag && up > left && diag > left)
                    {
                       T2[j,i] <- "V"
                    }
                   else if (diag == left && diag > up && left > up)
                   {
                      T2[j,i] <- "Z"
                   }
                   else
                   {
                      print(paste("ERROR: up",up,"left",left,"diag",diag))
                   }
                }
              }
           }
          print(T)
          print(T2)
          maxT <- max(T)
          print(paste("maxT=",maxT))

          # Print out T and T2 together:
          T3 <- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)
          colnames(T3) <- mycolnames
           rownames(T3) <- myrownames

          for (j in 2:(nrow(T))) # Go through each row (j) at a time
          {
              for (i in 2:(ncol(T))) # Go through each column (i) at a time
              {
                 value <- T[j,i]
                 value2 <- T2[j,i]
                 value3 <- paste(value,value2)
                 T3[j,i] <- value3
              }
           }
           print(T3)
}

printPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE)
{
     library(Biostrings)
     seq1aln <- pattern(alignment) # Get the alignment for the first sequence
     seq2aln <- subject(alignment) # Get the alignment for the second sequence
     alnlen  <- nchar(seq1aln)     # Find the number of columns in the alignment
     starts  <- seq(1, alnlen, by=chunksize)
     n       <- length(starts)
     seq1alnresidues <- 0
     seq2alnresidues <- 0
     for (i in 1:n) {
        chunkseq1aln <- substring(seq1aln, starts[i], starts[i]+chunksize-1)
        chunkseq2aln <- substring(seq2aln, starts[i], starts[i]+chunksize-1)
        # Find out how many gaps there are in chunkseq1aln:
        gaps1 <- countPattern("-",chunkseq1aln) # countPattern() is
from Biostrings library
        # Find out how many gaps there are in chunkseq2aln:
        gaps2 <- countPattern("-",chunkseq2aln) # countPattern() is
from Biostrings library
        # Calculate how many residues of the first sequence we have
printed so far in the alignment:
        seq1alnresidues <- seq1alnresidues + chunksize - gaps1
        # Calculate how many residues of the second sequence we have
printed so far in the alignment:
        seq2alnresidues <- seq2alnresidues + chunksize - gaps2
        if (returnlist == 'FALSE')
        {
           print(paste(chunkseq1aln,seq1alnresidues))
           print(paste(chunkseq2aln,seq2alnresidues))
           print(paste(' '))
        }
     }
     if (returnlist == 'TRUE')
     {
        vector1 <- s2c(substring(seq1aln, 1, nchar(seq1aln)))
        vector2 <- s2c(substring(seq2aln, 1, nchar(seq2aln)))
        mylist <- list(vector1, vector2)
        return(mylist)
     }
}

# This finds all ORFs in the forward strand of a sequence. If two ORFs
overlap and share a start/stop codon,
# it takes the longer one. Otherwise, overlapping ORFs are kept.
findORFsinSeq <- function(sequence)
{
   library(Biostrings)
   # Make vectors "positions" and "types" containing information on
the positions of ATGs in the sequence:
   mylist           <- findPotentialStartsAndStops(sequence)
   positions        <- mylist[[1]]
   types            <- mylist[[2]]
   # Make vectors "orfstarts" and "orfstops" to store the predicted
start and stop codons of ORFs
   orfstarts        <- numeric()
   orfstops         <- numeric()
   # Make a vector "orflengths" to store the lengths of the ORFs
   orflengths       <- numeric()
   # Print out the positions of ORFs in the sequence:
   numpositions     <- length(positions) # Find the length of vector "positions"
   if (numpositions >= 2)                # There must be at least one
start codon and one stop codon to have an ORF.
   {
      for (i in 1:(numpositions-1))
      {
         posi        <- positions[i]
         typei       <- types[i]
         found       <- 0
         while (found == 0)
         {
            for (j in (i+1):numpositions)
            {
               posj  <- positions[j]
               typej <- types[j]
               posdiff <- posj - posi
               posdiffmod3 <- posdiff %% 3
               orflength <- posj - posi + 3 # Add in the length of the
stop codon
               if (typei == "atg" && (typej == "taa" || typej == "tag"
|| typej == "tga") && posdiffmod3 == 0)
               {
                  # Check if we have already used the stop codon at
posj+2 in an ORF
                  numorfs <- length(orfstops)
                  usedstop <- -1
                  if (numorfs > 0)
                  {
                     for (k in 1:numorfs)
                     {
                        orfstopk <- orfstops[k]
                        if (orfstopk == (posj + 2)) { usedstop <- 1 }
                     }
                  }
                  if (usedstop == -1)
                  {
		     orfstarts <- append(orfstarts, posi, after=length(orfstarts))
                     orfstops <- append(orfstops, posj+2,
after=length(orfstops)) # Including the stop codon.
                     orflengths <- append(orflengths, orflength,
after=length(orflengths))
                  }
                  found <- 1
                  break
               }
	       if (j == numpositions) { found <- 1 }
            }
         }
      }
   }

   # Sort the final ORFs by start position:
   indices           <- order(orfstarts)
   orfstarts         <- orfstarts[indices]
   orfstops          <- orfstops[indices]

   # Find the lengths of the ORFs that we have
   orflengths        <- numeric()
   numorfs           <- length(orfstarts)
   for (i in 1:numorfs)
   {
      orfstart       <- orfstarts[i]
      orfstop        <- orfstops[i]
      orflength      <- orfstop - orfstart + 1
      orflengths     <- append(orflengths,orflength,after=length(orflengths))
   }

   mylist            <- list(orfstarts, orfstops, orflengths)
   return(mylist)
}

findPotentialStartsAndStops <- function(sequence)
{
   codons            <- c("atg", "taa", "tag", "tga") # Define a
vector with the sequences of potential start and stop codons
   # Find the number of occurrences of each type of potential start or
stop codon
   for (i in 1:4)
   {
      codon          <- codons[i]
      occurrences    <- matchPattern(codon, sequence) # Find all
occurrences of codon "codon" in sequence "sequence"
      codonpositions <- attr(occurrences,"start")     # Find the start
positions of all occurrences of "codon" in sequence "sequence"
      numoccurrences <- length(codonpositions)        # Find the total
number of potential start and stop codons in sequence "sequence"
      if (i == 1)
      {
         positions   <- codonpositions                # Make a copy of
vector "codonpositions" called "positions"
         types       <- rep(codon, numoccurrences)    # Make a vector
"types" containing "numoccurrences" copies of "codon"
      }
      else
      {
         # Add the vector "codonpositions" to the end of vector "positions":
         positions   <- append(positions, codonpositions,
after=length(positions))
         # Add the vector "rep(codon, numoccurrences)" to the end of
vector "types":
         types       <- append(types, rep(codon, numoccurrences),
after=length(types))
      }
   }
   # Sort the vectors "positions" and "types" in order of position
along the input sequence:
   indices           <- order(positions)
   positions         <- positions[indices]
   types             <- types[indices]
   # Return a list variable including vectors "positions" and "types":
   mylist            <- list(positions,types)
   return(mylist)
}

# Make a plot of the positions of the potential start and stop codons in
# each of the three reading frames in a sequence
plotPotentialStartsAndStops <- function(sequence)
{
   codons            <- c("atg", "taa", "tag", "tga") # Define a
vector with the sequences of potential start and stop codons
   # Find the number of occurrences of each type of potential start or
stop codon
   for (i in 1:4)
   {
      codon          <- codons[i]
      occurrences    <- matchPattern(codon, sequence) # Find all
occurrences of codon "codon" in sequence "sequence"
      codonpositions <- attr(occurrences,"start")     # Find the start
positions of all occurrences of "codon" in sequence "sequence"
      numoccurrences <- length(codonpositions)        # Find the total
number of potential start and stop codons in sequence "sequence"
      if (i == 1)
      {
         positions   <- codonpositions                # Make a copy of
vector "codonpositions" called "positions"
         types       <- rep(codon, numoccurrences)    # Make a vector
"types" containing "numoccurrences" copies of "codon"
      }
      else
      {
         # Add the vector "codonpositions" to the end of vector "positions":
         positions   <- append(positions, codonpositions,
after=length(positions))
         # Add the vector "rep(codon, numoccurrences)" to the end of
vector "types":
         types       <- append(types, rep(codon, numoccurrences),
after=length(types))
      }
   }
   # Sort the vectors "positions" and "types" in order of position
along the input sequence:
   indices           <- order(positions)
   positions         <- positions[indices]
   types             <- types[indices]
   # Make a plot showing the positions of the start and stop codons in
the input sequence:
   # Draw a line at y=0 from 1 to the length of the sequence:
   x                 <- c(1,nchar(sequence))
   y                 <- c(0,0)
   plot(x, y, ylim=c(0,3), type="l", axes=FALSE, xlab="Nucleotide",
ylab="Reading frame", main="Predicted start (red) and stop (blue)
codons")
   segments(1,1,nchar(sequence),1)
   segments(1,2,nchar(sequence),2)
   # Add the x-axis at y=0:
   axis(1, pos=0)
   # Add the y-axis labels:
   text(0.9,0.5,"+1")
   text(0.9,1.5,"+2")
   text(0.9,2.5,"+3")
   # Draw in each predicted start/stop codon:
   numcodons         <- length(positions)
   for (i in 1:numcodons)
   {
      position       <- positions[i]
      type           <- types[i]
      remainder      <- (position-1) %% 3
      if    (remainder == 0) # +1 reading frame
      {
         if (type == "atg") { segments(position,0,position,1,lwd=1,col="red") }
         else               { segments(position,0,position,1,lwd=1,col="blue")}
      }
      else if (remainder == 1)
      {
         if (type == "atg") { segments(position,1,position,2,lwd=1,col="red") }
         else               { segments(position,1,position,2,lwd=1,col="blue")}
      }
      else if (remainder == 2)
      {
         if (type == "atg") { segments(position,2,position,3,lwd=1,col="red") }
         else               { segments(position,2,position,3,lwd=1,col="blue")}
      }
   }
}

# This plots all ORFs in the forward strand of a sequence. If two ORFs
overlap and share a start/stop codon,
# it takes the longer one. Otherwise, overlapping ORFs are kept.
plotORFsinSeq <- function(sequence)
{
   # Make vectors "positions" and "types" containing information on
the positions of ATGs in the sequence:
   mylist           <- findPotentialStartsAndStops(sequence)
   positions        <- mylist[[1]]
   types            <- mylist[[2]]
   # Make vectors "orfstarts" and "orfstops" to store the predicted
start and stop codons of ORFs
   orfstarts        <- numeric()
   orfstops         <- numeric()
   # Make a vector "orflengths" to store the lengths of the ORFs
   orflengths       <- numeric()
   # Print out the positions of ORFs in the sequence:
   numpositions     <- length(positions) # Find the length of vector "positions"
   if (numpositions >= 2)                # There must be at least one
start codon and one stop codon to have an ORF.
   {
      for (i in 1:(numpositions-1))
      {
         posi        <- positions[i]
         typei       <- types[i]
         found       <- 0
         while (found == 0)
         {
            for (j in (i+1):numpositions)
            {
               posj  <- positions[j]
               typej <- types[j]
               posdiff <- posj - posi
               posdiffmod3 <- posdiff %% 3
               orflength <- posj - posi + 3 # Add in the length of the
stop codon
               if (typei == "atg" && (typej == "taa" || typej == "tag"
|| typej == "tga") && posdiffmod3 == 0)
               {
                  # Check if we have already used the stop codon at
posj+2 in an ORF
                  numorfs <- length(orfstops)
                  usedstop <- -1
                  if (numorfs > 0)
                  {
                     for (k in 1:numorfs)
                     {
                        orfstopk <- orfstops[k]
                        if (orfstopk == (posj + 2)) { usedstop <- 1 }
                     }
                  }
                  if (usedstop == -1)
                  {
		     orfstarts <- append(orfstarts, posi, after=length(orfstarts))
                     orfstops <- append(orfstops, posj+2,
after=length(orfstops)) # Including the stop codon.
                     orflengths <- append(orflengths, orflength,
after=length(orflengths))
                  }
                  found <- 1
                  break
               }
	       if (j == numpositions) { found <- 1 }
            }
         }
      }
   }

   # Sort the final ORFs by start position:
   indices           <- order(orfstarts)
   orfstarts         <- orfstarts[indices]
   orfstops          <- orfstops[indices]

   # Make a plot showing the positions of ORFs in the input sequence:
   # Draw a line at y=0 from 1 to the length of the sequence:
   x                 <- c(1,nchar(sequence))
   y                 <- c(0,0)
   plot(x, y, ylim=c(0,3), type="l", axes=FALSE, xlab="Nucleotide",
ylab="Reading frame", main="Predicted ORFs")
   segments(1,1,nchar(sequence),1)
   segments(1,2,nchar(sequence),2)
   # Add the x-axis at y=0:
   axis(1, pos=0)
   # Add the y-axis labels:
   text(0.9,0.5,"+1")
   text(0.9,1.5,"+2")
   text(0.9,2.5,"+3")
   # Make a plot of the ORFs in the sequence:
   numorfs           <- length(orfstarts)
   for (i in 1:numorfs)
   {
      orfstart       <- orfstarts[i]
      orfstop        <- orfstops[i]
      remainder      <- (orfstart-1) %% 3
      if    (remainder == 0) # +1 reading frame
      {
         rect(orfstart,0,orfstop,1,col="cyan",border="black")
      }
      else if (remainder == 1)
      {
         rect(orfstart,1,orfstop,2,col="cyan",border="black")
      }
      else if (remainder == 2)
      {
         rect(orfstart,2,orfstop,3,col="cyan",border="black")
      }
   }

}

# This carries out the Viterbi algorithm.
# Adapted from "Applied Statistics for Bioinformatics using R" by Wim
P. Krijnen, page 209
# ( cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf )
viterbi <- function(sequence, transitionmatrix, emissionmatrix)
{
   # Get the names of the states in the HMM:
   states <- rownames(theemissionmatrix)

   # Make the Viterbi matrix v:
   v <- makeViterbimat(sequence, transitionmatrix, emissionmatrix)

   # Go through each of the rows of the matrix v (where each row represents
   # a position in the DNA sequence), and find out which column has the
   # maximum value for that row (where each column represents one state of
   # the HMM):
   mostprobablestatepath <- apply(v, 1, function(x) which.max(x))

   # Print out the most probable state path:
   prevnucleotide <- sequence[1]
   prevmostprobablestate <- mostprobablestatepath[1]
   prevmostprobablestatename <- states[prevmostprobablestate]
   startpos <- 1
   for (i in 2:length(sequence))
   {
      nucleotide <- sequence[i]
      mostprobablestate <- mostprobablestatepath[i]
      mostprobablestatename <- states[mostprobablestate]
      if (mostprobablestatename != prevmostprobablestatename)
      {
         print(paste("Positions",startpos,"-",(i-1), "Most probable
state = ", prevmostprobablestatename))
         startpos <- i
      }
      prevnucleotide <- nucleotide
      prevmostprobablestatename <- mostprobablestatename
   }
   print(paste("Positions",startpos,"-",i, "Most probable state = ",
prevmostprobablestatename))
}

# This makes the matrix v using the Viterbi algorithm.
# Adapted from "Applied Statistics for Bioinformatics using R" by Wim
P. Krijnen, page 209
# ( cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf )
makeViterbimat <- function(sequence, transitionmatrix, emissionmatrix)
{
   # Change the sequence to uppercase
   sequence <- toupper(sequence)
   # Find out how many states are in the HMM
   numstates <- dim(transitionmatrix)[1]
   # Make a matrix with as many rows as positions in the sequence, and as many
   # columns as states in the HMM
   v <- matrix(NA, nrow = length(sequence), ncol = dim(transitionmatrix)[1])
   # Set the values in the first row of matrix v (representing the
first position of the sequence) to 0
   v[1, ] <- 0
   # Set the value in the first row of matrix v, first column to 1
   v[1,1] <- 1
   # Fill in the matrix v:
   for (i in 2:length(sequence)) # For each position in the DNA sequence:
   {
      for (l in 1:numstates) # For each of the states of in the HMM:
      {
         # Find the probabilility, if we are in state l, of choosing
the nucleotide at position in the sequence
         statelprobnucleotidei <- emissionmatrix[l,sequence[i]]

         # v[(i-1),] gives the values of v for the (i-1)th row of v,
ie. the (i-1)th position in the sequence.
         # In v[(i-1),] there are values of v at the (i-1)th row of
the sequence for each possible state k.
         # v[(i-1),k] gives the value of v at the (i-1)th row of the
sequence for a particular state k.

         # transitionmatrix[l,] gives the values in the lth row of the
transition matrix, xx should not be transitionmatrix[,l]?
         # probabilities of changing from a previous state k to a
current state l.

         # max(v[(i-1),] * transitionmatrix[l,]) is the maximum
probability for the nucleotide observed
         # at the previous position in the sequence in state k,
followed by a transition from previous
         # state k to current state l at the current nucleotide position.

         # Set the value in matrix v for row i (nucleotide position
i), column l (state l) to be:
         v[i,l] <-  statelprobnucleotidei * max(v[(i-1),] *
transitionmatrix[,l])
      }
   }

   return(v)
}

# Function to generate a DNA sequence, given a HMM and the length of
the sequence to be generated.
generatehmmseq <- function(transitionmatrix, emissionmatrix,
initialprobs, seqlength)
{
    nucleotides     <- c("A", "C", "G", "T")   # Define the alphabet
of nucleotides
    states          <- c("AT-rich", "GC-rich") # Define the names of the states
    mysequence      <- character()             # Create a vector for
storing the new sequence
    mystates        <- character()             # Create a vector for
storing the state that each position in the new sequence
                                               # was generated by
    # Choose the state for the first position in the sequence:
    firststate      <- sample(states, 1, rep=TRUE, prob=initialprobs)
    # Get the probabilities of the current nucleotide, given that we
are in the state "firststate":
    probabilities   <- emissionmatrix[firststate,]
    # Choose the nucleotide for the first position in the sequence:
    firstnucleotide <- sample(nucleotides, 1, rep=TRUE, prob=probabilities)
    mysequence[1]   <- firstnucleotide         # Store the nucleotide
for the first position of the sequence
    mystates[1]     <- firststate              # Store the state that
the first position in the sequence was generated by

    for (i in 2:seqlength)
    {
       prevstate    <- mystates[i-1]           # Get the state that
the previous nucleotide in the sequence was generated by
       # Get the probabilities of the current state, given that the
previous nucleotide was generated by state "prevstate"
       stateprobs   <- transitionmatrix[prevstate,]
       # Choose the state for the ith position in the sequence:
       state        <- sample(states, 1, rep=TRUE, prob=stateprobs)
       # Get the probabilities of the current nucleotide, given that
we are in the state "state":
       probabilities <- emissionmatrix[state,]
       # Choose the nucleotide for the ith position in the sequence:
       nucleotide   <- sample(nucleotides, 1, rep=TRUE, prob=probabilities)
       mysequence[i] <- nucleotide             # Store the nucleotide
for the current position of the sequence
       mystates[i]  <- state                   # Store the state that
the current position in the sequence was generated by
    }

    for (i in 1:length(mysequence))
    {
       nucleotide   <- mysequence[i]
       state        <- mystates[i]
       print(paste("Position", i, ", State", state, ", Nucleotide = ",
nucleotide))
    }
}

# Function to colour the edges in a path on a graph in colour "colour"
plotpath <- function(graphplot,path,colour)
{
   # Get the names of the vertices in the shortest path:
   vertices <- path$path_detail
   numvertices <- length(vertices)
   # Now make a list that contains the names of pairs of vertices that
are joined by edges:
   # There should be numvertices-1 pairs
   numedges <- numvertices - 1
   edges <- numeric(numedges)
   edges2 <- numeric(numedges)
   myvector <- vector()
   for (i in 1:numedges)
   {
      vertex1 <- vertices[i]
      vertex2 <- vertices[i+1]
      edge <- paste(vertex1,"~",vertex2,sep="")
      edge2 <- paste(vertex2,"~",vertex1,sep="")
      myvector[`edge`] <- colour   # Add named element to myvector
      myvector[`edge2`] <- colour  # Add named element to myvector
   }
   # Set the colour of the edges:
   edgeRenderInfo(graphplot) = list(col=myvector)
   renderGraph(graphplot)
}

# Function to make a graph based on protein-protein interaction data
in an input file
makeproteingraph <- function(myfile)
{
   library("graph")
   mytable <- read.table(file(myfile)) # Store the data in a data frame
   proteins1 <- mytable$V1
   proteins2 <- mytable$V2
   protnames <- c(levels(proteins1),levels(proteins2))
   # Find out how many pairs of proteins there are
   numpairs <- length(proteins1)
   # Find the unique protein names:
   uniquenames <-  unique(protnames)
   # Make a graph for these proteins with no edges:
   mygraph <- new("graphNEL", nodes = uniquenames)
   # Add edges to the graph:
   # See http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf for
more examples
   weights <- rep(1,numpairs)
   mygraph2 <- addEdge(as.vector(proteins1),as.vector(proteins2),mygraph,weights)

   return(mygraph2)
}

# Function to read in a regulatory network based on transcription
factor-target gene interaction
# data in an input file
readRegulatoryNetwork <- function(myfile)
{
   library("graph")
   mytable <- read.table(file(myfile),header=FALSE) # Store the data
in a data frame
   proteins <- mytable$V1
   # Change the protein names to gene names:
   proteins2 <- vector()
   for (i in 1:length(proteins))
   {
      protein <- as.character(proteins[i])
      substr(protein,1,1) <- tolower(substring(protein,1,1))
      proteins2[i] <- protein
   }
   genes <- mytable$V2
   types <- mytable$V3
   # Find out the number of edges:
   numedges <- length(proteins)
   # Find out how many transcription factor and target genes there are:
   nodenames <- union(proteins2,genes)
   nodenames <- unique(nodenames)
   # Make a graph for these transcription factors and target genes
with no edges:
   mygraph <- new("graphNEL", nodes=nodenames, edgemode="directed")
   # Add edges to the graph:
   # See http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf for
more examples
   weights <- rep(1,numedges)
   mygraph2 <- addEdge(as.vector(proteins2),as.vector(genes),mygraph,weights)
   # We can store the information on how to plot the graph:
   # Unfortunately, if we make a subgraph, this is lost
   # Set the types of the edges
   myvector <- vector()
   for (i in 1:numedges)
   {
      vertex1 <- proteins2[i]
      vertex2 <- genes[i]
      edge <- paste(vertex1,"~",vertex2,sep="")
      type <- types[i]
      if      (type == '+') { edgetype <- 'normal'  }
      else if (type == '-') { edgetype <- 'tee'     }
      else                  { edgetype <- 'normal'  } # both
repressive+activating or unknown xxx
      myvector[`edge`] <- edgetype # Add named element to myvector
   }
   # Set the type of the edges:
   edgeRenderInfo(mygraph2) = list(arrowhead=myvector)
   return(mygraph2)
}

# Function to make a regulatory network based on transcription
factor-target gene interaction
# data in an input file, and plot it
readAndPlotRegulatoryNetwork <- function(myfile)
{
   library("graph")
   mytable <- read.table(file(myfile),header=FALSE) # Store the data
in a data frame
   proteins <- mytable$V1
   # Change the protein names to gene names:
   proteins2 <- vector()
   for (i in 1:length(proteins))
   {
      protein <- as.character(proteins[i])
      substr(protein,1,1) <- tolower(substring(protein,1,1))
      proteins2[i] <- protein
   }
   genes <- mytable$V2
   types <- mytable$V3
   # Find out the number of edges:
   numedges <- length(proteins)
   # Find out how many transcription factor and target genes there are:
   nodenames <- union(proteins2,genes)
   nodenames <- unique(nodenames)
   # Make a graph for these transcription factors and target genes
with no edges:
   mygraph <- new("graphNEL", nodes=nodenames, edgemode="directed")
   # Add edges to the graph:
   # See http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf for
more examples
   weights <- rep(1,numedges)
   mygraph2 <- addEdge(as.vector(proteins2),as.vector(genes),mygraph,weights)
   # Set the types of the edges
   myvector <- vector()
   for (i in 1:numedges)
   {
      vertex1 <- proteins2[i]
      vertex2 <- genes[i]
      edge <- paste(vertex1,"~",vertex2,sep="")
      type <- types[i]
      if      (type == '+') { edgetype <- 'normal'  }
      else if (type == '-') { edgetype <- 'tee'     }
      else                  { edgetype <- 'normal'  } # xxx both
repressive and activating, or unknown
      myvector[`edge`] <- edgetype # Add named element to myvector
   }
   # Plot the graph:
   library("Rgraphviz")
   plot.new() # Make a new plot
   mygraphplot <- layoutGraph(mygraph2, layoutType="neato")
   # Set the colour of the edges:
   edgeRenderInfo(mygraphplot) = list(arrowhead=myvector)
   renderGraph(mygraphplot)
   return(mygraph2)
}

# Function to find network communities in a graph
findcommunities <- function(mygraph,minsize)
{
   # Load up the igraph library:
   library("igraph")
   # Set the counter for the number of communities:
   cnt <- 0
   # First find the connected components in the graph:
   myconnectedcomponents <- connectedComp(mygraph)
   # For each connected component, find the communities within that
connected component:
   numconnectedcomponents <- length(myconnectedcomponents)
   for (i in 1:numconnectedcomponents)
   {
      component <- myconnectedcomponents[[i]]
      # Find the number of nodes in this connected component:
      numnodes <- length(component)
      if (numnodes > 1) # We can only find communities if there is
more than one node
      {
         mysubgraph <- subGraph(component, mygraph)
         # Find the communities within this connected component:
         # print(component)
         myvector <- vector()
         mylist <- findcommunities2(mysubgraph,cnt,"FALSE",myvector,minsize)
         cnt <- mylist[[1]]
         myvector <- mylist[[2]]
      }
   }
   print(paste("There were",cnt,"communities in the input graph"))
}

# Function to find network communities in a connected component of a graph
findcommunities2 <- function(mygraph,cnt,plot,myvector,minsize)
{
   # Find the number of nodes in the input graph
   nodes <- nodes(mygraph)
   numnodes <- length(nodes)

   # Record the vertex number for each vertex name
   myvector <- vector()
   for (i in 1:numnodes)
   {
      node <- nodes[i] # "node" is the vertex name, i is the vertex number
      myvector[`node`] <- i  # Add named element to myvector
   }

   # Create a graph in the "igraph" library format, with numnodes nodes:
   newgraph <- graph.empty(n=numnodes,directed=FALSE)
   # First record which edges we have seen already in the "mymatrix" matrix,
   # so that we don't add any edge twice:
   mymatrix <- matrix(nrow=numnodes,ncol=numnodes)
   for (i in 1:numnodes)
   {
      for (j in 1:numnodes)
      {
         mymatrix[i,j] = 0
         mymatrix[j,i] = 0
      }
   }
   # Now add edges to the graph "newgraph":
   for (i in 1:numnodes)
   {
      node <- nodes[i] # "node" is the vertex name, i is the vertex number
      # Find the nodes that this node is joined to:
      neighbours <- adj(mygraph, node)
      neighbours <- neighbours[[1]] # Get the list of neighbours
      numneighbours <- length(neighbours)
      if (numneighbours >= 1) # If this node "node" has some edges to
other nodes
      {
         for (j in 1:numneighbours)
         {
            neighbour <- neighbours[j]
            # Get the vertex number
            neighbourindex <- myvector[neighbour]
            # Add a node in the new graph "newgraph" between vertices
(i-1) and (neighbourindex-1)
            # In graph "newgraph", the vertices are counted from 0 upwards.
            indexi <- i
            indexj <- neighbourindex
            # If we have not seen this edge already:
            if (mymatrix[indexi,indexj] == 0 && mymatrix[indexj,indexi] == 0)
            {
               mymatrix[indexi,indexj] <- 1
               mymatrix[indexj,indexi] <- 1
               # Add edges to the graph "newgraph"
               newgraph <- add.edges(newgraph, c(i-1, neighbourindex-1))
            }
         }
      }
   }
   # Set the names of the vertices in graph "newgraph":
   newgraph <- set.vertex.attribute(newgraph, "name", value=nodes)

   # Now find communities in the graph:
   communities <- spinglass.community(newgraph)
   # Find how many communities there are:
   sizecommunities <- communities$csize
   numcommunities <- length(sizecommunities)
   # Find which vertices belong to which communities:
   membership <- communities$membership
   # Get the names of vertices in the graph "newgraph":
   vertexnames <- get.vertex.attribute(newgraph, "name")

   # Print out the vertices belonging to each community:
   for (i in 1:numcommunities)
   {
      if (plot == TRUE) { cnt <- cnt + 1 }
      nummembers <- 0
      printout <- paste("Community",cnt,":")
      for (j in 1:length(membership))
      {
         community <- membership[j]
         community <- community + 1
         if (community == i) # If vertex j belongs to the ith community
         {
            vertexname <- vertexnames[j]
            if (plot == FALSE)
            {
               nummembers <- nummembers + 1
               # Print out the vertices belonging to the community
               printout <- paste(printout,vertexname)
            }
            else
            {
               # Colour in the vertices belonging to the community
               myvector[`vertexname`] <- cnt
            }
         }
      }
      if (plot == FALSE && nummembers >= minsize)
      {
         cnt <- cnt + 1
         print(printout)
      }
   }

   return(list(cnt,myvector))
}

# Function to plot network communities in a graph
plotcommunities <- function(mygraph)
{
   # Load the "igraph" library:
   library("igraph")
   # Make a plot of the graph
   graphplot <- layoutGraph(mygraph, layoutType="neato")
   renderGraph(graphplot)
   # Get the names of the nodes in the graph:
   vertices <- nodes(mygraph)
   numvertices <- length(vertices)
   # Now record the colour of each vertex in a vector "myvector":
   myvector <- vector()
   colour <- "red"
   for (i in 1:numvertices)
   {
      vertex <- vertices[i]
      myvector[`vertex`] <- colour   # Add named element to myvector
   }

   # Set the counter for the number of communities:
   cnt <- 0
   # First find the connected components in the graph:
   myconnectedcomponents <- connectedComp(mygraph)
   # For each connected component, find the communities within that
connected component:
   numconnectedcomponents <- length(myconnectedcomponents)
   for (i in 1:numconnectedcomponents)
   {
      component <- myconnectedcomponents[[i]]
      # Find the number of nodes in this connected component:
      numnodes <- length(component)
      if (numnodes > 1) # We can only find communities if there is
more than one node
      {
         mysubgraph <- subGraph(component, mygraph)
         # Find the communities within this connected component:
         # print(component)
         mylist <- findcommunities2(mysubgraph,cnt,"TRUE",myvector,0)
         cnt <- mylist[[1]]
         myvector <- mylist[[2]]
      }
   }

   # Get a set of cnt colours, where cnt is equal to the number of
communities found:
   mycolours <- rainbow(cnt)
   for (i in 1:length(mycolours))
   {
      # print(paste("Community",i,"is in",mycolours[i]))
   }
   # Set the colour of the vertices, so that vertices in each
community are of the same colour,
   # and vertices in different communities are different colours:
   myvector2 <- vector()
   for (i in 1:numvertices)
   {
      vertex <- vertices[i]
      community <- myvector[vertex]
      mycolour <- mycolours[community]
      # print(paste("Vertex",vertex,"is in community",community,"of
colour",mycolour))
      myvector2[`vertex`] <- mycolour
   }
   nodeRenderInfo(graphplot) = list(fill=myvector2)
   renderGraph(graphplot)
}

# Function to make a random graph
makerandomgraph <- function(numvertices,numedges)
{
   library("graph")
   # Make a vector with the names of the vertices
   mynames <- sapply(seq(1,numvertices),toString)
   myrandomgraph <- randomEGraph(mynames, edges = numedges)

   return(myrandomgraph)
}

# Function to find the connected component that contains a particular vertex
findcomponent <- function(graph,vertex)
{
   library("RBGL")
   found <- 0
   myconnectedcomponents <- connectedComp(graph)
   numconnectedcomponents <- length(myconnectedcomponents)
   for (i in 1:numconnectedcomponents)
   {
      componenti <- myconnectedcomponents[[i]]
      numvertices <- length(componenti)
      for (j in 1:numvertices)
      {
         vertexj <- componenti[j]
         if (vertexj == vertex)
         {
            found <- 1
            return(componenti)
         }
      }
   }
   print("ERROR: did not find vertex in the graph")
}

makeigraphgraph <- function(mygraph)
{
   library("igraph")
   # Find the number of nodes in the input graph
   nodes <- nodes(mygraph)
   numnodes <- length(nodes)

   # Record the vertex number for each vertex name
   myvector <- vector()
   for (i in 1:numnodes)
   {
      node <- nodes[i] # "node" is the vertex name, i is the vertex number
      myvector[`node`] <- i  # Add named element to myvector
   }

   # Create a graph in the "igraph" library format, with numnodes nodes:
   newgraph <- graph.empty(n=numnodes,directed=FALSE)
   # First record which edges we have seen already in the "mymatrix" matrix,
   # so that we don't add any edge twice:
   mymatrix <- matrix(nrow=numnodes,ncol=numnodes)
   for (i in 1:numnodes)
   {
      for (j in 1:numnodes)
      {
         mymatrix[i,j] = 0
         mymatrix[j,i] = 0
      }
   }
   # Now add edges to the graph "newgraph":
   for (i in 1:numnodes)
   {
      node <- nodes[i] # "node" is the vertex name, i is the vertex number
      # Find the nodes that this node is joined to:
      neighbours <- adj(mygraph, node)
      neighbours <- neighbours[[1]] # Get the list of neighbours
      numneighbours <- length(neighbours)
      if (numneighbours >= 1) # If this node "node" has some edges to
other nodes
      {
         for (j in 1:numneighbours)
         {
            neighbour <- neighbours[j]
            # Get the vertex number
            neighbourindex <- myvector[neighbour]
            # Add a node in the new graph "newgraph" between vertices
(i-1) and (neighbourindex-1)
            # In graph "newgraph", the vertices are counted from 0 upwards.
            indexi <- i
            indexj <- neighbourindex
            # If we have not seen this edge already:
            if (mymatrix[indexi,indexj] == 0 && mymatrix[indexj,indexi] == 0)
            {
               mymatrix[indexi,indexj] <- 1
               mymatrix[indexj,indexi] <- 1
               # Add edges to the graph "newgraph"
               newgraph <- add.edges(newgraph, c(i-1, neighbourindex-1))
            }
         }
      }
   }
   # Set the names of the vertices in graph "newgraph":
   newgraph <- set.vertex.attribute(newgraph, "name", value=nodes)

   return(newgraph)
}

# Function to print out the interactions of a particular protein, in a
particular experimental data set:
# exp is a vector containing the experimental data.
printproteininteractions <- function(exp, protein)
{
	for (i in 1:length(exp))
	{
		pair <- exp[i]
		thepair <- unlist(strsplit(pair, "~", fixed = TRUE))
		protein1 <- thepair[1]
		protein2 <- thepair[2]
		if (protein1 == protein)
		{
			print(thepair)
		}
		else if (protein2 == protein)
		{
			print(thepair)
		}
	}
}

# Function to use a Bayes network for protein-protein interaction data,
# to calculate likelihoods for particular protein protein interactions.
# "goldpos" is the set of gold standard positive protein-protein interactions
# "goldneg" is the set of gold standard negative protein-protein interactions
# "exp" is a list of vectors, each vector containing some experimental
data on protein-protein interactions
#   reported by a particular experiment
bayesforproteinpairs <- function(goldpos, goldneg, exp, protein1,
protein2, VALUEONLY=FALSE)
{
    # Make a vector containing the pair of interest:
    myvector <- vector()
    myvector[1] <- protein1
    myvector[2] <- protein2
    myvector <- sort(myvector)
    protpair2 <- paste(myvector[1],myvector[2],sep="~")

    # Calculate a likelihood ratio for protpair2:

    # Check how many gold standard positive pairs have the same
    # pattern of 'yes's & 'no's as protpair2 for the experimental data sets:
    numexpdatasets <- length(exp) # Number of experimental data sets
    numgoldpos <- length(goldpos) # Number of gold standard positives
    samepattern <- goldpos
    for (m in 1:numexpdatasets)
    {	
	expm <- exp[[m]]
	isfoundm <- is.element(protpair2,expm) # Says if protpair2 is in the
mth experimental data set
	if (isfoundm == "TRUE")
	{
	   samepattern <- intersect(samepattern,expm) # Find pairs that are
in goldpos and in expm
	}
	else
	{
	   samepattern <- setdiff(samepattern,expm) # Find pairs that are in
goldpos but not in expm
	}
    }
    #print(paste("Have same pattern in gold pos and expts:",samepattern))
    samepatterncnt <- length(samepattern)
    numeratoroflikelihoodratio <- samepatterncnt/numgoldpos
	
    # Check how many gold standard negative pairs have the same
    # pattern of 'yes's & 'no's for the experimental data sets:
    numgoldneg <- length(goldneg) # Number of gold standard negatives
    samepattern <- goldneg
    for (m in 1:numexpdatasets)
    {
 	expm <- exp[[m]]
	isfoundm <- is.element(protpair2,expm) # Says if protpair2 is in the
mth experimental data set
	if (isfoundm == "TRUE")
	{
    	    samepattern <- intersect(samepattern,expm) # Find pairs that
are in goldpos and in expm
        }
	else
	{
	    samepattern <- setdiff(samepattern,expm) # Find pairs that are in
goldpos but not in expm
	}
    }
    #print(paste("Have same pattern in gold neg and expts:",samepattern))
    samepatterncnt <- length(samepattern)
    denominatoroflikelihoodratio <- samepatterncnt/numgoldneg

    # Calculate the likelihood ratio:
    #print(paste("numerator:",numeratoroflikelihoodratio,"denominator:",denominatoroflikelihoodratio))
    likelihoodratio <- numeratoroflikelihoodratio/denominatoroflikelihoodratio
    # Print out the result for this protein-protein interaction:
    if (VALUEONLY == FALSE)
    {
       print(paste("Pair",protpair2,"confidence measure (likelihood
ratio) =",likelihoodratio))
    }
    else
    {
       return(likelihoodratio)
    }
}

# Function to calculate the rate of true positives and false positives, and
# false negatives in an experimental data set of protein-protein interactions:
calcproteinpairaccuracy <- function(goldpos, goldneg, exp)
{
   # Find the number of true interactions in the gold positive data:
   numgoldpos <- length(goldpos)

   # Find the number of true positive protein-protein interactions:
   truepos <- intersect(exp, goldpos)
   numtruepos <- length(truepos)

   # Find the number of false positive protein-protein interactions:
   falsepos <- intersect(exp, goldneg)
   numfalsepos <- length(falsepos)

   # Find the number of false negative protein-protein interactions:
   falseneg <- setdiff(goldpos, exp)
   numfalseneg <- length(falseneg)

   # Find the number of protein-protein interactions that overlap between the
   # experimental data set and gold standard data:
   numoverlap <- numtruepos + numfalsepos
   # Find the rate of true positives, and false positives:
   ratetruepos <- (numtruepos/numoverlap)*100
   ratefalsepos <- (numfalsepos/numoverlap)*100
   # Find the rate of false negatives:
   ratefalseneg <- (numfalseneg/numgoldpos)*100

   print(paste("Number of protein-protein pairs that overlap between
the gold standard and experimental data:",numoverlap))
   print(paste("Number of true positive pairs in the experimental
data:",numtruepos,"(",ratetruepos,"% of",numoverlap,"reported
interactions)"))
   print(paste("Number of false positive pairs in the experimental
data:",numfalsepos,"(",ratefalsepos,"% of",numoverlap,"reported
interactions)"))
   print(paste("Number of false negative pairs in the experimental
data:",numfalseneg,"(",ratefalseneg,"% of",numgoldpos,"true
interactions)"))
}

# Function to read in pairs of interacting proteins in an input file:
readproteinpairs <- function(file)
{
   MyData <- read.table(`file`,header=FALSE)
   mypairs <- paste(MyData$V1,MyData$V2,sep="~")

   return(mypairs)
}

# Use the Needleman-Wunsch algorithm to make an alignment of two strings
# S1 and S2
needlemanwunsch <- function(S1,S2,type="protein",gappenalty=-1)
{
   library(Biostrings)
   library(seqinr)
   S1 <- s2c(S1) # Convert the string to a vector of characters
   S2 <- s2c(S2) # Convert the string to a vector of characters
   lengthS1 <- length(S1)
   lengthS2 <- length(S2)
   # Check if it is a protein alignment, or DNA alignment:
   if (type == "protein")
   {
      # Read the BLOSUM matrix:
      data(BLOSUM45)
      s <- BLOSUM45[S2,S1]
      d <- -gappenalty # The gap penalty is -d
   }
   else
   {
      # Define the DNA scoring matrix:
      s1 <- nucleotideSubstitutionMatrix(match = 1, mismatch = -1,
baseOnly = TRUE)
      s <- s1[S2, S1]
      d <- -gappenalty # The gap penalty is -d
   }

   # Make the table T with lengthS1+1 columns and lengthS2+1 rows:
   # eg. if S1 = VIVADAVIS and S2 = VIVALASVEGAS, then we will
   # have a matrix T with 10 columns and 13 rows
   T <- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)

   # Make a matrix for storing the traceback:
   T2 <- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)

   # Initialise the first column and row of matrix T:
   T[1,] <- -seq(0,lengthS1*d,d) # Fill in the first row
   T[,1] <- -seq(0,lengthS2*d,d) # Fill in the first column
   # In R T[1,2] refers to row 1, column 2
   #      T[2,1] refers to row 2, column 1
   # print(T)

   # Initialise the first column and row of matrix T2:
   T2[1,] <- rep('-',lengthS1+1) # Fill in the first row
   T2[,1] <- rep('|',lengthS2+1) # Fill in the first column
   # print(T2)

   # Carry out the recurrence relation to fill in matrix T:
   for (j in 2:(nrow(T))) # Go through each row (j) at a time
   {
      for (i in 2:(ncol(T))) # Go through each column (i) at a time
      {
         print(paste("Calculating for col i=",i," row j=",j))
         # Set the value of T[i,j] ie. in column i and row j:
         T[j,i] <- max(c(T[j-1,i-1]+s[j-1,i-1],T[j,i-1]-d,T[j-1,i]-d))
         diag <- T[j-1,i-1]+s[j-1,i-1]
         up <- T[j-1,i]-d
         left <- T[j,i-1]-d
         if (diag > up && diag > left)
         {
            T2[j,i] <- ">"
         }
         else if (up > diag && up > left)
         {
            T2[j,i] <- "|"
         }
         else if (left > diag && left > up)
         {
            T2[j,i] <- "-"
         }
         else if (left == diag && up == diag)
         {
            T2[j,i] <- "*"
         }
         else if (up == left && up > diag && left > diag)
         {
            T2[j,i] <- "L"
         }
         else if (up == diag && up > left && diag > left)
         {
            T2[j,i] <- "V"
         }
         else if (diag == left && diag > up && left > up)
         {
            T2[j,i] <- "Z"
         }
         else
         {
            print(paste("ERROR: up",up,"left",left,"diag",diag))
         }
      }
   }
   print(T)
   print(T2)

   # Print out T and T2 together:
   T3 <- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)
   for (j in 2:(nrow(T))) # Go through each row (j) at a time
   {
      for (i in 2:(ncol(T))) # Go through each column (i) at a time
      {
         value <- T[j,i]
         value2 <- T2[j,i]
         value3 <- paste(value,value2)
         T3[j,i] <- value3
      }
   }
   print(T3)
}

# Function to retrieve Uniprot sequences using the SeqinR R library.
retrieveuniprotseqs <- function(seqnames)
{
	myseqs <- list() # Make a list to store the sequences
	library("seqinr")  # Load the SeqinR R library
	for (i in 1:length(seqnames))
	{
		seqname <- seqnames[i]
		print(paste("Retrieving sequence",seqname,"..."))
		choosebank("swissprot")
		queryname <- "query2"
		query <- paste("AC=",seqname,sep="")
		query(`queryname`,`query`)
		seq <- getSequence(query2$req) # Makes a vector "seq" containing the sequence
                if (is.list(seq) == TRUE)
                {
                   seq1b <- seq[[1]]
                   seq <- seq1b
                   if (is.list(seq1b) == TRUE)
                   {
                      seq1c <- seq[[1]]
                      seq <- seq1c
                   }
                }
		myseqs[[i]] <- seq
	}
	return(myseqs)
}

# Function to plot the attractors in the Boolean network
plotBooleanAttractors <- function(attractors)
{
   par(mfrow=c(1,length(attractors$attractors)))
   plotAttractors(attractors)
}

# Function to find network motifs in a gene regulatory network
findnetworkpatterns <- function(graph, patternsize, min=1)
{
   library(igraph)
   # Convert the graph into an "igraph" graph:
   mygraph2 <- igraph.from.graphNEL(graph)
   # Find network motifs in the graph "mygraph2":
   mygraphmotifs <- graph.motifs(mygraph2, patternsize)
   # Find how many motifs occur 1 or more times:
   numdifferentmotifs <- table(mygraphmotifs!=0)[[2]]
   # Make a plot for these motifs:
   #par(mfrow=c(3,4))
   #plot.new() # Make a new plot
   # Plot the 12 most common motifs:
   numplotted <- 0
   oo <- order(mygraphmotifs, decreasing=TRUE)
   # Find which motifs occur:
   nummotifs <- length(mygraphmotifs)
   for (i in 1:nummotifs)
   {
      motifi <- oo[i]
      motif <- mygraphmotifs[motifi] # The number of occurrences of this motif
      if (motif >= min) # There are some occurrences of this motif
      {
         # Find out what the motif looks like:
         motifgraph <- graph.isocreate(size=patternsize,
number=motifi-1, directed=TRUE)
         # Make a plot
         #if (numplotted < 12)
         #{
         #   numplotted <- numplotted + 1
         #   motifgraph2 <- igraph.to.graphNEL(motifgraph)
         #   plot(motifgraph2,attrs=list(node=list(width="2",height="2",fontsize="1"),graph=list(size="5")))
         #}
         edges <- E(motifgraph)
	 print(paste("This pattern occurs",motif,"times:"))
	 print(edges)
      }
   }
}

# Function to make a plot of a subgraph/graph in a regulatory network.
# By default, it makes a plot of the whole regulatory network.
# However, if you pass in a list of vertices in a subgraph, it will
plot just that subgraph.
plotRegulatoryNetwork <- function(graph, nodes=attr(graph,"nodes"),type="neato")
{
   # Make a subgraph of "graph" corresponding to the nodes in vector "nodes"
   mysubgraph <- subGraph(nodes, graph)
   # Copy the attributes of the edges from "graph" to "mysubgraph"
   myrenderinfo <- attr(graph, "renderInfo")
   myrenderinfo2 <- attr(myrenderinfo, "edges")
   myrenderinfo3 <- myrenderinfo2$arrowhead
   # Get all the edges in mysubgraph
   myedges <- listEdges(mysubgraph)
   numedges <- length(myedges)
   myvector <- vector()
   for (i in 1:numedges)
   {
      edge <- myedges[[i]]
      # Sometimes there can be an edge each way between two genes,
      # eg. FIS regulates CRP and CRP regulates FIS.
      numedges2 <- length(edge)
      for (j in 1:numedges2)
      {
         edge2 <- edge[[j]]
         # This seems to be necessary:
         if (j == 1 && numedges2 > 1) { edge2 <- edge2[[1]] }
         beginnode <- slot(edge2,"bNode")
         endnode <- slot(edge2,"eNode")
         edgename <- paste(beginnode,"~",endnode,sep="")
         # Find out the renderinfo for this edge in graph "graph"
         edgetype <- myrenderinfo3[`edgename`]
         edgetype <- edgetype[[1]]
         myvector[`edgename`] <- edgetype
      }
   }
   # Set the type of the edges:
   edgeRenderInfo(mysubgraph) = list(arrowhead=myvector)
   # Make a plot now:
   library("Rgraphviz")
   # Set the type of plot:
   if (type == "hierarchical") { type <- "dot" }
   mygraphplot <- layoutGraph(mysubgraph, layoutType=type)
   renderGraph(mygraphplot)	
}

# Function to get all the descendents of a particular node in a
directed graphNEL graph:
getdescendents <- function(graph, node)
{
   # First get all the nodes that are children of node "node":
   children <- adj(graph, node)[[1]]
   numchildren <- length(children)
   descendents <- vector() # Make a vector for storing all descendents
   # Record whether or not we have seen any of the children before
   allnodes <- nodes(graph)
   seen <- vector()
   for (i in 1:length(allnodes))
   {
      node <- allnodes[i]
      seen[`node`] <- 0
   }

   # For each of the children, check if it has its own children
   loopagain <- 1
   while(loopagain == 1 && numchildren > 0)
   {
      newchildren <- vector()
      foundnew <- 0
      # Record whether or not we have seen any of the children before
      for (i in 1:numchildren)
      {
         child <- children[i]
         seenchild <- seen[`child`]
         if (seenchild == 0)
         {
            seen[`child`] <- 1
            descendents <- append(descendents, child, after=length(descendents))
         }
      }
      for (i in 1:numchildren)
      {
         child <- children[i]
         grandchildren <- adj(graph, child)[[1]]
         numgrandchildren <- length(grandchildren)
	 if (numgrandchildren > 0 && child != node)
	 {
	    for (j in 1:numgrandchildren)
	    {
               grandchild <- grandchildren[j]
	       seengrandchild <- seen[`grandchild`]
	       if (grandchild != node && seengrandchild == 0)
	       {
	          newchildren <- append(newchildren, grandchild,
after=length(newchildren))
		  foundnew <- 1
	       }
            }
	 }
      }
      if (foundnew == 0)
      {
         loopagain <- 0
      }
      else
      {
         children <- newchildren
	 numchildren <- length(children)
      }
   }
   return(descendents)
}

# Function to retrieve GenBank sequences using the SeqinR R library.
retrievegenbankseqs <- function(seqnames, from="NA", to="NA")
{
	myseqs <- list()
	library("seqinr")
	for (i in 1:length(seqnames))
	{
		seqname <- seqnames[i]
		print(paste("Retrieving sequence",seqname,"..."))
                if (substring(seqname,1,3) == 'XM_' ||
substring(seqname,1,3) == 'NM_')
                {
                   choosebank("refseq")
                }
                else
                {
		   choosebank("genbank")
                }
		queryname <- "query2"
		query <- paste("AC=",seqname,sep="")
		query(`queryname`,`query`)
		seq <- getSequence(query2$req) # Makes a vector "seq" containing the sequence
                if (is.list(seq) == TRUE)
                {
                   seq1b <- seq[[1]]
                   seq <- seq1b
                   if (is.list(seq1b) == TRUE)
                   {
                      seq1c <- seq[[1]]
                      seq <- seq1c
                   }
                }
                if (from != 'NA' && to != 'NA')
                {
                   seq <- seq[from:to]
                }
	        myseqs[[i]] <- seq
	}
	return(myseqs)
}

# Function to retrieve GenBank virus sequences using the SeqinR R library.
retrievevirusseqs <- function(seqnames)
{
	myseqs <- list()
	library("seqinr")
	for (i in 1:length(seqnames))
	{
		seqname <- seqnames[i]
		print(paste("Retrieving sequence",seqname,"..."))
		choosebank("refseqViruses")
		queryname <- "query2"
		query <- paste("AC=",seqname,sep="")
		query(`queryname`,`query`)
		seq <- getSequence(query2$req) # Makes a vector "seq" containing the sequence
                if (is.list(seq) == TRUE)
                {
                   seq1b <- seq[[1]]
                   seq <- seq1b
                   if (is.list(seq1b) == TRUE)
                   {
                      seq1c <- seq[[1]]
                      seq <- seq1c
                   }
                }
		myseqs[[i]] <- seq
	}
	return(myseqs)
}


# Function to convert a graphNEL graph to an igraph graph.
# Copied 6-Feb-10 from:
# http://bazaar.launchpad.net/~igraph/igraph/0.6-main/annotate/head:/interfaces/R/igraph/R/conversion.R
igraph.from.graphNEL <- function(graphNEL, name=TRUE, weight=TRUE,
                                 unlist.attrs=TRUE)
{
  require(graph)
  if (!is(graphNEL, "graphNEL")) {
    stop("Not a graphNEL graph")
  }

  al <- lapply(edgeL(graphNEL), "[[", "edges")

  if (edgemode(graphNEL)=="undirected") {
    al <- mapply(SIMPLIFY=FALSE, seq_along(al), al, FUN=function(n, l) {
      c(l, rep(n, sum(l==n)))
    })
  }

  al <- lapply(al, function(x) x-1)

  mode <- if (edgemode(graphNEL)=="directed") "out" else "all"

  g <- graph.adjlist(al, directed=TRUE, duplicate=TRUE)
  if (name) {
    V(g)$name <- nodes(graphNEL)
  }

  ## Graph attributes

  g.n <- names(graphNEL@graphData)
  g.n <- g.n [ g.n != "edgemode" ]
  for (n in g.n) {
    g <- set.graph.attribute(g, n, graphNEL@graphData[[n]])
  }

  ## Vertex attributes
  v.n <- names(nodeDataDefaults(graphNEL))

  for (n in v.n) {
    val <- unname(nodeData(graphNEL, attr=n))
    if (unlist.attrs && all(sapply(val, length)==1)) { val <- unlist(val) }
    g <- set.vertex.attribute(g, n, value=val)
  }

  ## Edge attributes
  e.n <- names(edgeDataDefaults(graphNEL))
  if (!weight) { e.n <- e.n [ e.n != "weight" ] }

  if (length(e.n) > 0) {
    el <- get.edgelist(g)
    el <- paste(sep="|", el[,1], el[,2])
    for (n in e.n) {
      val <- unname(edgeData(graphNEL, attr=n)[el])
      if (unlist.attrs && all(sapply(val, length)==1)) { val <- unlist(val) }
      g <- set.edge.attribute(g, n, value=val)
    }
  }
  g
}

# Function to change an igraph graph into a graphNEL format graph.
# Copied from the development code for igraph
# http://bazaar.launchpad.net/~igraph/igraph/0.6-main/annotate/head:/interfaces/R/igraph/R/conversion.R
# on 26-Jan-2010
igraph.to.graphNEL <- function(graph) {
  if (!is.igraph(graph))
  {
    stop("Not an igraph graph")
  }

  require(graph)

  if ("name" %in% list.vertex.attributes(graph) &&
     is.character(V(graph)$name)) {
    name <- V(graph)$name
  } else {
    name <- as.character(seq(vcount(graph))-1)
  }

  edgemode <- if (is.directed(graph)) "directed" else "undirected"

  if ("weight" %in% list.edge.attributes(graph) &&
      is.numeric(E(graph)$weight)) {
    al <- get.adjedgelist(graph, "out")
    for (i in seq(along=al)) {
      edges <- get.edges(graph, al[[i]])
      edges <- ifelse( edges[,2]==i-1, edges[,1], edges[,2])
      weights <- E(graph)$weight[al[[i]]+1]
      al[[i]] <- list(edges=edges+1, weights=weights)
    }
  } else {
    al <- get.adjlist(graph, "out")
    al <- lapply(al, function(x) list(edges=x+1))
  }

  names(al) <- name

  res <- new("graphNEL", nodes=name, edgeL=al, edgemode=edgemode)

  ## Add graph attributes (other than 'directed')

  ## Are this "officially" supported at all?

  g.n <- list.graph.attributes(graph)

  if ("directed" %in% g.n) {
    warning("Cannot add graph attribute `directed'")
    g.n <- g.n[ g.n != "directed" ]
  }
  for (n in g.n) {
    res@graphData[[n]] <- get.graph.attribute(graph, n)
  }

  ## Add vertex attributes (other than 'name', that is already
  ## added as vertex names)

  v.n <- list.vertex.attributes(graph)

  v.n <- v.n[ v.n != "name" ]

  for (n in v.n) {
    nodeDataDefaults(res, attr=n) <- NA
    nodeData(res, attr=n) <- get.vertex.attribute(graph, n)
  }

  ## Add edge attributes (other than 'weight')

  e.n <- list.edge.attributes(graph)
  e.n <- e.n[ e.n != "weight" ]

  if (length(e.n) > 0) {
    el <- get.edgelist(graph)
    el <- paste(sep="|", el[,1], el[,2])
    for (n in e.n) {
      edgeDataDefaults(res, attr=n) <- NA
      res@edgeData@data[el] <- mapply(function(x,y) {
        xx <- c(x,y); names(xx)[length(xx)] <- n; xx },
                                      res@edgeData@data[el],
                                      get.edge.attribute(graph, n),
                                      SIMPLIFY=FALSE)
    }
  }

  res

}

# Function to generate X random sequences with a multinomial model, where the
# probabilities of the different letters are set equal to their frequencies
# in an input sequence "inputsequence". The input sequence comes in as a string
# of characters, eg. "ATSTWCWYSKLAV"
generateSeqsWithMultinomialModel <- function(inputsequence, X)
{
   # Change the input sequence into a vector of letters
   library("seqinr") # Load in the SeqinR library, so we can use function "s2c".
   inputsequencevector <- s2c(inputsequence)
   # Find the frequencies of the letters in the input sequence
"inputsequencevector":
   mylength <- length(inputsequencevector)
   mytable <- table(inputsequencevector)
   # Find the names of the letters in the sequence
   letters <- rownames(mytable)
   numletters <- length(letters)
   probabilities <- numeric() # Make a vector to store the
probabilities of letters
   for (i in 1:numletters)
   {
      letter <- letters[i]
      count <- mytable[[i]]
      probabilities[i] <- count/mylength
   }
   # Make X random sequences using the multinomial model with
probabilities "probabilities"
   seqs <- numeric(X)
   for (j in 1:X)
   {
      seq <- sample(letters, mylength, rep=TRUE, prob=probabilities) #
Sample with replacement
      seq <- c2s(seq)
      seqs[j] <- seq
   }
   # Return the vector of random sequences
   return(seqs)
}

# Function to convert a phylip format tree to an 'ape' format tree.
# In phylip, the bootstrap values are given for branches.
# In 'ape', it expects the bootstrap values to be given for nodes.
phylip.to.ape <- function(tree)
{
   # Get the labels of the tips (sequences)
   seqnames <- tree$tip.label
   numseqs <- length(seqnames)
   # Get the number of internal nodes in the tree:
   numinternalnodes <- tree$Nnode
   # Get the bootstrap values
   bootstraps <- tree$edge.length
   # Get the edges for the tree
   edges <- tree$edge
   # For each edge, find the node at the end of the edge (closer to the tips
   # of the tree):
   numedges <- nrow(edges)
   numnodes <- numinternalnodes + numseqs
   mybootstraps <- numeric(numinternalnodes)
   for (i in 1:numedges)
   {
      startnode <- edges[i,1]
      endnode <- edges[i,2]
      # Make sure the end node is not a tip:
      if (endnode > numseqs)
      {
         # Get the bootstrap value for this node:
         bootstrap <- bootstraps[i]
         mybootstraps[endnode] = bootstrap
      }
   }
   mybootstraps2 <- mybootstraps[(numseqs+1):(numnodes)]
   tree$node.label <- mybootstraps2

   return(tree)
}

plotproteinNJtreewithbootstraps <- function(alignment, theoutgroup)
{
   # define a function for making a tree:
   makemytree <- function(alignmentmat, outgroup=`theoutgroup`)
   {
      alignment <- as.alignment(alignmentmat)
      mydist <- dist.alignment(alignment)
      mytree <- nj(mydist)
      mytree <- makeLabel(mytree, space="") # get rid of spaces in tip names.
      myrootedtree <- root(mytree, outgroup, r=TRUE)
      return(myrootedtree)
   }
   # infer a tree
   mymat  <- as.matrix.alignment(alignment)
   myrootedtree <- makemytree(mymat, outgroup=theoutgroup)
   # bootstrap the tree
   mymat  <- as.matrix.alignment(alignment)
   myboot <- boot.phylo(myrootedtree, mymat, makemytree)
   # plot the tree
   plot.phylo(myrootedtree)
   nodelabels(myboot)
}

# Function to plot the range of fluxes, and optimal fluxes for reactions.
# model is the input model, eg. LIMEcoli
# Code copied from
http://cran.r-project.org/web/packages/LIM/vignettes/LIMecoli.pdf
plotFluxes <- function(model)
{
   library("LIM")
   LP <- Linp(model)
   xr <- Xranges(model)
   par(mfrow=c(1,2))
   nr <- model$NUnknowns
   ii <- 1:(nr/2)
   dotchart(LP$X[ii],xlim = range(xr),pch=16,cex=0.8)
   segments(xr[ii,1],1:nr,xr[ii,2],1:nr)
   ii <- (nr/2+1):nr
   dotchart(LP$X[ii],xlim = range(xr),pch=16,cex=0.8)
   segments(xr[ii,1],1:nr,xr[ii,2],1:nr)
   mtext(side= 3, cex=1.5, outer = TRUE, line=-1.5, "E coli Core
Metabolism, optimal solution and ranges")
}


# Function to identify periodically expressed genes in microarray time
series data.
# Code inspired by
http://www.cs.tut.fi/~ahdesmak/robustperiodic/doc/caulobacter.R, using
# the GeneCycle R library:
# The input data is a matrix with one row per gene and one column per
microarray.
findPeriodicGenes <- function(expndata)
{
   library("GeneCycle")
   # Put the data in the format required by fdrtool() and fisher.g.test(),
   # ie. one row per microarray and one column per gene,
   # and no missing values:
   expndata2 <- na.omit(expndata) # Discard genes that have missing
data ("NA" values)
   expndata3 <- t(expndata2)      # Find the transpose of matrix "expndata2"

   # Calculate p-values using Fisher's g test
   pval.mgenes <- fisher.g.test(expndata3)

   # Calculate q-values (ie. p-values that have been corrected for
multiple testing, as we
   # are carrying out tests for lots of different genes at once):
   fdr.out <- fdrtool(pval.mgenes, statistic="pvalue")

   # Print out the qvalues for genes that have qvalues of < 0.05:
   numgenes <- length(fdr.out$qval)
   genenames <- numeric()
   qvalues <- numeric()
   numgenestaken <- 0
   for (i in 1:numgenes)
   {
      qval <- fdr.out$qval[i]
      genename <- colnames(expndata3)[i]
      if (qval <= 0.05)
      {
         numgenestaken <- numgenestaken + 1
         genenames[numgenestaken] <- genename
         qvalues[numgenestaken] <- qval
      }
   }
   # Return a list variable that contains the names of periodically
expressed genes,
   # and the q-values:
   result <- list(genenames,qvalues)

   return(result)
}

# An R function to plot a heatmap for microarray data.
# The input data is the distance matrix "dist".
# Inspired by code in Hahne et al "Bioconductor Case Studies" page 140.
plotHeatmap <- function(dist)
{
   library("RColorBrewer")
   # Use red to correspond to high values and blue to low:
   hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
   hmcol <- rev(hmcol)
   # This uses hierarchical clustering:
   heatmap(as.matrix(dist), sym=TRUE, col=hmcol, distfun=function(x)
dist(x, method="manhattan"))
}

# An R function to plot the expression for all of the genes in a cluster.
# "clusters" is the result of k-means clustering, and num is the number of
# the cluster that we want to plot.
plotCluster <- function(clusters, num, expndata)
{
   # Find all the genes in this cluster:
   genes <- clusters$cluster[clusters$cluster == num]
   genes <- names(genes)
   numgenes <- length(genes)

   # Now make a plot with the expression of each gene across the samples:
   for (i in 1:numgenes)
   {
      gene <- genes[i]
      genedata <- expndata[`gene`,]
      if (i == 1)
      {
         ylab <- paste("Cluster",num)
         plot(genedata, type="l", col="blue", ylab=ylab)
      }
      else
      {
         points(genedata, type="l", col="blue")
      }
   }
}

# An R function that calculates the transition matrix for a Markov model,
# taking the matrix of observed dinucleotide relative frequencies as its
# input. Each row of the input matrix specifies a base, and each column
# specifies the following base:
# eg.
# D <- matrix(c(0.146, 0.052, 0.058, 0.089,
#               0.063, 0.029, 0.010, 0.056,
#               0.050, 0.030, 0.028, 0.051,
#               0.087, 0.047, 0.063, 0.140), byrow=TRUE, nrow=4)
makeTransitionMatrix <- function(D)
{
   # Initialise the transition matrix:
   P <- array(0,dim=c(4,4))

   # Calculate the transition matrix:
   for (i in 1:4) # For each row (i) of the transition matrix:
   {
      for (j in 1:4) # For each column (j) of the transition matrix:
      {
         value <- D[i,j]/sum(D[i,]) # sum(D[i,]) is the sum of the ith
row of the input matrix
         P[i,j] <- value
      }
   }

   return(P)
}

--047d7b86ef4c04daea04c84f04ed
Content-Type: text/html; charset=ISO-8859-1

<pre># Rfunctions.R
# Written by Avril Coghlan, <a href="mailto:a.coghlan@ucc.ie">a.coghlan@ucc.ie</a>

dnasmithwaterman &lt;- function(S1, S2, gapopen, gapextend, mymatch, mymismatch)
{
           library(Biostrings)
           library(seqinr)
           S1 &lt;- s2c(S1) # Convert the string to a vector of characters
           S2 &lt;- s2c(S2) # Convert the string to a vector of characters
           lengthS1 &lt;- length(S1)
           lengthS2 &lt;- length(S2)
           # Define the DNA scoring matrix:
           s1 &lt;- nucleotideSubstitutionMatrix(match = mymatch, mismatch = mymismatch, baseOnly = TRUE) 
           # match is usually +ve, mismatch -ve
           s &lt;- s1[S2, S1] 
           d &lt;- gapopen # The gap open penalty is normally -ve
           d2 &lt;- gapextend # The gap extension penalty is normally -ve

           # Make the table T with lengthS1+1 columns and lengthS2+1 rows:
           # eg. if S1 = VIVADAVIS and S2 = VIVALASVEGAS, then we will
           # have a matrix T with 10 columns and 13 rows
           T &lt;- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)
           mycolnames &lt;- c(&quot;NA&quot;,S1)
           myrownames &lt;- c(&quot;NA&quot;,S2)
           colnames(T) &lt;- mycolnames
           rownames(T) &lt;- myrownames

           # Make a matrix for storing the traceback:
           T2 &lt;- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)
           colnames(T2) &lt;- mycolnames
           rownames(T2) &lt;- myrownames

           # Initialise the first column and row of matrix T:
           rowlength &lt;- length(T[1,])
           for (k in 1:rowlength)
           {
              if          (k == 1) { T[1,k] &lt;- 0                    }
              else if (k == 2)     { T[1,k] &lt;- 0                    } 
              else                 { T[1,k] &lt;- 0                    }
           }
           collength &lt;- length(T[,1])
           for (k in 1:collength)
           {
              if          (k == 1) { T[k,1] &lt;- 0                    }
              else if (k == 2)     { T[k,1] &lt;- 0                    } 
              else                 { T[k,1] &lt;- 0                    }
           }
           
           # Initialise the first column and row of matrix T2:
           T2[1,] &lt;- rep(&#39;-&#39;,lengthS1+1) # Fill in the first row
           T2[,1] &lt;- rep(&#39;|&#39;,lengthS2+1) # Fill in the first column
           # print(T2)

           # Carry out the recurrence relation to fill in matrix T:
           for (j in 2:(nrow(T))) # Go through each row (j) at a time
           {
              for (i in 2:(ncol(T))) # Go through each column (i) at a time
              {
                 print(paste(&quot;Calculating for col i=&quot;,i,&quot; row j=&quot;,j))
                 # Set the value of T[i,j] ie. in column i and row j:
                 diag &lt;- T[j-1,i-1]+s[j-1,i-1]
                 if (T2[j-1,i] == &quot;|&quot;)      { up &lt;- T[j-1,i] + d2   } 
                 else                       { up &lt;- T[j-1,i] + d + d2    }
                 # Still use extension penalty for first position
                 if (T2[j,i-1] == &quot;-&quot;)      { left &lt;- T[j,i-1] + d2  }
                 else                       { left &lt;- T[j,i-1] + d  + d2 }
                 # Still use extension penalty for first position
                 T[j,i] &lt;- max(c(diag,up,left,0))
                 if (diag &lt; 0 &amp;&amp; up &lt; 0 &amp;&amp; left &lt; 0) 
                 {
                    T2[j,i] &lt;- &quot;+&quot;
                 }
                 else
                 {
                    if (diag &gt; up &amp;&amp; diag &gt; left)
                    {
                       T2[j,i] &lt;- &quot;&gt;&quot;
                    }
                    else if (up &gt; diag &amp;&amp; up &gt; left)
                    {
                       T2[j,i] &lt;- &quot;|&quot;
                    }
                    else if (left &gt; diag &amp;&amp; left &gt; up)
                    {
                        T2[j,i] &lt;- &quot;-&quot;
                    }
                    else if (left == diag &amp;&amp; up == diag)
                    {
                        T2[j,i] &lt;- &quot;*&quot;
                    }
                    else if (up == left &amp;&amp; up &gt; diag &amp;&amp; left &gt; diag)
                    {
                        T2[j,i] &lt;- &quot;L&quot;
                    }
                    else if (up == diag &amp;&amp; up &gt; left &amp;&amp; diag &gt; left)
                    {
                       T2[j,i] &lt;- &quot;V&quot;
                    }
                   else if (diag == left &amp;&amp; diag &gt; up &amp;&amp; left &gt; up)
                   {
                      T2[j,i] &lt;- &quot;Z&quot;
                   }
                   else
                   {
                      print(paste(&quot;ERROR: up&quot;,up,&quot;left&quot;,left,&quot;diag&quot;,diag))
                   }
                }
              }
           }
          print(T) 
          print(T2)
          maxT &lt;- max(T)
          print(paste(&quot;maxT=&quot;,maxT))

          # Print out T and T2 together:
          T3 &lt;- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)
          colnames(T3) &lt;- mycolnames
           rownames(T3) &lt;- myrownames

          for (j in 2:(nrow(T))) # Go through each row (j) at a time
          {
              for (i in 2:(ncol(T))) # Go through each column (i) at a time
              {
                 value &lt;- T[j,i]
                 value2 &lt;- T2[j,i]
                 value3 &lt;- paste(value,value2)
                 T3[j,i] &lt;- value3
              }
           }
           print(T3)
}

printPairwiseAlignment &lt;- function(alignment, chunksize=60, returnlist=FALSE)
{
     library(Biostrings)
     seq1aln &lt;- pattern(alignment) # Get the alignment for the first sequence
     seq2aln &lt;- subject(alignment) # Get the alignment for the second sequence
     alnlen  &lt;- nchar(seq1aln)     # Find the number of columns in the alignment
     starts  &lt;- seq(1, alnlen, by=chunksize)
     n       &lt;- length(starts)     
     seq1alnresidues &lt;- 0
     seq2alnresidues &lt;- 0
     for (i in 1:n) {
        chunkseq1aln &lt;- substring(seq1aln, starts[i], starts[i]+chunksize-1)
        chunkseq2aln &lt;- substring(seq2aln, starts[i], starts[i]+chunksize-1)
        # Find out how many gaps there are in chunkseq1aln:
        gaps1 &lt;- countPattern(&quot;-&quot;,chunkseq1aln) # countPattern() is from Biostrings library
        # Find out how many gaps there are in chunkseq2aln:                  
        gaps2 &lt;- countPattern(&quot;-&quot;,chunkseq2aln) # countPattern() is from Biostrings library
        # Calculate how many residues of the first sequence we have printed so far in the alignment:
        seq1alnresidues &lt;- seq1alnresidues + chunksize - gaps1
        # Calculate how many residues of the second sequence we have printed so far in the alignment:
        seq2alnresidues &lt;- seq2alnresidues + chunksize - gaps2
        if (returnlist == &#39;FALSE&#39;)
        {
           print(paste(chunkseq1aln,seq1alnresidues))
           print(paste(chunkseq2aln,seq2alnresidues))
           print(paste(&#39; &#39;))
        }
     }
     if (returnlist == &#39;TRUE&#39;)
     {
        vector1 &lt;- s2c(substring(seq1aln, 1, nchar(seq1aln)))
        vector2 &lt;- s2c(substring(seq2aln, 1, nchar(seq2aln)))
        mylist &lt;- list(vector1, vector2) 
        return(mylist)
     }
}

# This finds all ORFs in the forward strand of a sequence. If two ORFs overlap and share a start/stop codon,
# it takes the longer one. Otherwise, overlapping ORFs are kept.
findORFsinSeq &lt;- function(sequence)
{
   library(Biostrings)
   # Make vectors &quot;positions&quot; and &quot;types&quot; containing information on the positions of ATGs in the sequence:
   mylist           &lt;- findPotentialStartsAndStops(sequence)
   positions        &lt;- mylist[[1]]
   types            &lt;- mylist[[2]]
   # Make vectors &quot;orfstarts&quot; and &quot;orfstops&quot; to store the predicted start and stop codons of ORFs
   orfstarts        &lt;- numeric()
   orfstops         &lt;- numeric() 
   # Make a vector &quot;orflengths&quot; to store the lengths of the ORFs
   orflengths       &lt;- numeric()
   # Print out the positions of ORFs in the sequence:
   numpositions     &lt;- length(positions) # Find the length of vector &quot;positions&quot;
   if (numpositions &gt;= 2)                # There must be at least one start codon and one stop codon to have an ORF.
   {
      for (i in 1:(numpositions-1))
      {
         posi        &lt;- positions[i]
         typei       &lt;- types[i]
         found       &lt;- 0
         while (found == 0)
         {
            for (j in (i+1):numpositions)
            {
               posj  &lt;- positions[j]
               typej &lt;- types[j]
               posdiff &lt;- posj - posi
               posdiffmod3 &lt;- posdiff %% 3
               orflength &lt;- posj - posi + 3 # Add in the length of the stop codon
               if (typei == &quot;atg&quot; &amp;&amp; (typej == &quot;taa&quot; || typej == &quot;tag&quot; || typej == &quot;tga&quot;) &amp;&amp; posdiffmod3 == 0) 
               {
                  # Check if we have already used the stop codon at posj+2 in an ORF
                  numorfs &lt;- length(orfstops)
                  usedstop &lt;- -1    
                  if (numorfs &gt; 0)
                  {
                     for (k in 1:numorfs)
                     {
                        orfstopk &lt;- orfstops[k] 
                        if (orfstopk == (posj + 2)) { usedstop &lt;- 1 }
                     }
                  }
                  if (usedstop == -1)
                  {
		     orfstarts &lt;- append(orfstarts, posi, after=length(orfstarts))
                     orfstops &lt;- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
                     orflengths &lt;- append(orflengths, orflength, after=length(orflengths))
                  }
                  found &lt;- 1
                  break
               }
	       if (j == numpositions) { found &lt;- 1 }
            }
         }
      }
   }

   # Sort the final ORFs by start position:
   indices           &lt;- order(orfstarts)
   orfstarts         &lt;- orfstarts[indices] 
   orfstops          &lt;- orfstops[indices]

   # Find the lengths of the ORFs that we have
   orflengths        &lt;- numeric()
   numorfs           &lt;- length(orfstarts)
   for (i in 1:numorfs)
   {
      orfstart       &lt;- orfstarts[i]
      orfstop        &lt;- orfstops[i]
      orflength      &lt;- orfstop - orfstart + 1
      orflengths     &lt;- append(orflengths,orflength,after=length(orflengths))
   }

   mylist            &lt;- list(orfstarts, orfstops, orflengths)
   return(mylist)                    
}

findPotentialStartsAndStops &lt;- function(sequence)
{
   codons            &lt;- c(&quot;atg&quot;, &quot;taa&quot;, &quot;tag&quot;, &quot;tga&quot;) # Define a vector with the sequences of potential start and stop codons
   # Find the number of occurrences of each type of potential start or stop codon
   for (i in 1:4)                                     
   {
      codon          &lt;- codons[i]
      occurrences    &lt;- matchPattern(codon, sequence) # Find all occurrences of codon &quot;codon&quot; in sequence &quot;sequence&quot; 
      codonpositions &lt;- attr(occurrences,&quot;start&quot;)     # Find the start positions of all occurrences of &quot;codon&quot; in sequence &quot;sequence&quot;
      numoccurrences &lt;- length(codonpositions)        # Find the total number of potential start and stop codons in sequence &quot;sequence&quot;
      if (i == 1) 
      {
         positions   &lt;- codonpositions                # Make a copy of vector &quot;codonpositions&quot; called &quot;positions&quot;
         types       &lt;- rep(codon, numoccurrences)    # Make a vector &quot;types&quot; containing &quot;numoccurrences&quot; copies of &quot;codon&quot;   
      } 
      else
      {
         # Add the vector &quot;codonpositions&quot; to the end of vector &quot;positions&quot;:
         positions   &lt;- append(positions, codonpositions, after=length(positions))
         # Add the vector &quot;rep(codon, numoccurrences)&quot; to the end of vector &quot;types&quot;:
         types       &lt;- append(types, rep(codon, numoccurrences), after=length(types))
      }
   }
   # Sort the vectors &quot;positions&quot; and &quot;types&quot; in order of position along the input sequence:
   indices           &lt;- order(positions)
   positions         &lt;- positions[indices]
   types             &lt;- types[indices] 
   # Return a list variable including vectors &quot;positions&quot; and &quot;types&quot;:
   mylist            &lt;- list(positions,types)
   return(mylist)
}

# Make a plot of the positions of the potential start and stop codons in
# each of the three reading frames in a sequence
plotPotentialStartsAndStops &lt;- function(sequence)
{
   codons            &lt;- c(&quot;atg&quot;, &quot;taa&quot;, &quot;tag&quot;, &quot;tga&quot;) # Define a vector with the sequences of potential start and stop codons
   # Find the number of occurrences of each type of potential start or stop codon
   for (i in 1:4)                                     
   {
      codon          &lt;- codons[i]
      occurrences    &lt;- matchPattern(codon, sequence) # Find all occurrences of codon &quot;codon&quot; in sequence &quot;sequence&quot; 
      codonpositions &lt;- attr(occurrences,&quot;start&quot;)     # Find the start positions of all occurrences of &quot;codon&quot; in sequence &quot;sequence&quot;
      numoccurrences &lt;- length(codonpositions)        # Find the total number of potential start and stop codons in sequence &quot;sequence&quot;
      if (i == 1) 
      {
         positions   &lt;- codonpositions                # Make a copy of vector &quot;codonpositions&quot; called &quot;positions&quot;
         types       &lt;- rep(codon, numoccurrences)    # Make a vector &quot;types&quot; containing &quot;numoccurrences&quot; copies of &quot;codon&quot;   
      } 
      else
      {
         # Add the vector &quot;codonpositions&quot; to the end of vector &quot;positions&quot;:
         positions   &lt;- append(positions, codonpositions, after=length(positions))
         # Add the vector &quot;rep(codon, numoccurrences)&quot; to the end of vector &quot;types&quot;:
         types       &lt;- append(types, rep(codon, numoccurrences), after=length(types))
      }
   }
   # Sort the vectors &quot;positions&quot; and &quot;types&quot; in order of position along the input sequence:
   indices           &lt;- order(positions)
   positions         &lt;- positions[indices]
   types             &lt;- types[indices] 
   # Make a plot showing the positions of the start and stop codons in the input sequence:
   # Draw a line at y=0 from 1 to the length of the sequence:
   x                 &lt;- c(1,nchar(sequence))
   y                 &lt;- c(0,0)
   plot(x, y, ylim=c(0,3), type=&quot;l&quot;, axes=FALSE, xlab=&quot;Nucleotide&quot;, ylab=&quot;Reading frame&quot;, main=&quot;Predicted start (red) and stop (blue) codons&quot;)
   segments(1,1,nchar(sequence),1)
   segments(1,2,nchar(sequence),2)
   # Add the x-axis at y=0: 
   axis(1, pos=0) 
   # Add the y-axis labels:
   text(0.9,0.5,&quot;+1&quot;)
   text(0.9,1.5,&quot;+2&quot;)
   text(0.9,2.5,&quot;+3&quot;)
   # Draw in each predicted start/stop codon:
   numcodons         &lt;- length(positions)
   for (i in 1:numcodons)
   {
      position       &lt;- positions[i]
      type           &lt;- types[i]
      remainder      &lt;- (position-1) %% 3
      if    (remainder == 0) # +1 reading frame
      {
         if (type == &quot;atg&quot;) { segments(position,0,position,1,lwd=1,col=&quot;red&quot;) }
         else               { segments(position,0,position,1,lwd=1,col=&quot;blue&quot;)}
      }
      else if (remainder == 1)
      {
         if (type == &quot;atg&quot;) { segments(position,1,position,2,lwd=1,col=&quot;red&quot;) }
         else               { segments(position,1,position,2,lwd=1,col=&quot;blue&quot;)}
      }
      else if (remainder == 2)
      { 
         if (type == &quot;atg&quot;) { segments(position,2,position,3,lwd=1,col=&quot;red&quot;) }
         else               { segments(position,2,position,3,lwd=1,col=&quot;blue&quot;)}
      }
   }
}

# This plots all ORFs in the forward strand of a sequence. If two ORFs overlap and share a start/stop codon,
# it takes the longer one. Otherwise, overlapping ORFs are kept.
plotORFsinSeq &lt;- function(sequence)
{
   # Make vectors &quot;positions&quot; and &quot;types&quot; containing information on the positions of ATGs in the sequence:
   mylist           &lt;- findPotentialStartsAndStops(sequence)
   positions        &lt;- mylist[[1]]
   types            &lt;- mylist[[2]]
   # Make vectors &quot;orfstarts&quot; and &quot;orfstops&quot; to store the predicted start and stop codons of ORFs
   orfstarts        &lt;- numeric()
   orfstops         &lt;- numeric() 
   # Make a vector &quot;orflengths&quot; to store the lengths of the ORFs
   orflengths       &lt;- numeric()
   # Print out the positions of ORFs in the sequence:
   numpositions     &lt;- length(positions) # Find the length of vector &quot;positions&quot;
   if (numpositions &gt;= 2)                # There must be at least one start codon and one stop codon to have an ORF.
   {
      for (i in 1:(numpositions-1))
      {
         posi        &lt;- positions[i]
         typei       &lt;- types[i]
         found       &lt;- 0
         while (found == 0)
         {
            for (j in (i+1):numpositions)
            {
               posj  &lt;- positions[j]
               typej &lt;- types[j]
               posdiff &lt;- posj - posi
               posdiffmod3 &lt;- posdiff %% 3
               orflength &lt;- posj - posi + 3 # Add in the length of the stop codon
               if (typei == &quot;atg&quot; &amp;&amp; (typej == &quot;taa&quot; || typej == &quot;tag&quot; || typej == &quot;tga&quot;) &amp;&amp; posdiffmod3 == 0) 
               {
                  # Check if we have already used the stop codon at posj+2 in an ORF
                  numorfs &lt;- length(orfstops)
                  usedstop &lt;- -1    
                  if (numorfs &gt; 0)
                  {
                     for (k in 1:numorfs)
                     {
                        orfstopk &lt;- orfstops[k] 
                        if (orfstopk == (posj + 2)) { usedstop &lt;- 1 }
                     }
                  }
                  if (usedstop == -1)
                  {
		     orfstarts &lt;- append(orfstarts, posi, after=length(orfstarts))
                     orfstops &lt;- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
                     orflengths &lt;- append(orflengths, orflength, after=length(orflengths))
                  }
                  found &lt;- 1
                  break
               }
	       if (j == numpositions) { found &lt;- 1 }
            }
         }
      }
   }

   # Sort the final ORFs by start position:
   indices           &lt;- order(orfstarts)
   orfstarts         &lt;- orfstarts[indices] 
   orfstops          &lt;- orfstops[indices]

   # Make a plot showing the positions of ORFs in the input sequence:
   # Draw a line at y=0 from 1 to the length of the sequence:
   x                 &lt;- c(1,nchar(sequence))
   y                 &lt;- c(0,0)
   plot(x, y, ylim=c(0,3), type=&quot;l&quot;, axes=FALSE, xlab=&quot;Nucleotide&quot;, ylab=&quot;Reading frame&quot;, main=&quot;Predicted ORFs&quot;)
   segments(1,1,nchar(sequence),1)
   segments(1,2,nchar(sequence),2)
   # Add the x-axis at y=0: 
   axis(1, pos=0) 
   # Add the y-axis labels:
   text(0.9,0.5,&quot;+1&quot;)
   text(0.9,1.5,&quot;+2&quot;)
   text(0.9,2.5,&quot;+3&quot;)
   # Make a plot of the ORFs in the sequence:
   numorfs           &lt;- length(orfstarts)
   for (i in 1:numorfs)
   {
      orfstart       &lt;- orfstarts[i]
      orfstop        &lt;- orfstops[i]
      remainder      &lt;- (orfstart-1) %% 3
      if    (remainder == 0) # +1 reading frame
      {
         rect(orfstart,0,orfstop,1,col=&quot;cyan&quot;,border=&quot;black&quot;)
      }
      else if (remainder == 1)
      {
         rect(orfstart,1,orfstop,2,col=&quot;cyan&quot;,border=&quot;black&quot;)
      }
      else if (remainder == 2)
      { 
         rect(orfstart,2,orfstop,3,col=&quot;cyan&quot;,border=&quot;black&quot;)
      }
   }

}

# This carries out the Viterbi algorithm.
# Adapted from &quot;Applied Statistics for Bioinformatics using R&quot; by Wim P. Krijnen, page 209
# ( <a href="http://cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf">cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf</a> ) 
viterbi &lt;- function(sequence, transitionmatrix, emissionmatrix)
{
   # Get the names of the states in the HMM:
   states &lt;- rownames(theemissionmatrix)

   # Make the Viterbi matrix v:
   v &lt;- makeViterbimat(sequence, transitionmatrix, emissionmatrix)

   # Go through each of the rows of the matrix v (where each row represents
   # a position in the DNA sequence), and find out which column has the 
   # maximum value for that row (where each column represents one state of
   # the HMM):
   mostprobablestatepath &lt;- apply(v, 1, function(x) which.max(x))
   
   # Print out the most probable state path:
   prevnucleotide &lt;- sequence[1]
   prevmostprobablestate &lt;- mostprobablestatepath[1]
   prevmostprobablestatename &lt;- states[prevmostprobablestate]
   startpos &lt;- 1
   for (i in 2:length(sequence))
   {
      nucleotide &lt;- sequence[i]
      mostprobablestate &lt;- mostprobablestatepath[i]
      mostprobablestatename &lt;- states[mostprobablestate]
      if (mostprobablestatename != prevmostprobablestatename)
      {
         print(paste(&quot;Positions&quot;,startpos,&quot;-&quot;,(i-1), &quot;Most probable state = &quot;, prevmostprobablestatename))
         startpos &lt;- i
      }
      prevnucleotide &lt;- nucleotide
      prevmostprobablestatename &lt;- mostprobablestatename
   } 
   print(paste(&quot;Positions&quot;,startpos,&quot;-&quot;,i, &quot;Most probable state = &quot;, prevmostprobablestatename))
}

# This makes the matrix v using the Viterbi algorithm.
# Adapted from &quot;Applied Statistics for Bioinformatics using R&quot; by Wim P. Krijnen, page 209
# ( <a href="http://cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf">cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf</a> ) 
makeViterbimat &lt;- function(sequence, transitionmatrix, emissionmatrix)
{
   # Change the sequence to uppercase
   sequence &lt;- toupper(sequence)
   # Find out how many states are in the HMM
   numstates &lt;- dim(transitionmatrix)[1]
   # Make a matrix with as many rows as positions in the sequence, and as many
   # columns as states in the HMM
   v &lt;- matrix(NA, nrow = length(sequence), ncol = dim(transitionmatrix)[1])
   # Set the values in the first row of matrix v (representing the first position of the sequence) to 0
   v[1, ] &lt;- 0
   # Set the value in the first row of matrix v, first column to 1
   v[1,1] &lt;- 1
   # Fill in the matrix v:
   for (i in 2:length(sequence)) # For each position in the DNA sequence:
   {
      for (l in 1:numstates) # For each of the states of in the HMM:
      {
         # Find the probabilility, if we are in state l, of choosing the nucleotide at position in the sequence 
         statelprobnucleotidei &lt;- emissionmatrix[l,sequence[i]]

         # v[(i-1),] gives the values of v for the (i-1)th row of v, ie. the (i-1)th position in the sequence.
         # In v[(i-1),] there are values of v at the (i-1)th row of the sequence for each possible state k.
         # v[(i-1),k] gives the value of v at the (i-1)th row of the sequence for a particular state k.
         
         # transitionmatrix[l,] gives the values in the lth row of the transition matrix, xx should not be transitionmatrix[,l]?
         # probabilities of changing from a previous state k to a current state l.

         # max(v[(i-1),] * transitionmatrix[l,]) is the maximum probability for the nucleotide observed
         # at the previous position in the sequence in state k, followed by a transition from previous
         # state k to current state l at the current nucleotide position.

         # Set the value in matrix v for row i (nucleotide position i), column l (state l) to be:
         v[i,l] &lt;-  statelprobnucleotidei * max(v[(i-1),] * transitionmatrix[,l])  
      }
   }
   
   return(v)
}

# Function to generate a DNA sequence, given a HMM and the length of the sequence to be generated.
generatehmmseq &lt;- function(transitionmatrix, emissionmatrix, initialprobs, seqlength)
{
    nucleotides     &lt;- c(&quot;A&quot;, &quot;C&quot;, &quot;G&quot;, &quot;T&quot;)   # Define the alphabet of nucleotides
    states          &lt;- c(&quot;AT-rich&quot;, &quot;GC-rich&quot;) # Define the names of the states
    mysequence      &lt;- character()             # Create a vector for storing the new sequence
    mystates        &lt;- character()             # Create a vector for storing the state that each position in the new sequence
                                               # was generated by
    # Choose the state for the first position in the sequence:
    firststate      &lt;- sample(states, 1, rep=TRUE, prob=initialprobs)
    # Get the probabilities of the current nucleotide, given that we are in the state &quot;firststate&quot;:
    probabilities   &lt;- emissionmatrix[firststate,]
    # Choose the nucleotide for the first position in the sequence:
    firstnucleotide &lt;- sample(nucleotides, 1, rep=TRUE, prob=probabilities)
    mysequence[1]   &lt;- firstnucleotide         # Store the nucleotide for the first position of the sequence
    mystates[1]     &lt;- firststate              # Store the state that the first position in the sequence was generated by

    for (i in 2:seqlength)
    {
       prevstate    &lt;- mystates[i-1]           # Get the state that the previous nucleotide in the sequence was generated by
       # Get the probabilities of the current state, given that the previous nucleotide was generated by state &quot;prevstate&quot;
       stateprobs   &lt;- transitionmatrix[prevstate,]
       # Choose the state for the ith position in the sequence:
       state        &lt;- sample(states, 1, rep=TRUE, prob=stateprobs)
       # Get the probabilities of the current nucleotide, given that we are in the state &quot;state&quot;:
       probabilities &lt;- emissionmatrix[state,]
       # Choose the nucleotide for the ith position in the sequence:
       nucleotide   &lt;- sample(nucleotides, 1, rep=TRUE, prob=probabilities)
       mysequence[i] &lt;- nucleotide             # Store the nucleotide for the current position of the sequence
       mystates[i]  &lt;- state                   # Store the state that the current position in the sequence was generated by
    }

    for (i in 1:length(mysequence))
    {
       nucleotide   &lt;- mysequence[i]
       state        &lt;- mystates[i]
       print(paste(&quot;Position&quot;, i, &quot;, State&quot;, state, &quot;, Nucleotide = &quot;, nucleotide))
    } 
}

# Function to colour the edges in a path on a graph in colour &quot;colour&quot;       
plotpath &lt;- function(graphplot,path,colour)
{
   # Get the names of the vertices in the shortest path:
   vertices &lt;- path$path_detail
   numvertices &lt;- length(vertices)
   # Now make a list that contains the names of pairs of vertices that are joined by edges:
   # There should be numvertices-1 pairs 
   numedges &lt;- numvertices - 1
   edges &lt;- numeric(numedges)
   edges2 &lt;- numeric(numedges)
   myvector &lt;- vector()
   for (i in 1:numedges)
   {
      vertex1 &lt;- vertices[i]
      vertex2 &lt;- vertices[i+1]
      edge &lt;- paste(vertex1,&quot;~&quot;,vertex2,sep=&quot;&quot;)
      edge2 &lt;- paste(vertex2,&quot;~&quot;,vertex1,sep=&quot;&quot;)
      myvector[`edge`] &lt;- colour   # Add named element to myvector
      myvector[`edge2`] &lt;- colour  # Add named element to myvector   
   }
   # Set the colour of the edges:
   edgeRenderInfo(graphplot) = list(col=myvector)
   renderGraph(graphplot)
}

# Function to make a graph based on protein-protein interaction data in an input file
makeproteingraph &lt;- function(myfile)
{
   library(&quot;graph&quot;)
   mytable &lt;- read.table(file(myfile)) # Store the data in a data frame
   proteins1 &lt;- mytable$V1
   proteins2 &lt;- mytable$V2
   protnames &lt;- c(levels(proteins1),levels(proteins2))
   # Find out how many pairs of proteins there are
   numpairs &lt;- length(proteins1)
   # Find the unique protein names:
   uniquenames &lt;-  unique(protnames)
   # Make a graph for these proteins with no edges:
   mygraph &lt;- new(&quot;graphNEL&quot;, nodes = uniquenames)
   # Add edges to the graph:
   # See <a href="http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf">http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf</a> for more examples
   weights &lt;- rep(1,numpairs)
   mygraph2 &lt;- addEdge(as.vector(proteins1),as.vector(proteins2),mygraph,weights)
 
   return(mygraph2)      
}

# Function to read in a regulatory network based on transcription factor-target gene interaction
# data in an input file
readRegulatoryNetwork &lt;- function(myfile)
{
   library(&quot;graph&quot;)
   mytable &lt;- read.table(file(myfile),header=FALSE) # Store the data in a data frame
   proteins &lt;- mytable$V1
   # Change the protein names to gene names:
   proteins2 &lt;- vector()
   for (i in 1:length(proteins))
   {
      protein &lt;- as.character(proteins[i])
      substr(protein,1,1) &lt;- tolower(substring(protein,1,1))
      proteins2[i] &lt;- protein
   }
   genes &lt;- mytable$V2
   types &lt;- mytable$V3
   # Find out the number of edges:
   numedges &lt;- length(proteins)
   # Find out how many transcription factor and target genes there are:
   nodenames &lt;- union(proteins2,genes)
   nodenames &lt;- unique(nodenames) 
   # Make a graph for these transcription factors and target genes with no edges:
   mygraph &lt;- new(&quot;graphNEL&quot;, nodes=nodenames, edgemode=&quot;directed&quot;)
   # Add edges to the graph:
   # See <a href="http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf">http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf</a> for more examples
   weights &lt;- rep(1,numedges)
   mygraph2 &lt;- addEdge(as.vector(proteins2),as.vector(genes),mygraph,weights)
   # We can store the information on how to plot the graph: 
   # Unfortunately, if we make a subgraph, this is lost
   # Set the types of the edges
   myvector &lt;- vector()
   for (i in 1:numedges)
   {
      vertex1 &lt;- proteins2[i]
      vertex2 &lt;- genes[i]
      edge &lt;- paste(vertex1,&quot;~&quot;,vertex2,sep=&quot;&quot;)
      type &lt;- types[i]
      if      (type == &#39;+&#39;) { edgetype &lt;- &#39;normal&#39;  }
      else if (type == &#39;-&#39;) { edgetype &lt;- &#39;tee&#39;     }
      else                  { edgetype &lt;- &#39;normal&#39;  } # both repressive+activating or unknown xxx
      myvector[`edge`] &lt;- edgetype # Add named element to myvector
   }
   # Set the type of the edges:
   edgeRenderInfo(mygraph2) = list(arrowhead=myvector)
   return(mygraph2)
}

# Function to make a regulatory network based on transcription factor-target gene interaction
# data in an input file, and plot it
readAndPlotRegulatoryNetwork &lt;- function(myfile)
{
   library(&quot;graph&quot;)
   mytable &lt;- read.table(file(myfile),header=FALSE) # Store the data in a data frame
   proteins &lt;- mytable$V1
   # Change the protein names to gene names:
   proteins2 &lt;- vector()
   for (i in 1:length(proteins))
   {
      protein &lt;- as.character(proteins[i])
      substr(protein,1,1) &lt;- tolower(substring(protein,1,1))
      proteins2[i] &lt;- protein
   }
   genes &lt;- mytable$V2
   types &lt;- mytable$V3
   # Find out the number of edges:
   numedges &lt;- length(proteins)
   # Find out how many transcription factor and target genes there are:
   nodenames &lt;- union(proteins2,genes)
   nodenames &lt;- unique(nodenames) 
   # Make a graph for these transcription factors and target genes with no edges:
   mygraph &lt;- new(&quot;graphNEL&quot;, nodes=nodenames, edgemode=&quot;directed&quot;)
   # Add edges to the graph:
   # See <a href="http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf">http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf</a> for more examples
   weights &lt;- rep(1,numedges)
   mygraph2 &lt;- addEdge(as.vector(proteins2),as.vector(genes),mygraph,weights)
   # Set the types of the edges
   myvector &lt;- vector()
   for (i in 1:numedges)
   {
      vertex1 &lt;- proteins2[i]
      vertex2 &lt;- genes[i]
      edge &lt;- paste(vertex1,&quot;~&quot;,vertex2,sep=&quot;&quot;)
      type &lt;- types[i]
      if      (type == &#39;+&#39;) { edgetype &lt;- &#39;normal&#39;  }
      else if (type == &#39;-&#39;) { edgetype &lt;- &#39;tee&#39;     }
      else                  { edgetype &lt;- &#39;normal&#39;  } # xxx both repressive and activating, or unknown
      myvector[`edge`] &lt;- edgetype # Add named element to myvector
   }
   # Plot the graph:
   library(&quot;Rgraphviz&quot;)
   plot.new() # Make a new plot 
   mygraphplot &lt;- layoutGraph(mygraph2, layoutType=&quot;neato&quot;)  
   # Set the colour of the edges:
   edgeRenderInfo(mygraphplot) = list(arrowhead=myvector)
   renderGraph(mygraphplot)
   return(mygraph2)
}

# Function to find network communities in a graph
findcommunities &lt;- function(mygraph,minsize)
{
   # Load up the igraph library:
   library(&quot;igraph&quot;)
   # Set the counter for the number of communities:
   cnt &lt;- 0
   # First find the connected components in the graph:
   myconnectedcomponents &lt;- connectedComp(mygraph)
   # For each connected component, find the communities within that connected component:
   numconnectedcomponents &lt;- length(myconnectedcomponents)
   for (i in 1:numconnectedcomponents)
   {
      component &lt;- myconnectedcomponents[[i]]
      # Find the number of nodes in this connected component:
      numnodes &lt;- length(component)
      if (numnodes &gt; 1) # We can only find communities if there is more than one node
      {
         mysubgraph &lt;- subGraph(component, mygraph)
         # Find the communities within this connected component:
         # print(component)
         myvector &lt;- vector()
         mylist &lt;- findcommunities2(mysubgraph,cnt,&quot;FALSE&quot;,myvector,minsize)
         cnt &lt;- mylist[[1]]
         myvector &lt;- mylist[[2]]
      }
   }
   print(paste(&quot;There were&quot;,cnt,&quot;communities in the input graph&quot;))
} 

# Function to find network communities in a connected component of a graph    
findcommunities2 &lt;- function(mygraph,cnt,plot,myvector,minsize)
{
   # Find the number of nodes in the input graph
   nodes &lt;- nodes(mygraph)
   numnodes &lt;- length(nodes)

   # Record the vertex number for each vertex name
   myvector &lt;- vector()
   for (i in 1:numnodes)
   {
      node &lt;- nodes[i] # &quot;node&quot; is the vertex name, i is the vertex number
      myvector[`node`] &lt;- i  # Add named element to myvector
   }

   # Create a graph in the &quot;igraph&quot; library format, with numnodes nodes:
   newgraph &lt;- graph.empty(n=numnodes,directed=FALSE)
   # First record which edges we have seen already in the &quot;mymatrix&quot; matrix,
   # so that we don&#39;t add any edge twice:
   mymatrix &lt;- matrix(nrow=numnodes,ncol=numnodes)
   for (i in 1:numnodes)
   {
      for (j in 1:numnodes)
      {
         mymatrix[i,j] = 0
         mymatrix[j,i] = 0
      }
   }
   # Now add edges to the graph &quot;newgraph&quot;:
   for (i in 1:numnodes)
   {
      node &lt;- nodes[i] # &quot;node&quot; is the vertex name, i is the vertex number
      # Find the nodes that this node is joined to:
      neighbours &lt;- adj(mygraph, node)
      neighbours &lt;- neighbours[[1]] # Get the list of neighbours
      numneighbours &lt;- length(neighbours)
      if (numneighbours &gt;= 1) # If this node &quot;node&quot; has some edges to other nodes
      {
         for (j in 1:numneighbours)
         {
            neighbour &lt;- neighbours[j]
            # Get the vertex number
            neighbourindex &lt;- myvector[neighbour]
            # Add a node in the new graph &quot;newgraph&quot; between vertices (i-1) and (neighbourindex-1)
            # In graph &quot;newgraph&quot;, the vertices are counted from 0 upwards.
            indexi &lt;- i
            indexj &lt;- neighbourindex
            # If we have not seen this edge already:
            if (mymatrix[indexi,indexj] == 0 &amp;&amp; mymatrix[indexj,indexi] == 0)
            {
               mymatrix[indexi,indexj] &lt;- 1
               mymatrix[indexj,indexi] &lt;- 1 
               # Add edges to the graph &quot;newgraph&quot;
               newgraph &lt;- add.edges(newgraph, c(i-1, neighbourindex-1))
            }
         }   
      }
   }
   # Set the names of the vertices in graph &quot;newgraph&quot;:
   newgraph &lt;- set.vertex.attribute(newgraph, &quot;name&quot;, value=nodes)

   # Now find communities in the graph:
   communities &lt;- spinglass.community(newgraph)
   # Find how many communities there are:
   sizecommunities &lt;- communities$csize 
   numcommunities &lt;- length(sizecommunities)
   # Find which vertices belong to which communities:
   membership &lt;- communities$membership
   # Get the names of vertices in the graph &quot;newgraph&quot;:
   vertexnames &lt;- get.vertex.attribute(newgraph, &quot;name&quot;)

   # Print out the vertices belonging to each community:
   for (i in 1:numcommunities)
   {
      if (plot == TRUE) { cnt &lt;- cnt + 1 } 
      nummembers &lt;- 0
      printout &lt;- paste(&quot;Community&quot;,cnt,&quot;:&quot;) 
      for (j in 1:length(membership))
      {
         community &lt;- membership[j]
         community &lt;- community + 1
         if (community == i) # If vertex j belongs to the ith community
         {
            vertexname &lt;- vertexnames[j]
            if (plot == FALSE)
            {
               nummembers &lt;- nummembers + 1
               # Print out the vertices belonging to the community
               printout &lt;- paste(printout,vertexname)
            }
            else
            {
               # Colour in the vertices belonging to the community
               myvector[`vertexname`] &lt;- cnt   
            }
         } 
      }
      if (plot == FALSE &amp;&amp; nummembers &gt;= minsize) 
      { 
         cnt &lt;- cnt + 1
         print(printout) 
      } 
   }

   return(list(cnt,myvector))
}

# Function to plot network communities in a graph
plotcommunities &lt;- function(mygraph)
{
   # Load the &quot;igraph&quot; library:
   library(&quot;igraph&quot;)
   # Make a plot of the graph
   graphplot &lt;- layoutGraph(mygraph, layoutType=&quot;neato&quot;)  
   renderGraph(graphplot) 
   # Get the names of the nodes in the graph:        
   vertices &lt;- nodes(mygraph)
   numvertices &lt;- length(vertices)
   # Now record the colour of each vertex in a vector &quot;myvector&quot;:
   myvector &lt;- vector()
   colour &lt;- &quot;red&quot;
   for (i in 1:numvertices)
   {
      vertex &lt;- vertices[i]
      myvector[`vertex`] &lt;- colour   # Add named element to myvector
   }
   
   # Set the counter for the number of communities:
   cnt &lt;- 0
   # First find the connected components in the graph:
   myconnectedcomponents &lt;- connectedComp(mygraph)
   # For each connected component, find the communities within that connected component:
   numconnectedcomponents &lt;- length(myconnectedcomponents)
   for (i in 1:numconnectedcomponents)
   {
      component &lt;- myconnectedcomponents[[i]]
      # Find the number of nodes in this connected component:
      numnodes &lt;- length(component)
      if (numnodes &gt; 1) # We can only find communities if there is more than one node
      {
         mysubgraph &lt;- subGraph(component, mygraph)
         # Find the communities within this connected component:
         # print(component)
         mylist &lt;- findcommunities2(mysubgraph,cnt,&quot;TRUE&quot;,myvector,0)
         cnt &lt;- mylist[[1]]
         myvector &lt;- mylist[[2]]
      }
   }

   # Get a set of cnt colours, where cnt is equal to the number of communities found:
   mycolours &lt;- rainbow(cnt)
   for (i in 1:length(mycolours))
   {
      # print(paste(&quot;Community&quot;,i,&quot;is in&quot;,mycolours[i]))
   }
   # Set the colour of the vertices, so that vertices in each community are of the same colour,
   # and vertices in different communities are different colours:
   myvector2 &lt;- vector()
   for (i in 1:numvertices)
   {
      vertex &lt;- vertices[i]
      community &lt;- myvector[vertex]
      mycolour &lt;- mycolours[community] 
      # print(paste(&quot;Vertex&quot;,vertex,&quot;is in community&quot;,community,&quot;of colour&quot;,mycolour))
      myvector2[`vertex`] &lt;- mycolour
   }
   nodeRenderInfo(graphplot) = list(fill=myvector2)
   renderGraph(graphplot)
} 

# Function to make a random graph                  
makerandomgraph &lt;- function(numvertices,numedges)
{
   library(&quot;graph&quot;)
   # Make a vector with the names of the vertices
   mynames &lt;- sapply(seq(1,numvertices),toString) 
   myrandomgraph &lt;- randomEGraph(mynames, edges = numedges)

   return(myrandomgraph)
}

# Function to find the connected component that contains a particular vertex
findcomponent &lt;- function(graph,vertex)
{
   library(&quot;RBGL&quot;)
   found &lt;- 0
   myconnectedcomponents &lt;- connectedComp(graph)
   numconnectedcomponents &lt;- length(myconnectedcomponents)
   for (i in 1:numconnectedcomponents)
   {
      componenti &lt;- myconnectedcomponents[[i]] 
      numvertices &lt;- length(componenti)
      for (j in 1:numvertices)
      {
         vertexj &lt;- componenti[j]
         if (vertexj == vertex) 
         {
            found &lt;- 1
            return(componenti)
         } 
      }
   }
   print(&quot;ERROR: did not find vertex in the graph&quot;)
} 

makeigraphgraph &lt;- function(mygraph)
{
   library(&quot;igraph&quot;)
   # Find the number of nodes in the input graph
   nodes &lt;- nodes(mygraph)
   numnodes &lt;- length(nodes)

   # Record the vertex number for each vertex name
   myvector &lt;- vector()
   for (i in 1:numnodes)
   {
      node &lt;- nodes[i] # &quot;node&quot; is the vertex name, i is the vertex number
      myvector[`node`] &lt;- i  # Add named element to myvector
   }

   # Create a graph in the &quot;igraph&quot; library format, with numnodes nodes:
   newgraph &lt;- graph.empty(n=numnodes,directed=FALSE)
   # First record which edges we have seen already in the &quot;mymatrix&quot; matrix,
   # so that we don&#39;t add any edge twice:
   mymatrix &lt;- matrix(nrow=numnodes,ncol=numnodes)
   for (i in 1:numnodes)
   {
      for (j in 1:numnodes)
      {
         mymatrix[i,j] = 0
         mymatrix[j,i] = 0
      }
   }
   # Now add edges to the graph &quot;newgraph&quot;:
   for (i in 1:numnodes)
   {
      node &lt;- nodes[i] # &quot;node&quot; is the vertex name, i is the vertex number
      # Find the nodes that this node is joined to:
      neighbours &lt;- adj(mygraph, node)
      neighbours &lt;- neighbours[[1]] # Get the list of neighbours
      numneighbours &lt;- length(neighbours)
      if (numneighbours &gt;= 1) # If this node &quot;node&quot; has some edges to other nodes
      {
         for (j in 1:numneighbours)
         {
            neighbour &lt;- neighbours[j]
            # Get the vertex number
            neighbourindex &lt;- myvector[neighbour]
            # Add a node in the new graph &quot;newgraph&quot; between vertices (i-1) and (neighbourindex-1)
            # In graph &quot;newgraph&quot;, the vertices are counted from 0 upwards.
            indexi &lt;- i
            indexj &lt;- neighbourindex
            # If we have not seen this edge already:
            if (mymatrix[indexi,indexj] == 0 &amp;&amp; mymatrix[indexj,indexi] == 0)
            {
               mymatrix[indexi,indexj] &lt;- 1
               mymatrix[indexj,indexi] &lt;- 1 
               # Add edges to the graph &quot;newgraph&quot;
               newgraph &lt;- add.edges(newgraph, c(i-1, neighbourindex-1))
            }
         }   
      }
   }
   # Set the names of the vertices in graph &quot;newgraph&quot;:
   newgraph &lt;- set.vertex.attribute(newgraph, &quot;name&quot;, value=nodes)

   return(newgraph)
}

# Function to print out the interactions of a particular protein, in a particular experimental data set:
# exp is a vector containing the experimental data.
printproteininteractions &lt;- function(exp, protein)
{
	for (i in 1:length(exp))
	{
		pair &lt;- exp[i]
		thepair &lt;- unlist(strsplit(pair, &quot;~&quot;, fixed = TRUE))
		protein1 &lt;- thepair[1]
		protein2 &lt;- thepair[2]
		if (protein1 == protein)
		{
			print(thepair)
		}
		else if (protein2 == protein)
		{
			print(thepair)
		}
	}
}

# Function to use a Bayes network for protein-protein interaction data,
# to calculate likelihoods for particular protein protein interactions.
# &quot;goldpos&quot; is the set of gold standard positive protein-protein interactions
# &quot;goldneg&quot; is the set of gold standard negative protein-protein interactions
# &quot;exp&quot; is a list of vectors, each vector containing some experimental data on protein-protein interactions 
#   reported by a particular experiment
bayesforproteinpairs &lt;- function(goldpos, goldneg, exp, protein1, protein2, VALUEONLY=FALSE)
{
    # Make a vector containing the pair of interest:
    myvector &lt;- vector()
    myvector[1] &lt;- protein1
    myvector[2] &lt;- protein2
    myvector &lt;- sort(myvector)
    protpair2 &lt;- paste(myvector[1],myvector[2],sep=&quot;~&quot;)

    # Calculate a likelihood ratio for protpair2:

    # Check how many gold standard positive pairs have the same
    # pattern of &#39;yes&#39;s &amp; &#39;no&#39;s as protpair2 for the experimental data sets: 
    numexpdatasets &lt;- length(exp) # Number of experimental data sets
    numgoldpos &lt;- length(goldpos) # Number of gold standard positives
    samepattern &lt;- goldpos
    for (m in 1:numexpdatasets)
    {	
	expm &lt;- exp[[m]]
	isfoundm &lt;- is.element(protpair2,expm) # Says if protpair2 is in the mth experimental data set
	if (isfoundm == &quot;TRUE&quot;)
	{
	   samepattern &lt;- intersect(samepattern,expm) # Find pairs that are in goldpos and in expm
	}
	else
	{
	   samepattern &lt;- setdiff(samepattern,expm) # Find pairs that are in goldpos but not in expm
	}
    }
    #print(paste(&quot;Have same pattern in gold pos and expts:&quot;,samepattern))
    samepatterncnt &lt;- length(samepattern)
    numeratoroflikelihoodratio &lt;- samepatterncnt/numgoldpos
	
    # Check how many gold standard negative pairs have the same
    # pattern of &#39;yes&#39;s &amp; &#39;no&#39;s for the experimental data sets:
    numgoldneg &lt;- length(goldneg) # Number of gold standard negatives
    samepattern &lt;- goldneg
    for (m in 1:numexpdatasets)
    {
 	expm &lt;- exp[[m]]
	isfoundm &lt;- is.element(protpair2,expm) # Says if protpair2 is in the mth experimental data set
	if (isfoundm == &quot;TRUE&quot;)
	{
    	    samepattern &lt;- intersect(samepattern,expm) # Find pairs that are in goldpos and in expm
        }
	else
	{
	    samepattern &lt;- setdiff(samepattern,expm) # Find pairs that are in goldpos but not in expm
	}
    }
    #print(paste(&quot;Have same pattern in gold neg and expts:&quot;,samepattern))
    samepatterncnt &lt;- length(samepattern)
    denominatoroflikelihoodratio &lt;- samepatterncnt/numgoldneg

    # Calculate the likelihood ratio:
    #print(paste(&quot;numerator:&quot;,numeratoroflikelihoodratio,&quot;denominator:&quot;,denominatoroflikelihoodratio))
    likelihoodratio &lt;- numeratoroflikelihoodratio/denominatoroflikelihoodratio 
    # Print out the result for this protein-protein interaction:
    if (VALUEONLY == FALSE)
    {
       print(paste(&quot;Pair&quot;,protpair2,&quot;confidence measure (likelihood ratio) =&quot;,likelihoodratio))
    }
    else
    {
       return(likelihoodratio)
    }
}

# Function to calculate the rate of true positives and false positives, and
# false negatives in an experimental data set of protein-protein interactions:
calcproteinpairaccuracy &lt;- function(goldpos, goldneg, exp)
{
   # Find the number of true interactions in the gold positive data:
   numgoldpos &lt;- length(goldpos)

   # Find the number of true positive protein-protein interactions:
   truepos &lt;- intersect(exp, goldpos)
   numtruepos &lt;- length(truepos)
   
   # Find the number of false positive protein-protein interactions:
   falsepos &lt;- intersect(exp, goldneg)  
   numfalsepos &lt;- length(falsepos)
  
   # Find the number of false negative protein-protein interactions:
   falseneg &lt;- setdiff(goldpos, exp)    
   numfalseneg &lt;- length(falseneg)
 
   # Find the number of protein-protein interactions that overlap between the
   # experimental data set and gold standard data:
   numoverlap &lt;- numtruepos + numfalsepos
   # Find the rate of true positives, and false positives:
   ratetruepos &lt;- (numtruepos/numoverlap)*100
   ratefalsepos &lt;- (numfalsepos/numoverlap)*100
   # Find the rate of false negatives:
   ratefalseneg &lt;- (numfalseneg/numgoldpos)*100

   print(paste(&quot;Number of protein-protein pairs that overlap between the gold standard and experimental data:&quot;,numoverlap))
   print(paste(&quot;Number of true positive pairs in the experimental data:&quot;,numtruepos,&quot;(&quot;,ratetruepos,&quot;% of&quot;,numoverlap,&quot;reported interactions)&quot;))
   print(paste(&quot;Number of false positive pairs in the experimental data:&quot;,numfalsepos,&quot;(&quot;,ratefalsepos,&quot;% of&quot;,numoverlap,&quot;reported interactions)&quot;))
   print(paste(&quot;Number of false negative pairs in the experimental data:&quot;,numfalseneg,&quot;(&quot;,ratefalseneg,&quot;% of&quot;,numgoldpos,&quot;true interactions)&quot;))
}

# Function to read in pairs of interacting proteins in an input file:
readproteinpairs &lt;- function(file)
{
   MyData &lt;- read.table(`file`,header=FALSE)
   mypairs &lt;- paste(MyData$V1,MyData$V2,sep=&quot;~&quot;)

   return(mypairs) 
}

# Use the Needleman-Wunsch algorithm to make an alignment of two strings
# S1 and S2
needlemanwunsch &lt;- function(S1,S2,type=&quot;protein&quot;,gappenalty=-1) 
{
   library(Biostrings)
   library(seqinr)
   S1 &lt;- s2c(S1) # Convert the string to a vector of characters
   S2 &lt;- s2c(S2) # Convert the string to a vector of characters
   lengthS1 &lt;- length(S1)
   lengthS2 &lt;- length(S2)
   # Check if it is a protein alignment, or DNA alignment:
   if (type == &quot;protein&quot;)
   {
      # Read the BLOSUM matrix:
      data(BLOSUM45)
      s &lt;- BLOSUM45[S2,S1]
      d &lt;- -gappenalty # The gap penalty is -d
   }
   else
   {
      # Define the DNA scoring matrix:
      s1 &lt;- nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = TRUE)
      s &lt;- s1[S2, S1] 
      d &lt;- -gappenalty # The gap penalty is -d
   }

   # Make the table T with lengthS1+1 columns and lengthS2+1 rows:
   # eg. if S1 = VIVADAVIS and S2 = VIVALASVEGAS, then we will
   # have a matrix T with 10 columns and 13 rows
   T &lt;- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)

   # Make a matrix for storing the traceback:
   T2 &lt;- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)

   # Initialise the first column and row of matrix T:
   T[1,] &lt;- -seq(0,lengthS1*d,d) # Fill in the first row
   T[,1] &lt;- -seq(0,lengthS2*d,d) # Fill in the first column
   # In R T[1,2] refers to row 1, column 2
   #      T[2,1] refers to row 2, column 1
   # print(T) 

   # Initialise the first column and row of matrix T2:
   T2[1,] &lt;- rep(&#39;-&#39;,lengthS1+1) # Fill in the first row
   T2[,1] &lt;- rep(&#39;|&#39;,lengthS2+1) # Fill in the first column
   # print(T2)

   # Carry out the recurrence relation to fill in matrix T:
   for (j in 2:(nrow(T))) # Go through each row (j) at a time
   {
      for (i in 2:(ncol(T))) # Go through each column (i) at a time
      {
         print(paste(&quot;Calculating for col i=&quot;,i,&quot; row j=&quot;,j))
         # Set the value of T[i,j] ie. in column i and row j:
         T[j,i] &lt;- max(c(T[j-1,i-1]+s[j-1,i-1],T[j,i-1]-d,T[j-1,i]-d))
         diag &lt;- T[j-1,i-1]+s[j-1,i-1]
         up &lt;- T[j-1,i]-d
         left &lt;- T[j,i-1]-d
         if (diag &gt; up &amp;&amp; diag &gt; left)
         {
            T2[j,i] &lt;- &quot;&gt;&quot;
         }
         else if (up &gt; diag &amp;&amp; up &gt; left)
         {
            T2[j,i] &lt;- &quot;|&quot;
         }
         else if (left &gt; diag &amp;&amp; left &gt; up)
         {
            T2[j,i] &lt;- &quot;-&quot;
         }
         else if (left == diag &amp;&amp; up == diag)
         {
            T2[j,i] &lt;- &quot;*&quot;
         }
         else if (up == left &amp;&amp; up &gt; diag &amp;&amp; left &gt; diag)
         {
            T2[j,i] &lt;- &quot;L&quot;
         }
         else if (up == diag &amp;&amp; up &gt; left &amp;&amp; diag &gt; left)
         {
            T2[j,i] &lt;- &quot;V&quot;
         }
         else if (diag == left &amp;&amp; diag &gt; up &amp;&amp; left &gt; up)
         {
            T2[j,i] &lt;- &quot;Z&quot;
         }
         else
         {
            print(paste(&quot;ERROR: up&quot;,up,&quot;left&quot;,left,&quot;diag&quot;,diag))
         }
      }
   }
   print(T) 
   print(T2)

   # Print out T and T2 together:
   T3 &lt;- matrix(data=NA,ncol=lengthS1+1,nrow=lengthS2+1)
   for (j in 2:(nrow(T))) # Go through each row (j) at a time
   {
      for (i in 2:(ncol(T))) # Go through each column (i) at a time
      {
         value &lt;- T[j,i]
         value2 &lt;- T2[j,i]
         value3 &lt;- paste(value,value2)
         T3[j,i] &lt;- value3
      }
   }
   print(T3)
}

# Function to retrieve Uniprot sequences using the SeqinR R library.
retrieveuniprotseqs &lt;- function(seqnames)
{
	myseqs &lt;- list() # Make a list to store the sequences
	library(&quot;seqinr&quot;)  # Load the SeqinR R library
	for (i in 1:length(seqnames))
	{
		seqname &lt;- seqnames[i]
		print(paste(&quot;Retrieving sequence&quot;,seqname,&quot;...&quot;))
		choosebank(&quot;swissprot&quot;)
		queryname &lt;- &quot;query2&quot;
		query &lt;- paste(&quot;AC=&quot;,seqname,sep=&quot;&quot;)
		query(`queryname`,`query`)
		seq &lt;- getSequence(query2$req) # Makes a vector &quot;seq&quot; containing the sequence
                if (is.list(seq) == TRUE)
                {
                   seq1b &lt;- seq[[1]]
                   seq &lt;- seq1b
                   if (is.list(seq1b) == TRUE)
                   {
                      seq1c &lt;- seq[[1]]
                      seq &lt;- seq1c 
                   }
                }
		myseqs[[i]] &lt;- seq
	}
	return(myseqs)
}

# Function to plot the attractors in the Boolean network
plotBooleanAttractors &lt;- function(attractors)
{
   par(mfrow=c(1,length(attractors$attractors)))
   plotAttractors(attractors) 
}

# Function to find network motifs in a gene regulatory network
findnetworkpatterns &lt;- function(graph, patternsize, min=1)
{
   library(igraph)
   # Convert the graph into an &quot;igraph&quot; graph:
   mygraph2 &lt;- igraph.from.graphNEL(graph)
   # Find network motifs in the graph &quot;mygraph2&quot;:
   mygraphmotifs &lt;- graph.motifs(mygraph2, patternsize)
   # Find how many motifs occur 1 or more times:
   numdifferentmotifs &lt;- table(mygraphmotifs!=0)[[2]]
   # Make a plot for these motifs:
   #par(mfrow=c(3,4)) 
   #plot.new() # Make a new plot 
   # Plot the 12 most common motifs:
   numplotted &lt;- 0
   oo &lt;- order(mygraphmotifs, decreasing=TRUE)
   # Find which motifs occur:
   nummotifs &lt;- length(mygraphmotifs)
   for (i in 1:nummotifs)
   {  
      motifi &lt;- oo[i]
      motif &lt;- mygraphmotifs[motifi] # The number of occurrences of this motif
      if (motif &gt;= min) # There are some occurrences of this motif
      {
         # Find out what the motif looks like:
         motifgraph &lt;- graph.isocreate(size=patternsize, number=motifi-1, directed=TRUE)
         # Make a plot
         #if (numplotted &lt; 12)
         #{
         #   numplotted &lt;- numplotted + 1
         #   motifgraph2 &lt;- igraph.to.graphNEL(motifgraph)
         #   plot(motifgraph2,attrs=list(node=list(width=&quot;2&quot;,height=&quot;2&quot;,fontsize=&quot;1&quot;),graph=list(size=&quot;5&quot;)))
         #}
         edges &lt;- E(motifgraph)
	 print(paste(&quot;This pattern occurs&quot;,motif,&quot;times:&quot;))
	 print(edges)
      }
   }
}

# Function to make a plot of a subgraph/graph in a regulatory network.
# By default, it makes a plot of the whole regulatory network.
# However, if you pass in a list of vertices in a subgraph, it will plot just that subgraph. 
plotRegulatoryNetwork &lt;- function(graph, nodes=attr(graph,&quot;nodes&quot;),type=&quot;neato&quot;)
{
   # Make a subgraph of &quot;graph&quot; corresponding to the nodes in vector &quot;nodes&quot;
   mysubgraph &lt;- subGraph(nodes, graph)
   # Copy the attributes of the edges from &quot;graph&quot; to &quot;mysubgraph&quot;
   myrenderinfo &lt;- attr(graph, &quot;renderInfo&quot;)
   myrenderinfo2 &lt;- attr(myrenderinfo, &quot;edges&quot;)
   myrenderinfo3 &lt;- myrenderinfo2$arrowhead
   # Get all the edges in mysubgraph
   myedges &lt;- listEdges(mysubgraph)
   numedges &lt;- length(myedges)
   myvector &lt;- vector()
   for (i in 1:numedges)
   {
      edge &lt;- myedges[[i]]
      # Sometimes there can be an edge each way between two genes,
      # eg. FIS regulates CRP and CRP regulates FIS.
      numedges2 &lt;- length(edge)
      for (j in 1:numedges2)
      {
         edge2 &lt;- edge[[j]]
         # This seems to be necessary: 
         if (j == 1 &amp;&amp; numedges2 &gt; 1) { edge2 &lt;- edge2[[1]] }
         beginnode &lt;- slot(edge2,&quot;bNode&quot;)
         endnode &lt;- slot(edge2,&quot;eNode&quot;)
         edgename &lt;- paste(beginnode,&quot;~&quot;,endnode,sep=&quot;&quot;)
         # Find out the renderinfo for this edge in graph &quot;graph&quot;
         edgetype &lt;- myrenderinfo3[`edgename`]
         edgetype &lt;- edgetype[[1]]
         myvector[`edgename`] &lt;- edgetype
      }
   }
   # Set the type of the edges:
   edgeRenderInfo(mysubgraph) = list(arrowhead=myvector)
   # Make a plot now:
   library(&quot;Rgraphviz&quot;)
   # Set the type of plot:
   if (type == &quot;hierarchical&quot;) { type &lt;- &quot;dot&quot; }
   mygraphplot &lt;- layoutGraph(mysubgraph, layoutType=type)
   renderGraph(mygraphplot)	
}

# Function to get all the descendents of a particular node in a directed graphNEL graph:
getdescendents &lt;- function(graph, node)
{
   # First get all the nodes that are children of node &quot;node&quot;:
   children &lt;- adj(graph, node)[[1]]
   numchildren &lt;- length(children)
   descendents &lt;- vector() # Make a vector for storing all descendents
   # Record whether or not we have seen any of the children before
   allnodes &lt;- nodes(graph)
   seen &lt;- vector()
   for (i in 1:length(allnodes))
   {
      node &lt;- allnodes[i]
      seen[`node`] &lt;- 0
   }

   # For each of the children, check if it has its own children
   loopagain &lt;- 1
   while(loopagain == 1 &amp;&amp; numchildren &gt; 0)
   {
      newchildren &lt;- vector()
      foundnew &lt;- 0
      # Record whether or not we have seen any of the children before
      for (i in 1:numchildren)
      {
         child &lt;- children[i]
         seenchild &lt;- seen[`child`]
         if (seenchild == 0)
         {
            seen[`child`] &lt;- 1
            descendents &lt;- append(descendents, child, after=length(descendents))
         }
      }
      for (i in 1:numchildren)
      {
         child &lt;- children[i]
         grandchildren &lt;- adj(graph, child)[[1]]
         numgrandchildren &lt;- length(grandchildren)
	 if (numgrandchildren &gt; 0 &amp;&amp; child != node)
	 {
	    for (j in 1:numgrandchildren)
	    { 
               grandchild &lt;- grandchildren[j]
	       seengrandchild &lt;- seen[`grandchild`]
	       if (grandchild != node &amp;&amp; seengrandchild == 0)
	       {
	          newchildren &lt;- append(newchildren, grandchild, after=length(newchildren))
		  foundnew &lt;- 1
	       }
            }
	 }
      }
      if (foundnew == 0)
      { 
         loopagain &lt;- 0
      }
      else
      {
         children &lt;- newchildren
	 numchildren &lt;- length(children)
      }
   }
   return(descendents)
}

# Function to retrieve GenBank sequences using the SeqinR R library.
retrievegenbankseqs &lt;- function(seqnames, from=&quot;NA&quot;, to=&quot;NA&quot;)
{
	myseqs &lt;- list()
	library(&quot;seqinr&quot;)
	for (i in 1:length(seqnames))
	{
		seqname &lt;- seqnames[i]
		print(paste(&quot;Retrieving sequence&quot;,seqname,&quot;...&quot;))
                if (substring(seqname,1,3) == &#39;XM_&#39; || substring(seqname,1,3) == &#39;NM_&#39;)
                {
                   choosebank(&quot;refseq&quot;)
                }
                else
                {
		   choosebank(&quot;genbank&quot;)
                }
		queryname &lt;- &quot;query2&quot;
		query &lt;- paste(&quot;AC=&quot;,seqname,sep=&quot;&quot;)
		query(`queryname`,`query`)
		seq &lt;- getSequence(query2$req) # Makes a vector &quot;seq&quot; containing the sequence
                if (is.list(seq) == TRUE)
                {
                   seq1b &lt;- seq[[1]]
                   seq &lt;- seq1b
                   if (is.list(seq1b) == TRUE)
                   {
                      seq1c &lt;- seq[[1]]
                      seq &lt;- seq1c 
                   }
                }
                if (from != &#39;NA&#39; &amp;&amp; to != &#39;NA&#39;)
                {
                   seq &lt;- seq[from:to]
                }
	        myseqs[[i]] &lt;- seq
	}
	return(myseqs)
}

# Function to retrieve GenBank virus sequences using the SeqinR R library.
retrievevirusseqs &lt;- function(seqnames)
{
	myseqs &lt;- list()
	library(&quot;seqinr&quot;)
	for (i in 1:length(seqnames))
	{
		seqname &lt;- seqnames[i]
		print(paste(&quot;Retrieving sequence&quot;,seqname,&quot;...&quot;))
		choosebank(&quot;refseqViruses&quot;)
		queryname &lt;- &quot;query2&quot;
		query &lt;- paste(&quot;AC=&quot;,seqname,sep=&quot;&quot;)
		query(`queryname`,`query`)
		seq &lt;- getSequence(query2$req) # Makes a vector &quot;seq&quot; containing the sequence
                if (is.list(seq) == TRUE)
                {
                   seq1b &lt;- seq[[1]]
                   seq &lt;- seq1b
                   if (is.list(seq1b) == TRUE)
                   {
                      seq1c &lt;- seq[[1]]
                      seq &lt;- seq1c 
                   }
                }
		myseqs[[i]] &lt;- seq
	}
	return(myseqs)
}


# Function to convert a graphNEL graph to an igraph graph.
# Copied 6-Feb-10 from:
# <a href="http://bazaar.launchpad.net/~igraph/igraph/0.6-main/annotate/head:/interfaces/R/igraph/R/conversion.R">http://bazaar.launchpad.net/~igraph/igraph/0.6-main/annotate/head:/interfaces/R/igraph/R/conversion.R</a>
igraph.from.graphNEL &lt;- function(graphNEL, name=TRUE, weight=TRUE,
                                 unlist.attrs=TRUE) 
{
  require(graph)
  if (!is(graphNEL, &quot;graphNEL&quot;)) {
    stop(&quot;Not a graphNEL graph&quot;)
  }

  al &lt;- lapply(edgeL(graphNEL), &quot;[[&quot;, &quot;edges&quot;)

  if (edgemode(graphNEL)==&quot;undirected&quot;) {
    al &lt;- mapply(SIMPLIFY=FALSE, seq_along(al), al, FUN=function(n, l) {
      c(l, rep(n, sum(l==n)))
    })
  }

  al &lt;- lapply(al, function(x) x-1)

  mode &lt;- if (edgemode(graphNEL)==&quot;directed&quot;) &quot;out&quot; else &quot;all&quot;

  g &lt;- graph.adjlist(al, directed=TRUE, duplicate=TRUE) 
  if (name) {
    V(g)$name &lt;- nodes(graphNEL)
  }

  ## Graph attributes

  g.n &lt;- names(graphNEL@graphData)
  g.n &lt;- g.n [ g.n != &quot;edgemode&quot; ]
  for (n in g.n) {
    g &lt;- set.graph.attribute(g, n, graphNEL@graphData[[n]])
  }

  ## Vertex attributes
  v.n &lt;- names(nodeDataDefaults(graphNEL))

  for (n in v.n) {
    val &lt;- unname(nodeData(graphNEL, attr=n))
    if (unlist.attrs &amp;&amp; all(sapply(val, length)==1)) { val &lt;- unlist(val) }
    g &lt;- set.vertex.attribute(g, n, value=val)
  }

  ## Edge attributes
  e.n &lt;- names(edgeDataDefaults(graphNEL))
  if (!weight) { e.n &lt;- e.n [ e.n != &quot;weight&quot; ] }

  if (length(e.n) &gt; 0) {
    el &lt;- get.edgelist(g)
    el &lt;- paste(sep=&quot;|&quot;, el[,1], el[,2])
    for (n in e.n) {
      val &lt;- unname(edgeData(graphNEL, attr=n)[el])
      if (unlist.attrs &amp;&amp; all(sapply(val, length)==1)) { val &lt;- unlist(val) }
      g &lt;- set.edge.attribute(g, n, value=val)
    }
  }
  g 
}

# Function to change an igraph graph into a graphNEL format graph.
# Copied from the development code for igraph 
# <a href="http://bazaar.launchpad.net/~igraph/igraph/0.6-main/annotate/head:/interfaces/R/igraph/R/conversion.R">http://bazaar.launchpad.net/~igraph/igraph/0.6-main/annotate/head:/interfaces/R/igraph/R/conversion.R</a> 
# on 26-Jan-2010
igraph.to.graphNEL &lt;- function(graph) {
  if (!is.igraph(graph)) 
  {
    stop(&quot;Not an igraph graph&quot;)
  }

  require(graph)

  if (&quot;name&quot; %in% list.vertex.attributes(graph) &amp;&amp;
     is.character(V(graph)$name)) {
    name &lt;- V(graph)$name
  } else {
    name &lt;- as.character(seq(vcount(graph))-1)    
  }

  edgemode &lt;- if (is.directed(graph)) &quot;directed&quot; else &quot;undirected&quot;  

  if (&quot;weight&quot; %in% list.edge.attributes(graph) &amp;&amp;
      is.numeric(E(graph)$weight)) {
    al &lt;- get.adjedgelist(graph, &quot;out&quot;)
    for (i in seq(along=al)) {
      edges &lt;- get.edges(graph, al[[i]])
      edges &lt;- ifelse( edges[,2]==i-1, edges[,1], edges[,2])
      weights &lt;- E(graph)$weight[al[[i]]+1]
      al[[i]] &lt;- list(edges=edges+1, weights=weights)
    }
  } else {
    al &lt;- get.adjlist(graph, &quot;out&quot;)
    al &lt;- lapply(al, function(x) list(edges=x+1))
  }  
  
  names(al) &lt;- name

  res &lt;- new(&quot;graphNEL&quot;, nodes=name, edgeL=al, edgemode=edgemode)

  ## Add graph attributes (other than &#39;directed&#39;)

  ## Are this &quot;officially&quot; supported at all?

  g.n &lt;- list.graph.attributes(graph)

  if (&quot;directed&quot; %in% g.n) {
    warning(&quot;Cannot add graph attribute `directed&#39;&quot;)
    g.n &lt;- g.n[ g.n != &quot;directed&quot; ]
  }
  for (n in g.n) {
    res@graphData[[n]] &lt;- get.graph.attribute(graph, n)
  }

  ## Add vertex attributes (other than &#39;name&#39;, that is already
  ## added as vertex names)

  v.n &lt;- list.vertex.attributes(graph)

  v.n &lt;- v.n[ v.n != &quot;name&quot; ]

  for (n in v.n) {
    nodeDataDefaults(res, attr=n) &lt;- NA
    nodeData(res, attr=n) &lt;- get.vertex.attribute(graph, n)
  }

  ## Add edge attributes (other than &#39;weight&#39;)

  e.n &lt;- list.edge.attributes(graph)
  e.n &lt;- e.n[ e.n != &quot;weight&quot; ]

  if (length(e.n) &gt; 0) {
    el &lt;- get.edgelist(graph)
    el &lt;- paste(sep=&quot;|&quot;, el[,1], el[,2])
    for (n in e.n) {
      edgeDataDefaults(res, attr=n) &lt;- NA
      res@edgeData@data[el] &lt;- mapply(function(x,y) {
        xx &lt;- c(x,y); names(xx)[length(xx)] &lt;- n; xx },
                                      res@edgeData@data[el],
                                      get.edge.attribute(graph, n),
                                      SIMPLIFY=FALSE)
    }
  }
  
  res

}

# Function to generate X random sequences with a multinomial model, where the
# probabilities of the different letters are set equal to their frequencies
# in an input sequence &quot;inputsequence&quot;. The input sequence comes in as a string
# of characters, eg. &quot;ATSTWCWYSKLAV&quot; 
generateSeqsWithMultinomialModel &lt;- function(inputsequence, X)
{
   # Change the input sequence into a vector of letters
   library(&quot;seqinr&quot;) # Load in the SeqinR library, so we can use function &quot;s2c&quot;.
   inputsequencevector &lt;- s2c(inputsequence)
   # Find the frequencies of the letters in the input sequence &quot;inputsequencevector&quot;:
   mylength &lt;- length(inputsequencevector)
   mytable &lt;- table(inputsequencevector)
   # Find the names of the letters in the sequence
   letters &lt;- rownames(mytable)
   numletters &lt;- length(letters)
   probabilities &lt;- numeric() # Make a vector to store the probabilities of letters
   for (i in 1:numletters)
   {
      letter &lt;- letters[i]
      count &lt;- mytable[[i]]
      probabilities[i] &lt;- count/mylength
   } 
   # Make X random sequences using the multinomial model with probabilities &quot;probabilities&quot;
   seqs &lt;- numeric(X)
   for (j in 1:X)
   {
      seq &lt;- sample(letters, mylength, rep=TRUE, prob=probabilities) # Sample with replacement
      seq &lt;- c2s(seq)
      seqs[j] &lt;- seq
   }
   # Return the vector of random sequences
   return(seqs)
}

# Function to convert a phylip format tree to an &#39;ape&#39; format tree.
# In phylip, the bootstrap values are given for branches.
# In &#39;ape&#39;, it expects the bootstrap values to be given for nodes.
phylip.to.ape &lt;- function(tree) 
{
   # Get the labels of the tips (sequences)
   seqnames &lt;- tree$tip.label
   numseqs &lt;- length(seqnames)
   # Get the number of internal nodes in the tree:
   numinternalnodes &lt;- tree$Nnode
   # Get the bootstrap values
   bootstraps &lt;- tree$edge.length
   # Get the edges for the tree
   edges &lt;- tree$edge
   # For each edge, find the node at the end of the edge (closer to the tips
   # of the tree):
   numedges &lt;- nrow(edges)
   numnodes &lt;- numinternalnodes + numseqs
   mybootstraps &lt;- numeric(numinternalnodes) 
   for (i in 1:numedges)           
   {
      startnode &lt;- edges[i,1]
      endnode &lt;- edges[i,2]
      # Make sure the end node is not a tip:
      if (endnode &gt; numseqs)   
      {
         # Get the bootstrap value for this node:
         bootstrap &lt;- bootstraps[i]
         mybootstraps[endnode] = bootstrap
      }
   }
   mybootstraps2 &lt;- mybootstraps[(numseqs+1):(numnodes)]
   tree$node.label &lt;- mybootstraps2

   return(tree)
}

plotproteinNJtreewithbootstraps &lt;- function(alignment, theoutgroup)
{
   # define a function for making a tree:
   makemytree &lt;- function(alignmentmat, outgroup=`theoutgroup`)
   {
      alignment &lt;- as.alignment(alignmentmat)
      mydist &lt;- dist.alignment(alignment)
      mytree &lt;- nj(mydist)
      mytree &lt;- makeLabel(mytree, space=&quot;&quot;) # get rid of spaces in tip names.
      myrootedtree &lt;- root(mytree, outgroup, r=TRUE)
      return(myrootedtree)   
   }
   # infer a tree
   mymat  &lt;- as.matrix.alignment(alignment)
   myrootedtree &lt;- makemytree(mymat, outgroup=theoutgroup) 
   # bootstrap the tree
   mymat  &lt;- as.matrix.alignment(alignment)
   myboot &lt;- boot.phylo(myrootedtree, mymat, makemytree)
   # plot the tree
   plot.phylo(myrootedtree)
   nodelabels(myboot)
}

# Function to plot the range of fluxes, and optimal fluxes for reactions.
# model is the input model, eg. LIMEcoli
# Code copied from <a href="http://cran.r-project.org/web/packages/LIM/vignettes/LIMecoli.pdf">http://cran.r-project.org/web/packages/LIM/vignettes/LIMecoli.pdf</a>
plotFluxes &lt;- function(model)
{
   library(&quot;LIM&quot;)
   LP &lt;- Linp(model)
   xr &lt;- Xranges(model)
   par(mfrow=c(1,2))
   nr &lt;- model$NUnknowns
   ii &lt;- 1:(nr/2)
   dotchart(LP$X[ii],xlim = range(xr),pch=16,cex=0.8)
   segments(xr[ii,1],1:nr,xr[ii,2],1:nr)
   ii &lt;- (nr/2+1):nr
   dotchart(LP$X[ii],xlim = range(xr),pch=16,cex=0.8)
   segments(xr[ii,1],1:nr,xr[ii,2],1:nr)
   mtext(side= 3, cex=1.5, outer = TRUE, line=-1.5, &quot;E coli Core Metabolism, optimal solution and ranges&quot;)
}


# Function to identify periodically expressed genes in microarray time series data.
# Code inspired by <a href="http://www.cs.tut.fi/~ahdesmak/robustperiodic/doc/caulobacter.R">http://www.cs.tut.fi/~ahdesmak/robustperiodic/doc/caulobacter.R</a>, using
# the GeneCycle R library:
# The input data is a matrix with one row per gene and one column per microarray.
findPeriodicGenes &lt;- function(expndata)
{
   library(&quot;GeneCycle&quot;)
   # Put the data in the format required by fdrtool() and fisher.g.test(), 
   # ie. one row per microarray and one column per gene,
   # and no missing values:
   expndata2 &lt;- na.omit(expndata) # Discard genes that have missing data (&quot;NA&quot; values)
   expndata3 &lt;- t(expndata2)      # Find the transpose of matrix &quot;expndata2&quot;

   # Calculate p-values using Fisher&#39;s g test
   pval.mgenes &lt;- fisher.g.test(expndata3)

   # Calculate q-values (ie. p-values that have been corrected for multiple testing, as we
   # are carrying out tests for lots of different genes at once):
   fdr.out &lt;- fdrtool(pval.mgenes, statistic=&quot;pvalue&quot;)

   # Print out the qvalues for genes that have qvalues of &lt; 0.05:
   numgenes &lt;- length(fdr.out$qval)
   genenames &lt;- numeric()
   qvalues &lt;- numeric()
   numgenestaken &lt;- 0
   for (i in 1:numgenes)
   {
      qval &lt;- fdr.out$qval[i]
      genename &lt;- colnames(expndata3)[i]
      if (qval &lt;= 0.05)
      {
         numgenestaken &lt;- numgenestaken + 1
         genenames[numgenestaken] &lt;- genename
         qvalues[numgenestaken] &lt;- qval   
      }
   } 
   # Return a list variable that contains the names of periodically expressed genes,
   # and the q-values:
   result &lt;- list(genenames,qvalues)

   return(result) 
}

# An R function to plot a heatmap for microarray data.
# The input data is the distance matrix &quot;dist&quot;.
# Inspired by code in Hahne et al &quot;Bioconductor Case Studies&quot; page 140.
plotHeatmap &lt;- function(dist)
{
   library(&quot;RColorBrewer&quot;)
   # Use red to correspond to high values and blue to low: 
   hmcol &lt;- colorRampPalette(brewer.pal(10, &quot;RdBu&quot;))(256)
   hmcol &lt;- rev(hmcol)
   # This uses hierarchical clustering:
   heatmap(as.matrix(dist), sym=TRUE, col=hmcol, distfun=function(x) dist(x, method=&quot;manhattan&quot;))
}

# An R function to plot the expression for all of the genes in a cluster.
# &quot;clusters&quot; is the result of k-means clustering, and num is the number of
# the cluster that we want to plot.
plotCluster &lt;- function(clusters, num, expndata)
{
   # Find all the genes in this cluster:
   genes &lt;- clusters$cluster[clusters$cluster == num]
   genes &lt;- names(genes)
   numgenes &lt;- length(genes)

   # Now make a plot with the expression of each gene across the samples:
   for (i in 1:numgenes)
   {
      gene &lt;- genes[i]
      genedata &lt;- expndata[`gene`,]
      if (i == 1)
      {
         ylab &lt;- paste(&quot;Cluster&quot;,num)
         plot(genedata, type=&quot;l&quot;, col=&quot;blue&quot;, ylab=ylab)
      }
      else
      {
         points(genedata, type=&quot;l&quot;, col=&quot;blue&quot;)
      }
   }
}

# An R function that calculates the transition matrix for a Markov model,
# taking the matrix of observed dinucleotide relative frequencies as its
# input. Each row of the input matrix specifies a base, and each column
# specifies the following base:
# eg.
# D &lt;- matrix(c(0.146, 0.052, 0.058, 0.089,
#               0.063, 0.029, 0.010, 0.056,
#               0.050, 0.030, 0.028, 0.051,
#               0.087, 0.047, 0.063, 0.140), byrow=TRUE, nrow=4)
makeTransitionMatrix &lt;- function(D)
{
   # Initialise the transition matrix:
   P &lt;- array(0,dim=c(4,4))

   # Calculate the transition matrix: 
   for (i in 1:4) # For each row (i) of the transition matrix:
   {
      for (j in 1:4) # For each column (j) of the transition matrix:
      {
         value &lt;- D[i,j]/sum(D[i,]) # sum(D[i,]) is the sum of the ith row of the input matrix
         P[i,j] &lt;- value 
      }
   } 
   
   return(P)
}

