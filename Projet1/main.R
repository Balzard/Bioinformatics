source("./tuto/tuto.R")
library(ape)
library(ggplot2)

gc_content <- function(sequence, size, shift){
    solution = vector(length=nb_sliding_windows(length(sequence), size, shift))
    start = 1
    end = size

    for (i in 1:length(solution)){
        current_window = sequence[start:end]
        solution[i] = count(current_window, "g", "c") / length(current_window)
        start = start + shift
        end = if (i==(length(solution) - 1)) length(sequence) else end + shift
    }
    solution
}

# v = c("G", "T", "G", "A", "G", "C", "C", "G", "A", "G", "T", "G", "A", "C", "T", "C", "C", "A", "A", "T", "T", 
# "T", "G", "G", "A", "A", "A", "T", "A", "C", "T", "C", "C", "T", "C", "C", "G", "A")

# s = gc_content(v, 6, 3)
# s

mat = read.dna("./Projet1/ebola.fasta", format="fasta", as.character = TRUE)
l = as.vector(mat)
s = gc_content(l, 500, 250)
plot(s, type="l", col="red")

