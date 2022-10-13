source("./tuto/tuto.R")
library(ape)
library(ggplot2)

gc_content <- function(sequence, size, shift, nc1, nc2){
    solution = vector(length=nb_sliding_windows(length(sequence), size, shift))
    start = 1
    end = size

    for (i in 1:length(solution)){
        current_window = sequence[start:end]
        solution[i] = count(current_window, nc1, nc2) / length(current_window)
        start = start + shift
        end = if (i==(length(solution) - 1)) length(sequence) else end + shift
    }
    solution
}


plot_freq = function(file){
    mat = read.dna(file, format="fasta", as.character = TRUE) #./Projet1/ebola.fasta
    ebola_vector = as.vector(mat)
    gc_freq = gc_content(ebola_vector, 500, 250, "g", "c")
    at_freq = gc_content(ebola_vector, 500, 250, "a", "t")
    plot(at_freq, type="l", col="green", ylim=c(0.3,0.8), ylab="Freq")
    lines(gc_freq, type="l", col="red", lty=2)
    legend("topright", legend = c("A-T", "C-G"), col=c("green","red"), lty=1:2, cex=1)
}


ex_1.3.2 = function(){
    l1 = c(0.0975, 0.0589, 0.0489, 0.0863)
    l2 = c(0.0660, 0.0422, 0.0503, 0.0506)
    l3 = c(0.0626, 0.0447, 0.0430, 0.0588)
    l4 = c(0.0656, 0.0634, 0.0668, 0.0945)
    m = matrix(rbind(l1,l2,l3,l4), nrow = 4, ncol = 4)
    colnames(m) = c("*A", "*C", "*G", "*T")
    rownames(m) = c("A*", "C*", "G*", "T*")
    saveRDS(m, file="./Projet1/matrix_132.rds")
    m
}


dimer_freq = function(file){
    m = matrix(0, nrow = 4, ncol = 4)
    colnames(m) = c("*A", "*C", "*G", "*T")
    rownames(m) = c("A*", "C*", "G*", "T*")
    seq = read.dna(file, format="fasta", as.character = TRUE)
    seq = as.vector(seq)

    dic = c("a"=1, "c"=2, "g"=3, "t"=4)

    for(i in 1:length(seq) - 1){
        m[dic[seq[i]], dic[seq[i+1]]] = m[dic[seq[i]],dic[seq[i+1]]] + 1
    }

    computerFreq = function(x, l=length(seq)){
        return(x/l)
    }

    m = apply(m,2, computerFreq)
    saveRDS(m, file="./Projet1/matrix133.rds")
    m

}

ex_1.3.2()
dimer_freq("./Projet1/bovine.fasta")

