source("./tuto/tuto.R")
library(ape)
library(ggplot2)

openFile = function(file){
    mat = read.dna(file, format="fasta", as.character = TRUE) 
    seq = as.vector(mat)
    seq
}

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
    ebola_vector = openFile(file)
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
    seq = openFile(file)

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


odss_ratio = function(file){
    seq = aopenFile(file)
    seqLength = length(seq)

    getNucleotideProb = function(nc){
        length(which(seq==nc)) / seqLength
    }

    m_odds = matrix(0, nrow = 4, ncol = 4)
    colnames(m_odds) = c("a","c","g","t")
    rownames(m_odds) = c("a","c","g","t")
    nucProb = c("a"=0, "c"=0, "g"=0, "t"=0)

    for(i in names(nucProb)){
        nucProb[i] = getNucleotideProb(i)
    }

    m = dimer_freq(file)

    for(i in names(nucProb)){
        for(j in names(nucProb)){
            m_odds[i, j] = nucProb[i] * nucProb[j]
        }
    }

     m = m / m_odds
     saveRDS(m, file="./Projet1/odd_ratios.rds")
     m

}

openReadingFrames = function(file, k=c("0"=0, "10"=0, "50"=0, "100"=0, "300"=0, "500"=0),
                            startCodons = c("ttg", "ctg", "ata", "att", "atc", "atg", "gtg"),
                            stopCodons =c("tga", "taa", "tag")){
    seq = openFile(file)
    lengthSeq = length(seq)
    rest1 = lengthSeq %% 3
    rest2 = (lengthSeq - 1) %% 3
    rest3 = (lengthSeq - 2) %% 3

    reverse = c("a"="t", "t"="a", "g"="c", "c"="g")
    # readings frames
    rf1 = seq[1:(lengthSeq-rest1)]
    rf2 = seq[2:(lengthSeq-rest2)]
    rf3 = seq[3:(lengthSeq-rest3)]
    rf4 = replace(rev(seq), TRUE, reverse[rev(seq)])[1:(lengthSeq-rest1)]
    rf5 = replace(rev(seq), TRUE, reverse[rev(seq)])[2:(lengthSeq-rest2)]
    rf6 = replace(rev(seq), TRUE, reverse[rev(seq)])[3:(lengthSeq-rest3)]

    orf_starts = c(0,0,0,0,0,0)
    status = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)

    for(i in seq(1,lengthSeq,3)){
        if(i <= length(rf1)){
            codon1 = paste(rf1[i:(i+2)], collapse = "")
            codon4 = paste(rf4[i:(i+2)], collapse = "")
        }
        if(i <= length(rf2)){
            codon2 = paste(rf2[i:(i+2)], collapse = "")
            codon5 = paste(rf5[i:(i+2)], collapse = "") 
        }
        if(i <= length(rf3)){
            codon3 = paste(rf3[i:(i+2)], collapse = "")
            codon6 = paste(rf6[i:(i+2)], collapse = "") 

        codons = c(codon1, codon2, codon3, codon4, codon5, codon6)

        for(c in 1:length(codons)){

            if(codons[c] %in% startCodons && status[c] == FALSE){
                status[c] = TRUE
                orf_starts[c] = i
            }
            if(codons[c] %in% stopCodons && status[c] == TRUE && orf_starts[c] != i-1){
                status[c] = FALSE
                tmp = i - orf_starts[c]

                for(j in names(k)){
                    if(strtoi(j) <= tmp){
                        k[j] = k[j] + 1
                    }
                }
            }
        }
    }
    }
    saveRDS(k,"./Projet1/data/orf.rds")
    print(k)
}


#ORF()

openReadingFrames("./Projet1/data/bacterial_sequence.fasta")