library(ape)

nb_sliding_windows <- function(l, size, shift) {
    tmp <- size
    counter <- 1
    while (tmp < l) {
        tmp <- tmp + shift
        counter <- counter + 1
    }
    counter
}

count <- function (s, c1, c2) { # sequence, character
    length(which(s == c1 | s == c2))
}


frequency_in_sliding_window <- function(sequence, ntd, size, shift) {
    vec = vector(length = nb_sliding_windows(length(sequence), size, shift))
    start = 1
    end = size
    for(i in 1:length(vec)) {
        vec[i] = count(sequence[start:end], ntd) / length(sequence[start:end])
        start = start + shift
        end = if (i==(length(vec) - 1)) length(sequence) else end + shift
    
    }
    vec
}

v = c("G", "T", "G", "A", "G", "C", "C", "G") #"A", "G", "T", "G", "A", "C", "T", "C", "C", "A", "A", "T", "T", 
# "T", "G", "G", "A", "A", "A", "T", "A", "C", "T", "C", "C", "T", "C", "C", "G", "A")
x = count(v, "G", "C")
x

