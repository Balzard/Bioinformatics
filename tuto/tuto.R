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

count <- function (s, c) { # sequence, character
    length(which(s == c))
}


frequency_in_sliding_window <- function(sequence, ntd, size, shift) {
    vec = vector(length = nb_sliding_windows(length(sequence), size, shift))
    start = 1
    end = size
    for(i in 1:length(vec)) {
        if(i==length(vec)) {
            vec[i] = count(sequence[start:length(sequence)], ntd) / length(sequence[start:end])
        }
        else {
            vec[i] = count(sequence[start:end], ntd) / length(sequence[start:end])
        }
        start = start + shift
        end = end + shift
    
    }
    vec
}

v = c("C", "A", "T", "A", "C", "C", "T", "T", "C", "C", "A", "A", "A")
x = frequency_in_sliding_window(v, "T", 6, 6)
x

