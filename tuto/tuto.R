nb_sliding_windows <- function(l, size, shift){
    tmp = size
    counter = 1
    while(tmp < l){
        tmp = tmp + shift
        counter = counter + 1
    }
    counter
}

x = nb_sliding_windows(44,11,2)
x

