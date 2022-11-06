library(ape)

dlAndSaveData = function(){
    accessionsNumbers = readLines("./Projet2/AccessionNumbers.txt")
    for(i in accessionsNumbers){
        print(i)
        mat = read.GenBank(i,as.character = TRUE)
        seq = mat[[1]]
        # write.dna(mat, file = paste(c("./Projet2/data/", i, ".gb"), collapse = ""))
        print(length(seq))
    }
}

dlAndSaveData()

# for(i in filenames){
#     print(i)
#     mat = read.dna(paste(c("./Projet2/data/", i), collapse = ""), as.character = TRUE)
#     seq = mat[[1]]
#     print(length(mat))
# }
