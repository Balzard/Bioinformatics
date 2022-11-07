library(ape)

dlAndSaveData = function(){
    accessionsNumbers = c(readLines("./Projet2/AccessionNumbers.txt"))
    mat = read.GenBank(accessionsNumbers)
    write.dna(mat, file = "./Projet2/data.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10)
}

sequences = read.dna(file = "./Projet2/data.fasta", format = "fasta")
print(str(sequences))
