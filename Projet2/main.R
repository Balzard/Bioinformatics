library(ape)

dlAndSaveData = function(){
    accessionsNumbers = c(readLines("./AccessionNumbers.txt"))
    mat = read.GenBank(accessionsNumbers)
    write.dna(mat, file = "./data.fasta", format = "fasta", append = FALSE)
}

print(unlist(sequences[1]))


uniqchars <- function(x) unique(strsplit(x, "")[[1]]) 
matequal <- function(x, y)
  is.matrix(x) && is.matrix(y) && dim(x) == dim(y) && all(x == y)


sequences = read.dna(file = "data.fasta", format = "fasta")
x=0
for (i in 1:nrow(sequences)) {
  if(length(uniqchars(c2s(sequences[i,1:3822])))!=4){
    x=c(x,i)
    print(uniqchars(c2s(sequences[i,1:3822])))
  }
  a=i+1
  while(a<=100){
    if(matequal(sequences[i,1:3822],sequences[a,1:3822])){
      x=c(x,i)
      print("same")
    }
    a=a+1
    print(a)
  }
}
print(x)
sequences=sequences[-x,]
print(sequences)
