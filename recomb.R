chromInfo <- read.table("/Users/amit/data/RecombSims/chromInfo.bed", header=F)
names(chromInfo) <- c("chrom", "start", "end")
#1 Morgan is 100 Mbp
morgan <- 1e8

#crossovers are Poisson distributed, lambdas are the genetic length of the chromosome
lambdas <- chromInfo$end/morgan

#I learned sapply!
ncrossovers <- sapply(lambdas, function(x) rpois(1,x))
write("#crossover points", file="crossovers.txt")
for (i in 1:length(lambdas))
{

    points <- runif(ncrossovers[i],min=0, max=chromInfo$end[i])
    points <- sort(points)
    if (length(points) == 0) { 
          write ("None", file="crossovers.txt", sep="\t", append=T) } 
    else {
         write(points, file="crossovers.txt", sep="\t", append=T)
     }
    
}
