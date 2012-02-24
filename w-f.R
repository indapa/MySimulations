#simulate Wright Fisher model
#see this for more info: http://www.genetics.wustl.edu/bio5488/lecture_notes_2004/popgen_1.pdf


wright_fisher <- function(N,m,j) {

#Preliminaries
#number of chromosomes
N=N


#number of generations
m = m
x = numeric(m) 

#number of minor alleles in the first generation
x[1] = j


#Simulation
#binomial sampling at each generation - determines allele count in next generations
for (i in 2:m)
{
k=(x[i-1])/N
n=seq(0,N,1)
prob=dbinom(n,N,k)
x[i]=sample(0:N, 1, prob=prob)

}

plot( (x[1:m]/N), type="o", pch=19, xlab="Generation", ylim=c(0,1), ylab="allele frequency")
}


