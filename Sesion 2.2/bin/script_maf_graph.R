x <- readLines("biallelic_freq_output.frq")

x <- x[-1]

freq <- strsplit(x, "\\s+")

maf <- sapply(freq, function(x) {
    p1 <- as.numeric(sub(".*:", "", x[5]))
    p2 <- as.numeric(sub(".*:", "", x[6]))
    min(p1, p2)
})

pdf("MAF_histogram.pdf")

hist(
    maf,
    main = "MAF spectrum",
    xlab = "MAF"
)

dev.off()
