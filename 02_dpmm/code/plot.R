#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

dat <- read.table(args[1],header=F)
x <- dat[,1]
y <- dat[,2]

pdf(args[2])
  plot(x,y)
dev.off()
