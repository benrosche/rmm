rm(list=ls())

library(stringr)

f <- formula(outcome ~ var1 + var2 + mm(id, mmc(var3, var4), mmw(pupils^exp(teacher*b))) + var1 + var2) 



c(id, mmc, mmw)

sort(el(strsplit(as.character(f), " \\+ ")[1]))


f <- outcome ~ var1 + var2 + mm(id, mmc(var3, var4), mmw(pupils^exp(teacher*b)))

x <- sort(el(strsplit(as.character(f)[3], " \\+ ")))[1]

id  = gsub("^mm\\((.*)\\)$", "\\1", x)
mmc = gsub(".*mmc\\((.*?)\\).*", "\\1", f)
mmw = gsub(".*mmw\\((.*?\\))\\).*", "\\1", x)

c(id, mmc, mmw)