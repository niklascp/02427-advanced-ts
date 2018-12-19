library(ctsmr)
library(splines)
library(ggplot2)
library("GGally")
source("sdeTiTm.R")
load("Exercise3.RData")
fit1 <- sdeTiTm(AllDat,AllDat$yTi1,AllDat$Ph1)

summary(fit1,extended=TRUE)

Hour <- as.numeric(strftime(AllDat$date, format="%H"))

pred1 <- predict(fit1)
mean((pred1[[1]]$state$pred$Ti - AllDat$yTi1)^2)
# What is going on 10 AM?
# Try to fir a varying effective window area

df2 <- AllDat[, c('yTi1', 'yTi2', 'yTi3', 'yTi4', 'Ta', 'Gv', 'Ph1', 'Ph2')] 
p <- ggpairs(df2)
ggsave('latex/part2a-corr.png', plot = p, width = 7, height = 7)

p <- ggplot() +
    geom_point(aes(x = Hour, y = pred1[[1]]$state$pred$Ti - AllDat$yTi1)) +
    ylab('Error')
ggsave('latex/part2a-fit1-hour.png', plot = p, width = 7, height = 7)

idx <- (Hour>8 & Hour < 23) # It is impossible to fit a window area for the hours without any sun, so we limit the window area estimation to the hours with sun.
bs = bs(Hour[idx],df=5,intercept=TRUE) 
lengt(seq(9,22))
# What does the splines look like?
plot(bs[14:27,1],type='l')
lines(bs[ 14:27,2])
lines(bs[ 14:27,3])
lines(bs[ 14:27,4])
lines(bs[ 14:27,5])
nrow(bs)
bs1 <- bs2 <- bs3 <- bs4 <- bs5 <- bs6 <- numeric(dim(AllDat)[1])

# Idea: Use azimuth of the sun

bs1[idx] = bs[,1]
bs2[idx] = bs[,2]
bs3[idx] = bs[,3]
bs4[idx] = bs[,4]
bs5[idx] = bs[,5]

AllDat$bs1 = bs1
AllDat$bs2 = bs2
AllDat$bs3 = bs3
AllDat$bs4 = bs4
AllDat$bs5 = bs5


### You will have to implement sdeTITmAv ###
source("sdeTiTmAv.R")
fit2 <- sdeTiTmAv(AllDat,AllDat$yTi1,AllDat$Ph1)

plot(bs[14:27,1]*fit2$xm[3]+bs[14:27,2]*fit2$xm[4]+bs[14:27,3]*fit2$xm[5]+bs[14:27,4]*fit2$xm[6]+bs[14:27,5]*fit2$xm[7],type='l')

str(fit2$xm)
names(fit2$xm)

summary(fit2,extended=TRUE)

pred2 <- predict(fit2)
mean((pred2[[1]]$state$pred$Ti - AllDat$yTi1)^2)
p <- ggplot() +
    geom_point(aes(x = Hour, y = pred2[[1]]$state$pred$Ti - AllDat$yTi1)) +
    ylab('Error')
ggsave('latex/part2a-fit2-hour.png', plot = p, width = 7, height = 7)
