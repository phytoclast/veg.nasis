
#input data
#point intercepts
p.L <- 38+33
p.B <- 1
p.T <- 1
p.w0_1 <- 0
p.w1_2.5 <- 0
p.w2.5_7.5 <- 1
p.w7.5_25 <- 2
#line intercepts
l.w7.5_25 <- c(15,12,17,15)
l.w25 <- c(50)

#totals for point and line transects
#total count of points
p.total <- p.L+p.B+p.T+p.w0_1+p.w1_2.5+p.w2.5_7.5+p.w7.5_25
#total linear cm measured in line intercepts
l.total <- 2000+2000

#cover for larger woody debris in line intercept
l.w7.5_25 <- sum(l.w7.5_25)/l.total
l.w25 <- sum(l.w25)/l.total

#extra total due to line intercept
l.total2 <- l.w7.5_25+l.w25

#recalculate intercepts based on total points and discounting for the portion of ground covered by larger wood that wasn't accounted for in point intercepts.
x.L <- round(p.L/p.total/(1+l.total2)*100,1)
x.B <- round(p.B/p.total/(1+l.total2)*100,1)
x.T <- round(p.T/p.total/(1+l.total2)*100,1)
w0_1 <- round(p.w0_1/p.total/(1+l.total2)*100,1)
w1_2.5 <- round(p.w1_2.5/p.total/(1+l.total2)*100,1)
w2.5_7.5 <- round(p.w2.5_7.5/p.total/(1+l.total2)*100,1)

#add together estimates of this size class assessed by both point and line intercept because we opted to use 10 cm cutoff which is between the 7.5 and 25 cm thresholds
w7.5_25 <- round((p.w7.5_25/p.total/(1+l.total2)+l.w7.5_25)*100,1)
w25 <- round(l.w25*100,1)

#check totals
totals <- x.L+x.B+x.T+w0_1+w1_2.5+w2.5_7.5+w7.5_25+w25



