library(stringr)
library(dplyr)
library(ranger)


dogA=1
dogB=2
dogC=3
dogD=4
catA=1
catB=3
catC=2
catD=4
birdA=4
birdB=2
birdC=3
birdD=1

dogs=list(dogA=dogA,dogB=dogB,dogC=dogC,dogD=dogD)
cats=list(catA=catA,catB=catB,catC=catC,catD=catD)
birds=list(birdA=birdA,birdB=birdB,birdC=birdC,birdD=birdD)




dog.ind = sample(1:4,200,replace = T)
cat.ind = sample(1:4,200,replace = T)
bird.ind = sample(1:4,200,replace = T)

df <- cbind(dogs=names(dogs[dog.ind]),cats=names(cats[cat.ind]),birds=names(birds[bird.ind]),
            res = as.numeric(dogs[dog.ind])+as.numeric(cats[cat.ind])+as.numeric(birds[bird.ind])) |> as.data.frame()
df <- df |> mutate(dogs=as.factor(dogs),cats=as.factor(cats),birds=as.factor(birds), res=as.numeric(res))
df <- df |> mutate(res1 = ifelse(!birds %in% 'birdD', res, NA_real_))

rf <- ranger(res ~ dogs+cats+birds, data = df, sample.fraction = 1, mtry = 3, min.node.size = 1, respect.unordered.factors=T)
df$pred <-  predictions(predict(rf, df))

mean((df$res-df$pred)^2)

rf <- ranger(res1 ~ dogs+cats+birds, data = subset(df, !is.na(res1)), sample.fraction = 1, min.node.size = 1, respect.unordered.factors=T)
df$pred <-  predictions(predict(rf, df))

mean((df$res-df$pred)^2)


mod <- lm(res ~ dogs+cats+birds, data = df)
mod <- lm(res1 ~ dogs+cats+birds, data = df)


df$pred <-  predict(mod, df)

mean((df$res-df$pred)^2)

