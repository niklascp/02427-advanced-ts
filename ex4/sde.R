library('readr')
library('data.table')
library('ggplot2')
library('ctsmr')
library('forecast')

sde_data <- read_csv('sde_data.csv')
sde_data['Time'] <- sde_data['Time'] + 2*60*60
sde_data['Day'] <- parse_number(format(sde_data$Time, "%d"))
sde_data['Hour'] <- parse_number(format(sde_data$Time, "%H"))
sde_data['Minute'] <- parse_number(format(sde_data$Time, "%M"))
sde_data['Second'] <- parse_number(format(sde_data$Time, "%S"))
sde_data['TimeOfDay'] <- sde_data['Hour'] + sde_data['Minute'] / 60.0 + sde_data['Second'] / 60.0 / 60.0
sde_data['t'] <- sde_data$Day * 24 + sde_data$TimeOfDay
setDT(sde_data)

sde_data_day <- sde_data[4 <= Hour & Hour < 22, ]

# Add link number
lns = c(
    '448963789:1242123649:4459313610',
    '448963780:4459313610:1280331077',
    '448963780:1280331077:2599647834',
    '448963787:2599647834:4459313608')

sde_data_day$ln <- match(sde_data_day$LinkRef, lns)
sde_data_day$Trend_TodLoess10 <- 0
sde_data_day$Trend_TodLoess25 <- 0
sde_data_day$Trend_TodLoess50 <- 0

for (ln in lns) {
    ix = sde_data_day$LinkRef == ln

    # Fit LOESS Model from Time Of Day using only train data.
    loessMod10 <- loess(LinkRunTime ~ TimeOfDay, data=sde_data_day[ix][Day <= 10, ], span=0.10)
    loessMod25 <- loess(LinkRunTime ~ TimeOfDay, data=sde_data_day[ix][Day <= 10, ], span=0.25)
    loessMod50 <- loess(LinkRunTime ~ TimeOfDay, data=sde_data_day[ix][Day <= 10, ], span=0.50)
    sde_data_day[ix]$Trend_TodLoess10 <- predict(loessMod10, newdata = sde_data_day[ix]) 
    sde_data_day[ix]$Trend_TodLoess25 <- predict(loessMod25, newdata = sde_data_day[ix]) 
    sde_data_day[ix]$Trend_TodLoess50 <- predict(loessMod50, newdata = sde_data_day[ix]) 
}

# Remove any na
sde_data_day = sde_data_day[
    !is.na(sde_data_day$Trend_TodLoess10) &
    !is.na(sde_data_day$Trend_TodLoess25) &
    !is.na(sde_data_day$Trend_TodLoess50)]



ggplot(sde_data_day[ln == 3 & Day <= 2, ], aes(x=t)) +
    geom_point(aes(y=Speed)) +
    geom_line(aes(y=Trend_TodLoess10), color = 'red') +
    geom_line(aes(y=Trend_TodLoess25), color = 'green') +
    geom_line(aes(y=Trend_TodLoess50), color = 'blue') 

ggplot(sde_data_day[ln == 3, ], aes(x=TimeOfDay, y=LinkRunTime)) +
    geom_point() +
    geom_smooth() +
    ylim(0, 100)

ggplot(sde_data_day[Day == 1, ], aes(TimeOfDay, LinkRunTime, colour = LinkRef)) +
    geom_smooth()

ggplot(sde_data_day[ln == 3 & Day <= 5, ], aes(x=as.factor(Day), y=LinkRunTime)) +
    geom_boxplot()

ggPacf(sde_data_day[ln == 1, LinkRunTime], plot = FALSE)  %>% autoplot

data_train <- sde_data_day[ln == 3 & 8 <= Day & Day <= 8, ][order(t), ]
data_train2 <- sde_data_day[ln == 3 & 9 <= Day & Day <= 9, ][order(t), ]


ggplot(melt(data_train, id.vars="t", measure.vars=c("LinkRunTime", "Trend_TodLoess10", "Trend_TodLoess25", "Trend_TodLoess50")), aes(x = t, y = value, color = variable)) +
    geom_line()

min(data_train$LinkRunTime)
max(data_train$LinkRunTime)
mean(data_train$LinkRunTime)

#
# MODELS
#
model1 <- ctsm()
model1$addSystem(dU ~ theta * (mu - U) * dt + exp(sigma)*dw1 )
model1$addObs(y ~ U)
model1$setVariance(y ~ exp(e))
model1$setParameter(
    U = c(init = 14,1,150),
    theta = c(init = 1,-10,10),
    mu = c(init = 4,0,20),
    sigma = c(init = 1, lb = -30, ub = 10),
    e = c(init = 0, lb = -10, ub = 10))

# Run the parameter optimization
model1_data <- data.frame(t = data_train$t, y = data_train$LinkRunTime)
model1_data2 <- data.frame(t = data_train2$t, y = data_train2$LinkRunTime)
model1_fit <- model1$estimate(list(model1_data, model1_data2))
summary(model1_fit, extended = TRUE)

# Evaluate the model on train data
model1_pred_train_ <- predict(model1_fit)
model1_data$y_pred <- model1_pred_train_[[1]]$state$pred$U
sqrt(mean((model1_data$y_pred - model1_data$y)^2))
mean(abs(model1_data$y_pred - model1_data$y))


ggplot(melt(model1_data, id.vars="t", measure.vars=c("y", "y_pred")), aes(x = t, y = value, color = variable)) +
    geom_smooth()
    

model2 <- ctsm()
model2$addSystem(dU ~ theta * (mutod - U) * dt + exp(sigma)*dw1 )
model2$addInput(mutod)
model2$addObs(y ~ U)
model2$setVariance(y ~ exp(S))
model2$setParameter(
    U = c(init = 4,0,20),
    theta = c(init = 1,-10,10),
    sigma = c(init = 1, lb = -30, ub = 10),
    S = c(init= -30))

# Run the parameter optimization
model2_data <- data.frame(t = data_train$t, y = data_train$Speed, mutod = data_train$Trend_TodLoess25)
model2_fit <- model2$estimate(model2_data)
summary(model2_fit, extended = TRUE)

# Evaluate the model on train data
model2_pred_train_ <- predict(model2_fit)
model2_data$y_pred <- model2_pred_train_[[1]]$state$pred$U
mean((model2_data$y_pred - model2_data$y)^2)
mean(abs(model2_data$y_pred - model2_data$y))

ggplot(melt(model2_data, id.vars="t", measure.vars=c("y", "y_pred")), aes(x = t, y = value, color = variable)) +
    geom_smooth()
    

# Test
t_test = sde_data_day[ln == 3 & Day == 8, ][order(t), ]$t
y_test = sde_data_day[ln == 3 & Day == 8, ][order(t), ]$Speed
data_test <- data.frame(t = t_test, y = 0, true = y_test) 

pred_test_ <- predict(fit, newdata = data_test)
pred_test <- pred_test_[[1]]$pred$U
data_test$pred <- pred_test
mean(abs(pred_test - y_test))


sde_data_day[ln == 3 & Day == 8, ][order(t), ]

ix = t_test = sde_data_day[ln == 3 & Day == 8, ][order(t), ]$Hour == 5
data_test[ix]$y <- y_test[ix]

ggplot(data) +
    geom_smooth(aes_string(x = "t", y = "y"), color = 'red') +
    geom_smooth(aes_string(x = "t", y = "pred"), color = 'blue') 

ggplot(data_test) +
    geom_smooth(aes_string(x = "t", y = "true"), color = 'red') +
    geom_smooth(aes_string(x = "t", y = "pred"), color = 'blue') 

head(data, 1)