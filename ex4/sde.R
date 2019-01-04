library('readr')
library('data.table')
library('ggplot2')
library('ctsmr')
library('forecast')
library('zoo')
library('xts')

# Load data from CSV and add derived attributes
sde_data <- read_csv('sde_data.csv')
sde_data['Time'] <- sde_data['Time'] + 2*60*60
sde_data['Day'] <- parse_number(format(sde_data$Time, "%d"))
sde_data['Hour'] <- parse_number(format(sde_data$Time, "%H"))
sde_data['Minute'] <- parse_number(format(sde_data$Time, "%M"))
sde_data['Second'] <- parse_number(format(sde_data$Time, "%S"))
sde_data['TimeOfDay'] <- sde_data['Hour'] * 60 * 60 + sde_data['Minute'] * 60.0 + sde_data['Second']
sde_data['t'] <- (sde_data$Day - 1) * 24 * 60 * 60 + sde_data['Hour'] * 60 * 60 + sde_data['Minute'] * 60.0 + sde_data['Second']
setDT(sde_data)
sde_data <- sde_data[order(t), ]
sde_data[, dt := c(NA, diff(t)), by=LinkRef]

sde_data_day <- sde_data[8 <= Day & Day <= 19, ]

# Add link number
lns = c(
    '448963789:1242123649:4459313610',
    '448963780:4459313610:1280331077',
    '448963780:1280331077:2599647834',
    '448963787:2599647834:4459313608')

# Fit LOESS Model from Time Of Day for each Link using only train data (Day 1 - 10).
sde_data_day$ln <- match(sde_data_day$LinkRef, lns)
sde_data_day$Trend_TodLoess10 <- 0
sde_data_day$Trend_TodLoess25 <- 0
sde_data_day$Trend_TodLoess50 <- 0

for (ln in lns) {
    ix = sde_data_day$LinkRef == ln

    loessMod10 <- loess(LinkRunTime ~ TimeOfDay, data=sde_data_day[ix][Day <= 12, ], span=0.10)
    loessMod25 <- loess(LinkRunTime ~ TimeOfDay, data=sde_data_day[ix][Day <= 12, ], span=0.25)
    loessMod50 <- loess(LinkRunTime ~ TimeOfDay, data=sde_data_day[ix][Day <= 12, ], span=0.50)
    sde_data_day[ix]$Trend_TodLoess10 <- predict(loessMod10, newdata = sde_data_day[ix]) 
    sde_data_day[ix]$Trend_TodLoess25 <- predict(loessMod25, newdata = sde_data_day[ix]) 
    sde_data_day[ix]$Trend_TodLoess50 <- predict(loessMod50, newdata = sde_data_day[ix]) 
}

# Remove any na
sde_data_day = sde_data_day[
    !is.na(dt) &
    !is.na(Trend_TodLoess10) &
    !is.na(Trend_TodLoess25) &
    !is.na(Trend_TodLoess50), ]

round(sde_data_day[, list(mean(LinkRunTime), sd(LinkRunTime)), by = ln], 2)

ggplot(sde_data_day[ln == 3 & Day <= 12, ], aes(x=t)) +
    #geom_point(aes(y=LinkRunTime)) +
    geom_line(aes(y=Trend_TodLoess10), color = 'red') +
    geom_line(aes(y=Trend_TodLoess25), color = 'green') +
    geom_line(aes(y=Trend_TodLoess50), color = 'blue') 

ggplot(sde_data_day[ln == 3, ], aes(x=TimeOfDay, y=LinkRunTime)) +
    geom_point() +
    geom_smooth() +
    ylim(0, 100)

ggplot(sde_data_day[Day == 1, ], aes(TimeOfDay, LinkRunTime, colour = as.factor(ln))) +
    geom_smooth()

ggplot(sde_data_day[ln == 3 & Day <= 50, ], aes(x=as.factor(Day), y=LinkRunTime)) +
    geom_boxplot()

## TODO: Save this plot for paper
for (i in 1:4) {
    p <- ggPacf(sde_data_day[ln == i, LinkRunTime]) +
        ggtitle(paste0('Link ', i))

    ggsave(paste0(paste0('pacf-', i), '.pdf'), plot = p, width = 8, height = 3)
}

# Select days and links for train and test and remove measurements taken exactly 
# at the same time for the same link (yes, it happens!)
data_train <- sde_data_day[dt > 0 & ln == 1 & 8 <= Day & Day <= 12]
data_test <- sde_data_day[dt > 0 & ln == 1 & 15 <= Day & Day <= 19]

data_train_ts <- zoo(data_train$LinkRunTime, data_train$t)
data_test_ts <- zoo(data_test$LinkRunTime, data_test$t)

# Exponential smoothning (right alligned to avoid cheating)
tau <- 3600
es_train <- Vectorize(function(u) sum(ifelse(index(data_train_ts)<=u,coredata(data_train_ts) * exp((index(data_train_ts)-u)/tau),0)) / sum(ifelse(index(data_train_ts)<=u,exp((index(data_train_ts)-u)/tau),0)))
es_test <- Vectorize(function(u) sum(ifelse(index(data_test_ts)<=u,coredata(data_test_ts) * exp((index(data_test_ts)-u)/tau),0)) / sum(ifelse(index(data_test_ts)<=u,exp((index(data_test_ts)-u)/tau),0)))

data_train$LinkRunTime_ <- es_train(data_train$t)
data_test$LinkRunTime_ <- es_test(data_test$t)

# Inspect ES model
ggplot(melt(data_train, id.vars="t", measure.vars=c("LinkRunTime_", "Trend_TodLoess10", "Trend_TodLoess25", "Trend_TodLoess50")), aes(x = t, y = value, color = variable)) +
    geom_line()

ggplot(melt(data_test, id.vars="t", measure.vars=c("LinkRunTime_", "Trend_TodLoess10", "Trend_TodLoess25", "Trend_TodLoess50")), aes(x = t, y = value, color = variable)) +
    geom_line()

#
# MODELS
#
# Simple CAR Model
model1 <- ctsm()
model1$addSystem(dU ~ theta * ((mu - U)) * dt + exp(sigma)*dw1 )
model1$addObs(y ~ U)
model1$setVariance(y ~ exp(e))
model1$setParameter(
    U = c(init = 2,1,1000),
    theta = c(init = 1, lb = 1E-5, ub = 1E5),
    mu = c(init = 14,1,1000),
    sigma = c(init = 1, lb = 1E-5, ub = 1E5),
    e = c(init = 1, lb = 1E-5, ub = 1E5))

# Run the parameter optimization
model1_data <- data.frame(data_train[, list(t, y = LinkRunTime)])
model1_fit <- model1$estimate(model1_data)
summary(model1_fit, extended = TRUE)

# Evaluate the model on train data
model1_pred_train_ <- predict(model1_fit)
model1_data$y_pred <- model1_pred_train_[[1]]$state$pred$U
sqrt(mean((model1_data$y_pred - data_train$LinkRunTime)^2))
mean(abs(model1_data$y_pred - data_train$LinkRunTime))

# Evaluate the model on test data using sliding windows
results <- data.frame();

for (i in 5:(32*5)) {
    print(paste0('Predicting Window #', i))

    w <- 45 * 60 # Window size in seconds
    t_min <- min((data_test$Day - 1) * 24 * 60 * 60) + (i - 1) * w
    t_max <- min((data_test$Day - 1) * 24 * 60 * 60) + i * w

    # Slice data into bootstrap (known y) and evaluation (unknown y) windows, and combine.
    model1_test_boostrap <- data_test[t < t_min, list(t, y = LinkRunTime)]
    model1_test_window <- data_test[t_min <= t & t < t_max, list(t, y = NA)]
    model1_test_data <- rbind(model1_test_boostrap, model1_test_window)

    model1_test_pred_ <- predict(model1_fit, newdata = model1_test_data)
    model1_test_data$y_pred <- model1_test_pred_[[1]]$pred$U
    model1_test_data$t_ahead <- model1_test_data$t - t_min
    model1_test_data$window <- i
    model1_test_data[model1_test_data$t_ahead >= 0]$y <- data_test[t_min <= t & t < t_max, LinkRunTime]
    results <- rbind(results, model1_test_data[model1_test_data$t_ahead >= 0])
}
round(sqrt(mean((results$y_pred - results$y)^2)), 2)
round(mean(abs(results$y_pred - results$y)), 2)

ggplot(melt(results, id.vars="t", measure.vars=c("y", "y_pred")), aes(x = t, y = value, color = variable)) +
    geom_line()

# ====================================================

# Sesonal CAR Model
data_train[, dtod := c(0, diff(Trend_TodLoess50) / diff(t))]
data_test[, dtod := c(0, diff(Trend_TodLoess50) / diff(t))]

ggplot(melt(data_train, id.vars="t", measure.vars=c("LinkRunTime_", "Trend_TodLoess50", "dtod")), aes(x = t, y = value, color = variable)) +
    geom_line() +
    facet_grid(rows = vars(variable), scales = "free_y") 

# Remove/reduce disconituity between days
data_train[dtod > .001, ]$dtod <- 0 
data_train[dtod < -.001, ]$dtod <- 0 
data_test[dtod > .001, ]$dtod <- 0 
data_test[dtod < -.001, ]$dtod <- 0 

ggplot(melt(data_train, id.vars="t", measure.vars=c("LinkRunTime_", "Trend_TodLoess50", "dtod")), aes(x = t, y = value, color = variable)) +
    geom_line() +
    facet_grid(rows = vars(variable), scales = "free_y") 

ggsave('dtod-1.pdf', width = 10, height = 12)

model2 <- ctsm()
model2$addSystem(dU ~ (theta1 * (mu - U) + theta2 * dtod) * dt + exp(sigma)*dw1 )
model2$addInput(dtod)
model2$addObs(y ~ U)
model2$setVariance(y ~ exp(e))
model2$setParameter(
    U = c(init = 2,1,1000),
    theta1 = c(init = 1, lb = 1E-5, ub = 1E5),
    theta2 = c(init = 1, lb = 1E-5, ub = 1E5),
    mu = c(init = 14,1,1000),
    sigma = c(init = 1, lb = 1E-5, ub = 1E5),
    e = c(init = 1, lb = 1E-5, ub = 1E5))

# Run the parameter optimization
model2_data <- data.frame(data_train[, list(t, dtod, y = LinkRunTime)])
model2_fit <- model2$estimate(model2_data)
summary(model2_fit, extended = TRUE)

# Evaluate the model on train data
model2_pred_train_ <- predict(model2_fit)
model2_data$y_pred <- model2_pred_train_[[1]]$state$pred$U
sqrt(mean((model2_data$y_pred - model2_data$y)^2))
mean(abs(model2_data$y_pred - model2_data$y))

# Plot evalution on train data
ggplot(melt(model2_data, id.vars="t", measure.vars=c("y", "y_pred")), aes(x = t, y = value, color = variable)) +
    geom_line()

# Evaluate the model on test data using sliding windows
results <- data.frame();

for (i in 5:(32*5)) {
    print(paste0('Predicting Window #', i))

    w <- 45 * 60 # Window size in seconds
    t_min <- min((data_test$Day - 1) * 24 * 60 * 60) + (i - 1) * w
    t_max <- min((data_test$Day - 1) * 24 * 60 * 60) + i * w

    # Slice data into bootstrap (known y) and evaluation (unknown y) windows, and combine.
    model2_test_boostrap <- data_test[t < t_min, list(t, dtod = dtod, y = LinkRunTime)]
    model2_test_window <- data_test[t_min <= t & t < t_max, list(t, dtod = dtod, y = NA)]
    model2_test_data <- rbind(model2_test_boostrap, model2_test_window)

    model2_test_pred_ <- predict(model2_fit, newdata = model2_test_data)
    model2_test_data$y_pred <- model2_test_pred_[[1]]$pred$U
    model2_test_data$t_ahead <- model2_test_data$t - t_min
    model2_test_data$window <- i
    model2_test_data[model2_test_data$t_ahead >= 0]$y <- data_test[t_min <= t & t < t_max, LinkRunTime]
    results <- rbind(results, model2_test_data[model2_test_data$t_ahead >= 0])
}
round(sqrt(mean((results$y_pred - results$y)^2)), 2)
round(mean(abs(results$y_pred - results$y)), 2)

ggplot(melt(results, id.vars="t", measure.vars=c("y", "y_pred")), aes(x = t, y = value, color = variable)) +
    geom_line()

# ====================================================

# Test baselines
round(sqrt(mean((data_test$Trend_TodLoess10 - data_test$LinkRunTime)^2)), 2)
round(sqrt(mean((data_test$Trend_TodLoess25 - data_test$LinkRunTime)^2)), 2)
round(sqrt(mean((data_test$Trend_TodLoess50 - data_test$LinkRunTime)^2)), 2)

round(mean(abs((data_test$Trend_TodLoess10 - data_test$LinkRunTime))), 2)
round(mean(abs((data_test$Trend_TodLoess25 - data_test$LinkRunTime))), 2)
round(mean(abs((data_test$Trend_TodLoess50 - data_test$LinkRunTime))), 2)

# Evaluate the model on test data using sliding windows
results_es <- data.frame();

for (i in 5:(32*5)) {
    print(paste0('Predicting Window #', i))

    w <- 45 * 60 # Window size in seconds
    t_min <- min((data_test$Day - 1) * 24 * 60 * 60) + (i - 1) * w
    t_max <- min((data_test$Day - 1) * 24 * 60 * 60) + i * w

    # Slice data into bootstrap (known y) and evaluation (unknown y) windows, and combine.
    model_es_test_boostrap <- data_test[t < t_min, list(t, y = LinkRunTime_)]
    model_es_test_window <- data_test[t_min <= t & t < t_max, list(t, y = LinkRunTime)]
    model_es_test_window$y_pred <- tail(model_es_test_boostrap$y, n=1)
    model_es_test_window$t_ahead <- model_es_test_window$t - t_min
    model_es_test_window$window <- i
    results_es <- rbind(results_es, model_es_test_window)
}
round(sqrt(mean((results_es$y_pred - results_es$y)^2)), 2)
round(mean(abs(results_es$y_pred - results_es$y)), 2)

# Evaluate the HA model on test data using train for estimating the HA over time of day in 15 minute blocks.
ha_train_ts <- zoo(data_train$LinkRunTime, as.POSIXct(data_train$Time))
ha_train_xts <- as.xts(ha_train_ts)
ha_train_15min <- align.time(ha_train_xts, n = 15 * 60)

model_ha_ <- data.table(t = format(index(ha_train_15min), '%H:%M'), y = coredata(ha_train_15min)[,1]);
model_ha <- model_ha_[, list(y = mean(y)), by = t]
setkey(model_ha, t)

ha_test_ts <- zoo(data_test$LinkRunTime, as.POSIXct(data_test$Time))
ha_test_xts <- as.xts(ha_test_ts)
ha_test_15min <- align.time(ha_test_xts, n = 15 * 60)

ha_test <- data.table(t = data_test$t, t_ha = format(index(ha_test_15min), '%H:%M'), y = coredata(ha_test_15min)[,1]);
ha_test$y_pred <- model_ha[ha_test$t_ha]$y

round(sqrt(mean((ha_test$y_pred - ha_test$y)^2, na.rm = TRUE)), 2)
round(mean(abs(ha_test$y_pred - ha_test$y), na.rm = TRUE), 2)