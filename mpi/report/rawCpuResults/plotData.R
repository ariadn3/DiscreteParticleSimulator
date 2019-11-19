library(ggplot2)
library(dplyr)
library(tidyr)

SAVE = TRUE

# setwd("E:/Brownie/mpi/report/rawCpuResults")

# Extrapolation
seqBenchmarkData = read.csv("seqBenchmarkData.csv", stringsAsFactors = FALSE)
parBenchmarkData = read.csv("parBenchmarkData.csv", stringsAsFactors = FALSE)
I7seqExtrapolate = lm(wall.clock.time ~ poly(N, 2), data = seqBenchmarkData[seqBenchmarkData$configType=="PureI7",])
# summary(I7seqExtrapolate)$r.squared
XeSseqExtrapolate = lm(wall.clock.time ~ poly(N, 2), data = seqBenchmarkData[seqBenchmarkData$configType=="PureXeS",])
# summary(XeSseqExtrapolate)$r.squared
I7parExtrapolate = lm(wall.clock.time ~ poly(N, 2), data = parBenchmarkData[parBenchmarkData$configType=="PureI7",])
# summary(I7parExtrapolate)$r.squared
XeSparExtrapolate = lm(wall.clock.time ~ poly(N, 2), data = parBenchmarkData[parBenchmarkData$configType=="PureXeS",])
# summary(XeSparExtrapolate)$r.squared

I7seqExtrapolateFunc = function(x) { predict.lm(I7seqExtrapolate, data.frame(N=c(x))) }
I7parExtrapolateFunc = function(x) { predict.lm(I7parExtrapolate, data.frame(N=c(x))) }
XeSseqExtrapolateFunc = function(x) { predict.lm(XeSseqExtrapolate, data.frame(N=c(x))) }
XeSparExtrapolateFunc = function(x) { predict.lm(XeSparExtrapolate, data.frame(N=c(x))) }
I7speedupExtrapolateFunc = function(x) { I7seqExtrapolateFunc(x) / I7parExtrapolateFunc(x) }
XeSspeedupExtrapolateFunc = function(x) { XeSseqExtrapolateFunc(x) / XeSparExtrapolateFunc(x) }

parData = read.csv("parData.csv", stringsAsFactors = FALSE)
parData$processes = as.factor(parData$processes)
parData$seqEstRuntime = 0
parData$parEstRuntime = 0
parData[parData$configType == "PureI7",]$seqEstRuntime = predict.lm(I7seqExtrapolate, data.frame(N = parData[parData$configType == "PureI7",]$N))
parData[parData$configType == "PureI7",]$parEstRuntime = predict.lm(I7parExtrapolate, data.frame(N = parData[parData$configType == "PureI7",]$N))
parData[parData$configType == "PureXeS",]$seqEstRuntime = predict.lm(XeSseqExtrapolate, data.frame(N = parData[parData$configType == "PureXeS",]$N))
parData[parData$configType == "PureXeS",]$parEstRuntime = predict.lm(XeSparExtrapolate, data.frame(N = parData[parData$configType == "PureXeS",]$N))
parData$configType = as.factor(parData$configType)
parData$seqSpeedup = parData$seqEstRuntime/parData$wall.clock.time
parData$parSpeedup = parData$parEstRuntime/parData$wall.clock.time

SAVE_LOCATION = "../processedCpuResults/"

defaultTheme = theme(text=element_text(size=14),
                     axis.text=element_text(size=12),
                     legend.text=element_text(size=11),
                     axis.title=element_text(size=14))
machineLabs = labs(linetype = 'Run config', color = 'Processes')


extraPoints = data.frame(N = c(3000, 4000, 6000, 8000, 12000, 16000, 24000))
extraPoints$I7Seqest = I7seqExtrapolateFunc(extraPoints$N)
extraPoints$XeSSeqest = XeSseqExtrapolateFunc(extraPoints$N)
ggplot(seqBenchmarkData) +
    geom_point(aes(N, wall.clock.time)) +
    stat_function(fun = I7seqExtrapolateFunc, aes(color = 'i7 OpenMP')) +
    stat_function(fun = XeSseqExtrapolateFunc, aes(color = 'XeS OpenMP')) +
    geom_point(aes(N, I7Seqest), data = extraPoints, color = 'red') +
    geom_point(aes(N, XeSSeqest), data = extraPoints, color = 'red') +
    labs(x = "N, Particle count",
         y = "Execution time (s)")
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "seqExtrapolate.png"), dpi = 300), NA)

extraPoints = data.frame(N = c(8000, 12000, 16000, 24000))
extraPoints$I7Seqest = I7seqExtrapolateFunc(extraPoints$N)
extraPoints$XeSSeqest = XeSseqExtrapolateFunc(extraPoints$N)
extraPoints$I7Parest = I7parExtrapolateFunc(extraPoints$N)
extraPoints$XeSParest = XeSparExtrapolateFunc(extraPoints$N)
ggplot(parBenchmarkData) +
    geom_point(aes(N, wall.clock.time)) +
    stat_function(fun = I7parExtrapolateFunc, aes(color = 'i7 OpenMP')) +
    stat_function(fun = XeSparExtrapolateFunc, aes(color = 'XeS OpenMP')) +
    geom_point(aes(N, I7Parest), data = extraPoints, color = 'red') +
    geom_point(aes(N, XeSParest), data = extraPoints, color = 'red') +
    labs(x = "N, Particle count",
         y = "Execution time (s)",
         color = 'Implementation')
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "parExtrapolate.png"), dpi = 300), NA)

parBenchmarkData$estSpeedup = 0
parBenchmarkData[parBenchmarkData$configType == 'PureI7',]$estSpeedup = I7seqExtrapolateFunc(parBenchmarkData[parBenchmarkData$configType == 'PureI7',]$N) / parBenchmarkData[parBenchmarkData$configType == 'PureI7',]$wall.clock.time
parBenchmarkData[parBenchmarkData$configType == 'PureXeS',]$estSpeedup = XeSseqExtrapolateFunc(parBenchmarkData[parBenchmarkData$configType == 'PureXeS',]$N) / parBenchmarkData[parBenchmarkData$configType == 'PureXeS',]$wall.clock.time
extraPoints$I7estSpeedup = extraPoints$I7Seqest / extraPoints$I7Parest
extraPoints$XeSestSpeedup = extraPoints$XeSSeqest / extraPoints$XeSParest
ggplot(parBenchmarkData) +
    geom_point(aes(N, estSpeedup)) +
    stat_function(fun = I7speedupExtrapolateFunc, aes(color = 'i7 OpenMP')) +
    stat_function(fun = XeSspeedupExtrapolateFunc, aes(color = 'XeS OpenMP')) +
    geom_point(aes(N, I7estSpeedup), data = extraPoints, color = 'red') +
    geom_point(aes(N, XeSestSpeedup), data = extraPoints, color = 'red') +
    labs(x = "N, Particle count",
         y = "Speedup",
         color = 'Implementation')
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "parSpeedupExtrapolate.png"), dpi = 300), NA)

### Standard runs ###

filter(parData, N < 32000, configType == 'PureXeS') %>%
    ggplot() +
    geom_point(aes(N, wall.clock.time, color=processes)) +
    geom_line(aes(N, wall.clock.time, color=processes)) +
    labs(x = "N, Particle count",
         y = "Execution time (s)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "PureXeS-varyNandProcesses.png"), dpi = 300), NA)

filter(parData, N < 32000, configType == 'PureXeS') %>%
    ggplot() +
    geom_point(aes(N, seqSpeedup, color=processes)) +
    geom_line(aes(N, seqSpeedup, color=processes)) +
    geom_hline(yintercept = 1, linetype='dashed') +
    labs(x = "N, Particle count",
         y = "Speedup (extrapolated sequential)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "PureXeS-seqSpeedup.png"), dpi = 300), NA)

filter(parData, N < 32000, configType == 'PureI7') %>%
    ggplot() +
    geom_point(aes(N, wall.clock.time, color=processes)) +
    geom_line(aes(N, wall.clock.time, color=processes)) +
    labs(x = "N, Particle count",
         y = "Execution time (s)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "PureI7-varyNandProcesses.png"), dpi = 300), NA)

filter(parData, N < 32000, configType == 'PureI7') %>%
    ggplot() +
    geom_point(aes(N, seqSpeedup, color=processes)) +
    geom_line(aes(N, seqSpeedup, color=processes)) +
    geom_hline(yintercept = 1, linetype='dashed') +
    labs(x = "N, Particle count",
         y = "Speedup (extrapolated sequential)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "PureI7-seqSpeedup.png"), dpi = 300), NA)
