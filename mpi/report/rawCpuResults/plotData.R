library(ggplot2)
library(dplyr)
library(tidyr)

SAVE = FALSE

seqBenchmarkData = read.csv("seqBenchmarkData.csv", stringsAsFactors = FALSE)
parBenchmarkData = read.csv("parBenchmarkData.csv", stringsAsFactors = FALSE)
I7seqExtrapolate = lm(wall.clock.time ~ poly(N, 2), data = seqBenchmarkData[seqBenchmarkData$configType=="PureI7",])
summary(I7seqExtrapolate)$r.squared
XeSseqExtrapolate = lm(wall.clock.time ~ poly(N, 2), data = seqBenchmarkData[seqBenchmarkData$configType=="PureXeS",])
summary(XeSseqExtrapolate)$r.squared
I7parExtrapolate = lm(wall.clock.time ~ poly(N, 2), data = parBenchmarkData[parBenchmarkData$configType=="PureI7",])
summary(I7parExtrapolate)$r.squared
XeSparExtrapolate = lm(wall.clock.time ~ poly(N, 2), data = parBenchmarkData[parBenchmarkData$configType=="PureXeS",])
summary(XeSparExtrapolate)$r.squared

predict.lm(I7seqExtrapolate, data.frame(N = 2000))

# setwd("E:/Brownie/cuda/report/rawResults/")
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

filter(parData, N < 32000, configType == 'PureXeS') %>%
    ggplot() +
    geom_point(aes(N, parSpeedup, color=processes)) +
    geom_line(aes(N, parSpeedup, color=processes)) +
    geom_hline(yintercept = 1, linetype='dashed') +
    labs(x = "N, Particle count",
         y = "Speedup (extrapolated OpenMP)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "PureXeS-parSpeedup.png"), dpi = 300), NA)

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

filter(parData, N < 32000, configType == 'PureI7') %>%
    ggplot() +
    geom_point(aes(N, parSpeedup, color=processes)) +
    geom_line(aes(N, parSpeedup, color=processes)) +
    geom_hline(yintercept = 1, linetype='dashed') +
    labs(x = "N, Particle count",
         y = "Speedup (extrapolated OpenMP)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "PureI7-parSpeedup.png"), dpi = 300), NA)
