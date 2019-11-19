library(ggplot2)
library(dplyr)
library(tidyr)

SAVE = TRUE

# setwd("E:/Brownie/cuda/report/rawResults/")
gpuData = read.csv("gpuData.csv")

SAVE_LOCATION = "../processedGpuResults/"

defaultTheme = theme(text=element_text(size=14),
                     axis.text=element_text(size=12),
                     legend.text=element_text(size=11),
                     axis.title=element_text(size=14))
machineLabs = labs(linetype = 'Run config', color = 'Machine')

# Extrapolation
seqBenchmarkData = read.csv("../rawCpuResults/seqBenchmarkData.csv", stringsAsFactors = FALSE)
parBenchmarkData = read.csv("../rawCpuResults/parBenchmarkData.csv", stringsAsFactors = FALSE)
I7seqExtrapolate = lm(wall.clock.time ~ poly(N, 2), data = seqBenchmarkData[seqBenchmarkData$configType=="PureI7",])
# summary(I7seqExtrapolate)$r.squared
XeSseqExtrapolate = lm(wall.clock.time ~ poly(N, 2), data = seqBenchmarkData[seqBenchmarkData$configType=="PureXeS",])
# summary(XeSseqExtrapolate)$r.squared

I7seqExtrapolateFunc = function(x) { predict.lm(I7seqExtrapolate, data.frame(N=c(x))) }
XeSseqExtrapolateFunc = function(x) { predict.lm(XeSseqExtrapolate, data.frame(N=c(x))) }

### Standard runs ###
tempData = gpuData
tempData$i7SeqEst = I7seqExtrapolateFunc(tempData$N)
tempData$XeSSeqEst = XeSseqExtrapolateFunc(tempData$N)
tempData$i7Speedup = tempData$i7SeqEst / tempData$time
tempData$XeSSpeedup = tempData$XeSSeqEst / tempData$time

ggplot(tempData) +
    geom_point(aes(N, time, color=machine)) +
    geom_line(aes(N, time, color=machine, linetype=config)) +
    labs(x = "N, Particle count",
         y = "Execution time (s)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "gpu-varyN.png"), dpi = 300), NA)

tempData %>%
    filter(machine != 'TitanV') %>%
    ggplot() +
    geom_point(aes(N, i7Speedup, color=machine)) +
    geom_line(aes(N, i7Speedup, color=machine, linetype=config)) +
    labs(x = "N, Particle count",
         y = "Speedup (extrapolated sequential)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "gpu-speedupI7.png"), dpi = 300), NA)

tempData %>%
    filter(machine != 'TitanV') %>%
    ggplot() +
    geom_point(aes(N, XeSSpeedup, color=machine)) +
    geom_line(aes(N, XeSSpeedup, color=machine, linetype=config)) +
    labs(x = "N, Particle count",
         y = "Speedup (extrapolated sequential)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "gpu-speedupXeS.png"), dpi = 300), NA)

### Float runs ###

floatData = read.csv("floatData.csv")
tempFloatData = floatData
tempFloatData$type = as.factor("float")
tempGpuData = gpuData
tempGpuData = tempGpuData[,-c(8:11)]
colnames(tempGpuData)[7] = "time"
tempGpuData$type = as.factor("double")
varyTypeData = rbind(tempGpuData, tempFloatData)
rm(tempGpuData, tempFloatData)
colnames
# varyTypeData = left_join(gpuData, select(floatData, N, L, r, steps, config, machine), by = c('N', 'L', 'r', 'steps', 'config', 'machine'))

ggplot(varyTypeData[varyTypeData$config == "fmad",]) +
    geom_point(aes(N, time, colour=machine)) +
    geom_line(aes(N, time, color=machine, linetype=type)) +
    labs(x = "N, Particle count",
         y = "Execution time (s)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "double-float-comparison-fmad.png"), dpi = 300), NA)

ggplot(varyTypeData[varyTypeData$config == "nofmad",]) +
    geom_point(aes(N, time, colour=machine)) +
    geom_line(aes(N, time, color=machine, linetype=type)) +
    labs(x = "N, Particle count",
         y = "Execution time (s)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "double-float-comparison-nofmad.png"), dpi = 300), NA)

rm(varyTypeData)
