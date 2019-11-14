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

### Standard runs ###

ggplot(gpuData) +
    geom_point(aes(N, time, color=machine)) +
    geom_line(aes(N, time, color=machine, linetype=config)) +
    labs(x = "N, Particle count",
         y = "Execution time (s)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "gpu-varyN.png"), dpi = 300), NA)

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
