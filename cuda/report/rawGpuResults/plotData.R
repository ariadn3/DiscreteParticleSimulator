library(ggplot2)
library(dplyr)

setwd("E:/Brownie/cuda/report/rawResults/")
cpuData = read.csv("cpuData.csv")
gpuData = read.csv("gpuData.csv")

SAVE_LOCATION = "../gpuProcessedResults/"

defaultTheme = theme(text=element_text(size=14),
                     axis.text=element_text(size=12),
                     legend.text=element_text(size=12),
                     axis.title=element_text(size=14))
machineLabs = labs(linetype = 'Run config', color = 'Machine')

ggplot(gpuData) +
    geom_point(aes(N, time, color=machine)) +
    geom_line(aes(N, time, color=machine, linetype=config)) +
    labs(title = "Plot of time (s) against particle count",
         x = "N, Particle count",
         y = "Execution time (s)") +
    defaultTheme + machineLabs
ggsave(paste0(SAVE_LOCATION, "gpu-varyN.png"), dpi = 300)

ggplot(cpuData) + 
    geom_point(aes(N, time)) +
    geom_line(aes(N, time)) +
    labs(title = "Plot of time (s) against particle count",
         x = "N, Particle count",
         y = "Execution time (s)") +
    defaultTheme + machineLabs
ggsave(paste0(SAVE_LOCATION, "cpu-varyN.png"), dpi = 300)

colnames(cpuData)[which(colnames(cpuData) == "time")] = "cpuTime"
colnames(gpuData)[which(colnames(gpuData) == "time")] = "gpuTime"
varyNData = left_join(gpuData,select(cpuData, N, L, r, steps, cpuTime), by = c('N', 'L', 'r', 'steps'))
varyNData$speedup = varyNData$cpuTime / varyNData$gpuTime

ggplot(varyNData) +
    geom_point(aes(N, speedup, color=machine)) +
    geom_line(aes(N, speedup, color=machine, linetype=config)) +
    labs(title = "Plot of speedup against particle count",
         x = "N, Particle count",
         y = "Speedup") +
    defaultTheme + machineLabs
ggsave(paste0(SAVE_LOCATION, "gpu-speedup.png"), dpi=300)