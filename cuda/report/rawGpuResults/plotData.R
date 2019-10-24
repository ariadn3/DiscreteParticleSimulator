library(ggplot2)
library(dplyr)
library(tidyr)

# setwd("E:/Brownie/cuda/report/rawResults/")
cpuData = read.csv("cpuData.csv")
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
ggsave(paste0(SAVE_LOCATION, "gpu-varyN-withJetson.png"), dpi = 300)

gpuData = gpuData[!(gpuData$machine == "JetsonTX2"),]

ggplot(gpuData) +
    geom_point(aes(N, time, color=machine)) +
    geom_line(aes(N, time, color=machine, linetype=config)) +
    labs(x = "N, Particle count",
         y = "Execution time (s)") +
    defaultTheme + machineLabs
ggsave(paste0(SAVE_LOCATION, "gpu-varyN.png"), dpi = 300)

ggplot(cpuData) + 
    geom_point(aes(N, time)) +
    geom_line(aes(N, time)) +
    labs(x = "N, Particle count",
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

### Analysis of GPU Kernel runtimes ###

propData = gather(gpuData, "functType", "funcTime", checkWallCollision:updateParticles)
propData = propData[propData$machine == "TitanRTX" & propData$config == "fmad",]
propData = propData %>% 
    group_by(N, L, r, steps, config, machine) %>%
    mutate(gpuTotalTime = sum(funcTime)) %>%
    mutate(gpuProp = funcTime/gpuTotalTime, totalProp = funcTime/gpuTime)

ggplot(propData) +
    geom_point(aes(N, gpuProp, color=functType)) +
    geom_line(aes(N, gpuProp, color=functType)) +
    geom_hline(yintercept = 1) +
    theme(legend.position="bottom") + 
    defaultTheme
ggsave(paste0(SAVE_LOCATION, "titan-fmad-proportion-GPU.png"), dpi = 300)

ggplot(propData) +
    geom_point(aes(N, totalProp, color=functType)) +
    geom_line(aes(N, totalProp, color=functType)) +
    geom_hline(yintercept = 1) +
    theme(legend.position="bottom") + 
    defaultTheme
ggsave(paste0(SAVE_LOCATION, "titan-fmad-proportion-Total.png"), dpi = 300)

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
    defaultTheme + machineLabs
ggsave(paste0(SAVE_LOCATION, "double-float-comparison-fmad.png"), dpi = 300)

ggplot(varyTypeData[varyTypeData$config == "nofmad",]) +
    geom_point(aes(N, time, colour=machine)) +
    geom_line(aes(N, time, color=machine, linetype=type)) +
    defaultTheme + machineLabs
ggsave(paste0(SAVE_LOCATION, "double-float-comparison-nofmad.png"), dpi = 300)

### ChunkSize runs ###

chunkSizeData = read.csv("chunkSizeData.csv")
chunkSizeData = chunkSizeData[chunkSizeData$config == "fmad",]

n2000Data = chunkSizeData[chunkSizeData$N == 2000,]
n32000Data = chunkSizeData[chunkSizeData$N == 32000,]

# Filter times

filterData = read.csv("filterTimes.csv")
ggplot(filterData, ) +
    geom_point(aes(N, avgTime)) +
    geom_line(aes(N, avgTime)) +
    labs(x = "N", y = "Average time (per step) to qsort") +
    defaultTheme
ggsave(paste0(SAVE_LOCATION, "avgSortTime.png"), dpi = 300)
    
# Divergence analysis

diverged = read.csv("divergedParticles.csv")
diverged = gather(diverged, "type", "particles", fmad:nofmad)

ggplot(diverged) +
    geom_hline(yintercept=1000) +
    geom_line(aes(step, particles, color=type)) +
    labs(x = "Step number", y = "Particles diverged") +
    xlim(0, 175) + 
    defaultTheme
ggsave(paste0(SAVE_LOCATION, "divergedParticles.png"), dpi = 300)
