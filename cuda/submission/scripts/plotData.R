library(ggplot2)
library(dplyr)
library(tidyr)

SAVE = FALSE

setwd("E:/Brownie/cuda/report/rawGpuResults/")
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
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "gpu-varyN-withJetson.png"), dpi = 300), NA)

gpuData = gpuData[!(gpuData$machine == "JetsonTX2"),]

ggplot(gpuData) +
    geom_point(aes(N, time, color=machine)) +
    geom_line(aes(N, time, color=machine, linetype=config)) +
    labs(x = "N, Particle count",
         y = "Execution time (s)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "gpu-varyN.png"), dpi = 300), NA)

ggplot(cpuData) + 
    geom_point(aes(N, time)) +
    geom_line(aes(N, time)) +
    labs(x = "N, Particle count",
         y = "Execution time (s)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "cpu-varyN.png"), dpi = 300), NA)

colnames(cpuData)[which(colnames(cpuData) == "time")] = "cpuTime"
colnames(gpuData)[which(colnames(gpuData) == "time")] = "gpuTime"
varyNData = left_join(gpuData,select(cpuData, N, L, r, steps, cpuTime), by = c('N', 'L', 'r', 'steps'))
varyNData$speedup = varyNData$cpuTime / varyNData$gpuTime

ggplot(varyNData) +
    geom_point(aes(N, speedup, color=machine)) +
    geom_line(aes(N, speedup, color=machine, linetype=config)) +
    labs(x = "N, Particle count",
         y = "Speedup") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "gpu-speedup.png"), dpi=300))

rm(varyNData)

### Analysis of GPU Kernel runtimes ###

propData = gather(gpuData, "functType", "funcTime", checkWallCollision:updateParticles)
propData = propData[propData$machine == "TitanRTX" & propData$config == "fmad",]
propData = propData %>% 
    group_by(N, L, r, steps, config, machine) %>%
    mutate(gpuTotalTime = sum(funcTime)) %>%
    mutate(gpuProp = funcTime/gpuTotalTime)

ggplot(propData) +
    geom_point(aes(N, gpuProp, color=functType)) +
    geom_line(aes(N, gpuProp, color=functType)) +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    theme(legend.position="bottom") + 
    labs(y = "Proportion of GPU kernel time",
         color = "Kernel") +
    defaultTheme
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "titan-fmad-proportion-GPU.png"), dpi = 300))

propData = gather(gpuData, "functType", "funcTime", checkWallCollision:updateParticles)
propData = propData[propData$machine == "TitanRTX",]

ggplot(propData) +
    geom_point(aes(N, funcTime, color=functType)) +
    geom_line(aes(N, funcTime, color=functType, linetype=config)) +
    labs(y = "GPU kernel time (s)",
         linetype = "Run config",
         color = "Kernel") + 
    defaultTheme 
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "titan-fmad-kernelTime-all.png"), dpi = 300))
ggplot(filter(propData, functType == 'checkCollision')) +
    geom_point(aes(N, funcTime, color=functType)) +
    geom_line(aes(N, funcTime, color=functType, linetype=config)) +
    labs(y = "GPU kernel time (s)",
         linetype = "Run config",
         color = "Kernel") +
    defaultTheme 
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "titan-fmad-kernelTime-checkCollision.png"), dpi = 300))
ggplot(filter(propData, functType != 'checkCollision')) +
    geom_point(aes(N, funcTime, color=functType)) +
    geom_line(aes(N, funcTime, color=functType, linetype=config)) +
    labs(y = "GPU kernel time (s)",
         linetype = "Run config",
         color = "Kernel") + 
    defaultTheme
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "titan-fmad-kernelTime-noCheckCollision.png"), dpi = 300))

rm(propData)

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
    labs(y = "Execution time (s)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "double-float-comparison-fmad.png"), dpi = 300), NA)

ggplot(varyTypeData[varyTypeData$config == "nofmad",]) +
    geom_point(aes(N, time, colour=machine)) +
    geom_line(aes(N, time, color=machine, linetype=type)) +
    labs(y = "Execution time (s)") +
    defaultTheme + machineLabs
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "double-float-comparison-nofmad.png"), dpi = 300), NA)


rm(varyTypeData)

### ChunkSize runs ###

chunkSizeData = read.csv("chunkSizeData.csv")
chunkSizeData = chunkSizeData[chunkSizeData$config == "fmad",]

n32000Data = chunkSizeData[chunkSizeData$N == 32000,]

n32000Data %>%
    filter(X1DchunkSize == 64) %>%
    ggplot() +
    geom_point(aes(X2DchunkSize, time)) +
    geom_line(aes(X2DchunkSize, time)) +
    labs(x = 'Chunk size of 2D kernels',
         y = 'Execution time (s)') +
    defaultTheme
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "chunkSize-2D-32000-fmad.png"), dpi = 300), NA)
n32000Data %>%
    filter(X2DchunkSize == 64) %>%
    ggplot() +
    geom_point(aes(X1DchunkSize, time)) +
    geom_line(aes(X1DchunkSize, time)) +
    labs(x = 'Chunk size of 1D kernels',
         y = 'Execution time (s)') +
    defaultTheme
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "chunkSize-1D-32000-fmad.png"), dpi = 300), NA)

rm(chunkSizeData)
rm(n32000Data)

# Filter times

filterData = read.csv("filterTimes.csv")
ggplot(filterData) +
    geom_point(aes(N, avgTime)) +
    geom_line(aes(N, avgTime)) +
    labs(x = "N", y = "Average time (per step) to qsort (s)") +
    defaultTheme
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "avgSortTime.png"), dpi = 300), NA)

rm(filterData)

# Divergence analysis

diverged = read.csv("divergedParticles.csv")
diverged = gather(diverged, "type", "particles", fmad:nofmad)

ggplot(diverged) +
    geom_hline(yintercept = 1000, linetype = "dashed") +
    geom_line(aes(step, particles, color = type)) +
    labs(x = "Step number",
         y = "Particles diverged",
         color = "Run config") +
    xlim(0, 175) + 
    defaultTheme
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "divergedParticles.png"), dpi = 300), NA)

rm(diverged)

# Collision type analysis

collisions = read.csv("collisionType.csv")
ggplot(collisions) +
    geom_point(aes(N, ppCollisions/pwCollisions, color=type), alpha = 0.5) +
    geom_line(aes(N, ppCollisions/pwCollisions, color=type), alpha = 0.5) +
    labs(x = "N", y = "P-P/P-W collision ratio") + 
    defaultTheme
ifelse(SAVE, ggsave(paste0(SAVE_LOCATION, "collisionType.png"), dpi = 300), NA)
