library(ggplot2)
library(dplyr)

setwd("E:/Brownie/openmp/report/rawResults/")
seqData = read.csv("seqData.csv")
parData = read.csv("parData.csv")

DEFAULT_N = 1000
DEFAULT_L = 20000
DEFAULT_r = 1
DEFAULT_STEPS = 1000

SAVE_LOCATION = "../processedResults/"

defaultTheme = theme(text=element_text(size=14),
                     axis.text=element_text(size=12),
                     legend.text=element_text(size=12),
                     axis.title=element_text(size=14))
machineLabs = labs(color = 'Machine')

varyNData = with(seqData, seqData[which(L == DEFAULT_L & r == DEFAULT_r & steps == DEFAULT_STEPS),])
ggplot(varyNData) +
    geom_point(aes(N, wall.clock.time, color=machine)) + 
    geom_line(aes(N, wall.clock.time, color=machine)) +
    labs(title = "Plot of wall clock time (s) against particle count",
         x  = "N, Particle count",
         y = "Wall clock time (s)") +
    defaultTheme + machineLabs
ggsave(paste0(SAVE_LOCATION, "seq-varyN.png"), dpi = 300)
write.csv(varyNData, paste0(SAVE_LOCATION, "seq-varyN.csv"), row.names = FALSE)
rm(varyNData)

varyLData = with(seqData, seqData[which(N == DEFAULT_N & r == DEFAULT_r & steps == DEFAULT_STEPS),])
ggplot(varyLData) +
    geom_point(aes(L, wall.clock.time, color=machine)) +
    geom_line(aes(L, wall.clock.time, color=machine)) +
    labs(title = "Plot of wall clock time (s) against box length (units)",
         x  = "L, Box length (units)",
         y = "Wall clock time (s)") +
    defaultTheme + machineLabs
ggsave(paste0(SAVE_LOCATION, "seq-varyL.png"), dpi = 300)
write.csv(varyLData, paste0(SAVE_LOCATION, "seq-varyL.csv"), row.names = FALSE)
rm(varyLData)

varyRData = with(seqData, seqData[which(N == DEFAULT_N & L == DEFAULT_L & steps == DEFAULT_STEPS),])
ggplot(varyRData) +
    geom_point(aes(r, wall.clock.time, color=machine)) +
    geom_line(aes(r, wall.clock.time, color=machine)) +
    labs(title = "Plot of wall clock time (s) against particle radius length (units)",
         x  = "r, Particle radius length (units)",
         y = "Wall clock time (s)") +
    defaultTheme + machineLabs
ggsave(paste0(SAVE_LOCATION, "seq-varyR.png"), dpi = 300)
write.csv(varyRData, paste0(SAVE_LOCATION, "seq-varyR.csv"), row.names = FALSE)
rm(varyRData)

varyStepsData = with(seqData, seqData[which(N == DEFAULT_N & L == DEFAULT_L & r == DEFAULT_r),])
ggplot(varyStepsData) +
    geom_point(aes(steps, wall.clock.time, color=machine)) +
    geom_line(aes(steps, wall.clock.time, color=machine)) +
    labs(title = "Plot of wall clock time (s) against number of steps",
         x  = "Steps",
         y = "Wall clock time (s)") +
    defaultTheme + machineLabs
ggsave(paste0(SAVE_LOCATION, "seq-varySteps.png"), dpi = 300)
write.csv(varyStepsData, paste0(SAVE_LOCATION, "seq-varySteps.csv"), row.names = FALSE)
rm(varyStepsData)

for (n in sort(unique(parData$N))) {
    varyThreadsData = with(parData, parData[which(n == N & L == DEFAULT_L & r == DEFAULT_r & steps == DEFAULT_STEPS),])
    print(ggplot(varyThreadsData) +
              geom_point(aes(threads, wall.clock.time, color=machine)) + 
              geom_line(aes(threads, wall.clock.time, color=machine)) +
              labs(title = paste("Plot of wall clock time (s) against number of threads for N =", n),
                   x  = "Thread count",
                   y = "Wall clock time (s)",
                   legend = "Machine",
                   caption = "Axis for wall clock time is in log10 scale") +
              scale_y_log10() +
              defaultTheme +
              machineLabs)
    ggsave(paste0(SAVE_LOCATION, "par-", n, "N-varyThreads.png"), dpi = 300)
    
    print(ggplot(varyThreadsData) +
              geom_point(aes(threads, context.switches, color=machine)) + 
              geom_line(aes(threads, context.switches, color=machine)) +
              labs(title = paste("Plot of context switches against number of threads for N =", n),
                   x  = "Thread count",
                   y = "Context switches",
                   legend = "Machine") +
              defaultTheme +
              machineLabs)
    ggsave(paste0(SAVE_LOCATION, "par-", n, "N-contextSwitches.png"), dpi = 300)
    
    seqPart = with(seqData, seqData[which(n == N & L == DEFAULT_L & r == DEFAULT_r & steps == DEFAULT_STEPS),])
    colnames(seqPart)[colnames(seqPart)=="wall.clock.time"] = "seqTime"
    varyThreadsData = left_join(varyThreadsData, select(seqPart, machine, seqTime), by = 'machine')
    varyThreadsData$speedup = varyThreadsData$seqTime / varyThreadsData$wall.clock.time
    
    print(ggplot(varyThreadsData) +
              geom_point(aes(threads, speedup, color=machine)) + 
              geom_line(aes(threads, speedup, color=machine)) +
              labs(title = paste("Plot of speedup against number of threads for N =", n),
                   x  = "Thread count",
                   y = "Speedup") +
              defaultTheme +
              machineLabs)
    ggsave(paste0(SAVE_LOCATION, "par-", n, "N-speedup.png"), dpi = 300)
    write.csv(varyThreadsData, paste0(SAVE_LOCATION, "par-", n, "N-varyThreads.csv"), row.names = FALSE)
}

rm(varyThreadsData, seqPart)
rm(list=ls())
