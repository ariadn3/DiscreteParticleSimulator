library(ggplot2)
library(dplyr)

setwd("E:/Brownie/openmp/report/rawResults/")
seqData = read.csv("seqData.csv")
parData = read.csv("parData.csv")
par2Data = read.csv("parData2.csv")

DEFAULT_N = 2000
DEFAULT_L = 20000
DEFAULT_r = 1
DEFAULT_STEPS = 1000

SAVE_LOCATION = "../extraProcessedResults/"

defaultTheme = theme(text=element_text(size=12),
                     axis.text=element_text(size=12),
                     legend.text=element_text(size=12),
                     axis.title=element_text(size=14))
machineLabs = labs(color = 'Machine')

for (g in sort(unique(par2Data$gridSize))) {
    varyThreadsData = with(par2Data, par2Data[which(N == DEFAULT_N & L == DEFAULT_L & r == DEFAULT_r & steps == DEFAULT_STEPS & g == gridSize),])
    
    parPart = with(parData, parData[which(N == DEFAULT_N & L == DEFAULT_L & r == DEFAULT_r & steps == DEFAULT_STEPS),])
    colnames(parPart)[colnames(parPart)=="wall.clock.time"] = "parTime"
    parPart = select(parPart, machine, threads, parTime)
    varyThreadsData = left_join(varyThreadsData, parPart, by = c('machine', 'threads'))
    
    seqPart = with(seqData, seqData[which(N == DEFAULT_N & L == DEFAULT_L & r == DEFAULT_r & steps == DEFAULT_STEPS),])
    colnames(seqPart)[colnames(seqPart)=="wall.clock.time"] = "seqTime"
    seqPart = select(seqPart, machine, seqTime)
    varyThreadsData = left_join(varyThreadsData, seqPart, by = 'machine')
    varyThreadsData$speedup = varyThreadsData$seqTime / varyThreadsData$wall.clock.time
    
    print(ggplot(varyThreadsData) +
              geom_point(aes(threads, wall.clock.time, color=machine)) + 
              geom_line(aes(threads, wall.clock.time, color=machine)) +
              geom_point(aes(threads, parTime, color=machine)) +
              geom_line(aes(threads, parTime, color=machine), linetype='dashed') +
              labs(title = paste("Plot of wall clock time (s) against number of threads for grid size = ", g),
                   x  = "Thread count",
                   y = "Wall clock time (s)",
                   legend = "Machine",
                   caption = "Axis for wall clock time is in log10 scale") +
              scale_y_log10() +
              defaultTheme +
              machineLabs)
    ggsave(paste0(SAVE_LOCATION, paste0("optPar-gridSize", g, "-varyThreads.png")), dpi = 300)
    
    print(ggplot(varyThreadsData) +
              geom_point(aes(threads, speedup, color=machine)) + 
              geom_line(aes(threads, speedup, color=machine)) +
              geom_point(aes(threads, seqTime/parTime, color=machine)) +
              geom_line(aes(threads, seqTime/parTime, color=machine), linetype='dashed') +
              labs(title = paste("Plot of speedup against number of threads for grid size =", g),
                   x  = "Thread count",
                   y = "Speedup") +
              defaultTheme +
              machineLabs)
    ggsave(paste0(SAVE_LOCATION, paste0("optPar-gridSize", g, "-speedup.png"))), dpi = 300)
}

rm(list=ls())
