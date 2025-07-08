library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(magrittr)

sampleFile = "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_files/Can_amplicons_NRF1KD_QuasR_input.txt"
samples = "amplicon_DE_data"
RegionOfInterest = GRanges("chr7", IRanges(18990973, 18991389))
RegionOfInterest_ext = GRanges("chr7", IRanges(18990931, 18991430))
nrf1_motifs = GRanges(
  rep("chr7", 5),
  IRanges(
    start = c(18991023, 18991044, 18991098, 18991148, 18991224),
    end = c(18991033, 18991054, 18991108, 18991158, 18991234)
  )
)
nrf1_motifs$TF = "NRF1"

CallContextMethylation(
  sampleFile = sampleFile, samples = samples, 
  genome = BSgenome.Mmusculus.UCSC.mm10, RegionOfInterest = RegionOfInterest_ext, 
  coverage = 20, ConvRate.thr = NULL, returnSM = TRUE, clObj = NULL, verbose = FALSE
) -> Methylation

FootprintCharter(
  MethSM = Methylation[[2]], 
  RegionOfInterest = RegionOfInterest, 
  RegionOfInterest_ext = RegionOfInterest_ext, 
  TFBSs = nrf1_motifs, 
  coverage = 30, 
  verbose = FALSE
) -> FC_results

# A
x.axis.breaks = c(start(RegionOfInterest_ext), end(RegionOfInterest_ext))

PlotAvgSMF(
  MethGR = Methylation[[1]],
  RegionOfInterest = RegionOfInterest_ext,
  TFBSs = nrf1_motifs
) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1.25) +
  scale_x_continuous(limits = x.axis.breaks, breaks = x.axis.breaks, labels = format(x.axis.breaks, nsmall=1, big.mark=",")) +
  scale_y_continuous(limits = c(-.1,1), breaks = c(0,1), labels = function(x) format(100*x,digits = 2)) +
  ggtitle(NULL) +
  xlab(as.character(seqnames(RegionOfInterest_ext))) +
  ylab("SMF (1 - meth %)") +
  theme_bw() +
  theme(
    text = element_text(size = 18), axis.title.y = element_text(vjust = -5), axis.title.x = element_text(vjust = 6), legend.position = "none",
    panel.border = element_rect(linewidth = 1, fill = "transparent")
  ) -> pl
ggsave("fig1a.pdf", pl, width = 10, height = 10, units = "cm")

# B
PlotSM(
  MethSM = Methylation[[2]],
  RegionOfInterest = RegionOfInterest_ext, 
  sorting.strategy = "None", 
  SortedReads = NULL
) -> sm.plot
x.axis.breaks = as.integer(c(18990931, 18991430))
sm.plot +
  facet_null() +
  ylab(paste0("2329 molecules")) +
  scale_x_continuous(limits = x.axis.breaks, breaks = x.axis.breaks, labels = format(x.axis.breaks, nsmall=1, big.mark=",")) +
  scale_discrete_manual(aesthetics = "fill", values = c("black", "grey")) + 
  theme_classic() + 
  xlab("chr7") + 
  theme(
    text = element_text(size = 18), axis.title.y = element_text(vjust = 0), axis.title.x = element_text(vjust = 5), 
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_rect(linewidth = 1, fill = "transparent"),
    legend.position = "none"
  ) -> pl
png("fig1b.png", width = 10, height = 10, units = "cm", res = 300)
pl
dev.off()

# C
partitions.order = c(1,9,4,10,13,5,6,2,14,15,12,3,16,11,7,8)
ordered.molecules = lapply(FC_results$partitioned.molecules, function(x){x[rev(partitions.order)]})

PlotSM(
  MethSM = Methylation[[2]], 
  RegionOfInterest = RegionOfInterest_ext, 
  SortedReads = ordered.molecules, 
  sorting.strategy = "custom"
) +
  facet_null() +
  ylab(paste0("895 molecules")) +
  scale_x_continuous(limits = x.axis.breaks, breaks = x.axis.breaks, labels = format(x.axis.breaks, nsmall=1, big.mark=",")) +
  scale_discrete_manual(aesthetics = "fill", values = c("black", "grey")) + 
  theme_classic() + 
  xlab("chr7") + 
  theme(
    text = element_text(size = 18), axis.title.y = element_text(vjust = 0), axis.title.x = element_text(vjust = 5), 
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_rect(linewidth = 1, fill = "transparent"),
    legend.position = "none"
  ) -> pl
png("fig1c.png", width = 10, height = 10, units = "cm", res = 300)
pl
dev.off()

# D
PlotFootprints(
  MethSM = Methylation[[2]], 
  partitioned.molecules = FC_results$partitioned.molecules, 
  footprints.df = FC_results$footprints.df, 
  TFBSs = nrf1_motifs
) -> footprint.plot

footprint.plot.2 = footprint.plot
footprint.plot.2$layers[[1]] = footprint.plot$layers[[5]]
footprint.plot.2$layers[[2]] = footprint.plot$layers[[4]]
footprint.plot.2$layers[[3]] = footprint.plot$layers[[2]]
footprint.plot.2$layers[[4]] = footprint.plot$layers[[3]]
footprint.plot.2$layers[[6]] = NULL
footprint.plot.2$layers[[5]] = NULL
footprint.plot.2$layers[[1]]$aes_params$alpha = 1
footprint.plot.2$data %<>%
  mutate(partition.nr = factor(partition.nr, levels = c(1,9,4,10,13,5,6,2,14,15,12,3,16,11,7,8))) %>%
  filter(partition.nr %in% c(1,6,16))
footprint.plot.2 +
  geom_point(size = .75) +
  geom_line(linewidth = .75) +
  xlab(as.character(seqnames(RegionOfInterest))) +
  ylab("SMF (1 - meth %)") +
  facet_wrap(~partition.nr, ncol = 1) +
  ggtitle(NULL, subtitle = NULL) +
  scale_x_continuous(limits = x.axis.breaks, breaks = x.axis.breaks, labels = format(x.axis.breaks, nsmall=1, big.mark=",")) +
  scale_y_continuous(limits = c(-.5,1), breaks = c(0,1), labels = function(x) format(100*x,digits = 2)) +
  theme_bw() +
  theme(text = element_text(size = 18), axis.title.x = element_text(vjust = 6), 
        panel.border = element_rect(linewidth = 1, fill = "transparent"), 
        strip.text = element_blank(), strip.background = element_blank()) -> pl
ggsave("fig1d.pdf", pl, width = 10, height = 10, units = "cm")

# E
FC_results$footprints.df %>%
  mutate(biological.state = case_when(
    biological.state == "noise" ~ "accessible",
    biological.state == "unrecognized" ~ "nucleosome",
    !biological.state %in% c("noise", "unrecognized") ~ biological.state
  )) %>%
  Plot_FootprintCharter_SM(
    footprints.df = ., 
    RegionOfInterest = RegionOfInterest_ext, 
    partitions.order = partitions.order
  ) +
  facet_null() +
  ylab(paste0("895 molecules")) +
  scale_x_continuous(limits = x.axis.breaks, breaks = x.axis.breaks, labels = format(x.axis.breaks, nsmall=1, big.mark=",")) +
  theme_classic() + 
  xlab("chr7") + 
  theme(
    text = element_text(size = 18), axis.title.y = element_text(vjust = 0), axis.title.x = element_text(vjust = 5), 
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_rect(linewidth = 1, fill = "transparent"),
    legend.position = "none"
  ) -> pl
png("fig1e.png", width = 10, height = 10, units = "cm", res = 300)
pl
dev.off()
