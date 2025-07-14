library(SingleMoleculeFootprinting)
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(magrittr)
library(parallel)
library(ComplexHeatmap)
library(peakRAM)

# A
Methylation = qs::qread(system.file("extdata", "Methylation_4.qs", package="SingleMoleculeFootprinting"))
MethSM = Methylation[[2]]
RegionOfInterest = GRanges("chr6", IRanges(88106000, 88106500))
RegionOfInterest = IRanges::resize(RegionOfInterest, 80, "center")

iteration_df = expand_grid(
  k = c(2,4,8,16,32, 64, 128), 
  n = c(1,3,5,10,20,50,100)
  )

Reduce(rbind, mclapply(seq(nrow(iteration_df)), function(i){
  
  k_value = iteration_df$k[i]
  n_value = iteration_df$n[i]
  
  FootprintCharter(
    MethSM = MethSM,
    RegionOfInterest = RegionOfInterest,
    k = k_value,
    n = n_value, 
    cytosine.coverage.thr = min(5, n_value)
  ) -> FC_results
  
  data.frame(
    k = k_value,
    n = n_value,
    computed_clusters = max(as.integer(FC_results$footprints.df$partition.nr))
  )
  
}, mc.preschedule = TRUE, mc.cores = 7)) -> results

results %>%
  spread(k, computed_clusters) %>%
  mutate(n = factor(n, levels = rev(n))) %>%
  arrange(n) %>%
  column_to_rownames("n") %>%
  as.matrix() %>%
  Heatmap(
    cluster_columns = FALSE, cluster_rows = FALSE, 
    row_title = "n", column_title = "k", 
    row_title_rot = 0, row_names_side = "left", row_names_centered = TRUE,
    column_title_side = "bottom", column_names_rot = 0, column_names_centered = TRUE,
    col = circlize::colorRamp2(colors = jcolors::jcolors_contin("pal12")(max(.)), breaks = seq(max(.))), heatmap_legend_param = list(at = c(2,4,8,16,32,64,128)), name = "output clusters"
    ) -> pl
pdf("/g/krebs/barzaghi/analyses/15.02.24_bioinformatics_application_note_figures/supplementary_fig1a.pdf", width = 5, height = 4)
pl
dev.off()

# B
Reduce(rbind, mclapply(seq(nrow(iteration_df)), function(i){
  
  k_value = iteration_df$k[i]
  n_value = iteration_df$n[i]
  
  FootprintCharter(
    MethSM = MethSM,
    RegionOfInterest = RegionOfInterest,
    k = k_value,
    n = n_value,
    cytosine.coverage.thr = min(5, n_value)
  ) -> FC_results
  
  FC_results$footprints.df %>%
    filter(biological.state == "noise") %>%
    dplyr::select(footprint.idx) %>%
    mutate(
      k = k_value,
      n = n_value,
    ) -> return.df
  
  if(nrow(return.df) == 0){
    return.df = data.frame(
      footprint.idx = NA, k = k_value,n = n_value
    )
  }
  
  return(return.df)
  
}, mc.preschedule = TRUE, mc.cores = 7)) -> results

results %>%
  distinct(footprint.idx, k, n) %>%
  group_by(k,n) %>%
  summarise(nr_noise_footprints = n(), .groups = "drop") %>%
  spread(k, nr_noise_footprints, fill = 0) %>%
  mutate(n = factor(n, levels = rev(n))) %>%
  arrange(n) %>%
  column_to_rownames("n") %>%
  as.matrix() %>%
  Heatmap(
    cluster_columns = FALSE, cluster_rows = FALSE, 
    row_title = "n", column_title = "k", 
    row_title_rot = 0, row_names_side = "left", row_names_centered = TRUE,
    column_title_side = "bottom", column_names_rot = 0, column_names_centered = TRUE,
    col = circlize::colorRamp2(colors = jcolors::jcolors_contin("pal12")(max(.)), breaks = seq(max(.))), heatmap_legend_param = list(at = c(1,3,10,15,20)), name = "nr noise footprints"
  ) -> pl
pdf("/g/krebs/barzaghi/analyses/15.02.24_bioinformatics_application_note_figures/supplementary_fig1b.pdf", width = 5, height = 4)
pl
dev.off()

# C
sampleFile = "/g/krebs/barzaghi/HTS/SMF/MM/QuasR_input_files/Can_amplicons_NRF1KD_QuasR_input.txt"
samples = "amplicon_DE_data"
RegionOfInterest = GRanges("chr7", IRanges(18990973, 18991389))
RegionOfInterest_ext = GRanges("chr7", IRanges(18990931, 18991430))

CallContextMethylation(
  sampleFile = sampleFile, samples = samples, 
  genome = BSgenome.Mmusculus.UCSC.mm10, RegionOfInterest = RegionOfInterest_ext, 
  coverage = 20, ConvRate.thr = NULL, returnSM = TRUE, clObj = NULL, verbose = FALSE
) -> Methylation
MethSM = Methylation[[2]]

iteration_df = expand_grid(
  nr_molecules = c(100, 250, 500, 1000, 1500, 2000, 5000, 10000),
  nr_cytosines = c(50, 100, 300, 500, 1000, 2000, 3000, 10000)
)

Reduce(rbind, mclapply(seq(nrow(iteration_df)), function(i){
  
  print(i)
  
  nr_molecules = iteration_df$nr_molecules[i]
  nr_cytosines = iteration_df$nr_cytosines[i]
  
  MethSM_subset = list()
  MethSM_subset$amplicon_DE_data = MethSM$amplicon_DE_data[
    sample(rownames(MethSM$amplicon_DE_data), size = nr_molecules, replace = TRUE),
    sort(sample(colnames(MethSM$amplicon_DE_data), size = nr_cytosines, replace = TRUE))
    ]
  
  start = Sys.time()
  FootprintCharter(
    MethSM = MethSM_subset,
    RegionOfInterest = RegionOfInterest,
  ) -> FC_results
  end = Sys.time()
  
  data.frame(
    nr_molecules = nr_molecules,
    nr_cytosines = nr_cytosines,
    running_time = end - start
  )
  
}, mc.preschedule = TRUE, mc.cores = 1)) -> results

results %>%
  mutate(running_time = ceiling(as.numeric(running_time))) %>%
  spread(nr_cytosines, running_time) %>%
  mutate(nr_molecules = factor(nr_molecules, levels = rev(nr_molecules))) %>%
  arrange(nr_molecules) %>%
  column_to_rownames("nr_molecules") %>%
  as.matrix() %>%
  Heatmap(
    cluster_columns = FALSE, cluster_rows = FALSE, 
    row_title = "Nr molecules", column_title = "Nr cytosines", 
    row_title_rot = 90, row_names_side = "left", row_names_centered = TRUE,
    column_title_side = "bottom", column_names_rot = 0, column_names_centered = TRUE,
    col = circlize::colorRamp2(colors = jcolors::jcolors_contin("pal12")(max(.)), breaks = seq(max(.))), heatmap_legend_param = list(at = c(4,8,16,32,64)), name = "running time (sec)"
  ) -> pl
pdf("/g/krebs/barzaghi/analyses/15.02.24_bioinformatics_application_note_figures/supplementary_fig1c.pdf", width = 5.5, height = 4.4)
pl
dev.off()

# D
Reduce(rbind, mclapply(seq(nrow(iteration_df)), function(i){
  
  print(i)
  
  nr_molecules = iteration_df$nr_molecules[i]
  nr_cytosines = iteration_df$nr_cytosines[i]
  
  memory_profiling = NULL
  
  while(is.null(memory_profiling)){
    
    print("memory profiling still null")
    
    MethSM_subset = list()
    MethSM_subset$amplicon_DE_data = MethSM$amplicon_DE_data[
      sample(rownames(MethSM$amplicon_DE_data), size = nr_molecules, replace = TRUE),
      sort(sample(colnames(MethSM$amplicon_DE_data), size = nr_cytosines, replace = TRUE))
    ]
    
    tryCatch({
      peakRAM({
        FootprintCharter(
          MethSM = MethSM_subset,
          RegionOfInterest = RegionOfInterest,
          verbose = FALSE
        ) -> FC_results
      })
    }, error = function(err){NULL}) -> memory_profiling
    
  }
  
  data.frame(
    nr_molecules = nr_molecules,
    nr_cytosines = nr_cytosines,
    peakRAM = memory_profiling$Peak_RAM_Used_MiB
  )
  
}, mc.preschedule = TRUE, mc.cores = 8)) -> results

results %>%
  mutate(peakRAM = ceiling(as.numeric(peakRAM*1.048576))) %>%
  spread(nr_cytosines, peakRAM) %>%
  mutate(nr_molecules = factor(nr_molecules, levels = rev(nr_molecules))) %>%
  arrange(nr_molecules) %>%
  column_to_rownames("nr_molecules") %>%
  as.matrix() %>%
  Heatmap(
    cluster_columns = FALSE, cluster_rows = FALSE, 
    row_title = "Nr molecules", column_title = "Nr cytosines", 
    row_title_rot = 90, row_names_side = "left", row_names_centered = TRUE,
    column_title_side = "bottom", column_names_rot = 0, column_names_centered = TRUE,
    col = circlize::colorRamp2(colors = jcolors::jcolors_contin("pal12")(max(.)), breaks = seq(max(.))), heatmap_legend_param = list(at = c(50,100,500,1000,4000)), name = "peak RAM (MB)"
  ) -> pl
pdf("/g/krebs/barzaghi/analyses/15.02.24_bioinformatics_application_note_figures/supplementary_fig1d.pdf", width = 5.5, height = 4.4)
pl
dev.off()