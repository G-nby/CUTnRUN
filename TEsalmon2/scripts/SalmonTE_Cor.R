message("Step 1: Loading required libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(stringr))  #Bo modified

write.results <- function(dat) {
  res <- dat$res
  summary <- dat$summary
  sheet.fmt <- dat$sheet.fmt
  path <- dat$path
  norm_counts_df <- dat$norm_counts # gby modified
  norm_counts_res <- left_join(norm_counts_df, res, by = "name") # gby modified
  res_filtered <- res[res$class != "Simple_repeat", ] # gby modified0924
  norm_counts_res_filtered <- norm_counts_res[norm_counts_res$class != "Simple_repeat", ] # gby modified0924
  if(sheet.fmt == "tsv") {
    write.table(norm_counts_res_filtered, file=file.path(path,"normalized_counts.tsv"), sep="\t", row.names = F) # gby modified
    write.table(res_filtered, file=file.path(path,"results.tsv"), sep="\t",  row.names = F)
    write.table(summary$clade$table, file=file.path(path,"results.clade.tsv"), sep="\t", row.names = F)
    write.table(summary$class$table, file=file.path(path,"results.class.tsv"), sep="\t", row.names = F)
  } else if (sheet.fmt == "csv") {
    write.csv(norm_counts_res_filtered, file=file.path(path,"normalized_counts.csv"), row.names = F)
    write.csv(res_filtered, file=file.path(path,"results.csv"), row.names = F)
    write.csv(summary$clade$table, file=file.path(path,"results.clade.csv"), row.names = F)
    write.csv(summary$class$table, file=file.path(path,"results.class.csv"), row.names = F)
  } else {
    WriteXLS(norm_counts_res_filtered, ExcelFileName=file.path(path,"normalized_counts.xls"), row.names=F)
    WriteXLS(res_filtered, ExcelFileName=file.path(path,"results.xls"), row.names=F)
    WriteXLS(summary$clade$table, ExcelFileName=file.path(path,"results.clade.xls"), row.names=F)
    WriteXLS(summary$class$table, ExcelFileName=file.path(path,"results.class.xls"), row.names=F)
  }
}

write.figures <- function(dat) {
  path <- dat$path
  ext <- dat$fig.fmt
  ggsave(file.path(path, paste0("ma.plot.", ext)), dat$summary$ma.plot)
  ggsave(file.path(path, paste0("clade.", ext)), dat$summary$clade$figure)
  ggsave(file.path(path, paste0("class.", ext)), dat$summary$class$figure)
  ggsave(file.path(path, paste0("volcano.plot.", ext)), dat$summary$volcano.plot)
  ggsave(file.path(path, paste0("pdist.plot.", ext)), dat$summary$pdist.plot)
  ggsave(file.path(path, paste0("pval_volcano.plot.", ext)), dat$summary$pval_volcano.plot)
  ggsave(file.path(path, paste0("pval_ma.plot.", ext)), dat$summary$pval_ma.plot)
}


#do.deseq2 <- function(dat) {
#  count <- dat$count
#  col_data <- dat$col_data
#  dds <- DESeqDataSetFromMatrix(countData = round(count),
#                                colData = col_data,
#                                design = ~condition)
#  dds <- dds[ rowSums(counts(dds)) > 10, ]
#  dds <- DESeq(dds)
#  res <- results(dds)
#  df_res <- data.frame(res)
#  dat$res <- df_res %>% rownames_to_column("name")
#  dat$dds <- dds
#  
#  message(sprintf("herehereherehereherehereherehere:"))
#  message("Available resultsNames: ", paste(resultsNames(dds), collapse = ", "))
#  
#  res$padj < 0.01 & abs(res$log2FoldChange) > 0.5
#  dat
#}

# gby modified,for valid with chushai dispersion error
do.deseq2 <- function(dat) {
  count <- dat$count
  samples <- dat$samples
  dds <- DESeqDataSetFromMatrix(countData = round(count),
                                colData = col_data,
                                design = ~condition)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  
  dds <- tryCatch(
    {
      DESeq(dds)
    },
    error = function(e) {
      message("Caught error message:")
      message(e$message)
      
      if (grepl("invalid 'x'", e$message)) {
      message("!! Loess failed, retrying with fitType='mean' ...")
      dds <- DESeq(dds, fitType="mean")
      return(dds)
      }
      
      if (grepl("newsplit: out of vertex space", e$message)) {
      message("!! I guess Locfit failed, retrying with fitType='mean' ...")
      dds <- DESeq(dds, fitType="mean")
      return(dds)
      }
      
      if (grepl("all gene-wise dispersion estimates", e$message)) {
        message("!!Dispersion fit failed, using gene-wise estimates instead...")
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        dds <- nbinomWaldTest(dds)
        return(dds)
      } 
      
      else {
        stop(e)
      }
    }
  )
  
  res <- results(dds)
  df_res <- data.frame(res)
  
  norm_counts <- counts(dds, normalized = TRUE)
  norm_counts_df <- data.frame(norm_counts)
  
  
  dat$res <- df_res %>% rownames_to_column("name")
  
  dat$norm_counts <- norm_counts_df %>% rownames_to_column("name")
  
  dat$dds <- dds
  
  message("herehereherehereherehereherehere:")
  message("Available resultsNames: ", paste(resultsNames(dds), collapse = ", "))
  
  res$padj < 0.01 & abs(res$log2FoldChange) > 0.5
  dat
}



do.summary <- function(dat) {
  calc.stat <- function(group) {
    tmp <- data.frame(group=dat$res[,group], 
                      value=dat$res[,dat$y.name]) 
    colnames(tmp) <- c("group", "value")
    tmp <- tmp %>%
      filter(!is.na(value))
    tb <- tmp %>%
      group_by(group) %>% 
      summarise( n = n(), 
                 mean = mean(value, na.rm=T), 
                 sd = sd(value, na.rm=T),
                 p.value = tryCatch({t.test(value)$p.value}, 
                                    error = function(e) NA ) )
    fig <- tmp %>% 
      ggplot(aes(x=group, y=value)) +
      geom_boxplot(width=0.5) + 
      geom_jitter(width=0.1, alpha=0.5) + theme(legend.position = "none") +
      geom_hline(yintercept = 0, color='red') +
      xlab(group) + ylab(dat$y.name) + theme_minimal() + 
      theme(text = element_text(size = 12)) +
      theme(legend.position = "none") + 
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    return(list(table = tb, figure = fig))
  }
  dat$summary <- list(class = calc.stat("class"),
                      clade = calc.stat("clade"))
  return(dat)
}


# a code to draw MA-plot
#draw.MAplot <- function(dat) {
#  return(dat$res %>%
#    filter( !is.na(padj) ) %>%
#    mutate( sig = ifelse(padj < 0.05, "DE", "NC" )) %>%
#    ggplot( aes_string(x="baseMean", y=dat$y.name))+
#    geom_hline(yintercept = 0, col = "red", alpha=0.5) +
#    geom_point( aes(colour=sig) ) + 
#    scale_colour_manual(values = c("red", "black"), limits = c("DE", "NC")) +
#    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                  labels = trans_format("log10", math_format(10^.x)) ) +
#    theme(legend.position = "none", text = element_text(size = 18))  + theme_minimal())
#}

draw.MAplot <- function(dat) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Please install package 'ggrepel' for better label display.")
  }
  
# prepare data
  plot_data <- dat$res %>%
    filter(!is.na(padj) & class != "Simple_repeat") %>%  
    mutate(
      sig = ifelse(padj < 0.05 & abs(!!sym(dat$y.name)) > 0.5, "DE", "NC"),  
      group = ifelse(sig == "DE", as.character(clade), "Non-DE") 
    )
  
  p <- ggplot(plot_data, aes_string(x = "baseMean", y = dat$y.name)) +
    geom_hline(yintercept = 0, col = "red", alpha = 0.5) +  
    

    geom_point(aes(color = group), alpha = 0.7) +
    
    ggrepel::geom_text_repel(
      aes(
        label = ifelse(sig == "DE", name, NA),
        colour = group, 
      ),
      size = 3,
      max.overlaps = 50,      # max overlap tag
      min.segment.length = 0.2,
      bg.color = "white", 
      bg.r = 0.15, 
    ) +
    
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    
    theme_minimal() +
    theme(
      text = element_text(size = 18),             
      legend.position = "bottom",                 
      legend.title = element_blank()              
    ) +
    labs(x = "Mean Expression", y = dat$y.name)   
  

  clades <- unique(plot_data$group[plot_data$group != "Non-DE"])
  
  color_values <- c("Non-DE" = "grey70") 
  
  # generate color for clades
  if(length(clades) > 0) {
    clade_colors <- scales::hue_pal()(length(clades))
    names(clade_colors) <- clades
    color_values <- c(color_values, clade_colors)
  }
  
  p + scale_color_manual(values = color_values)
}

draw.pval_MAplot <- function(dat) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Please install package 'ggrepel' for better label display.")
  }
  
# prepare data
  plot_data <- dat$res %>%
    filter(!is.na(padj) & class != "Simple_repeat") %>%  
    mutate(
      sig = ifelse(pvalue < 0.05 & abs(!!sym(dat$y.name)) > 0.5, "DE", "NC"),  
      group = ifelse(sig == "DE", as.character(clade), "Non-DE") 
    )
  
  p <- ggplot(plot_data, aes_string(x = "baseMean", y = dat$y.name)) +
    geom_hline(yintercept = 0, col = "red", alpha = 0.5) +  
    

    geom_point(aes(color = group), alpha = 0.7) +
    
    ggrepel::geom_text_repel(
      aes(
        label = ifelse(sig == "DE", name, NA),
        colour = group, 
      ),
      size = 3,
      max.overlaps = 50,      # max overlap tag
      min.segment.length = 0.2,
      bg.color = "white", 
      bg.r = 0.15, 
    ) +
    
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    
    theme_minimal() +
    theme(
      text = element_text(size = 18),             
      legend.position = "bottom",                 
      legend.title = element_blank()              
    ) +
    labs(x = "Mean Expression", y = dat$y.name)   
  

  clades <- unique(plot_data$group[plot_data$group != "Non-DE"])
  
  color_values <- c("Non-DE" = "grey70") 
  
  # generate color for clades
  if(length(clades) > 0) {
    clade_colors <- scales::hue_pal()(length(clades))
    names(clade_colors) <- clades
    color_values <- c(color_values, clade_colors)
  }
  
  p + scale_color_manual(values = color_values)
}

# volcanoplot function
draw.volcanoplot <- function(dat) {
  plot_data <- dat$res %>%
    filter(!is.na(padj) & class != "Simple_repeat") %>%
    mutate(
      sig = ifelse(padj < 0.05 & abs(!!sym(dat$y.name)) > 0.5, "DE", "NC"),
      negLog10Padj = -log10(padj),
      group = ifelse(sig == "DE", as.character(clade), "Non-DE")
    )

  p <- ggplot(plot_data, aes_string(x = dat$y.name, y = "negLog10Padj")) +
    geom_point(aes(color = group), alpha = 0.7) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
    ggrepel::geom_text_repel(
      aes(label = ifelse(sig == "DE", name, NA), color = group),
      size = 3, max.overlaps = 50,
      min.segment.length = 0.2,
      bg.color = "white", bg.r = 0.15
    ) +
    theme_minimal() +
    labs(x = dat$y.name, y = expression(-log[10](padj))) +
    theme(text = element_text(size = 18),
          legend.position = "bottom",
          legend.title = element_blank())

  clades <- unique(plot_data$group[plot_data$group != "Non-DE"])
  color_values <- c("Non-DE" = "grey70")
  if(length(clades) > 0) {
    clade_colors <- scales::hue_pal()(length(clades))
    names(clade_colors) <- clades
    color_values <- c(color_values, clade_colors)
  }

  p + scale_color_manual(values = color_values)
}

draw.pval_volcanoplot <- function(dat) {
  plot_data <- dat$res %>%
    filter(!is.na(pvalue) & class != "Simple_repeat") %>%
    mutate(
      sig = ifelse(pvalue < 0.05 & abs(!!sym(dat$y.name)) > 0.5, "DE", "NC"),
      negLog10Pval = -log10(pvalue),
      group = ifelse(sig == "DE", as.character(clade), "Non-DE")
    )

  p <- ggplot(plot_data, aes_string(x = dat$y.name, y = "negLog10Pval")) +
    geom_point(aes(color = group), alpha = 0.7) +
    geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
    ggrepel::geom_text_repel(
      aes(label = ifelse(sig == "DE", name, NA), color = group),
      size = 3, max.overlaps = 50,
      min.segment.length = 0.2,
      bg.color = "white", bg.r = 0.15
    ) +
    theme_minimal() +
    labs(x = dat$y.name, y = expression(-log[10](pvalue))) +
    theme(text = element_text(size = 18),
          legend.position = "bottom",
          legend.title = element_blank())

  clades <- unique(plot_data$group[plot_data$group != "Non-DE"])
  color_values <- c("Non-DE" = "grey70")
  if(length(clades) > 0) {
    clade_colors <- scales::hue_pal()(length(clades))
    names(clade_colors) <- clades
    color_values <- c(color_values, clade_colors)
  }

  p + scale_color_manual(values = color_values)
}

# draw pdistplot to check the p distribution sothat we can know is the padj is suitable to apply
draw.pdistplot <- function(dat) {
  plot_data <- dat$res %>% 
    filter(!is.na(pvalue))  

  ggplot(plot_data, aes(x = pvalue)) +
    geom_histogram(binwidth = 0.02, fill = "steelblue", color = "white", boundary = 0) +
    geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = "P-value Distribution",
         x = "Raw P-value",
         y = "Frequency") +
    theme(text = element_text(size = 18))
}



SalmonTE <- function(count, samples, annotation,
                     sheet.fmt = "csv", fig.fmt = "pdf", path = ".") {
  dat <- list(
    count = count, 
    samples = samples
  )
  if(!is.null(samples)) {
    dat$SampleID <- factor(dat$SampleID, level = samples)
  }
  

  if( analysis == "DE" ) {
    if(!is.null(condition_level)) {
      dat$SampleID <- factor(dat$SampleID, level = samples)
    }
    dat <- do.deseq2(dat)
    dat$y.name <- "log2FoldChange"
  } else {
    dat$res <- do.lm(dat)
    dat$y.name <- "b.value"
  }
  dat$annotation <- annotation
  dat$sheet.fmt <- sheet.fmt
  dat$path <- path
  dat$fig.fmt <- fig.fmt
  dat$res <- suppressWarnings(left_join(dat$res, dat$annotation, by="name"))
  nc <- ncol(dat$res)
  dat$res <- dat$res[, c(1, nc-1, nc, which(!(1:nc %in% c(1,nc-1,nc))))]
  
  dat <- do.summary(dat) # get summary in clade and class, and got fig here, later draw&save in write_figures
  dat$summary$ma.plot <- draw.MAplot(dat) # get ma.plot here, and later draw&save in write_figures
  dat$summary$volcano.plot <- draw.volcanoplot(dat) #gby added
  dat$summary$pdist.plot <- draw.pdistplot(dat) # gby added
  dat$summary$pval_volcano.plot <- draw.pval_volcanoplot(dat) #gby added
  dat$summary$pval_ma.plot <- draw.pval_MAplot(dat) # gby added
  dat
}

GenerateOutput <- function(dat) {
  write.results(dat)
  write.figures(dat)
  save(dat, file = file.path(dat$path, "data.Rdata"))
}

message("Step 2: Loading input data...")
args <- commandArgs(T)

#analysis <- args[5]
samples <- NULL
if(!is.na(args[6])) {
  samples <- str_split(args[6], ",", simplify = T)
}
count <- read.csv(file.path(args[1], "EXPR.csv"), row.names="TE", check.names = FALSE, header = TRUE)  #header = TRUE
col_data <- read.csv(file.path(args[1], "condition.csv"), row.names = "SampleID")
keep_samples <- rownames(col_data)[col_data$SampleID %in% samples] #filter wanted cols  # gby added
message(keep_samples)
count <- count[, keep_samples, drop = FALSE]  # gby added
col_data <- col_data[keep_samples, , drop = FALSE]  # gby added
annotation <- read.csv(file.path(args[1], "clades.csv"))
#message(sprintf("Step 3: Running the %s analysis...", analysis))
dat <- SalmonTE(count, samples, annotation, args[2], args[3], args[4])

message(sprintf("Step 4: Generating output..."))
suppressMessages(GenerateOutput(dat))
message(sprintf("Step 5: The statistical analysis has been completed. Please check '%s' directory to see the analysis result!", args[4]))
