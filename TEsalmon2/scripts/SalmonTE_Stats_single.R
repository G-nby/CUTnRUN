message("Step 1: Loading required libraries...")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(stringr))  #Bo modified
suppressPackageStartupMessages(library(edgeR))  #gby modified

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
  ggsave(file.path(path, paste0("corr.name.", ext)), dat$summary$corr.name)
  ggsave(file.path(path, paste0("corr.name.log2fc.", ext)), dat$summary$corr.log2fc)
  ggsave(file.path(path, paste0("corr.clade.", ext)), dat$summary$corr.clade)
  ggsave(file.path(path, paste0("corr.class.", ext)), dat$summary$corr.class)
}


do.deseq2 <- function(dat, BCV = 0.1) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Please install package 'edgeR'")
  }
  count <- dat$count
  #message(head(count))
  #message(str(dat))
  #message(fixed_dispersion)
  if(ncol(count) != 2) stop("do.deseq2 expects count matrix with exactly 2 samples (ctrl,treat).")
  sample_names <- colnames(count)
  ctrl_name <- sample_names[1]
  treat_name <- sample_names[2]
  message(ctrl_name)
  message(treat_name)

  # build DGEList and normalize
  d <- edgeR::DGEList(counts = round(as.matrix(count)))
  #d$samples$group <- factor(c("ctrl","treat"))
  d$samples$group <- factor(c(ctrl_name, treat_name),
                          levels = c(ctrl_name, treat_name))
  d <- edgeR::calcNormFactors(d, method = "TMM")

  # normalized CPM
  norm_counts <- edgeR::cpm(d, normalized.lib.sizes = TRUE)
  norm_counts_df <- data.frame(norm_counts, check.names = FALSE)
  norm_counts_df <- tibble::rownames_to_column(norm_counts_df, var = "name")

  # compute log2FC (treat / ctrl) , pseudo for 0 cannot bu divisor
  pseudo <- 1
  log2FC_vec <- log2((norm_counts[, treat_name] + pseudo) / (norm_counts[, ctrl_name] + pseudo))

  # produce fake p-values via exactTest with fixed dispersion
  # BCV 0.4 for human , 0.1 for other model animals
  #message("BCV class = ", class(BCV))
  #message("BCV value = ", BCV)
  BCV <- as.numeric(BCV)  
  disp <- BCV^2
  et <- tryCatch({
    edgeR::exactTest(d, dispersion = disp)   # BCV^2 ²ÅÊÇ dispersion
  }, error = function(e) {
    message("edgeR exactTest failed: ", e$message)
    NULL
  })

  if(!is.null(et)) {
    tbl <- as.data.frame(edgeR::topTags(et, n = nrow(count), sort.by = "none")) # preserve order
    # topTags returns logFC as log2(treat/ctrl) if pair set correctly
    pval <- tbl$PValue
  } else {
    pval <- rep(NA, nrow(count))
    tbl <- data.frame(logFC = log2FC_vec, PValue = pval)
  }
  
  #message(str(et))
  #message(str(tbl))
  #message(str(pval))
  # assemble result dataframe
  df_res <- data.frame(
    name = rownames(count),
    log2FoldChange = as.numeric(log2FC_vec),
    pvalue = as.numeric(pval),
    baseMean = rowMeans(norm_counts)
  , stringsAsFactors = FALSE)
  #message(str(df_res))
  

  df_res$padj <- stats::p.adjust(df_res$pvalue, method = "fdr")

  dat$res <- df_res
  dat$norm_counts <- norm_counts_df
  dat$dds <- NULL

  message("do.deseq2 (edgeR fallback) done. Using BCV = ", BCV)
  dat
  #message(str(dat))
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


draw.corr_name <- function(dat, top_n = 30) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Please install package 'ggrepel' for better label display.")
  }
  # dat$norm_counts: data.frame name, sample1, sample2
  nc <- dat$norm_counts
  # assume columns: name, <ctrl>, <treat>
  sample_cols <- colnames(nc)[-1]
  ctrl <- sample_cols[1]; treat <- sample_cols[2]
  df <- nc %>% dplyr::rename(x = !!rlang::sym(ctrl), y = !!rlang::sym(treat)) %>%
    left_join(dat$res %>% dplyr::select(name, log2FoldChange), by = "name") %>%
    left_join(dat$annotation, by = "name") # ensure class/clade exist

  # correlation stats
  ct <- tryCatch(cor.test(df$x, df$y, method = "pearson"), error = function(e) NULL)
  rtxt <- if(!is.null(ct)) sprintf("r = %.3f\np = %.2g", ct$estimate, ct$p.value) else "r = NA\np = NA"

  # label topN by abs(log2FC)
  lab <- df %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% head(top_n)

  p <- ggplot(df, aes(x = x, y = y, color = class)) +
    geom_point(alpha = 0.7, size = 1.5) +
    ggrepel::geom_text_repel(data = lab, aes(label = name), size = 2.5, max.overlaps = 30) +
    theme_minimal() +
    labs(x = ctrl, y = treat, title = paste0("Correlation (name-level)")) +
    theme(legend.position = "bottom", text = element_text(size = 12)) +
    annotate("text", x = Inf, y = -Inf, label = rtxt, hjust = 1.1, vjust = -0.1, size = 3)

  p
}

draw.corr_log2fc <- function(dat, top_n = 30) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Please install package 'ggrepel'")
  }
  nc <- dat$norm_counts
  sample_cols <- colnames(nc)[-1]
  ctrl <- sample_cols[1]; treat <- sample_cols[2]

  df <- nc %>%
    dplyr::rename(raw_x = !!rlang::sym(ctrl),
                  raw_y = !!rlang::sym(treat)) %>%
    mutate(x = log2(raw_x + 1),
           y = log2(raw_y + 1)) %>% 
    left_join(dat$res %>% dplyr::select(name, log2FoldChange), by = "name") %>%
    left_join(dat$annotation, by = "name")

  # correlation on log scale
  ct <- tryCatch(cor.test(df$x, df$y, method = "pearson"), error = function(e) NULL)
  rtxt <- if(!is.null(ct)) sprintf("r = %.3f\np = %.2g", ct$estimate, ct$p.value) else "r = NA\np = NA"

  lab <- df %>%
    arrange(desc(abs(log2FoldChange))) %>%
    head(top_n)

  ggplot(df, aes(x = x, y = y, color = class)) +
    geom_point(alpha = 0.7, size = 1.5) +
    ggrepel::geom_text_repel(data = lab, aes(label = name), size = 2.5, max.overlaps = 30) +
    theme_minimal() +
    labs(x = paste0("log2(",ctrl,"+1)"),
         y = paste0("log2(",treat,"+1)"),
         title = "Correlation (log2 CPM)") +
    theme(legend.position = "bottom", text = element_text(size = 12)) +
    annotate("text", x = Inf, y = -Inf, label = rtxt,
             hjust = 1.1, vjust = -0.1, size = 3)
}

# draw aggregated correlation: level = "clade" or "class"
draw.corr_agg <- function(dat, level = c("clade","class")) {
  level <- match.arg(level)
  nc <- dat$norm_counts
  sample_cols <- colnames(nc)[-1]
  ctrl <- sample_cols[1]; treat <- sample_cols[2]
  df <- nc %>% dplyr::rename(x = !!rlang::sym(ctrl), y = !!rlang::sym(treat)) %>%
    left_join(dat$annotation, by = "name")

  agg <- df %>% group_by_at(level) %>%
    summarise(x = mean(x, na.rm = TRUE), y = mean(y, na.rm = TRUE)) %>% ungroup()

  ct <- tryCatch(cor.test(agg$x, agg$y, method = "pearson"), error = function(e) NULL)
  rtxt <- if(!is.null(ct)) sprintf("r = %.3f\np = %.2g", ct$estimate, ct$p.value) else "r = NA\np = NA"

  p <- ggplot(agg, aes(x = x, y = y, color = .data[[level]])) +
    geom_point(size = 2) +
    ggrepel::geom_text_repel(aes(label = .data[[level]]), size = 3) +
    theme_minimal() +
    labs(x = ctrl, y = treat, title = paste0("Correlation (", level, "-level)")) +
    theme(legend.position = "none", text = element_text(size = 12)) +
    annotate("text", x = Inf, y = -Inf, label = rtxt, hjust = 1.1, vjust = -0.1, size = 3)

  p
}

draw.MAplot <- function(dat) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Please install package 'ggrepel' for better label display.")
  }
  
# prepare data
  plot_data <- dat$res %>%
    filter(!is.na(padj) & class != "Simple_repeat") %>%  
    mutate(
      #sig = ifelse(padj < 0.05 & abs(!!sym(dat$y.name)) > 0.5, "DE", "NC"),
      sig = ifelse(padj < dat$pvalmax & abs(!!sym(dat$y.name)) > dat$log2fc_min, "DE", "NC"),  
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
      #sig = ifelse(pvalue < 0.05 & abs(!!sym(dat$y.name)) > 0.5, "DE", "NC"),
      sig = ifelse(pvalue < dat$pvalmax & abs(!!sym(dat$y.name)) > dat$log2fc_min, "DE", "NC"),  
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
      #sig = ifelse(padj < 0.05 & abs(!!sym(dat$y.name)) > 0.5, "DE", "NC"),
      sig = ifelse(padj < dat$pvalmax & abs(!!sym(dat$y.name)) > dat$log2fc_min, "DE", "NC"),
      negLog10Padj = -log10(padj),
      group = ifelse(sig == "DE", as.character(clade), "Non-DE")
    )

  p <- ggplot(plot_data, aes_string(x = dat$y.name, y = "negLog10Padj")) +
    geom_point(aes(color = group), alpha = 0.7) +
    #geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = "dashed") +
    #geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
    geom_vline(xintercept = c(-dat$log2fc_min, dat$log2fc_min), col = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(dat$pvalmax), col = "red", linetype = "dashed") +
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
      #sig = ifelse(pvalue < 0.05 & abs(!!sym(dat$y.name)) > 0.5, "DE", "NC"),
      sig = ifelse(pvalue < dat$pvalmax & abs(!!sym(dat$y.name)) > dat$log2fc_min, "DE", "NC"),
      negLog10Pval = -log10(pvalue),
      group = ifelse(sig == "DE", as.character(clade), "Non-DE")
    )

  p <- ggplot(plot_data, aes_string(x = dat$y.name, y = "negLog10Pval")) +
    geom_point(aes(color = group), alpha = 0.7) +
    #geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = "dashed") +
    #geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
    geom_vline(xintercept = c(-dat$log2fc_min, dat$log2fc_min), col = "red", linetype = "dashed") +
    geom_hline(yintercept = -log10(dat$pvalmax), col = "red", linetype = "dashed") +
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



SalmonTE <- function(count, samples, annotation, col_data,
                     BCV,
                     sheet.fmt = "csv", fig.fmt = "pdf", path = ".", log2fc_min = 0.5,
                     pvalmax = 0.05) {
  col_data <- col_data[samples, , drop = FALSE]
  count <- count[, samples, drop = FALSE]
  dat <- list(
    count = count, 
    col_data = col_data
  )
  dat$log2fc_min <- log2fc_min
  dat$pvalmax    <- pvalmax

  if(!is.null(samples)) {
    rownames(dat$col_data) <- factor(rownames(dat$col_data), level = samples)
  }
  dat <- do.deseq2(dat,BCV)
  dat$y.name <- "log2FoldChange"
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
  dat$summary$corr.name  <- draw.corr_name(dat, top_n = 30)
  dat$summary$corr.log2fc  <- draw.corr_log2fc(dat, top_n = 30)
  dat$summary$corr.clade <- draw.corr_agg(dat, level = "clade")
  dat$summary$corr.class <- draw.corr_agg(dat, level = "class")
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
if(!is.na(args[5])) {
  samples <- str_split(args[5], ",", simplify = T)
}
BCV = args[6]
log2fc_min = as.numeric(args[7])
pvalmax       = as.numeric(args[8])
message(log2fc_min)
message(pvalmax)
count <- read.csv(file.path(args[1], "EXPR.csv"), row.names="TE", check.names = FALSE, header = TRUE)  #header = TRUE
col_data <- read.csv(file.path(args[1], "condition.csv"), row.names = "SampleID")
keep_samples <- rownames(col_data)[rownames(col_data) %in% samples] #filter wanted cols  # gby added
message(keep_samples)
count <- count[, keep_samples, drop = FALSE]  # gby added
col_data <- col_data[keep_samples, , drop = FALSE]  # gby added
annotation <- read.csv(file.path(args[1], "clades.csv"))
#message(sprintf("Step 3: Running the %s analysis...", analysis))
dat <- SalmonTE(count, samples, annotation, col_data, BCV, args[2], args[3], args[4],log2fc_min, pvalmax)

message(sprintf("Step 4: Generating output..."))
suppressMessages(GenerateOutput(dat))
message(sprintf("Step 5: The statistical analysis has been completed. Please check '%s' directory to see the analysis result!", args[4]))
