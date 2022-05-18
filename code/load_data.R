source(here('code/efficiency_nanoseq_functions.R'))

load_data <- function(fdir, pattern, samples, read_func=read.delim) {
    df <- list.files(
        fdir,
        full.names = TRUE,
        recursive = TRUE,
        pattern = pattern) %>%
        lapply(., read_func)

    for (i in seq(length(samples))) {
        df[[i]]$Sample <- samples[i]
    }
    df <- rbindlist(df)

    # add nuclease + protocol info
    df$protocol <- 'NanoSeq'
    df$protocol[grep('Nux', df$Sample)] <- 'xGen'
    df$nuclease <- str_split(df$Sample, regex("N(uxg|an)|_")) %>%
        lapply(., tail, 2) %>% lapply(., dplyr::first) %>% unlist()

    return(df)
}


load_variants <- function(fdir, sample_names) {
    vars <- list.files(
        fdir,
        full.names = TRUE,
        recursive = TRUE,
        pattern = '.vcf') %>%
        lapply(., read.vcfR, verbose = FALSE)

    var_df <- NULL
    for (i in 1:length(vars)) {
        svars <- vars[[i]]@fix %>% data.frame()
        if (nrow(svars) > 0) {
            svars$sample <- sample_names[i]
            var_df <- rbind(var_df, svars)
        }
    }
    var_df$id <- paste(var_df$CHROM, var_df$POS, sep = '_')

    return(var_df)
}

extract_std <- function(genome_results) {
    std <- genome_results[grep('std', genome_results$BamQC.report),] %>%
        strsplit(., split='=') %>%
        last() %>% last() %>%
        gsub(' |X', '', .) %>% as.numeric()
    return(std)
}

load_cov_stats <- function(cov, qualimap_dir, samples) {
    cov_stats <- list.files(
        qualimap_dir,
        full.names = TRUE,
        recursive = TRUE,
        pattern = 'genome_results.txt') %>%
        lapply(., read.delim) %>%
        lapply(., extract_std) %>%
        unlist()

    cov_stats <- data.frame(cov_std=cov_stats, Sample=samples)
    cov_stats <- data.table(cov)[,mean(Coverage), by=Sample] %>%
        data.frame() %>% inner_join(., cov_stats, by='Sample')
    colnames(cov_stats)[2] <- 'cov_mean'

    return(cov_stats)
}

load_nanoseq_stats <- function(nanoseq_dir) {
    tsvs <- list.files(nanoseq_dir,
                       pattern = '.tsv',
                       full.names = TRUE) %>%
        grep('GC', ., value = TRUE, invert = TRUE) %>%
        lapply(., read.delim)

    for (i in seq(length(samples))) {
        tsvs[[i]]$Sample <- samples[i]
        tsvs[[i]]$metric <- rownames(tsvs[[i]])
    }
    tsvs <- rbindlist(tsvs)
    colnames(tsvs)[1] <- 'value'
    tsvs <- tsvs[!is.na(tsvs$value)]

    # add GC deviation
    tmp <- tsvs[grep('GC', tsvs$metric)]
    gc <- tmp[tmp$metric == 'GC_SINGLE']
    gc$value <- abs(tmp[tmp$metric == 'GC_BOTH']$value - gc$value)
    gc$metric <- 'GC_DEVIATION'
    tsvs <- rbind(tsvs, gc)

    # extract protocol and nuclease labels
    tsvs$protocol <- 'NanoSeq'
    tsvs$protocol[grep('Nux', tsvs$Sample)] <- 'xGen'
    tsvs$nuclease <- str_split(tsvs$Sample, regex("N(uxg|an)|_")) %>%
        lapply(., tail, 2) %>% lapply(., dplyr::first) %>% unlist()

    return(tsvs)
}

load_rbs_data <- function(rinfo_dir) {
    rbs <- list.files(
            rinfo_dir,
            full.names = TRUE,
            recursive = TRUE,
            pattern = '\\.txt.gz') %>%
        lapply(., fread) %>%
        lapply(., get_rbs)

    return(rbs)
}

get_rbs <- function(x) {
    rb_plus <- x[x$strand == '+',]
    rb_minus <- x[x$strand == '-',]

    rbs <- full_join(rb_plus, rb_minus, by = c('chrom', 'pos', 'mpos', 'umi')) %>% distinct()
    rbs <- rbs[,c('chrom', 'pos', 'mpos', 'umi', '0.x', '0.y')]
    colnames(rbs)[5:6] <- c('x', 'y')
    rbs$x[is.na(rbs$x)]<- 0; rbs$y[is.na(rbs$y)]<- 0;

    return(rbs)
}

load_markdup_data <- function(markdup_dir, sample_names) {
    mdup <- list.files(
        markdup_dir,
        full.names = TRUE,
        recursive = TRUE,
        pattern = 'txt') %>%
        paste('grep -E "Library|LIBRARY"', .) %>%
        lapply(., fread) %>%
        suppressMessages()

    names(mdup) <- sample_names
    mdup <- mdup %>% lapply(., function(x){x$PERCENT_DUPLICATION}) %>% unlist()

    return(mdup)
}

get_qmap_coverage <- function(qualimap_dir, sample_names) {
    qmap <- load_data(qualimap_dir, 'genome_results.txt', sample_names)
    qmap <- qmap[qmap$BamQC.report %like% 'mean cov',]
    qmap$coverage <-
        str_split(qmap$BamQC.report, ' = ') %>%
        lapply(., last) %>%
        unlist() %>%
        gsub('X|,', '', .) %>%
        lapply(., as.numeric) %>%
        unlist()
    qmap <- qmap[,c('Sample', 'coverage')]
    return(qmap)
}
