# Efficiency metric calculations
calculate_singletons <- function(rbs) {
    total_reads <- sum(rbs$x, rbs$y)
    singletons <- sum(rbs$x == 1 & rbs$y == 0 | rbs$x == 0 & rbs$y == 1)
    frac_singletons <- singletons / total_reads
    return(frac_singletons)
}

calculate_family_stats <- function(rbs) {
    rbs$size <- rbs$x + rbs$y
    return(c(total_families = nrow(rbs),
             family_mean = mean(rbs$size),
             family_median = median(rbs$size),
             family_max = max(rbs$size),
             families_gt1 = sum(rbs$x > 1 | rbs$y > 1),
             single_families = sum(rbs$x == 1 & rbs$y == 0 | rbs$x == 0 & rbs$y == 1),
             paired_families = sum(rbs$x > 0 & rbs$y > 0),
             paired_and_gt1 = sum(rbs$x > 1 & rbs$y > 1)))
}

calculate_metrics <- function(rbs) {
    metrics <- data.frame(sample = names(rbs))
    metrics$frac_singletons <- lapply(rbs, calculate_singletons) %>% unlist()
    metrics$efficiency <- lapply(rbs, calculate_efficiency) %>% unlist()
    metrics$drop_out_rate <- lapply(rbs, calculate_missed_fraction) %>% unlist()
    metrics <- lapply(rbs, calculate_gc) %>%
                    data.frame %>% t() %>%
                    cbind(metrics, .)
    metrics <- lapply(rbs, calculate_family_stats) %>%
                    data.frame %>% t() %>%
                    cbind(metrics, .)
    return(metrics)
}

calculate_metrics_single <- function(rbs) {
    metrics <- NULL
    metrics$frac_singletons <- calculate_singletons(rbs)
    metrics$efficiency <- calculate_efficiency(rbs)
    metrics$drop_out_rate <- calculate_missed_fraction(rbs)

    gc_stats <- calculate_gc(rbs)
    family_stats <- calculate_family_stats(rbs)
    metrics <- data.frame(metrics, t(gc_stats))
    metrics <- data.frame(metrics, t(family_stats))

    return(metrics)
}

calc_metrics_new_rbs <- function(rinfo_dir, pattern="\\.txt.gz", cores=8) {
    metrics <-
        list.files(
            rinfo_dir,
            full.names = TRUE,
            recursive = TRUE,
            pattern = pattern
        ) %>%
        mclapply(., fread, mc.cores=cores) %>%
        mclapply(., calculate_metrics_single, mc.cores=cores)

    return(metrics)
}

# functions below are adapted from
# R/efficiency_nanoseq.R and perl/efficiency_nanoseq.pl
# from https://github.com/cancerit/NanoSeq

# from cancerit/NanoSeq documentation:
# "This is the number of duplex bases divided by the number of sequenced bases."
calculate_efficiency <- function(rbs) {
    bases_ok_rbs <- nrow(rbs[rbs$x > 1 & rbs$y > 1,]) * ((rlen - skips) * 2)
    total_reads <- sum(c(rbs$x, rbs$y))
    bases_sequenced <- total_reads * rlen * 2
    eff <- bases_ok_rbs / bases_sequenced
    return(eff)
}

# from cancerit/NanoSeq documentation:
# "This shows the fraction of read bundles missing one of the two
# original strands beyond what would be expected under random sampling
# (assuming a binomial process).
calculate_missed_fraction <- function(rbs) {
    rbs <- data.frame(rbs)
    rbs$size <- rbs$x + rbs$y
    rbs$size <- pmin(rbs$size, 10)
    total_missed <- 0
    for(size in c(4:10)) {
        exp_orphan <- (0.5 ** size) * 2
        total_this_size = nrow(rbs[which(rbs$size == size),])
        if(total_this_size > 0) {
            with_both_strands <- nrow(rbs[which(rbs$size == size & rbs$x > 0 & rbs$y > 0),])
            obs_orphan <- 1 - with_both_strands / total_this_size
            missed <- (obs_orphan - exp_orphan) * total_this_size
            total_missed <- total_missed + missed
        }
    }
    total_missed_fraction = total_missed / nrow(rbs[which(rbs$size >= 4),])
    return(total_missed_fraction)
}

# from cancerit/NanoSeq documentation:
# The GC content of RBs with both strands and with just one strand.
# I return the difference between the two values.
calculate_gc <- function(rbs, sample_n = 10000, max_gap = 100000) {
    rbs <- data.frame(rbs)
    colnames(rbs)[5:6] = c('plus', 'minus')

    # remove any chroms not in the sizes vector
    # rbs <- rbs[rbs$chrom %in% names(genome_max),]

    # remove records with mate positions exceeding genome max
    rbs$end <- rbs$mpos + rlen - skips
    if(length(genome_max) > 1) {
        for(i in 1:length(genome_max)) {
            chrom <- names(genome_max[i])
            chrom_max <- genome_max[i]
            rbs <- rbs[!(rbs$chrom == chrom & rbs$end > chrom_max),]
        }
    } else {
        rbs <- rbs[!rbs$end > genome_max,]
    }

    # remove records with large distances between mates
    rbs <- rbs[rbs$end - rbs$pos < max_gap,]

    # remove zero-records
    rbs <- rbs[rbs$pos != 0,]

    rbs_both <- rbs[which(rbs$minus + rbs$plus >= 4 & rbs$minus >= 2 & rbs$plus >= 2),]
    rbs_both <- rbs_both[sample(1:nrow(rbs_both), min(sample_n, nrow(rbs_both))),]

    rbs_single <- rbs[which(rbs$minus + rbs$plus > 4 & (rbs$minus == 0 | rbs$plus == 0)),]
    rbs_single <- rbs_single[sample(1:nrow(rbs_single), min(sample_n, nrow(rbs_single))),]

    seqs_both <- GRanges(rbs_both$chrom,
                         IRanges(start=rbs_both$pos,
                                 end=rbs_both$end)) %>%
        scanFa(genomeFile, .) %>%
        as.vector()

    seqs_single <- GRanges(rbs_single$chrom,
                           IRanges(start=rbs_single$pos,
                                   end=rbs_single$end)) %>%
        scanFa(genomeFile, .) %>%
        as.vector()

    seqs_both_collapsed <- paste(seqs_both, collapse = '')
    seqs_single_collapsed <- paste(seqs_single, collapse = '')

    tri_both <- DNAString(seqs_both_collapsed) %>% trinucleotideFrequency(.)
    tri_single <- DNAString(seqs_single_collapsed) %>% trinucleotideFrequency(.)

    tri_both_freqs <- tri_both  / sum(tri_both)
    tri_single_freqs <- tri_single / sum(tri_single)

    gc_both <- s2c(seqs_both_collapsed) %>% GC()
    gc_single <- s2c(seqs_single_collapsed) %>% GC()

    return(c(gc_single=gc_single, gc_both=gc_both, gc_deviation=abs(gc_single - gc_both)))
}
