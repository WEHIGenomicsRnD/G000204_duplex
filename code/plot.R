plot_metric <- function(tsvs, metric, title) {
    tmp <- data.frame(tsvs)[grep(metric, tsvs$metric),]
    p <- ggplot(tmp, aes(Sample, value, fill = metric)) +
        geom_histogram(stat = 'identity', position = 'dodge') +
        theme_bw() +
        coord_flip() +
        ggtitle(title)
    return(p)
}

plot_metric_boxplot <- function(tsvs, x, metric, title) {
    tmp <- data.frame(tsvs)[grep(metric, tsvs$metric),]
    p <- ggplot(tmp, aes_string(x, 'value')) +
        geom_boxplot() +
        theme_bw() +
        ggtitle(title)
    return(p)
}
