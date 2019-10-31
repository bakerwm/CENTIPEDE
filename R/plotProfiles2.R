
#' Create plots 
#' @param Profile matrix of probability
#' @param Mlen length of the motif
#' @param xlab character, 
#' @param ylab character,
#' @param mainTitle character,
#' @param legPos where the legend
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export 
plotProfile2 <- function(Profile,
                         Mlen=0,
                         xlab="Dist. to motif (bp)",
                         ylab="Cut-site probability",
                         mainTitle="",
                         legPos="topright"){
  # #using ggplot2
  # library(ggplot2)
  # library(dplyr)
  df <- as.data.frame(Profile)
  colnames(df) <- "score"
  
  # calculate: flank, width
  flank <- nrow(df) / 4 - Mlen / 2
  
  ## labels on x-axis
  # 1, 50, 101, 119, 169, 219
  breaks1 <- seq(1, flank + 1, length.out = 3)
  breaks2 <- breaks1 + flank + Mlen - 1
  breaks  <- c(breaks1, breaks2)
  # -100, -50, 0, 0, 50, 100
  tag <- seq(0, flank, length.out = 3)
  tags  <- c(-rev(tag), tag)
  
  ## prepare data
  coords <- seq_len(nrow(df)/2)
  df <- df %>%
    dplyr::mutate(position = c(coords, rev(coords)),
                  strand   = rep(c("fwd", "rev"), each = nrow(df) / 2))
  
  ## plot
  p <- ggplot(df, aes(position, score, color = strand)) +
    geom_line(size = 0.7) +
    geom_vline(xintercept = c(flank + 1, flank + Mlen), linetype = 2,
               color = "grey40") +
    annotate("segment", 
             x    = c((flank +1):(flank + Mlen)), 
             xend = c((flank +1):(flank + Mlen)),
             y    = 0, 
             yend = max(df$score) * 0.01,
             colour = "blue") +
    scale_x_continuous(breaks = breaks,
                       labels = tags, 
                       minor_breaks = c(101:119)) +
    scale_color_manual(values = c("darkblue", "darkred")) +
    xlab("Distance to motif (bp)") +
    ylab("Cut-site Probability") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(
      legend.position = c(0.9, 0.9),
      panel.border = element_blank(),
      panel.grid   = element_blank(),
      legend.title = element_blank(),
      axis.line    = element_line(color = "black", size = 0.5),
      axis.ticks   = element_line(color = "black", size = 0.5),
      axis.text    = element_text(color = "black", size = 10)
    )
  
  ## add title
  if (! missing(mainTitle)) {
    p <- p + ggtitle(mainTitle)
  }
  
  return(p)
}

