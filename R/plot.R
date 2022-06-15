#' @include class.R
NULL

#' Plot Weight Functions
#'
#' @template n
#' @param r1 Number of responses in one basket
#' @template weights
#'
#' @return A plot.
#' @export
#'
#' @examples
#' plot_weights(n = 20, r1 = 10, weight_fun = weights_prob)
plot_weights <- function(n, r1, weight_fun, weight_params = list()) {
  design <- setupOneStageBasket(k = 2, theta0 = 0.2)

  if (!any(lapply(weight_params, length) > 1)) {
    weight_mat <- do.call(weight_fun, args = c(weight_params, design = design,
      n = n))
    plotdf <- data.frame(r = 0:n, weight = weight_mat[r1 + 1, ])

    ggplot2::ggplot(plotdf, ggplot2::aes(x = r, y = weight)) +
      ggplot2::geom_line() +
      ggplot2::xlim(0, n) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_bw()
  } else if (sum(lapply(weight_params, length) > 1) == 1) {
    index <- which(lapply(weight_params, length) > 1)
    l <- length(weight_params[[index]])
    plotdf <- data.frame()

    for (i in 1:l) {
      weight_loop <- c(weight_params[-index], list(weight_params[[index]][i]))
      names(weight_loop)[length(weight_loop)] <- names(index)
      weight_mat <- do.call(weight_fun, args = c(weight_loop, design = design,
        n = n))
      weight_mat_temp <- cbind(rep(weight_params[[index]][i], n + 1),
        0:n, weight_mat[r1 + 1, ])
      plotdf <- rbind(plotdf, weight_mat_temp)
    }

    names(plotdf) <- c("param", "r", "weight")
    plotdf$param <- factor(plotdf$param)
    ggplot2::ggplot(plotdf, ggplot2::aes(x = r, y = weight)) +
      ggplot2::geom_line(ggplot2::aes(group = param, colour = param)) +
      ggplot2::xlim(0, n) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_bw() +
      ggplot2::scale_color_discrete(name = names(index))
  } else {
    stop("only one weight parameter can have length > 1")
  }
}




