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
#' plot_weights(n = 20, r1 = 10, weight_fun = weights_jsd)
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
  } else if (sum(lapply(weight_params, length) > 1) == 2) {
    index <- which(lapply(weight_params, length) > 1)
    param_grid <- expand.grid(weight_params)
    plotdf <- data.frame()

    for (i in 1:nrow(param_grid)) {
      weight_mat <- do.call(weight_fun, args = c(as.list(param_grid[i, ]),
        design = design, n = n))
      weight_mat_temp <- cbind(
        rep(param_grid[i, index[1]], n + 1),
        rep(param_grid[i, index[2]], n + 1),
        0:n,
        weight_mat[r1 + 1, ]
      )
      plotdf <- rbind(plotdf, weight_mat_temp)
    }
    names(plotdf) <- c("param1", "param2", "r", "weight")
    plotdf$param1 <- factor(plotdf$param1)
    plotdf$param2 <- factor(plotdf$param2)
    ggplot2::ggplot(plotdf, ggplot2::aes(x = r, y = weight)) +
      ggplot2::geom_line(ggplot2::aes(group = param1, colour = param1)) +
      ggplot2::facet_grid(
        cols = ggplot2::vars(param2),
        labeller = ggplot2::as_labeller(function(x)
          paste(names(index)[2], "=", x))
      ) +
      ggplot2::xlim(0, n) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_bw() +
      ggplot2::scale_color_discrete(name = names(index)[1])
  } else {
    stop("at most two weight parameter can have length > 1")
  }
}




