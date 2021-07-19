library(ggplot2)
ggmiami2 <-
  function (data,
            split_by,
            split_at,
            chr = "chr",
            pos = "pos",
            p = "p",
            chr_colors = c("black", "grey"),
            upper_chr_colors = NULL,
            lower_chr_colors = NULL,
            upper_ylab = "-log10(p)",
            lower_ylab = "-log10(p)",
            genome_line = 5e-08,
            genome_line_color = "red",
            suggestive_line_upper = 1e-05,
            suggestive_line_bottom = 1e-05,
            suggestive_line_color = "blue",
            hits_label_col = NULL,
            hits_label = NULL,
            top_n_hits = 5,
            upper_labels_df = NULL,
            lower_labels_df = NULL,
            upper_highlight = NULL,
            upper_highlight_col = NULL,
            upper_highlight_color = "green",
            lower_highlight = NULL,
            lower_highlight_col = NULL,
            lower_highlight_color = "green"
  )
  {
    plot_data <- prep_miami_data(
      data = data,
      split_by = split_by,
      split_at = split_at,
      chr = chr,
      pos = pos,
      p = p
    )
    if (all(!is.null(chr_colors), any(!is.null(upper_chr_colors),
                                      !is.null(lower_chr_colors)))) {
      stop(
        "You have specified both chr_colors and upper_chr_colors and/or\n         lower_chr_colors. This package does not know how to use both\n         information simultaneously. Please only use one method for coloring:\n         either chr_colors, for making upper and lower plot have the same\n         colors, or upper_chr_colors + lower_chr_colors for specifying different\n         colors for upper and lower plot."
      )
    }
    if (upper_ylab == "-log10(p)") {
      upper_ylab <- expression("-log"[10] * "(p)")
    }
    else {
      upper_ylab <- bquote(atop(.(upper_ylab), "-log"[10] *
                                  "(p)"))
    }
    if (lower_ylab == "-log10(p)") {
      lower_ylab <- expression("-log"[10] * "(p)")
    }
    else {
      lower_ylab <- bquote(atop(.(lower_ylab), "-log"[10] *
                                  "(p)"))
    }
    upper_plot <-
      ggplot2::ggplot() + ggplot2::geom_point(
        data = plot_data$upper,
        aes(
          x = .data$rel_pos,
          y = .data$logged_p,
          color = as.factor(.data$chr)
        ),
        size = 0.25
      ) + ggplot2::scale_x_continuous(
        labels = plot_data$axis$chr,
        breaks = plot_data$axis$chr_center,
        expand = ggplot2::expansion(mult = 0.01),
        guide = ggplot2::guide_axis(check.overlap = TRUE)
      ) +
      ggplot2::scale_y_continuous(limits = c(0, plot_data$maxp),
                                  expand = ggplot2::expansion(mult = c(0.02, 0))) +
      ggplot2::labs(x = "", y = upper_ylab) + ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(b = 0, l = 10)
      )
    lower_plot <-
      ggplot2::ggplot() + ggplot2::geom_point(
        data = plot_data$lower,
        aes(
          x = .data$rel_pos,
          y = .data$logged_p,
          color = as.factor(.data$chr)
        ),
        size = 0.25
      ) + ggplot2::scale_x_continuous(
        breaks = plot_data$axis$chr_center,
        position = "top",
        expand = ggplot2::expansion(mult = 0.01)
      ) +
      ggplot2::scale_y_reverse(limits = c(plot_data$maxp, 0),
                               expand = ggplot2::expansion(mult = c(0, 0.02))) +
      ggplot2::labs(x = "", y = lower_ylab) + ggplot2::theme_classic() +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 0,
                                      l = 10)
      )
    if (all(!is.null(chr_colors),
            is.null(upper_chr_colors),
            is.null(lower_chr_colors))) {
      if (length(chr_colors) == 2) {
        chr_colors <- rep(chr_colors, length.out = nrow(plot_data$axis))
      }
      else if (length(chr_colors) == nrow(plot_data$axis)) {
        chr_colors <- chr_colors
      }
      else {
        stop(
          "The number of colors specified in {chr_colors} does not match the\n         number of chromosomes to be displayed."
        )
      }
      upper_plot <-
        upper_plot + ggplot2::scale_color_manual(values = chr_colors)
      lower_plot <-
        lower_plot + ggplot2::scale_color_manual(values = chr_colors)
    }
    else if (all(is.null(chr_colors),
                 !is.null(upper_chr_colors),
                 !is.null(lower_chr_colors))) {
      if (length(upper_chr_colors) == 2) {
        upper_chr_colors <-
          rep(upper_chr_colors, length.out = nrow(plot_data$axis))
      }
      else if (length(upper_chr_colors) == nrow(plot_data$axis)) {
        upper_chr_colors <- upper_chr_colors
      }
      else {
        stop(
          "The number of colors specified in {upper_chr_colors} does not match\n           the number of chromosomes to be displayed."
        )
      }
      if (length(lower_chr_colors) == 2) {
        lower_chr_colors <-
          rep(lower_chr_colors, length.out = nrow(plot_data$axis))
      }
      else if (length(lower_chr_colors) == nrow(plot_data$axis)) {
        lower_chr_colors <- lower_chr_colors
      }
      else {
        stop(
          "The number of colors specified in {lower_chr_colors} does not match\n           the number of chromosomes to be displayed."
        )
      }
      upper_plot <-
        upper_plot + ggplot2::scale_color_manual(values = upper_chr_colors)
      lower_plot <-
        lower_plot + ggplot2::scale_color_manual(values = lower_chr_colors)
    }
    else if (all(is.null(chr_colors), any(is.null(upper_chr_colors),
                                          is.null(lower_chr_colors)))) {
      stop(
        "It looks like you've specified one of upper or lower chr colors\n         without specifying the other. This package needs both colors, unless\n         you want the upper and lower plot to have the same colors, which is\n         done using {chr_colors}."
      )
    }
    if (!is.null(suggestive_line_upper)) {
      upper_plot <-
        upper_plot + ggplot2::geom_hline(
          yintercept = -log10(suggestive_line_upper),
          color = suggestive_line_color,
          linetype = "solid",
          size = 0.4
        )
    }
    if (!is.null(suggestive_line_bottom)) {
      lower_plot <-
        lower_plot + ggplot2::geom_hline(
          yintercept = -log10(suggestive_line_bottom),
          color = suggestive_line_color,
          linetype = "solid",
          size = 0.4
        ) + annotate("segment", x = 10, xend = 120212500, y = 16, yend = 16,  colour = "red",size = 0.8) + 
        annotate('text', x = 252276783, y = 16, label = 'FDR < 0.05', color = 'black',size = 2.5)
    }
    if (!is.null(genome_line)) {
      upper_plot <-
        upper_plot + ggplot2::geom_hline(
          yintercept = -log10(genome_line),
          color = genome_line_color,
          linetype = "dashed",
          size = 0.4
        )
      lower_plot <-
        lower_plot + ggplot2::geom_hline(
          yintercept = -log10(genome_line),
          color = genome_line_color,
          linetype = "dashed",
          size = 0.4
        )
    }
    if (all(
      !is.null(hits_label_col), 
      any(!is.null(upper_labels_df),
          !is.null(lower_labels_df)))
    ) {
      stop(
        "You have specified both hits_label_col and a *_labels_df. This\n         package does not know how to use both information simultaneously.\n         Please only use one method for labelling: either hits_label_col (with\n         or without hits_label), or *_labels_df."
      )
    }
    if (all(!is.null(hits_label_col),
            is.null(upper_labels_df),
            is.null(lower_labels_df))) {
      upper_labels_df <- make_miami_labels(
        data = plot_data$upper,
        hits_label_col = hits_label_col,
        hits_label = hits_label,
        top_n_hits = top_n_hits
      )
      lower_labels_df <- make_miami_labels(
        data = plot_data$lower,
        hits_label_col = hits_label_col,
        hits_label = hits_label,
        top_n_hits = top_n_hits
      )
      upper_plot <-
        upper_plot + ggrepel::geom_label_repel(
          data = upper_labels_df,
          aes(
            x = .data$rel_pos,
            y = .data$logged_p,
            label = .data$label
          ),
          size = 2,
          segment.size = 0.2,
          point.padding = 0.3,
          ylim = c(plot_data$maxp /
                     2, NA),
          min.segment.length = 0,
          force = 2,
          max.overlaps  = Inf,
          box.padding = 0.5
        )
      lower_plot <-
        lower_plot + ggrepel::geom_label_repel(
          data = lower_labels_df,
          aes(
            x = .data$rel_pos,
            y = .data$logged_p,
            label = .data$label
          ),
          size = 2,
          max.overlaps  = Inf,
          segment.size = 0.2,
          point.padding = 0.3,
          ylim = c(NA, -(plot_data$maxp /
                           2)),
          min.segment.length = 0,
          force = 2,
          box.padding = 0.5
        )
    }
    if (all(is.null(hits_label_col), !is.null(upper_labels_df))) {
      checkmate::assertNames(colnames(upper_labels_df),
                             identical.to = c("rel_pos",
                                              "logged_p", "label"))
      upper_plot <-
        upper_plot + ggrepel::geom_text_repel(
          data = upper_labels_df,
          aes(
            x = .data$rel_pos,
            y = .data$logged_p,
            label = .data$label
          ),
          size = 2,
          segment.size = 0.2,
          max.overlaps  = Inf,
          point.padding = 0.3,
          ylim = c(plot_data$maxp /
                     2, NA),
          min.segment.length = 0,
          force = 2,
          box.padding = 0.5
        )
    }
    if (all(is.null(hits_label_col), !is.null(lower_labels_df))) {
      checkmate::assertNames(colnames(lower_labels_df),
                             identical.to = c("rel_pos",
                                              "logged_p", "label"))
      lower_plot <-
        lower_plot + ggrepel::geom_text_repel(
          data = lower_labels_df,
          aes(
            x = .data$rel_pos,
            y = .data$logged_p,
            label = .data$label
          ),
          size = 2,
          max.overlaps  = Inf,
          segment.size = 0.2,
          point.padding = 0.3,
          ylim = c(NA, -(plot_data$maxp /
                           2)),
          min.segment.length = 0,
          force = 2,
          box.padding = 0.5
        )
    }
    if (all(!is.null(upper_highlight), !is.null(upper_highlight_col))) {
      upper_highlight_df <- highlight_miami(
        data = plot_data$upper,
        highlight = upper_highlight,
        highlight_col = upper_highlight_col,
        highlight_color = upper_highlight_color
      )
      upper_plot <-
        upper_plot + ggplot2::geom_point(
          data = upper_highlight_df,
          aes(x = .data$rel_pos, y = .data$logged_p),
          color = upper_highlight_df$color,
          size = 0.25
        )
    }
    if (all(!is.null(lower_highlight), !is.null(lower_highlight_col))) {
      lower_highlight_df <- highlight_miami(
        data = plot_data$lower,
        highlight = lower_highlight,
        highlight_col = lower_highlight_col,
        highlight_color = lower_highlight_color
      )
      lower_plot <-
        lower_plot + ggplot2::geom_point(
          data = lower_highlight_df,
          aes(x = .data$rel_pos, y = .data$logged_p),
          color = lower_highlight_df$color,
          size = 0.25
        )
    }
    gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
  }