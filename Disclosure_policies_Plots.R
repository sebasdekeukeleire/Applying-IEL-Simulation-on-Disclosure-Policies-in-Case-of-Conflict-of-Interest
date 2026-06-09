install.packages(c("ggplot2", "patchwork"))

library(ggplot2)
library(patchwork)

data_dir   <- getwd()
output_dir <- getwd()

# nash benchmarks
NASH_ALPHA_B  <- 0
NASH_TAU_B    <- 1
NASH_DELTA_B  <- 0
NASH_MSG_BIAS <- 4.5

# random play references
s_vals   <- rep(1:10, each = 10)
a_vals   <- rep(1:10, times = 10)
RAND_CLI <- 100 - mean((s_vals - a_vals)^2)   # 83.5
RAND_ADV <- 100 - mean((10 - 1:10)^2)         # 71.5

COL_NO  <- "#2166ac"
COL_MAN <- "#d6604d"

theme_thesis <- function() {
  theme_classic(base_size = 11) +
    theme(
      legend.position    = "bottom",
      legend.title       = element_blank(),
      legend.key.width   = unit(1.4, "cm"),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold"),
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4),
      plot.title         = element_text(face = "bold", size = 12),
      plot.subtitle      = element_text(size = 9.5, colour = "grey40",
                                        margin = margin(b = 8))
    )
}

# --- load data ---
no_disc  <- read.csv(file.path(data_dir, "IEL_v2_no_disclosure_Npairs5_J100_T50_REP100.csv"))
man_disc <- read.csv(file.path(data_dir, "IEL_v2_mandatory_disclosure_Npairs5_J100_T50_REP100.csv"))

no_disc$treatment  <- "No Disclosure"
man_disc$treatment <- "Mandatory Disclosure"

df <- rbind(no_disc, man_disc)
df$treatment <- factor(df$treatment, levels = c("No Disclosure", "Mandatory Disclosure"))
df$msg_bias  <- abs(df$message - df$state)

df_B <- df[df$type == 0, ]

# --- two-stage aggregation ---
two_stage_avg <- function(data, var) {
  s1 <- aggregate(data[[var]] ~ rep + period + treatment, data = data,
                  FUN = function(x) mean(x, na.rm = TRUE))
  names(s1)[4] <- "value"
  
  mn  <- aggregate(value ~ period + treatment, data = s1, FUN = mean)
  sd  <- aggregate(value ~ period + treatment, data = s1, FUN = sd)
  n   <- aggregate(value ~ period + treatment, data = s1, FUN = length)
  
  out       <- mn
  names(out)[3] <- "mean"
  out$se    <- sd$value / sqrt(n$value)
  out$lo    <- out$mean - 1.96 * out$se
  out$hi    <- out$mean + 1.96 * out$se
  return(out)
}

agg_alpha_B <- two_stage_avg(df_B, "alpha_B")
agg_tau_B   <- two_stage_avg(df_B[df_B$treatment == "Mandatory Disclosure", ], "tau_B")
agg_delta_B <- two_stage_avg(df_B[df_B$treatment == "Mandatory Disclosure", ], "delta_B")
agg_msgbias <- two_stage_avg(df_B, "msg_bias")
agg_cli_pay <- two_stage_avg(df,   "payoff_cli")
agg_adv_pay <- two_stage_avg(df_B, "payoff_adv")

# --- figure helpers ---
col_scale_c <- scale_colour_manual(values = c("No Disclosure" = COL_NO, "Mandatory Disclosure" = COL_MAN))
col_scale_f <- scale_fill_manual(  values = c("No Disclosure" = COL_NO, "Mandatory Disclosure" = COL_MAN))

make_conv_panel <- function(agg_data, y_label, nash_val, nash_text,
                            y_limits = NULL, show_legend = TRUE) {
  p <- ggplot(agg_data, aes(x = period, y = mean, colour = treatment, fill = treatment)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = nash_val, linetype = "dashed", colour = "grey40", linewidth = 0.55) +
    annotate("text", x = 48, y = nash_val, label = nash_text,
             hjust = 1, vjust = -0.5, size = 2.9, colour = "grey35") +
    col_scale_c + col_scale_f +
    labs(x = "Period", y = y_label) +
    theme_thesis()
  
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  if (!show_legend)        p <- p + theme(legend.position = "none")
  return(p)
}

make_payoff_panel <- function(agg_data, y_label, rand_ref, rand_label) {
  ggplot(agg_data, aes(x = period, y = mean, colour = treatment, fill = treatment)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = rand_ref, linetype = "dotted", colour = "grey50", linewidth = 0.55) +
    annotate("text", x = 48, y = rand_ref, label = rand_label,
             hjust = 1, vjust = -0.5, size = 2.9, colour = "grey40") +
    col_scale_c + col_scale_f +
    labs(x = "Period", y = y_label) +
    coord_cartesian(ylim = c(60, 100)) +
    theme_thesis()
}

# --- figure 1: strategy convergence ---
fig1 <- (p1a | p1b | p1c) +
  plot_annotation(
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

# --- figure 2: message inflation ---
fig2 <- ggplot(agg_msgbias, aes(x = period, y = mean, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = NASH_MSG_BIAS, linetype = "dashed", colour = "grey40", linewidth = 0.55) +
  annotate("text", x = 2, y = NASH_MSG_BIAS, label = "E[|m-s|] | Nash: 4.5",
           hjust = 0, vjust = -0.5, size = 2.9, colour = "grey35") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.45) +
  annotate("text", x = 2, y = 0, label = "Fully truthful: 0",
           hjust = 0, vjust = -0.5, size = 2.9, colour = "grey50") +
  col_scale_c + col_scale_f +
  labs(
    x = "Period", y = "|m - s|  (message inflation)"
  ) +
  coord_cartesian(ylim = c(0, 6)) +
  theme_thesis()

# --- figure 3: payoff trajectories ---
fig3 <- (p3a | p3b) +
  plot_annotation(
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# --- save ---
cat("\nSaving figures to:", output_dir, "\n")

ggsave(file.path(output_dir, "fig1_strategy_convergence.pdf"),  plot = fig1, width = 10, height = 4.5)
ggsave(file.path(output_dir, "fig2_message_inflation.pdf"),     plot = fig2, width = 7,  height = 4.5)
ggsave(file.path(output_dir, "fig3_payoff_trajectories.pdf"),   plot = fig3, width = 8,  height = 4)
ggsave(file.path(output_dir, "fig1_strategy_convergence.png"),  plot = fig1, width = 10, height = 4.5, dpi = 300)
ggsave(file.path(output_dir, "fig2_message_inflation.png"),     plot = fig2, width = 7,  height = 4.5, dpi = 300)
ggsave(file.path(output_dir, "fig3_payoff_trajectories.png"),   plot = fig3, width = 8,  height = 4,   dpi = 300)

cat("Saved:\n  fig1_strategy_convergence.pdf / .png\n  fig2_message_inflation.pdf / .png\n  fig3_payoff_trajectories.pdf / .png\nDone!\n")