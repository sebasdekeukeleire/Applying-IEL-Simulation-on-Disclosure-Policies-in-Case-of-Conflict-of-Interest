# Summary statistics table (periods 41-50) with standard errors

library(dplyr)

no_disc  <- read.csv("IEL_v2_no_disclosure_Npairs5_J100_T50_REP100.csv")
man_disc <- read.csv("IEL_v2_mandatory_disclosure_Npairs5_J100_T50_REP100.csv")

no_disc$treatment  <- "No Disclosure"
man_disc$treatment <- "Mandatory Disclosure"

df <- rbind(no_disc, man_disc)
df$msg_bias <- abs(df$message - df$state)
df_B        <- df[df$type == 0, ]

late   <- df[df$period >= 41, ]
late_B <- late[late$type == 0, ]

# --- two-stage mean and SE ---
stage1_mean_se <- function(data, var) {
  s1 <- aggregate(data[[var]] ~ rep + treatment, data = data,
                  FUN = function(x) mean(x, na.rm = TRUE))
  names(s1)[3] <- "value"
  
  mn  <- aggregate(value ~ treatment, data = s1, FUN = mean)
  sd  <- aggregate(value ~ treatment, data = s1, FUN = sd)
  n   <- aggregate(value ~ treatment, data = s1, FUN = length)
  
  out           <- mn
  names(out)[2] <- "mean"
  out$se        <- sd$value / sqrt(n$value)
  out$var       <- var
  return(out)
}

alpha_B_tbl  <- stage1_mean_se(late_B, "alpha_B")
msg_bias_tbl <- stage1_mean_se(late_B, "msg_bias")
tau_B_tbl    <- stage1_mean_se(
  late_B[late_B$treatment == "Mandatory Disclosure", ], "tau_B")
delta_B_tbl  <- stage1_mean_se(
  late_B[late_B$treatment == "Mandatory Disclosure", ], "delta_B")
cli_pay_tbl  <- stage1_mean_se(late,   "payoff_cli")
adv_pay_tbl  <- stage1_mean_se(late_B, "payoff_adv")
tau_no_tbl   <- stage1_mean_se(late_B, "tau_no")
delta_no_tbl <- stage1_mean_se(late_B, "delta_no")

# --- helper to extract mean (SE) ---
fmt <- function(tbl, tr) {
  row <- tbl[tbl$treatment == tr, ]
  if (nrow(row) == 0) return(NA)
  sprintf("%.3f (%.3f)", row$mean, row$se)
}

# --- assemble table ---
table_rows <- data.frame(
  Variable = c("alpha_B", "|m-s|", "tau_no", "delta_no", "tau_B", "delta_B",
               "Client payoff", "Biased adv. payoff"),
  No_Disclosure = c(
    fmt(alpha_B_tbl,  "No Disclosure"),
    fmt(msg_bias_tbl, "No Disclosure"),
    fmt(tau_no_tbl,   "No Disclosure"),
    fmt(delta_no_tbl, "No Disclosure"),
    NA, NA,
    fmt(cli_pay_tbl,  "No Disclosure"),
    fmt(adv_pay_tbl,  "No Disclosure")
  ),
  Mandatory = c(
    fmt(alpha_B_tbl,  "Mandatory Disclosure"),
    fmt(msg_bias_tbl, "Mandatory Disclosure"),
    fmt(tau_no_tbl,   "Mandatory Disclosure"),
    fmt(delta_no_tbl, "Mandatory Disclosure"),
    fmt(tau_B_tbl,    "Mandatory Disclosure"),
    fmt(delta_B_tbl,  "Mandatory Disclosure"),
    fmt(cli_pay_tbl,  "Mandatory Disclosure"),
    fmt(adv_pay_tbl,  "Mandatory Disclosure")
  ),
  Nash = c("0", "4.5", "1", "0", "1", "0", "-", "-")
)

print(table_rows)

# --- total welfare ---
late$total_welfare <- late$payoff_cli + late$payoff_adv
tw_tbl <- stage1_mean_se(late, "total_welfare")
print(tw_tbl)

late_A <- late[late$type == 1, ]
late_A$total_welfare <- late_A$payoff_cli + late_A$payoff_adv
tw_A_tbl <- stage1_mean_se(late_A, "total_welfare")
print(tw_A_tbl)

late_B$total_welfare <- late_B$payoff_cli + late_B$payoff_adv
tw_B_tbl <- stage1_mean_se(late_B, "total_welfare")
print(tw_B_tbl)

# --- convergence speed ---
conv_stage1 <- aggregate(alpha_B ~ rep + period + treatment,
                         data = df_B, FUN = mean, na.rm = TRUE)
conv_stage2 <- aggregate(alpha_B ~ period + treatment,
                         data = conv_stage1, FUN = mean)

first_below <- function(tr) {
  sub <- conv_stage2[conv_stage2$treatment == tr, ]
  sub <- sub[order(sub$period), ]
  idx <- which(sub$alpha_B < 0.1)
  if (length(idx) == 0) return(NA)
  sub$period[idx[1]]
}

cat("\nConvergence speed - first period (rep-averaged) alpha_B < 0.1:\n")
cat("  No Disclosure:        period", first_below("No Disclosure"),        "\n")
cat("  Mandatory Disclosure: period", first_below("Mandatory Disclosure"), "\n")

# --- Welch t-tests: payoffs ---
late_cli_rep <- aggregate(payoff_cli ~ rep + treatment, data = late,   FUN = mean)
late_adv_rep <- aggregate(payoff_adv ~ rep + treatment, data = late_B, FUN = mean)

no_cli_vals  <- late_cli_rep$payoff_cli[late_cli_rep$treatment == "No Disclosure"]
man_cli_vals <- late_cli_rep$payoff_cli[late_cli_rep$treatment == "Mandatory Disclosure"]
no_adv_vals  <- late_adv_rep$payoff_adv[late_adv_rep$treatment == "No Disclosure"]
man_adv_vals <- late_adv_rep$payoff_adv[late_adv_rep$treatment == "Mandatory Disclosure"]

tt_cli <- t.test(man_cli_vals, no_cli_vals, paired = FALSE, var.equal = FALSE)
tt_adv <- t.test(man_adv_vals, no_adv_vals, paired = FALSE, var.equal = FALSE)

cat("\nWelch t-test: CLIENT payoff (periods 41-50)\n")
cat("  Mean difference:", round(mean(man_cli_vals) - mean(no_cli_vals), 3), "\n")
cat("  t =", round(tt_cli$statistic, 3), " df =", round(tt_cli$parameter, 1),
    " p =", round(tt_cli$p.value, 4), "\n")

cat("\nWelch t-test: BIASED ADVISOR payoff (periods 41-50)\n")
cat("  Mean difference:", round(mean(man_adv_vals) - mean(no_adv_vals), 3), "\n")
cat("  t =", round(tt_adv$statistic, 3), " df =", round(tt_adv$parameter, 1),
    " p =", round(tt_adv$p.value, 4), "\n")

# --- Welch t-test: total welfare ---
late_tw_rep <- aggregate(total_welfare ~ rep + treatment, data = late, FUN = mean)

no_tw_vals  <- late_tw_rep$total_welfare[late_tw_rep$treatment == "No Disclosure"]
man_tw_vals <- late_tw_rep$total_welfare[late_tw_rep$treatment == "Mandatory Disclosure"]

tt_tw <- t.test(man_tw_vals, no_tw_vals, paired = FALSE, var.equal = FALSE)

cat("\nWelch t-test: TOTAL WELFARE (periods 41-50)\n")
cat("  Mean difference:", round(mean(man_tw_vals) - mean(no_tw_vals), 3), "\n")
cat("  t =", round(tt_tw$statistic, 3), " df =", round(tt_tw$parameter, 1),
    " p =", round(tt_tw$p.value, 4), "\n")

# --- Welch t-test: total welfare Type B only ---
late_tw_B_rep <- aggregate(total_welfare ~ rep + treatment, data = late_B, FUN = mean)

no_tw_B_vals  <- late_tw_B_rep$total_welfare[late_tw_B_rep$treatment == "No Disclosure"]
man_tw_B_vals <- late_tw_B_rep$total_welfare[late_tw_B_rep$treatment == "Mandatory Disclosure"]

tt_tw_B <- t.test(man_tw_B_vals, no_tw_B_vals, paired = FALSE, var.equal = FALSE)

cat("\nWelch t-test: TOTAL WELFARE Type B only (periods 41-50)\n")
cat("  Mean difference:", round(mean(man_tw_B_vals) - mean(no_tw_B_vals), 3), "\n")
cat("  t =", round(tt_tw_B$statistic, 3), " df =", round(tt_tw_B$parameter, 1),
    " p =", round(tt_tw_B$p.value, 4), "\n")


# --- robustness table ---

data_dir <- getwd()

# Safe rep-level mean: extracts column directly, no formula interface
s1_mean_safe <- function(dat, var) {
  col <- dat[[var]]
  if (is.null(col) || all(is.na(col))) return(NA_real_)
  rep_ids <- sort(unique(dat$rep))
  rep_avgs <- sapply(rep_ids, function(r) {
    mean(col[dat$rep == r], na.rm = TRUE)
  })
  mean(rep_avgs, na.rm = TRUE)
}

summarise_run <- function(filepath) {
  d          <- read.csv(filepath)
  d$msg_bias <- abs(d$message - d$state)
  max_period <- max(d$period)
  late_loc   <- d[d$period > max_period - 10, ]
  late_B_loc <- late_loc[late_loc$type == 0, ]
  list(
    alpha_B    = s1_mean_safe(late_B_loc, "alpha_B"),
    msg_bias   = s1_mean_safe(late_B_loc, "msg_bias"),
    tau_no     = s1_mean_safe(late_B_loc, "tau_no"),
    delta_no   = s1_mean_safe(late_B_loc, "delta_no"),
    tau_B      = s1_mean_safe(late_B_loc, "tau_B"),
    delta_B    = s1_mean_safe(late_B_loc, "delta_B"),
    cli_payoff = s1_mean_safe(late_loc,   "payoff_cli"),
    adv_payoff = s1_mean_safe(late_B_loc, "payoff_adv")
  )
}

robustness_runs <- list(
  list(
    label   = "Baseline (T=50, J=100, P_ex=0.033, sigma=0.10)",
    no_file = "IEL_v2_no_disclosure_Npairs5_J100_T50_REP100.csv",
    ma_file = "IEL_v2_mandatory_disclosure_Npairs5_J100_T50_REP100.csv"
  ),
  list(
    label   = "T=100",
    no_file = "IEL_v3_no_disclosure_Npairs5_J100_T100_REP100_Pex0033_sig010.csv",
    ma_file = "IEL_v3_mandatory_disclosure_Npairs5_J100_T100_REP100_Pex0033_sig010.csv"
  ),
  list(
    label   = "T=200",
    no_file = "IEL_v3_no_disclosure_Npairs5_J100_T200_REP100_Pex0033_sig010.csv",
    ma_file = "IEL_v3_mandatory_disclosure_Npairs5_J100_T200_REP100_Pex0033_sig010.csv"
  ),
  list(
    label   = "P_ex=0.01",
    no_file = "IEL_v3_no_disclosure_Npairs5_J100_T50_REP100_Pex0010_sig010.csv",
    ma_file = "IEL_v3_mandatory_disclosure_Npairs5_J100_T50_REP100_Pex0010_sig010.csv"
  ),
  list(
    label   = "P_ex=0.10",
    no_file = "IEL_v3_no_disclosure_Npairs5_J100_T50_REP100_Pex0100_sig010.csv",
    ma_file = "IEL_v3_mandatory_disclosure_Npairs5_J100_T50_REP100_Pex0100_sig010.csv"
  ),
  list(
    label   = "sigma=0.05",
    no_file = "IEL_v3_no_disclosure_Npairs5_J100_T50_REP100_Pex0033_sig005.csv",
    ma_file = "IEL_v3_mandatory_disclosure_Npairs5_J100_T50_REP100_Pex0033_sig005.csv"
  ),
  list(
    label   = "sigma=0.20",
    no_file = "IEL_v3_no_disclosure_Npairs5_J100_T50_REP100_Pex0033_sig020.csv",
    ma_file = "IEL_v3_mandatory_disclosure_Npairs5_J100_T50_REP100_Pex0033_sig020.csv"
  ),
  list(
    label   = "J=50",
    no_file = "IEL_v3_no_disclosure_Npairs5_J50_T50_REP100_Pex0033_sig010.csv",
    ma_file = "IEL_v3_mandatory_disclosure_Npairs5_J50_T50_REP100_Pex0033_sig010.csv"
  ),
  list(
    label   = "J=200",
    no_file = "IEL_v3_no_disclosure_Npairs5_J200_T50_REP100_Pex0033_sig010.csv",
    ma_file = "IEL_v3_mandatory_disclosure_Npairs5_J200_T50_REP100_Pex0033_sig010.csv"
  )
)

rob_rows <- lapply(robustness_runs, function(run) {
  no <- summarise_run(file.path(data_dir, run$no_file))
  ma <- summarise_run(file.path(data_dir, run$ma_file))
  data.frame(
    Check  = run$label,
    aB_no  = round(no$alpha_B,    3),
    aB_ma  = round(ma$alpha_B,    3),
    msg_no = round(no$msg_bias,   3),
    msg_ma = round(ma$msg_bias,   3),
    tno_no = round(no$tau_no,     3),
    tno_ma = round(ma$tau_no,     3),
    dno_no = round(no$delta_no,   3),
    dno_ma = round(ma$delta_no,   3),
    tB_ma  = round(ma$tau_B,      3),
    dB_ma  = round(ma$delta_B,    3),
    cli_no = round(no$cli_payoff, 3),
    cli_ma = round(ma$cli_payoff, 3),
    adv_no = round(no$adv_payoff, 3),
    adv_ma = round(ma$adv_payoff, 3),
    stringsAsFactors = FALSE
  )
})

rob_table <- do.call(rbind, rob_rows)
names(rob_table) <- c(
  "Check",
  "aB (No)", "aB (Mand)",
  "|m-s| (No)", "|m-s| (Mand)",
  "tno (No)", "tno (Mand)",
  "dno (No)", "dno (Mand)",
  "tB (Mand)", "dB (Mand)",
  "Cli (No)", "Cli (Mand)",
  "Adv (No)", "Adv (Mand)"
)

print(rob_table, row.names = FALSE)


# --- WELCH T-TESTS: advisor behaviour across treatments --- 


get_rep_means <- function(data, var, tr) {
  sub <- data[data$treatment == tr, ]
  col <- sub[[var]]
  if (is.null(col)) stop(paste("Column not found:", var))
  rep_ids <- sort(unique(sub$rep))
  sapply(rep_ids, function(r) mean(col[sub$rep == r], na.rm = TRUE))
}

alpha_rep_no  <- get_rep_means(late_B, "alpha_B",  "No Disclosure")
alpha_rep_man <- get_rep_means(late_B, "alpha_B",  "Mandatory Disclosure")
msg_rep_no    <- get_rep_means(late_B, "msg_bias", "No Disclosure")
msg_rep_man   <- get_rep_means(late_B, "msg_bias", "Mandatory Disclosure")

tt_alpha <- t.test(alpha_rep_man, alpha_rep_no, paired = FALSE, var.equal = FALSE)
tt_msg   <- t.test(msg_rep_man,   msg_rep_no,   paired = FALSE, var.equal = FALSE)

cat("\nWelch t-test: alpha_B across treatments (periods 41-50)\n")
cat("  Mean difference:", round(mean(alpha_rep_man) - mean(alpha_rep_no), 4), "\n")
cat("  t =", round(tt_alpha$statistic, 3), " df =", round(tt_alpha$parameter, 1),
    " p =", round(tt_alpha$p.value, 4), "\n")

cat("\nWelch t-test: |m-s| across treatments (periods 41-50)\n")
cat("  Mean difference:", round(mean(msg_rep_man) - mean(msg_rep_no), 4), "\n")
cat("  t =", round(tt_msg$statistic, 3), " df =", round(tt_msg$parameter, 1),
    " p =", round(tt_msg$p.value, 4), "\n")