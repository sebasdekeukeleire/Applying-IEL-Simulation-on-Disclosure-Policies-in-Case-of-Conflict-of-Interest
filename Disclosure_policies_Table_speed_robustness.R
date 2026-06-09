# Summary statistics table (periods 41–50) with standard errors 

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
  
  out        <- mn
  names(out)[2] <- "mean"
  out$se     <- sd$value / sqrt(n$value)
  out$var    <- var
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

# --- helper to extract mean (SE) ---
fmt <- function(tbl, tr) {
  row <- tbl[tbl$treatment == tr, ]
  if (nrow(row) == 0) return(NA)
  sprintf("%.3f (%.3f)", row$mean, row$se)
}

# --- assemble table ---
table_rows <- data.frame(
  Variable = c("alpha_B", "|m-s|", "tau_B", "delta_B",
               "Client payoff", "Biased adv. payoff"),
  No_Disclosure = c(
    fmt(alpha_B_tbl,  "No Disclosure"),
    fmt(msg_bias_tbl, "No Disclosure"),
    NA,
    NA,
    fmt(cli_pay_tbl,  "No Disclosure"),
    fmt(adv_pay_tbl,  "No Disclosure")
  ),
  Mandatory = c(
    fmt(alpha_B_tbl,  "Mandatory Disclosure"),
    fmt(msg_bias_tbl, "Mandatory Disclosure"),
    fmt(tau_B_tbl,    "Mandatory Disclosure"),
    fmt(delta_B_tbl,  "Mandatory Disclosure"),
    fmt(cli_pay_tbl,  "Mandatory Disclosure"),
    fmt(adv_pay_tbl,  "Mandatory Disclosure")
  ),
  Nash = c("0", "4.5", "1", "0", "—", "—")
)

print(table_rows)

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

cat("\nConvergence speed — first period (rep-averaged) alpha_B < 0.1:\n")
cat("  No Disclosure:        period", first_below("No Disclosure"),        "\n")
cat("  Mandatory Disclosure: period", first_below("Mandatory Disclosure"), "\n")

# --- Welch t-tests ---
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



# ROBUSTNESS TABLE 
# Late-period means for key variables across all parameter combinations.
# Late is always the final 10 periods of each run.

# ---  define all robustness runs ---

data_dir <- getwd()

robustness_runs <- list(
  list(
    label   = "Baseline (T=50, J=100, P_ex=0.033, σ=0.10)",
    no_file = "IEL_v2_no_disclosure_Npairs5_J100_T50_REP100.csv",
    ma_file = "IEL_v2_mandatory_disclosure_Npairs5_J100_T50_REP100.csv"
  ),
  # Time horizon
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
  # Experimentation rate
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
  # Perturbation size
  list(
    label   = "σ=0.05",
    no_file = "IEL_v3_no_disclosure_Npairs5_J100_T50_REP100_Pex0033_sig005.csv",
    ma_file = "IEL_v3_mandatory_disclosure_Npairs5_J100_T50_REP100_Pex0033_sig005.csv"
  ),
  list(
    label   = "σ=0.20",
    no_file = "IEL_v3_no_disclosure_Npairs5_J100_T50_REP100_Pex0033_sig020.csv",
    ma_file = "IEL_v3_mandatory_disclosure_Npairs5_J100_T50_REP100_Pex0033_sig020.csv"
  ),
  # Portfolio size
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

# --- function to compute late-period means for one treatment file ---
summarise_run <- function(filepath) {
  d            <- read.csv(filepath)
  d$msg_bias   <- abs(d$message - d$state)
  max_period   <- max(d$period)
  late         <- d[d$period > max_period - 10, ]
  late_B       <- late[late$type == 0, ]
  
  # Two-stage mean: average within rep first, then across reps
  s1_mean <- function(dat, var) {
    s1 <- aggregate(dat[[var]] ~ rep, data = dat,
                    FUN = function(x) mean(x, na.rm = TRUE))
    mean(s1[[2]], na.rm = TRUE)
  }
  
  list(
    alpha_B    = s1_mean(late_B, "alpha_B"),
    msg_bias   = s1_mean(late_B, "msg_bias"),
    cli_payoff = s1_mean(late,   "payoff_cli"),
    adv_payoff = s1_mean(late_B, "payoff_adv")
  )
}

# --- build the table ---
rob_rows <- lapply(robustness_runs, function(run) {
  no <- summarise_run(file.path(data_dir, run$no_file))
  ma <- summarise_run(file.path(data_dir, run$ma_file))
  data.frame(
    Check           = run$label,
    aB_no           = round(no$alpha_B,    3),
    aB_ma           = round(ma$alpha_B,    3),
    msg_no          = round(no$msg_bias,   3),
    msg_ma          = round(ma$msg_bias,   3),
    cli_no          = round(no$cli_payoff, 3),
    cli_ma          = round(ma$cli_payoff, 3),
    adv_no          = round(no$adv_payoff, 3),
    adv_ma          = round(ma$adv_payoff, 3),
    stringsAsFactors = FALSE
  )
})

rob_table <- do.call(rbind, rob_rows)
names(rob_table) <- c(
  "Check",
  "αB (No Disc)", "αB (Mand)",
  "|m-s| (No Disc)", "|m-s| (Mand)",
  "Cli. payoff (No)", "Cli. payoff (Mand)",
  "Adv. payoff (No)", "Adv. payoff (Mand)"
)

print(rob_table, row.names = FALSE)