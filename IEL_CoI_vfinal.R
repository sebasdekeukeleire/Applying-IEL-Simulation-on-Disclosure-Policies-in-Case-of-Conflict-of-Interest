# IEL Simulation — Conflict of Interest Cheap Talk
#
# WHAT THIS CODE DOES:
#
# We simulate a game between an ADVISOR and a CLIENT, repeated many times.
#
# Each round:
#   1. A true "state" s is drawn randomly from 1 to 10
#   2. The advisor is randomly assigned a type:
#        - Type A (unbiased): wants to give correct advice
#        - Type B (biased):   always wants the client to pick action 10
#   3. The advisor sends a message m (1-10) to the client
#   4. The client sees the message (and in mandatory disclosure, the type)
#      and picks an action a (1-10)
#   5. Payoffs are calculated: client wants a = s, biased advisor wants a = 10
#
# HOW STRATEGIES EVOLVE (Individual Evolutionary Learning):
#   Each agent carries J candidate strategies (numbers between 0 and 1).
#   Every round:
#     - Each strategy is randomly nudged slightly (experimentation)
#     - We calculate: "what would each strategy have earned this round?"
#     - Better strategies replace worse ones (tournament)
#     - Next round's strategy is picked by lottery (better = more likely)
#
# TWO TREATMENTS:
#   - No disclosure:        client never knows advisor type
#   - Mandatory disclosure: client always knows advisor type
#
#   N_pairs controls the within-repetition sample size. In each of the REP
#   independent repetitions, N_pairs advisor-client pairs play the game
#   simultaneously. All pairs are independent draws of the same game, they
#   share no information and are used purely to generate more observations
#   per repetition. Summary statistics are first averaged across pairs within
#   a repetition, then across repetitions. Increasing N_pairs reduces
#   within-repetition noise; increasing REP reduces across-repetition noise.
#
# ROBUSTNESS PARAMETERS:
#   To run robustness checks, change any of the parameters below and re-run.
#   Output filenames are generated automatically from the parameter values,
#   so each combination produces a uniquely named CSV



rm(list = ls())
library(truncnorm)

treatments <- c("no_disclosure", "mandatory_disclosure")


N_pairs   <- 5
T         <- 50
REP       <- 100
J         <- 100
P_ex      <- 0.033
dev       <- 0.20    
seed_base <- 100

clamp_int <- function(x) as.integer(max(1, min(10, round(x))))
client_payoff  <- function(state, action)  100 - (state - action)^2
typeA_payoff   <- function(state, action)  100 - (state - action)^2
typeB_payoff   <- function(action)         100 - (10 - action)^2

# alpha_B in [0,1]: alpha_B = 1 -> truthful message (m = s),
#                   alpha_B = 0 -> maximally inflated message (m = 10)
get_advisor_message <- function(state, type, alpha_B) {
  if (type == 1) clamp_int(state)
  else clamp_int(alpha_B * state + (1 - alpha_B) * 10)
}

# tau_raw in [0,1] is rescaled to tau in [1,10].
# If m <= tau: client trusts the message and sets action = m.
# If m >  tau: client discounts toward the prior mean (5.5) using weight delta.
get_client_action <- function(m, tau_raw, delta) {
  tau <- 1 + 9 * tau_raw
  if (m <= tau) clamp_int(m)
  else clamp_int(delta * m + (1 - delta) * 5.5)
}

# CLIENT STRATEGY PARAMETERIZATION
#
# The client's optimal strategy depends on what information is given.
# Under no_disclosure, he/she cannot observe advisor type, so he/she applies
# a single (tau, delta) pair regardless. Under mandatory_disclosure, he/she
# knows the advisor's type and can condition his/her response:
#
#   tau_no  (raw, rescaled to [1,10]):
#     Trust threshold when type is unobserved, or when type = A is revealed.
#     Under no_disclosure this is used always; under mandatory_disclosure it
#     is used only for type-A advisors.
#
#   delta_no:
#     Discounting weight toward prior mean when threshold is exceeded,
#     applied in the same cases as tau_no.
#
#   tau_B  (raw, rescaled to [1,10]):
#     Trust threshold applied only under mandatory_disclosure when type = B
#     is revealed. Nash prediction: tau_B -> 1 (client ignores all messages
#     from known-biased advisors).
#
#   delta_B:
#     Discounting weight toward prior mean for known-biased advisors.
#     Nash prediction: delta_B -> 0 (client fully discounts biased messages).
#
# In total, the client strategy space is four-dimensional. Under no_disclosure
# only components 1-2 are active; components 3-4 are unused but still updated
# by the IEL algorithm (they evolve freely without selection pressure).

get_client_components <- function(tau_no, delta_no, tau_B, delta_B, treatment, type) {
  if (treatment == "no_disclosure" || type == 1)
    list(tau = tau_no, delta = delta_no)
  else
    list(tau = tau_B, delta = delta_B)
}

run_simulation <- function(treatment) {

  cat("Running treatment:", treatment, "\n")
  all_results <- list()

  for (rep in 1:REP) {

    set.seed(seed_base + rep)

    # Each pair has one advisor and one client, each carrying J candidate strategies.
    # Advisor strategies: one-dimensional (alpha_B in [0,1]).
    # Client strategies:  four-dimensional (tau_no, delta_no, tau_B, delta_B), all in [0,1].
    adv_strategies <- lapply(1:N_pairs, function(p) matrix(runif(J), nrow = J, ncol = 1))
    cli_strategies <- lapply(1:N_pairs, function(p) matrix(runif(J * 4), nrow = J, ncol = 4))
    adv_weights    <- lapply(1:N_pairs, function(p) rep(1, J))
    cli_weights    <- lapply(1:N_pairs, function(p) rep(1, J))
    adv_active     <- sample(1:J, N_pairs, replace = TRUE)
    cli_active     <- sample(1:J, N_pairs, replace = TRUE)

    rep_results <- list()

    for (t in 1:T) {

      state_vec <- sample(1:10, N_pairs, replace = TRUE)

      # Type draw: 50/50 unbiased/biased, fixed by the experiment design.
      type_vec  <- sample(c(1, 0), N_pairs, replace = TRUE)

      message_vec    <- numeric(N_pairs)
      action_vec     <- numeric(N_pairs)
      payoff_adv_vec <- numeric(N_pairs)
      payoff_cli_vec <- numeric(N_pairs)
      alpha_B_used   <- rep(NA_real_, N_pairs)
      tau_no_used    <- numeric(N_pairs)
      delta_no_used  <- numeric(N_pairs)
      tau_B_used     <- numeric(N_pairs)
      delta_B_used   <- numeric(N_pairs)
      fixed_vec      <- logical(N_pairs)

      # --- play the game ---
      for (p in 1:N_pairs) {
        state <- state_vec[p]
        type  <- type_vec[p]

        # FIXED-RULE (mandatory disclosure, type A):
        # When the advisor is known to be unbiased and the treatment is mandatory
        # disclosure, the dominant strategy is for the advisor to report truthfully
        # (m = state) and for the client to follow the message (action = m).
        # Technically, the client observes only the message m, not the state directly.
        # Since the unbiased advisor always sends m = state, the client sets
        # action = m = state. This equilibrium is hard-coded rather than
        # leaving it to the IEL algorithm, shortcutting a trivial convergence
        # problem and keeping the simulation focused on the strategically
        # interesting case (type-B interactions).
        
        if (treatment == "mandatory_disclosure" && type == 1) {
          message_vec[p]    <- state                   # advisor reports truthfully
          action_vec[p]     <- message_vec[p]          # client follows the message
          payoff_cli_vec[p] <- client_payoff(state, action_vec[p])
          payoff_adv_vec[p] <- typeA_payoff(state, action_vec[p])
          fixed_vec[p]      <- TRUE
          next
        }

        fixed_vec[p] <- FALSE

        # Retrieve active strategy components for this pair.
        # Renamed from generic 's_*' to avoid clash with 'state' (the game state variable).
        alpha_B_act  <- adv_strategies[[p]][adv_active[p], 1]
        tau_no_act   <- cli_strategies[[p]][cli_active[p], 1]
        delta_no_act <- cli_strategies[[p]][cli_active[p], 2]
        tau_B_act    <- cli_strategies[[p]][cli_active[p], 3]
        delta_B_act  <- cli_strategies[[p]][cli_active[p], 4]

        m    <- get_advisor_message(state, type, alpha_B_act)
        comp <- get_client_components(tau_no_act, delta_no_act, tau_B_act, delta_B_act, treatment, type)
        a    <- get_client_action(m, comp$tau, comp$delta)

        message_vec[p]    <- m
        action_vec[p]     <- a
        payoff_cli_vec[p] <- client_payoff(state, a)
        payoff_adv_vec[p] <- if (type == 1) typeA_payoff(state, a) else typeB_payoff(a)
        alpha_B_used[p]   <- alpha_B_act
        tau_no_used[p]    <- tau_no_act
        delta_no_used[p]  <- delta_no_act
        tau_B_used[p]     <- tau_B_act
        delta_B_used[p]   <- delta_B_act
      }

      # --- experimentation ---
      # Each strategy component is independently nudged with probability P_ex,
      # drawing from a truncated normal centred on its current value.
      adv_exp <- adv_strategies
      cli_exp <- cli_strategies

      for (p in 1:N_pairs) {
        for (j in 1:J) {
          if (runif(1) <= P_ex)
            adv_exp[[p]][j, 1] <- rtruncnorm(1, a = 0, b = 1, mean = adv_exp[[p]][j, 1], sd = dev)
          for (k in 1:4) {
            if (runif(1) <= P_ex)
              cli_exp[[p]][j, k] <- rtruncnorm(1, a = 0, b = 1, mean = cli_exp[[p]][j, k], sd = dev)
          }
        }
      }

      # --- foregone payoffs ---
      # For each candidate strategy j, we compute what payoff it would have earned
      # this round given the realised state and the opponent's actual action.
      #
      # ADVISOR foregone payoff:
      #   Strategy alpha_j implies a counterfactual message m_j.
      #   The client's actual action is then adjusted by the message difference:
      #   a_j = clamp(a_actual + (m_j - m_actual))
      #   This is a linear additive shift and is a standard approximation in IEL
      #   foregone payoff calculations (cf. Arifovic et al., 2024), justified by
      #   the linearity of the client's best-response in the neighbourhood of
      #   the equilibrium.
      #
      # CLIENT foregone payoff:
      #   Strategy (tau_j, delta_j) implies a counterfactual action a_j given the
      #   actual message m_actual that was observed this round.
      adv_foregone <- list()
      cli_foregone <- list()

      for (p in 1:N_pairs) {
        if (fixed_vec[p]) next

        state  <- state_vec[p]
        type   <- type_vec[p]
        m_act  <- message_vec[p]
        a_act  <- action_vec[p]

        adv_foregone[[p]] <- numeric(J)
        cli_foregone[[p]] <- numeric(J)

        for (j in 1:J) {
          m_j <- get_advisor_message(state, type, adv_exp[[p]][j, 1])
          a_j <- clamp_int(a_act + (m_j - m_act))
          adv_foregone[[p]][j] <- if (type == 1) typeA_payoff(state, a_j) else typeB_payoff(a_j)

          comp_j <- get_client_components(cli_exp[[p]][j, 1], cli_exp[[p]][j, 2],
                                          cli_exp[[p]][j, 3], cli_exp[[p]][j, 4],
                                          treatment, type)
          cli_foregone[[p]][j] <- client_payoff(state, get_client_action(m_act, comp_j$tau, comp_j$delta))
        }
      }

      # --- tournament ---
      adv_new    <- adv_exp
      cli_new    <- cli_exp
      adv_wt_new <- list()
      cli_wt_new <- list()

      for (p in 1:N_pairs) {
        if (fixed_vec[p]) {
          adv_wt_new[[p]] <- adv_weights[[p]]
          cli_wt_new[[p]] <- cli_weights[[p]]
          next
        }

        adv_wt_new[[p]] <- numeric(J)
        cli_wt_new[[p]] <- numeric(J)

        for (j in 1:J) {
          c1 <- sample(1:J, 1); c2 <- sample(1:J, 1)
          winner  <- if (adv_foregone[[p]][c1] >= adv_foregone[[p]][c2]) c1 else c2
          adv_new[[p]][j, ]  <- adv_exp[[p]][winner, ]
          adv_wt_new[[p]][j] <- adv_foregone[[p]][winner]

          c1 <- sample(1:J, 1); c2 <- sample(1:J, 1)
          winner  <- if (cli_foregone[[p]][c1] >= cli_foregone[[p]][c2]) c1 else c2
          cli_new[[p]][j, ]  <- cli_exp[[p]][winner, ]
          cli_wt_new[[p]][j] <- cli_foregone[[p]][winner]
        }

        if (sum(adv_wt_new[[p]]) <= 0) adv_wt_new[[p]] <- rep(1, J)
        if (sum(cli_wt_new[[p]]) <= 0) cli_wt_new[[p]] <- rep(1, J)
      }

      adv_strategies <- adv_new
      cli_strategies <- cli_new
      adv_weights    <- adv_wt_new
      cli_weights    <- cli_wt_new

      # --- selection ---
      for (p in 1:N_pairs) {
        adv_active[p] <- sample(1:J, 1, prob = adv_weights[[p]])
        cli_active[p] <- sample(1:J, 1, prob = cli_weights[[p]])
      }

      rep_results[[t]] <- data.frame(
        rep        = rep,
        period     = t,
        pair       = 1:N_pairs,
        state      = state_vec,
        type       = type_vec,
        message    = message_vec,
        action     = action_vec,
        payoff_adv = payoff_adv_vec,
        payoff_cli = payoff_cli_vec,
        alpha_B    = alpha_B_used,
        tau_no     = 1 + 9 * tau_no_used,
        delta_no   = delta_no_used,
        tau_B      = 1 + 9 * tau_B_used,
        delta_B    = delta_B_used
      )
    }

    all_results[[rep]] <- do.call(rbind, rep_results)
    if (rep %% 10 == 0 || rep == 1) cat("  Completed rep", rep, "of", REP, "\n")
  }

  return(do.call(rbind, all_results))
}

results <- list()
for (tr in treatments) results[[tr]] <- run_simulation(tr)

print_summary <- function(df, treatment) {

  cat("\n============================================================\n")
  cat("TREATMENT:", treatment, "\n")
  cat("============================================================\n")

  cat("\nOverall averages (all pairs):\n")
  cat("  Avg client payoff:   ", round(mean(df$payoff_cli), 2), "\n")
  cat("  Avg advisor payoff:  ", round(mean(df$payoff_adv), 2), "\n")
  cat("  Avg message - state: ", round(mean(df$message - df$state), 3), "\n")
  cat("  Avg action  - state: ", round(mean(df$action  - df$state), 3), "\n")

  cat("\nAverages by advisor type (1 = A unbiased, 0 = B biased):\n")
  print(round(aggregate(cbind(message, action, payoff_adv, payoff_cli) ~ type,
                        data = df, FUN = mean), 3))

  df_B   <- df[df$type == 0, ]
  first5 <- df_B[df_B$period <= 5, ]
  last5  <- df_B[df_B$period >= (T - 4), ]

  cat("\nStrategy convergence — biased advisor (Type B) periods only:\n")
  cat("  alpha_B:  ", round(mean(first5$alpha_B, na.rm = TRUE), 3),
      "->", round(mean(last5$alpha_B, na.rm = TRUE), 3), "  (Nash: 0)\n")
  cat("  tau_no:   ", round(mean(first5$tau_no), 2),
      "->", round(mean(last5$tau_no), 2), "  (threshold on [1,10])\n")
  cat("  delta_no: ", round(mean(first5$delta_no), 3),
      "->", round(mean(last5$delta_no), 3), "  (Nash: 0)\n")

  if (treatment == "mandatory_disclosure") {
    cat("  tau_B:    ", round(mean(first5$tau_B), 2),
        "->", round(mean(last5$tau_B), 2), "  (Nash: ~1)\n")
    cat("  delta_B:  ", round(mean(first5$delta_B), 3),
        "->", round(mean(last5$delta_B), 3), "  (Nash: 0)\n")
  }
}

for (tr in treatments) print_summary(results[[tr]], tr)


output_dir <- getwd()
cat("\nSaving results to:", output_dir, "\n")


sig_tag <- paste0("sig", gsub("\\.", "", formatC(dev,  format = "f", digits = 2)))
Pex_tag <- paste0("Pex", gsub("\\.", "", formatC(P_ex, format = "f", digits = 3)))

for (tr in treatments) {
  filename <- paste0(
    "IEL_v3_", tr,
    "_Npairs", N_pairs,
    "_J",      J,
    "_T",      T,
    "_REP",    REP,
    "_",       Pex_tag,
    "_",       sig_tag,
    ".csv"
  )
  write.csv(results[[tr]], file = file.path(output_dir, filename), row.names = FALSE)
  cat("Saved:", filename, "\n")
}

cat("\nDone!\n")
