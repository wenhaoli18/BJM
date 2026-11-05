## Summaries and figures for the transportBJM() simulation study
suppressMessages(library(ggplot2))
files <- list.files("sim/full", "^rep_[0-9]+\\.rds$", full.names = TRUE)
files <- files[!grepl("rep_000", files)]          # pilot replicate excluded
x <- do.call(rbind, lapply(files, readRDS))
x$method <- sub("200$", "_n200", sub("100$", "_n100", sub("50$", "_n50", x$method)))
x$method <- sub("^citl", "citl", x$method)
sc_order <- c("A k0.5 n100", "A k0.5 n400", "A k1 n100", "A k1 n400", "A k2 n100", "A k2 n400",
              "B assay k1", "B assay k2", "B link k1", "B link k2", "C sparse k2", "C dense k2", "D CR k2")
x$scenario <- factor(x$scenario, sc_order)
cat("replicates:", length(unique(x$rep)), "\n")
dir.create("sim/out", showWarnings = FALSE)

## ---- 1. estimation of the log baseline-hazard shift alpha (once per rep) -----
a <- x[x$landmark == 0 & !is.na(x$alpha_hat), ]
a$alpha_true <- log(a$k)
a$cover <- abs(a$alpha_hat - a$alpha_true) <= 1.96 * a$alpha_se
tab_alpha <- do.call(rbind, lapply(split(a, list(a$scenario, a$method), drop = TRUE), function(d) data.frame(
  scenario = d$scenario[1], method = d$method[1], true_HR = d$k[1],
  mean_HR_hat = mean(exp(d$alpha_hat)), bias_alpha = mean(d$alpha_hat - d$alpha_true),
  sd_alpha = sd(d$alpha_hat), rmse_alpha = sqrt(mean((d$alpha_hat - d$alpha_true)^2)),
  mean_se = mean(d$alpha_se), coverage = mean(d$cover))))
tab_alpha <- tab_alpha[order(tab_alpha$scenario, tab_alpha$method), ]
write.csv(tab_alpha, "sim/out/table_alpha.csv", row.names = FALSE)

## ---- 2. prediction accuracy ---------------------------------------------------
x$calib_ratio <- x$predicted / x$oracle        # mean predicted / mean true risk
tab_pred <- do.call(rbind, lapply(split(x, list(x$scenario, x$landmark, x$method), drop = TRUE), function(d) data.frame(
  scenario = d$scenario[1], landmark = d$landmark[1], method = d$method[1],
  rmse_vs_true_risk = mean(d$rmse_oracle), calib_ratio = mean(d$calib_ratio),
  abs_log_calib = mean(abs(log(d$calib_ratio))),
  brier_excess = mean(d$brier - d$brier_oracle), auc = mean(d$auc, na.rm = TRUE))))
tab_pred <- tab_pred[order(tab_pred$scenario, tab_pred$landmark, tab_pred$method), ]
write.csv(tab_pred, "sim/out/table_prediction.csv", row.names = FALSE)

## ---- figures (palette validated with the dataviz validator, light mode) -------
pal <- c("source_only" = "#2a78d6", "unlabeled" = "#eb6834", "unlabeled_noref" = "#1baf7a",
         "citl_n50" = "#eda100", "baseline_n50" = "#e87ba4", "local_refit" = "#008300")
shp <- c("source_only" = 16, "unlabeled" = 17, "unlabeled_noref" = 2, "citl_n50" = 15,
         "baseline_n50" = 18, "local_refit" = 4)
theme_viz <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e5e0", linewidth = 0.3),
        axis.text = element_text(colour = "#5f5e5a"), axis.title = element_text(colour = "#3b3a37"),
        strip.text = element_text(colour = "#1f1e1c", face = "bold", hjust = 0),
        legend.position = "top", legend.title = element_blank(),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA),
        plot.title = element_text(face = "bold", colour = "#1f1e1c"),
        plot.subtitle = element_text(colour = "#5f5e5a"))

## Figure 1: estimated / true hazard ratio
f1 <- a[a$method %in% c("unlabeled", "unlabeled_noref", "citl_n50"), ]
f1$ratio <- exp(f1$alpha_hat - f1$alpha_true)
s1 <- aggregate(ratio ~ scenario + method, f1, function(v) c(m = exp(mean(log(v))),
                lo = exp(quantile(log(v), 0.1)), hi = exp(quantile(log(v), 0.9))))
s1 <- cbind(s1[1:2], as.data.frame(s1$ratio)); names(s1)[3:5] <- c("m", "lo", "hi")
p1 <- ggplot(s1, aes(scenario, m, colour = method, shape = method)) +
  geom_hline(yintercept = 1, colour = "#9a9992", linewidth = 0.4) +
  geom_linerange(aes(ymin = lo, ymax = hi), position = position_dodge(0.6), linewidth = 0.6) +
  geom_point(position = position_dodge(0.6), size = 2.6) +
  scale_colour_manual(values = pal, breaks = names(pal)) + scale_shape_manual(values = shp, breaks = names(pal)) +
  scale_y_log10(breaks = c(0.5, 0.7, 1, 1.4, 2)) +
  labs(x = NULL, y = "Estimated / true hazard ratio (log scale)",
       title = "Estimating the baseline-hazard shift",
       subtitle = "Geometric mean and 10th-90th percentile over replicates; 1 = unbiased") +
  theme_viz + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave("sim/out/fig1_alpha.png", p1, width = 9, height = 4.8, dpi = 200, bg = "#fcfcfb")

## Figure 2: RMSE of predicted risk vs the true risk
f2 <- tab_pred[tab_pred$method %in% names(pal), ]
f2$landmark <- factor(paste0("Landmark ", f2$landmark, " y, 2-year risk"))
p2 <- ggplot(f2, aes(scenario, rmse_vs_true_risk, colour = method, shape = method)) +
  geom_point(position = position_dodge(0.7), size = 2.4) +
  facet_wrap(~ landmark, ncol = 1) +
  scale_colour_manual(values = pal, breaks = names(pal)) + scale_shape_manual(values = shp, breaks = names(pal)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "RMSE vs. true risk",
       title = "Accuracy of the predicted risks in the target population",
       subtitle = "Root mean squared difference between each patient's predicted and true risk, averaged over replicates") +
  theme_viz + theme(panel.spacing.y = unit(1.2, "lines")) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave("sim/out/fig2_rmse.png", p2, width = 9, height = 6.5, dpi = 200, bg = "#fcfcfb")

options(width = 160)
print(tab_alpha, digits = 3, row.names = FALSE)
print(tab_pred[tab_pred$method %in% names(pal), ], digits = 3, row.names = FALSE)
