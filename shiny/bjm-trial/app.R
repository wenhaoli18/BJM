# BJM clinical trial toolkit (Shiny)
#
# Run from R:  shiny::runApp("shiny/bjm-trial")
#
# Two uses of BJM in clinical trials, illustrated on pbc3. Each "example
# trial" splits pbc3 at random into a historical half, which the models are
# fit on (as a sponsor would use a natural-history cohort or earlier
# trials), and a trial half, whose patients get random entry dates over a
# 3-year accrual period.
#
#   Event forecast  At a data cut, forecast the calendar date of the D-th
#                   event, from each at-risk patient's biomarker history
#                   (simulateTrajectory()), against two Weibull forecasts.
#   Enrichment      Screen candidates for S years and enroll those with the
#                   highest predicted risk (predictRisk()), against Cox
#                   models on the biomarkers at entry or at S.
#
# Research demo only -- not for clinical or trial-design use.

library(shiny)
library(bslib)
library(ggplot2)
suppressMessages({library(BJM); library(survival)})

data(pbc3)
pbc3$id <- as.character(pbc3$id)

## ---------------------------------------------------------------
## Translations (English is the default language)
## ---------------------------------------------------------------
# Chinese strings are only ever vector *values*, never names: names given
# inside c() get mangled to <U+XXXX> when R runs in a non-UTF-8 locale.
# Entries with %1$s etc. are sprintf() templates.
i18n <- list(
  title        = c(en = "BJM Clinical Trial Toolkit", zh = "BJM 临床试验工具箱"),
  trial        = c(en = "Example trial", zh = "示例试验"),
  trial_help   = c(en = "Each example trial splits pbc3 at random: the models are fit on one half (the historical data), and the other half is the trial, with entry dates spread over 3 years.",
                   zh = "每个示例试验把 pbc3 随机分成两半：一半作为历史数据用来拟合模型，另一半作为试验，入组日期分布在 3 年内。"),
  disclaimer   = c(en = "An illustration on the PBC data, not a tool for designing or running real trials.",
                   zh = "这是基于 PBC 数据的演示，不能用于设计或执行真实试验。"),
  tab_fc       = c(en = "Event forecast", zh = "事件数预测"),
  tab_en       = c(en = "Enrichment", zh = "富集设计"),
  progress     = c(en = "Computing...", zh = "正在计算…"),

  fc_help      = c(en = "An event-driven trial is analysed when D events have occurred. At the data cut, each method forecasts the calendar date of the D-th event; the dotted line shows when it actually happened.",
                   zh = "事件驱动的试验在累计发生 D 个事件时进行分析。在数据截止时，每种方法预测第 D 个事件发生的日期；虚线是它实际发生的时间。"),
  fc_D         = c(en = "Target events D (% of the trial's events)", zh = "目标事件数 D（占试验全部事件的 %）"),
  fc_cut       = c(en = "Data cut: share of D observed (%)", zh = "数据截止：已观察到的 D 的比例 (%)"),
  fc_status    = c(en = "Data cut at year %1$s of the trial: %2$s of D = %3$s events observed, %4$s patients still at risk.",
                   zh = "数据截止于试验第 %1$s 年：已观察到 %2$s 个事件（D = %3$s），仍有 %4$s 名患者处于风险中。"),
  h_fc_plot    = c(en = "Cumulative events: observed and forecast", zh = "累计事件数：已观察与预测"),
  cap_fc_plot  = c(en = "Black: events observed by the data cut. Band and line: 5th-95th percentile and median of the forecast. Dashed: the target D; dotted: the actual date of the D-th event.",
                   zh = "黑线：截止前已观察到的事件。带与线：预测的第 5-95 百分位数与中位数。虚线（水平）：目标 D；点线（竖直）：第 D 个事件的实际日期。"),
  h_fc_tab     = c(en = "Forecast date of the D-th event", zh = "第 D 个事件的预测日期"),
  col_method   = c(en = "Method", zh = "方法"),
  col_fc       = c(en = "Forecast (year)", zh = "预测（年）"),
  col_pi       = c(en = "90% interval", zh = "90% 区间"),
  col_err      = c(en = "Error (months)", zh = "误差（月）"),
  actual_row   = c(en = "Actual", zh = "实际"),
  m_bjm        = c(en = "BJM (biomarker histories)", zh = "BJM（生物标志物轨迹）"),
  m_wb_hist    = c(en = "Weibull, historical data", zh = "Weibull，历史数据"),
  m_wb_blind   = c(en = "Weibull, blinded trial data", zh = "Weibull，盲态试验数据"),
  x_cal        = c(en = "Trial calendar time (years)", zh = "试验日历时间 (年)"),
  y_events     = c(en = "Cumulative events", zh = "累计事件数"),
  h_fc_bt      = c(en = "Backtest over 50 example trials", zh = "50 个示例试验的回测"),
  fc_bt_note   = c(en = "In this backtest BJM forecasts the D-th event too early, by about 8 to 13 months, with intervals that are too narrow. The cause is calibration, not ranking: BJM ranks patients' risks well (C-index 0.83 at 6 years), but overestimates the absolute risk of patients in the middle of the risk range, so it predicts too many events too soon. A forecast of event dates needs well-calibrated absolute risks; ranking alone is enough for enrichment (next tab).",
                   zh = "在这次回测中，BJM 预测的第 D 个事件日期偏早约 8 到 13 个月，区间也偏窄。原因在于校准而不是排序：BJM 对患者风险的排序很好（6 年 C-index 0.83），但高估了中等风险患者的绝对风险，因此预测出过多、过早的事件。事件日期预测需要校准良好的绝对风险；而富集设计只依赖排序（见下一个标签页）。"),
  bt_at        = c(en = "Data cut at", zh = "数据截止于"),
  bt_mae       = c(en = "Mean abs. error (months)", zh = "平均绝对误差（月）"),
  bt_bias      = c(en = "Bias (months)", zh = "偏差（月）"),
  bt_cover     = c(en = "90% interval coverage", zh = "90% 区间覆盖率"),
  bt_width     = c(en = "Interval width (months)", zh = "区间宽度（月）"),
  of_D         = c(en = "%s%% of D", zh = "D 的 %s%%"),

  en_help      = c(en = "Candidates are screened for S years. Those still event-free are ranked by their predicted risk of an event in the next 3 years, and only the highest-risk share is enrolled. A higher event rate means fewer patients are needed for the same number of events. Event rates are what actually happened to the enrolled candidates (Kaplan-Meier).",
                   zh = "候选患者先经过 S 年的筛选期。仍未发生事件的患者按其未来 3 年的预测风险排序，只入组风险最高的一部分。事件率越高，达到同样事件数所需的患者越少。事件率是入组候选患者实际发生的情况（Kaplan-Meier 估计）。"),
  en_S         = c(en = "Screening period S (years)", zh = "筛选期 S (年)"),
  en_share     = c(en = "Share of candidates enrolled (%)", zh = "入组的候选患者比例 (%)"),
  en_hr        = c(en = "Hazard ratio to detect", zh = "要检测的风险比 HR"),
  en_power     = c(en = "Power", zh = "检验效能"),
  vb_p         = c(en = "BJM: 3-year event rate of the enrolled", zh = "BJM：入组者的 3 年事件率"),
  vb_n         = c(en = "BJM: patients needed", zh = "BJM：所需患者数"),
  vb_n_all     = c(en = "Without enrichment: %s", zh = "不富集：%s"),
  vb_scr       = c(en = "BJM: patients to screen", zh = "BJM：需筛选的患者数"),
  h_en_tab     = c(en = "Risk models compared", zh = "风险模型比较"),
  col_c        = c(en = "C-index", zh = "C 指数"),
  col_p        = c(en = "Event rate", zh = "事件率"),
  col_n        = c(en = "Patients needed", zh = "所需患者"),
  col_scr      = c(en = "Patients to screen", zh = "需筛选患者"),
  m_baseline   = c(en = "Cox, biomarkers at entry", zh = "Cox，入组时的生物标志物"),
  m_landmark   = c(en = "Cox, last biomarkers at S", zh = "Cox，S 时最近的生物标志物"),
  m_bjm_en     = c(en = "BJM, biomarker history to S", zh = "BJM，截至 S 的生物标志物轨迹"),
  h_en_plot    = c(en = "Patients needed and screened, by share enrolled", zh = "按入组比例：所需与需筛选的患者数"),
  cap_en_plot  = c(en = "Solid: patients needed. Dashed: patients to screen to enrol them. Patients needed come from Schoenfeld's number of events for the hazard ratio (two-sided 5% level, 1:1), divided by the pooled event rate of the two arms.",
                   zh = "实线：所需患者数。虚线：为入组这些患者需要筛选的人数。所需患者数由 Schoenfeld 公式按所设 HR 算出的事件数（双侧 5% 显著性水平，1:1 分配）除以两组合并事件率得到。"),
  x_share      = c(en = "Share of candidates enrolled", zh = "入组的候选患者比例"),
  y_patients   = c(en = "Patients", zh = "患者数"),
  h_en_bt      = c(en = "Backtest over 50 example trials", zh = "50 个示例试验的回测"),
  en_bt_note   = c(en = "With 2 years of screening, BJM ranks candidates best of the three and needs the fewest patients; with 1 year (about three visits) its history adds little over the last values.",
                   zh = "筛选期为 2 年时，BJM 在三种模型中排序最好、所需患者最少；筛选期只有 1 年（约 3 次随访）时，轨迹相比最近一次测量值的增益很小。"),
  bt_S         = c(en = "Screening", zh = "筛选期"),
  bt_share     = c(en = "Enrolled", zh = "入组比例"),
  yr           = c(en = "%s y", zh = "%s 年")
)
tr <- function(key, lang) unname(i18n[[key]][lang])

## ---------------------------------------------------------------
## Data and models
## ---------------------------------------------------------------
bio_vars <- c("serBilir", "albumin", "prothrombin")
base <- pbc3[!duplicated(pbc3$id), c("id", "years", "status3", "age", "sex")]
visits <- pbc3[, c("id", "year", "age", "sex", bio_vars)]
trans <- survivalTrans(c(1, 3, 5, 7))
long_fixed <- setNames(lapply(bio_vars, function(b)
  as.formula(paste(b, "~ year + age + sex + years + years * year"))), bio_vars)
long_random <- setNames(rep(list(~ year | id), length(bio_vars)), bio_vars)
accrual <- 3
H_en <- 3      # enrichment: trial follow-up after screening (years)
fc_sims <- 300

# an example trial: the split, entry dates, and the models fit on the
# historical half
make_trial <- function(seed) {
  set.seed(1000 + seed)
  hist_ids <- sample(base$id, nrow(base) / 2)
  hist <- base[base$id %in% hist_ids, ]
  trial <- base[!base$id %in% hist_ids, ]
  trial$entry <- runif(nrow(trial), 0, accrual)
  hv <- pbc3[pbc3$id %in% hist_ids & pbc3$status3 == 1, ]
  list(hist = hist, trial = trial,
       surv = survivalSub(hist, Surv(years, status3) ~ age + sex, NULL),
       long = suppressWarnings(longitudinalSub(hv, long_fixed, long_random)),
       weibull = survreg(Surv(years, status3) ~ 1, data = hist, dist = "weibull"))
}

# Weibull sampler, conditional on being event-free at t0
weibull_sampler <- function(fit) {
  k <- 1 / fit$scale; lambda <- exp(coef(fit))
  function(t0, n) lambda * ((t0 / lambda)^k - log(runif(n)))^(1 / k)
}

## --- Event forecast ------------------------------------------------
# For each method, a list of n_sim vectors: the calendar dates of all
# events (observed and simulated) in one simulated future of the trial
forecast <- function(tm, D, cut) {
  set.seed(1)
  trial <- tm$trial
  enrolled <- trial[trial$entry < cut, ]
  n_future <- sum(trial$entry >= cut)
  enrolled$t_cut <- cut - enrolled$entry
  enrolled$t_obs <- pmin(enrolled$years, enrolled$t_cut)
  enrolled$d_obs <- as.integer(enrolled$status3 == 1 & enrolled$years <= enrolled$t_cut)
  enrolled$lost <- enrolled$status3 == 0 & enrolled$years < enrolled$t_cut
  at_risk <- enrolled[enrolled$d_obs == 0 & !enrolled$lost, ]
  observed <- with(enrolled[enrolled$d_obs == 1, ], entry + t_obs)
  drop_rate <- max(sum(enrolled$lost), 0.5) / sum(enrolled$t_obs)
  rate_in <- nrow(enrolled) / min(cut, accrual)

  draw_bjm <- function(history, t0) {
    s <- simulateTrajectory(history, tm$long, tm$surv, prediction_time = t0, times = t0,
                            time_variable = "year", trans$survival_variable_all,
                            trans$survival_trans_function, n_sim = fc_sims,
                            seed = sample.int(1e6, 1))
    ifelse(s$event == 1, s$event_time, Inf)
  }
  bjm_risk <- vapply(seq_len(nrow(at_risk)), function(i) {
    p <- at_risk[i, ]
    h <- visits[visits$id == p$id & visits$year <= p$t_cut, ]
    h$years <- NA
    draw_bjm(h, p$t_cut)
  }, numeric(fc_sims))
  bjm_future <- vapply(seq_len(n_future), function(j) {
    donor <- enrolled[sample.int(nrow(enrolled), 1), ]   # covariates not known yet
    h <- data.frame(id = paste0("new", j), year = 0, age = donor$age, sex = donor$sex,
                    serBilir = NA, albumin = NA, prothrombin = NA, years = NA)
    draw_bjm(h, 0)
  }, numeric(fc_sims))
  rw_hist <- weibull_sampler(tm$weibull)
  rw_blind <- weibull_sampler(survreg(Surv(t_obs, d_obs) ~ 1, data = enrolled, dist = "weibull"))

  sims <- lapply(c(bjm = "bjm", wb_hist = "wb_hist", wb_blind = "wb_blind"), function(m)
    lapply(seq_len(fc_sims), function(s) {
      t_risk <- switch(m, bjm = bjm_risk[s, ],
                       wb_hist = rw_hist(at_risk$t_cut, nrow(at_risk)),
                       wb_blind = rw_blind(at_risk$t_cut, nrow(at_risk)))
      ev <- at_risk$entry + t_risk
      ev[ev > cut + rexp(nrow(at_risk), drop_rate)] <- Inf
      t_new <- switch(m, bjm = bjm_future[s, ], wb_hist = rw_hist(0, n_future),
                      wb_blind = rw_blind(0, n_future))
      t_new[t_new > rexp(n_future, drop_rate)] <- Inf
      ev_new <- cut + cumsum(rexp(n_future, rate_in)) + t_new
      sort(c(observed, ev, ev_new)[is.finite(c(observed, ev, ev_new))])
    }))
  list(sims = sims, observed = sort(observed), cut = cut, D = D,
       n_at_risk = nrow(at_risk))
}

## --- Enrichment ----------------------------------------------------
first_vals <- function(ids) {
  v <- pbc3[pbc3$id %in% ids & pbc3$year == 0, c("id", bio_vars)]
  names(v)[-1] <- paste0(bio_vars, "_0"); v
}
last_vals <- function(ids, S) {
  v <- pbc3[pbc3$id %in% ids & pbc3$year <= S, c("id", "year", bio_vars)]
  do.call(rbind, lapply(split(v, v$id), function(d) {
    d <- d[order(d$year), ]
    data.frame(id = d$id[1], t(sapply(bio_vars, function(b) {
      x <- d[[b]][!is.na(d[[b]])]; if (length(x)) x[length(x)] else NA })))
  }))
}
cox_risk <- function(fit, newdata, from, H) {
  sf <- survfit(fit, newdata = newdata)
  s_at <- function(t) as.numeric(summary(sf, times = t, extend = TRUE)$surv)
  1 - s_at(from + H) / (if (from > 0) s_at(from) else 1)
}

# predicted risk of an event in (S, S + H] for each candidate still
# event-free and followed at S, under the three models
enrichment_risks <- function(tm, S) {
  hist <- tm$hist
  cand <- tm$trial[tm$trial$years > S, ]
  fit_b <- coxph(Surv(years, status3) ~ age + sex + serBilir_0 + albumin_0 + prothrombin_0,
                 data = merge(hist, first_vals(hist$id)))
  hl <- merge(hist[hist$years > S, ], last_vals(hist$id[hist$years > S], S))
  hl$t_lm <- hl$years - S
  fit_l <- coxph(Surv(t_lm, status3) ~ age + sex + serBilir + albumin + prothrombin, data = hl)
  cb <- merge(cand, first_vals(cand$id)); cl <- merge(cand, last_vals(cand$id, S))
  risk <- cand[, c("id", "years", "status3")]
  risk$baseline <- cox_risk(fit_b, cb[match(risk$id, cb$id), ], S, H_en)
  risk$landmark <- cox_risk(fit_l, cl[match(risk$id, cl$id), ], 0, H_en)
  risk$bjm <- vapply(risk$id, function(i) {
    h <- visits[visits$id == i & visits$year <= S, ]
    h$years <- NA
    unname(predictRisk(h, tm$long, tm$surv, prediction_time = S, horizon = H_en,
                       time_variable = "year", trans$survival_variable_all,
                       trans$survival_trans_function, bandcount1 = 10,
                       bandcount2 = 10)$risk_prob_1)
  }, numeric(1))
  risk
}

# event probability over (S, S + H] of the enrolled, by Kaplan-Meier
km_prob <- function(d, S) {
  t <- pmin(d$years - S, H_en); e <- d$status3 == 1 & d$years - S <= H_en
  1 - summary(survfit(Surv(t, e) ~ 1), times = H_en, extend = TRUE)$surv
}
# Schoenfeld events for the hazard ratio, over the pooled event rate (1:1)
patients_needed <- function(p, hr, power) {
  d <- 4 * (qnorm(0.975) + qnorm(power))^2 / log(hr)^2
  d / ((p + 1 - (1 - p)^hr) / 2)
}
enrol_top <- function(risk, m, q, S) {
  ok <- !is.na(risk[[m]])
  cutoff <- quantile(risk[[m]][ok], 1 - q, names = FALSE)
  km_prob(risk[ok & risk[[m]] >= cutoff, ], S)
}
c_index <- function(risk, m, S) {
  ok <- !is.na(risk[[m]])
  d <- risk[ok, ]
  concordance(Surv(pmin(d$years - S, H_en), d$status3 == 1 & d$years - S <= H_en) ~
                I(-d[[m]]))$concordance
}

## --- Backtests (50 example trials; the scripts are in backtests/) ----
fc_backtest <- data.frame(
  at = rep(c(30, 50, 70), each = 3),
  method = rep(c("bjm", "wb_hist", "wb_blind"), 3),
  mae = c(13.4, 8.4, 13.1, 11.8, 7.4, 10.2, 8.7, 5.0, 6.5),
  bias = c(-13.0, -1.5, -7.3, -11.2, -2.8, -8.5, -8.2, -2.1, -4.4),
  cover = c(0.32, 0.80, 0.58, 0.34, 0.78, 0.58, 0.46, 0.86, 0.72),
  width = c(14.5, 27.2, 24.2, 13.8, 24.0, 20.5, 12.8, 20.5, 19.0))
en_backtest <- data.frame(
  S = rep(c(1, 2), each = 6), share = rep(rep(c(50, 25), each = 3), 2),
  method = rep(c("baseline", "landmark", "bjm"), 4),
  c = c(0.821, 0.848, 0.847, 0.821, 0.848, 0.847, 0.834, 0.853, 0.869, 0.834, 0.853, 0.869),
  p = c(0.406, 0.415, 0.416, 0.598, 0.643, 0.648, 0.459, 0.447, 0.471, 0.638, 0.684, 0.716),
  n = c(703, 686, 684, 467, 431, 427, 620, 638, 605, 438, 404, 383))

fc_methods <- c("bjm", "wb_hist", "wb_blind")
en_methods <- c("baseline", "landmark", "bjm")
method_label <- function(m, lang)
  tr(c(bjm = "m_bjm", wb_hist = "m_wb_hist", wb_blind = "m_wb_blind",
       baseline = "m_baseline", landmark = "m_landmark")[[m]], lang)
en_label <- function(m, lang) if (m == "bjm") tr("m_bjm_en", lang) else method_label(m, lang)
method_colours <- c(bjm = "#2980b9", wb_hist = "#7f8c8d", wb_blind = "#d35400",
                    baseline = "#7f8c8d", landmark = "#d35400")

fmt <- function(x, d = 1) formatC(x, format = "f", digits = d)
pct <- function(p) sprintf("%.0f%%", 100 * p)
html_table <- function(head, rows)
  tags$div(class = "table-responsive",
           tags$table(class = "table table-sm mb-0",
                      tags$thead(tags$tr(lapply(head, tags$th))),
                      tags$tbody(lapply(rows, function(r) tags$tr(lapply(r, tags$td))))))

## ---------------------------------------------------------------
## UI
## ---------------------------------------------------------------
# every static piece of text is an output, so it follows input$lang
# (renderUI rather than renderText: renderText goes through cat(), which
# escapes Chinese to <U+XXXX> in a non-UTF-8 locale)
txt <- function(key) uiOutput(paste0("txt_", key), inline = TRUE)
caption <- function(key) card_footer(tags$small(class = "text-muted", txt(key)))

fc_panel <- nav_panel(
  txt("tab_fc"),
  tags$p(class = "text-muted", txt("fc_help")),
  layout_columns(
    col_widths = breakpoints(sm = 12, lg = c(6, 6)),
    sliderInput("fc_D", txt("fc_D"), min = 50, max = 90, value = 80, step = 10),
    sliderInput("fc_cut", txt("fc_cut"), min = 20, max = 80, value = 50, step = 10)),
  tags$p(uiOutput("fc_status")),
  card(card_header(txt("h_fc_plot")), plotOutput("fc_plot", height = 520), caption("cap_fc_plot")),
  card(card_header(txt("h_fc_tab")), uiOutput("fc_table")),
  card(card_header(txt("h_fc_bt")), uiOutput("fc_bt"),
       card_footer(tags$small(txt("fc_bt_note"))))
)

en_panel <- nav_panel(
  txt("tab_en"),
  tags$p(class = "text-muted", txt("en_help")),
  layout_columns(
    col_widths = breakpoints(sm = 12, lg = c(3, 3, 3, 3)),
    selectInput("en_S", txt("en_S"), c(1, 2), selected = 2),
    sliderInput("en_share", txt("en_share"), min = 20, max = 100, value = 50, step = 5),
    sliderInput("en_hr", txt("en_hr"), min = 0.5, max = 0.9, value = 0.7, step = 0.05),
    selectInput("en_power", txt("en_power"), c("80%" = 0.8, "90%" = 0.9))),
  layout_columns(
    fill = FALSE,
    value_box(txt("vb_p"), uiOutput("vb_p"), theme = "primary"),
    value_box(txt("vb_n"), uiOutput("vb_n"), theme = "success"),
    value_box(txt("vb_scr"), uiOutput("vb_scr"), theme = "secondary")),
  card(card_header(txt("h_en_tab")), uiOutput("en_table")),
  card(card_header(txt("h_en_plot")), plotOutput("en_plot", height = 360), caption("cap_en_plot")),
  card(card_header(txt("h_en_bt")), uiOutput("en_bt"),
       card_footer(tags$small(txt("en_bt_note"))))
)

ui <- page_sidebar(
  title = txt("title"),
  window_title = "BJM Clinical Trial Toolkit",
  fillable = FALSE,  # let the page scroll instead of squashing the plots
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  # bslib keeps the main area at least 576px wide while the sidebar is open,
  # which overflows a narrow window or the iframe on the pkgdown site
  tags$head(tags$style(".bslib-sidebar-layout { --bslib-page-main-min-width: 0px; }")),
  sidebar = sidebar(
    width = 280,
    radioButtons("lang", "Language / 语言",
                 setNames(c("en", "zh"), c("English", "中文")),
                 selected = "en", inline = TRUE),
    selectInput("seed", txt("trial"), 1:5, selected = 1),
    tags$small(class = "text-muted", txt("trial_help")),
    tags$small(class = "text-muted", txt("disclaimer"))
  ),
  navset_card_underline(fc_panel, en_panel)
)

## ---------------------------------------------------------------
## Server
## ---------------------------------------------------------------
server <- function(input, output, session) {

  lang <- reactive(if (identical(input$lang, "zh")) "zh" else "en")
  for (key in names(i18n)) local({
    k <- key
    output[[paste0("txt_", k)]] <- renderUI(tr(k, lang()))
  })

  tm <- reactive(withProgress(message = tr("progress", isolate(lang())),
                              make_trial(as.integer(input$seed))))

  ## --- Event forecast ----------------------------------------------
  fc_setup <- reactive({
    t <- tm()$trial
    ev <- sort(with(t[t$status3 == 1, ], entry + years))
    D <- round(input$fc_D / 100 * length(ev))
    list(D = D, actual = ev[D], cut = ev[max(1, round(input$fc_cut / 100 * D))] + 1e-6)
  })
  fc <- reactive({
    s <- fc_setup()
    withProgress(message = tr("progress", isolate(lang())), forecast(tm(), s$D, s$cut))
  })
  fc_dates <- reactive({
    f <- fc()
    sapply(f$sims, function(m) vapply(m, function(ev) if (length(ev) >= f$D) ev[f$D] else Inf,
                                      numeric(1)))
  })

  output$fc_status <- renderUI({
    f <- fc(); lg <- lang()
    tags$span(class = "text-muted",
              sprintf(tr("fc_status", lg), fmt(f$cut), length(f$observed), f$D, f$n_at_risk))
  })

  output$fc_plot <- renderPlot({
    f <- fc(); s <- fc_setup(); lg <- lang()
    grid <- seq(f$cut, max(s$actual, quantile(fc_dates()[is.finite(fc_dates())], 0.99)) + 1,
                length.out = 120)
    bands <- do.call(rbind, lapply(fc_methods, function(m) {
      counts <- vapply(f$sims[[m]], function(ev) findInterval(grid, ev), numeric(length(grid)))
      data.frame(method = m, time = grid,
                 lo = apply(counts, 1, quantile, 0.05), med = apply(counts, 1, median),
                 hi = apply(counts, 1, quantile, 0.95))
    }))
    labels <- setNames(vapply(fc_methods, method_label, "", lang = lg), fc_methods)
    bands$method <- factor(labels[bands$method], levels = labels)
    obs <- data.frame(time = c(0, f$observed, f$cut), n = c(0, seq_along(f$observed), length(f$observed)))
    ggplot(bands, aes(time)) +
      geom_ribbon(aes(ymin = lo, ymax = hi, fill = method), alpha = 0.25) +
      geom_line(aes(y = med, colour = method), linewidth = 1) +
      geom_step(data = obs, aes(time, n), linewidth = 0.8) +
      geom_hline(yintercept = f$D, linetype = 2) +
      geom_vline(xintercept = s$actual, linetype = 3, linewidth = 0.8) +
      facet_wrap(~ method, ncol = 1) +
      scale_fill_manual(values = setNames(method_colours[fc_methods], labels), guide = "none") +
      scale_colour_manual(values = setNames(method_colours[fc_methods], labels), guide = "none") +
      labs(x = tr("x_cal", lg), y = tr("y_events", lg)) +
      theme_bw(base_size = 13)
  })

  output$fc_table <- renderUI({
    d <- fc_dates(); a <- fc_setup()$actual; lg <- lang()
    rows <- c(list(list(tags$strong(tr("actual_row", lg)), tags$strong(fmt(a, 2)), "", "")),
              lapply(fc_methods, function(m) {
                q <- quantile(d[, m], c(0.05, 0.5, 0.95), names = FALSE)
                list(method_label(m, lg), tags$strong(fmt(q[2], 2)),
                     sprintf("(%s, %s)", fmt(q[1], 2), fmt(q[3], 2)), fmt(12 * (q[2] - a)))
              }))
    html_table(c(tr("col_method", lg), tr("col_fc", lg), tr("col_pi", lg), tr("col_err", lg)), rows)
  })

  output$fc_bt <- renderUI({
    lg <- lang()
    rows <- lapply(seq_len(nrow(fc_backtest)), function(i) {
      b <- fc_backtest[i, ]
      list(sprintf(tr("of_D", lg), b$at), method_label(b$method, lg), fmt(b$mae), fmt(b$bias),
           fmt(b$cover, 2), fmt(b$width))
    })
    html_table(c(tr("bt_at", lg), tr("col_method", lg), tr("bt_mae", lg), tr("bt_bias", lg),
                 tr("bt_cover", lg), tr("bt_width", lg)), rows)
  })

  ## --- Enrichment --------------------------------------------------
  en_risk <- reactive({
    S <- as.numeric(input$en_S)
    withProgress(message = tr("progress", isolate(lang())), enrichment_risks(tm(), S))
  })
  en_at <- function(m, q) {
    S <- as.numeric(input$en_S)
    p <- enrol_top(en_risk(), m, q, S)
    n <- patients_needed(p, input$en_hr, as.numeric(input$en_power))
    c(p = p, n = n, screened = n / q)
  }
  en_all <- reactive(patients_needed(km_prob(en_risk(), as.numeric(input$en_S)),
                                     input$en_hr, as.numeric(input$en_power)))

  output$vb_p <- renderUI(pct(en_at("bjm", input$en_share / 100)[["p"]]))
  output$vb_n <- renderUI({
    x <- en_at("bjm", input$en_share / 100)
    tagList(round(x[["n"]]),
            tags$p(class = "fs-6 mb-0", sprintf(tr("vb_n_all", lang()), round(en_all()))))
  })
  output$vb_scr <- renderUI(round(en_at("bjm", input$en_share / 100)[["screened"]]))

  output$en_table <- renderUI({
    lg <- lang(); S <- as.numeric(input$en_S); q <- input$en_share / 100
    rows <- lapply(en_methods, function(m) {
      x <- en_at(m, q)
      list(en_label(m, lg), fmt(c_index(en_risk(), m, S), 3), pct(x[["p"]]),
           tags$strong(round(x[["n"]])), round(x[["screened"]]))
    })
    html_table(c(tr("col_method", lg), tr("col_c", lg), tr("col_p", lg), tr("col_n", lg),
                 tr("col_scr", lg)), rows)
  })

  output$en_plot <- renderPlot({
    lg <- lang()
    shares <- seq(0.2, 1, by = 0.05)
    curves <- do.call(rbind, lapply(en_methods, function(m) do.call(rbind, lapply(shares, function(q) {
      x <- en_at(m, q)
      data.frame(method = m, share = q, needed = x[["n"]], screened = x[["screened"]])
    }))))
    labels <- setNames(vapply(en_methods, en_label, "", lang = lg), en_methods)
    curves$method <- factor(labels[curves$method], levels = labels)
    ggplot(curves, aes(share, colour = method)) +
      geom_line(aes(y = needed), linewidth = 1) +
      geom_line(aes(y = screened), linewidth = 0.8, linetype = 2) +
      geom_vline(xintercept = input$en_share / 100, colour = "grey50") +
      scale_x_continuous(tr("x_share", lg), labels = scales::percent) +
      scale_colour_manual(values = setNames(method_colours[en_methods], labels), name = NULL) +
      labs(y = tr("y_patients", lg)) +
      theme_bw(base_size = 13) + theme(legend.position = "bottom") +
      guides(colour = guide_legend(ncol = 1))
  })

  output$en_bt <- renderUI({
    lg <- lang()
    rows <- lapply(seq_len(nrow(en_backtest)), function(i) {
      b <- en_backtest[i, ]
      list(sprintf(tr("yr", lg), b$S), paste0(b$share, "%"), en_label(b$method, lg),
           fmt(b$c, 3), pct(b$p), b$n)
    })
    html_table(c(tr("bt_S", lg), tr("bt_share", lg), tr("col_method", lg), tr("col_c", lg),
                 tr("col_p", lg), tr("col_n", lg)), rows)
  })
}

shinyApp(ui, server)
