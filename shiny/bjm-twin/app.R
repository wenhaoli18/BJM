# BJM digital twin demo (Shiny)
#
# Run from R:  shiny::runApp("shiny/bjm-twin")
#
# A patient's twin is a set of futures drawn by simulateTrajectory(),
# conditional on the measurements up to the year the twin is built at.
# "Replay" rebuilds the twin as that year moves through follow-up;
# "Ask the twin" answers questions from the drawn futures; "Trial design"
# uses every pbc3 patient's baseline twin as a prognostic covariate in a
# randomized trial (PROCOVA) and shows the sample size it saves.
# Research demo only -- not for clinical use.

library(shiny)
library(bslib)
library(ggplot2)
suppressMessages(library(BJM))

data(pbc3)

## ---------------------------------------------------------------
## Translations (English is the default language)
## ---------------------------------------------------------------
# Chinese strings are only ever vector *values*, never names: names given
# inside c() get mangled to <U+XXXX> when R runs in a non-UTF-8 locale.
# Entries with %1$s etc. are sprintf() templates.
i18n <- list(
  title        = c(en = "BJM Digital Twin Demo", zh = "BJM 数字孪生演示"),
  source       = c(en = "Data source", zh = "数据来源"),
  src_demo     = c(en = "Example patient (pbc3)", zh = "示例患者 (pbc3)"),
  src_upload   = c(en = "Upload CSV", zh = "上传 CSV"),
  pid          = c(en = "Patient ID", zh = "患者 ID"),
  file         = c(en = "Follow-up data for one patient (.csv)", zh = "单个患者的随访数据 (.csv)"),
  browse       = c(en = "Browse...", zh = "选择文件..."),
  template     = c(en = "Download CSV template", zh = "下载 CSV 模板"),
  landmark     = c(en = "Twin built at year", zh = "孪生体构建时点 (年)"),
  H            = c(en = "Forecast horizon (years)", zh = "预测跨度 (年)"),
  n_sim        = c(en = "Number of simulated futures", zh = "模拟未来的数量"),
  bio          = c(en = "Biomarker to plot", zh = "绘图的生物标志物"),
  disclaimer   = c(en = "The twin is a statistical predictive model fit on the PBC data, not a mechanistic or causal model. Research demo only; not for clinical decisions.",
                   zh = "孪生体是基于 PBC 数据拟合的统计预测模型，不是机理模型或因果模型。仅用于研究演示，不可用于临床决策。"),
  tab_replay   = c(en = "Replay the twin", zh = "回放孪生体"),
  tab_ask      = c(en = "Ask the twin", zh = "向孪生体提问"),
  replay_help  = c(en = "Drag the slider, or press its play button: the twin is rebuilt from the measurements up to that year. Filled points were used by the twin; open points were measured later and show how its forecast compares with what happened.",
                   zh = "拖动滑块，或点击滑块下的播放按钮：孪生体会用截至该年的测量值重新构建。实心点是孪生体使用的数据；空心点是之后才测得的，可用来对比预测与实际。"),
  vb_any       = c(en = "Death or transplant within window", zh = "窗口内死亡或移植"),
  vb_death     = c(en = "Death within window", zh = "窗口内死亡"),
  vb_tx        = c(en = "Transplant within window", zh = "窗口内移植"),
  vb_actual    = c(en = "What actually happened", zh = "实际结局"),
  actual_death = c(en = "Died at year %.1f", zh = "第 %.1f 年死亡"),
  actual_tx    = c(en = "Transplanted at year %.1f", zh = "第 %.1f 年接受移植"),
  actual_cens  = c(en = "Event-free at last follow-up (year %.1f)", zh = "末次随访时未发生事件 (第 %.1f 年)"),
  actual_na    = c(en = "Unknown (uploaded data)", zh = "未知（上传数据）"),
  h_bio_fc     = c(en = "Forecast biomarker", zh = "生物标志物预测"),
  h_cif        = c(en = "Forecast events", zh = "事件预测"),
  cap_bio_fc   = c(en = "Line and band: median and 5th-95th percentile among the simulated futures still event-free at each time.",
                   zh = "线与带：每个时点仍未发生事件的模拟未来中的中位数与第 5-95 百分位数。"),
  cap_cif      = c(en = "Share of the simulated futures with each event by each time. The dotted line marks what actually happened, if it falls in the window.",
                   zh = "各时点前发生各类事件的模拟未来所占比例。虚线标出实际发生的结局（如落在窗口内）。"),
  used         = c(en = "used by the twin", zh = "孪生体已使用"),
  later        = c(en = "measured later", zh = "之后测得"),
  x_time       = c(en = "Follow-up time (years)", zh = "随访时间 (年)"),
  y_prob       = c(en = "Cumulative probability", zh = "累积概率"),
  ev_any       = c(en = "Death or transplant", zh = "死亡或移植"),
  ev_death     = c(en = "Death", zh = "死亡"),
  ev_tx        = c(en = "Transplant", zh = "移植"),
  ev_free      = c(en = "Event-free", zh = "未发生事件"),
  h_spag       = c(en = "A sample of the simulated futures", zh = "模拟未来的抽样"),
  cap_spag     = c(en = "Each grey line is one simulated future of the biomarker, stopping at its event; its end point shows how that future ends within the window. Blue: median of all futures still event-free.",
                   zh = "每条灰线是一个模拟的未来指标轨迹，在事件发生时停止；终点颜色表示该未来在窗口内的结局。蓝线：所有仍未发生事件的未来的中位数。"),
  outcome      = c(en = "End of the future", zh = "未来的结局"),
  h_q1         = c(en = "How likely is an event?", zh = "发生事件的可能性有多大？"),
  q_within     = c(en = "Within (years)", zh = "时间范围 (年)"),
  a_q1         = c(en = "Probability within %1$s years (by year %2$s)", zh = "%1$s 年内（第 %2$s 年前）的概率"),
  h_q2         = c(en = "Will the biomarker cross a threshold?", zh = "生物标志物会越过某个阈值吗？"),
  q_bio        = c(en = "Biomarker", zh = "生物标志物"),
  q_dir        = c(en = "Direction", zh = "方向"),
  dir_above    = c(en = "above", zh = "高于"),
  dir_below    = c(en = "below", zh = "低于"),
  q_thr        = c(en = "Threshold", zh = "阈值"),
  a_cross      = c(en = "Goes %1$s %2$s at some time within %3$s years, before any event",
                   zh = "%3$s 年内、在发生事件之前，某个时点%1$s %2$s"),
  a_end        = c(en = "Event-free at year %1$s, with the biomarker %2$s %3$s then",
                   zh = "第 %1$s 年时仍未发生事件，且当时指标%2$s %3$s"),
  h_q3         = c(en = "When would the event happen?", zh = "事件会在什么时候发生？"),
  a_median     = c(en = "Median year of death or transplant, among the futures with one before year %s",
                   zh = "在第 %s 年前发生死亡或移植的未来中，事件发生年份的中位数"),
  a_free       = c(en = "Probability of being event-free through year %s", zh = "到第 %s 年仍未发生事件的概率"),
  yr_fmt       = c(en = "year %.1f", zh = "第 %.1f 年"),
  none         = c(en = "none of the futures", zh = "没有任何未来"),
  progress     = c(en = "Simulating futures...", zh = "正在模拟未来…"),
  err_missing  = c(en = "Missing columns: ", zh = "缺少列: "),
  err_nodata   = c(en = "No visits up to this year", zh = "截至该年没有随访记录"),
  tab_trial    = c(en = "Trial design", zh = "试验设计"),
  sidebar_trial = c(en = "The patient settings apply to the other two tabs.",
                    zh = "患者相关设置仅用于另外两个标签页。"),
  trial_help   = c(en = "In a randomized trial, each participant's twin can be built from their baseline data alone and its prediction of the endpoint used as a covariate in the analysis (the PROCOVA idea). Because the prediction uses only data from before randomization, the adjustment does not bias the treatment effect, but it shrinks its variance by a factor of 1 - ρ², where ρ is the correlation between the prediction and the outcome. Here ρ is estimated from the 312 pbc3 participants' baseline twins (200 futures each).",
                   zh = "在随机对照试验中，可以只用每位受试者的基线数据构建其孪生体，并把孪生体对终点的预测作为分析中的协变量（即 PROCOVA 的思路）。由于预测只用到随机化之前的数据，这种调整不会使治疗效应产生偏倚，但能把其方差缩小为原来的 1 - ρ² 倍，其中 ρ 是预测值与结局之间的相关系数。这里的 ρ 由 pbc3 中 312 名受试者的基线孪生体（每人 200 个模拟未来）估计。"),
  trial_warn   = c(en = "An illustration on the pbc3 data, not a trial-design tool for real use. The twins were fit on the same pbc3 data, outcomes included, so ρ is optimistic; in a real trial the model is fit on historical data and frozen before the trial starts, and the sample size comes from a prespecified, validated calculation.",
                   zh = "这只是基于 pbc3 数据的演示，不是可用于实际试验设计的工具。孪生体是在同一份 pbc3 数据（包括结局）上拟合的，因此 ρ 偏乐观；在真实试验中，模型应基于历史数据拟合并在试验开始前固定，样本量也应由预先设定并经过验证的计算得出。"),
  h_design     = c(en = "Trial settings", zh = "试验设置"),
  endpoint     = c(en = "Endpoint", zh = "终点"),
  ep_bio       = c(en = "Biomarker at year X", zh = "第 X 年的生物标志物"),
  ep_event     = c(en = "Death or transplant within X years", zh = "X 年内死亡或移植"),
  ep_year      = c(en = "Year X", zh = "第 X 年"),
  effect_d     = c(en = "Treatment effect (difference in means, in SDs of the outcome)",
                   zh = "治疗效应（均值差，以结局的标准差为单位）"),
  effect_rrr   = c(en = "Treatment effect (relative risk reduction, %)",
                   zh = "治疗效应（相对风险降低，%）"),
  alpha        = c(en = "Significance level (two-sided)", zh = "显著性水平（双侧）"),
  power        = c(en = "Power", zh = "检验效能"),
  vb_rho       = c(en = "Correlation ρ of twin and outcome", zh = "孪生体与结局的相关系数 ρ"),
  vb_n0        = c(en = "Sample size without the twin", zh = "不使用孪生体的样本量"),
  vb_n1        = c(en = "Sample size with the twin", zh = "使用孪生体的样本量"),
  vb_saved     = c(en = "Fewer participants", zh = "受试者减少"),
  h_trial_tab  = c(en = "Analysis compared", zh = "分析方法比较"),
  col_analysis = c(en = "Covariate", zh = "协变量"),
  col_rho      = c(en = "ρ", zh = "ρ"),
  col_n        = c(en = "Total sample size", zh = "总样本量"),
  an_none      = c(en = "None (unadjusted)", zh = "无（未调整）"),
  an_base      = c(en = "Baseline value of the biomarker", zh = "该生物标志物的基线值"),
  an_base_bili = c(en = "Baseline bilirubin", zh = "基线胆红素"),
  an_twin      = c(en = "Twin prediction", zh = "孪生体预测"),
  trial_base_note = c(en = "Adjusting for the baseline value alone is shown for comparison: a twin is worth it when it combines several biomarkers and covariates into a better prediction than any one of them.",
                      zh = "仅调整基线值的结果作为对照列出：只有当孪生体把多个生物标志物和协变量综合成比其中任何一个都更好的预测时，它才真正有价值。"),
  trial_n_note = c(en = "%1$s pbc3 participants with an outcome; control-arm risk %2$s. 1:1 randomization; normal-approximation sample size, multiplied by 1 - ρ² for an adjusted analysis.",
                   zh = "pbc3 中有结局数据的受试者 %1$s 名；对照组风险 %2$s。1:1 随机化；样本量采用正态近似，调整分析时乘以 1 - ρ²。"),
  trial_n_note_bio = c(en = "%s pbc3 participants with a visit near year X. 1:1 randomization; normal-approximation sample size, multiplied by 1 - ρ² for an adjusted analysis.",
                       zh = "pbc3 中在第 X 年附近有随访的受试者 %s 名。1:1 随机化；样本量采用正态近似，调整分析时乘以 1 - ρ²。"),
  h_scatter    = c(en = "Twin prediction and outcome in pbc3", zh = "pbc3 中孪生体预测与结局"),
  cap_scatter  = c(en = "One point per participant: the twin's prediction from baseline data against the outcome observed.",
                   zh = "每个点代表一名受试者：基于基线数据的孪生体预测值与实际观察到的结局。"),
  h_ncurve     = c(en = "Sample size by treatment effect", zh = "不同治疗效应下的样本量"),
  cap_ncurve   = c(en = "Total sample size for the chosen significance level and power; the dashed line marks the chosen effect.",
                   zh = "在所选显著性水平和检验效能下的总样本量；虚线标出所选治疗效应。"),
  x_twin_pred  = c(en = "Twin prediction", zh = "孪生体预测"),
  x_twin_risk  = c(en = "Twin's predicted risk", zh = "孪生体预测的风险"),
  y_observed   = c(en = "Observed at year %s", zh = "第 %s 年观察值"),
  y_n          = c(en = "Total sample size", zh = "总样本量"),
  obs_event    = c(en = "Event within %s years", zh = "%s 年内发生事件"),
  obs_free     = c(en = "Event-free at year %s", zh = "第 %s 年未发生事件")
)
tr <- function(key, lang) unname(i18n[[key]][lang])

## ---------------------------------------------------------------
## Model fitting (once, at start-up), as in the Digital twins article:
## death or transplantation (status3), which of the two (status4:
## 0 death, 1 transplantation), and three biomarkers
## ---------------------------------------------------------------
bio_vars <- c("serBilir", "albumin", "prothrombin")
bio_labels <- list(
  en = c("log(serum bilirubin)", "Albumin", "log(prothrombin time)"),
  zh = c("log(血清胆红素)", "白蛋白", "log(凝血酶原时间)"))

long_fixed <- lapply(bio_vars, function(b)
  as.formula(paste(b, "~ year + age + sex + years + years * year")))
long_random <- rep(list(~ year | id), length(bio_vars))
names(long_fixed) <- names(long_random) <- bio_vars

surv_fit <- survivalSub(pbc3[!duplicated(pbc3$id), ],
                        form_marginal_surv = Surv(years, status3) ~ age + sex,
                        form_conditional_cr = status4 ~ years + age + sex)
long_fit <- longitudinalSub(pbc3[pbc3$status3 == 1, ], long_fixed, long_random)
trans    <- survivalTrans(c(1, 3, 5, 7))

# simulateTrajectory() reports events after the last follow-up as
# event-free through it, so forecasts stop there
max_follow_up <- max(pbc3$years)

## ---------------------------------------------------------------
## Trial design (once, at start-up), as in "Twins in a randomized trial"
## in the Digital twins article: each participant's twin is built from
## their baseline measurements alone, and its prediction of the endpoint
## is the mean over 200 futures. With truncate = FALSE the biomarker draws
## are not cut off at the drawn event time.
## ---------------------------------------------------------------
trial_years <- 1:5
baseline <- pbc3[pbc3$year == 0, ]
baseline_twins <- baseline
baseline_twins$years <- NA
baseline_twins$status3 <- NA
baseline_twins$status4 <- NA
trial_sims <- as.data.frame(simulateTrajectory(
  baseline_twins, long_fit, surv_fit, prediction_time = 0, times = trial_years,
  time_variable = "year", trans$survival_variable_all,
  trans$survival_trans_function, n_sim = 200, truncate = FALSE, seed = 1))
trial_sims$event_by <- trial_sims$event == 1 & trial_sims$event_time <= trial_sims$year
# one row per participant and year: the twin's mean biomarkers and risk
twin_pred <- aggregate(cbind(serBilir, albumin, prothrombin, event_by) ~ id + year,
                       data = trial_sims, FUN = mean)
rm(trial_sims, baseline_twins)

# participants' twin prediction, observed outcome and baseline covariate
# for an endpoint; censored participants without an outcome are left out
trial_data <- function(endpoint, x, b = "serBilir") {
  pred <- twin_pred[twin_pred$year == x, ]
  if (endpoint == "bio") {
    visit <- pbc3[abs(pbc3$year - x) < 0.5 & !is.na(pbc3[[b]]), ]
    visit <- visit[!duplicated(visit$id), ]
    d <- data.frame(id = visit$id, outcome = visit[[b]])
    d$twin <- pred[[b]][match(d$id, pred$id)]
    d$base <- baseline[[b]][match(d$id, baseline$id)]
  } else {
    known <- baseline[baseline$years >= x | baseline$status3 == 1, ]
    d <- data.frame(id = known$id,
                    outcome = as.numeric(known$status3 == 1 & known$years <= x),
                    base = known$serBilir)
    d$twin <- pred$event_by[match(d$id, pred$id)]
  }
  d[complete.cases(d), ]
}

# total sample size of a 1:1 trial, two-sided test, normal approximation;
# an analysis adjusted for a covariate with correlation rho needs 1 - rho^2
# times as many (PROCOVA)
n_continuous <- function(d, alpha, power, rho = 0)
  2 * ceiling(2 * (qnorm(1 - alpha / 2) + qnorm(power))^2 / d^2 * (1 - rho^2))
n_binary <- function(p0, rrr, alpha, power, rho = 0) {
  p1 <- p0 * (1 - rrr)
  2 * ceiling((qnorm(1 - alpha / 2) + qnorm(power))^2 *
                (p0 * (1 - p0) + p1 * (1 - p1)) / (p0 - p1)^2 * (1 - rho^2))
}

required_cols <- c("year", "age", "sex", bio_vars)
template <- pbc3[pbc3$id == 2, required_cols]

prepare_patient <- function(dat, lang = "en") {
  missing <- setdiff(required_cols, names(dat))
  if (length(missing))
    stop(tr("err_missing", lang), paste(missing, collapse = ", "))
  dat <- dat[order(dat$year), required_cols]
  dat$id <- "patient"
  dat
}

## ---------------------------------------------------------------
## Twin helpers
## ---------------------------------------------------------------
# n_sim futures from the measurements up to `landmark`
build_twin <- function(patient, landmark, horizon, n_sim) {
  history <- patient[patient$year <= landmark + 1e-8, ]
  history$years <- NA    # event time and type are not known yet
  history$status4 <- NA
  times <- seq(landmark, min(landmark + horizon, max_follow_up), by = 0.25)
  sims <- as.data.frame(simulateTrajectory(
    history, long_fit, surv_fit, prediction_time = landmark, times = times,
    time_variable = "year", trans$survival_variable_all,
    trans$survival_trans_function, n_sim = n_sim, seed = 1))
  futures <- sims[!duplicated(sims$sim), ]
  list(sims = sims, futures = futures, landmark = landmark, times = times,
       end = max(times))
}

# P(event of each kind by time t), from one row per future
event_prob <- function(futures, t) {
  by_t <- futures$event == 1 & futures$event_time <= t + 1e-8
  c(any = mean(by_t),
    death = mean(by_t & futures$status4 == 0),
    tx = mean(by_t & futures$status4 == 1))
}

pct <- function(p) sprintf("%.0f%%", 100 * p)
fmt_year <- function(x) format(round(x, 2))

event_colours <- function(lang)
  setNames(c("#c0392b", "#8e44ad", "#7f8c8d"),
           c(tr("ev_death", lang), tr("ev_tx", lang), tr("ev_free", lang)))

## ---------------------------------------------------------------
## UI
## ---------------------------------------------------------------
# every static piece of text is an output, so it follows input$lang
# (renderUI rather than renderText: renderText goes through cat(), which
# escapes Chinese to <U+XXXX> in a non-UTF-8 locale)
txt <- function(key) uiOutput(paste0("txt_", key), inline = TRUE)

source_choices <- function(lang)
  setNames(c("demo", "upload"), c(tr("src_demo", lang), tr("src_upload", lang)))
bio_choices <- function(lang) setNames(bio_vars, bio_labels[[lang]])
dir_choices <- function(lang)
  setNames(c("above", "below"), c(tr("dir_above", lang), tr("dir_below", lang)))
endpoint_choices <- function(lang)
  setNames(c("bio", "event"), c(tr("ep_bio", lang), tr("ep_event", lang)))

replay_panel <- nav_panel(
  txt("tab_replay"), value = "replay",
  tags$p(class = "text-muted", txt("replay_help")),
  layout_columns(
    fill = FALSE,
    value_box(txt("vb_any"), uiOutput("vb_any"), theme = "danger"),
    value_box(txt("vb_death"), uiOutput("vb_death"), theme = "secondary"),
    value_box(txt("vb_tx"), uiOutput("vb_tx"), theme = "secondary"),
    value_box(txt("vb_actual"), uiOutput("vb_actual"), theme = "primary")
  ),
  layout_columns(
    col_widths = breakpoints(sm = 12, xl = c(6, 6)),  # stacked unless wide
    card(card_header(txt("h_bio_fc")), plotOutput("replay_bio", height = 330),
         card_footer(tags$small(class = "text-muted", txt("cap_bio_fc")))),
    card(card_header(txt("h_cif")), plotOutput("replay_cif", height = 330),
         card_footer(tags$small(class = "text-muted", txt("cap_cif"))))
  )
)

ask_panel <- nav_panel(
  txt("tab_ask"), value = "ask",
  card(card_header(txt("h_spag")), plotOutput("spaghetti", height = 330),
       card_footer(tags$small(class = "text-muted", txt("cap_spag")))),
  layout_columns(
    col_widths = breakpoints(sm = 12, lg = c(6, 6, 12), xl = c(4, 4, 4)),
    card(card_header(txt("h_q1")),
         sliderInput("q1_x", txt("q_within"), min = 0.5, max = 5, value = 3, step = 0.5),
         uiOutput("a_q1")),
    card(card_header(txt("h_q2")),
         selectInput("q_bio", txt("q_bio"), bio_choices("en")),
         layout_columns(
           radioButtons("q_dir", txt("q_dir"), dir_choices("en"), inline = TRUE),
           numericInput("q_thr", txt("q_thr"), value = 1, step = 0.1)),
         sliderInput("q2_x", txt("q_within"), min = 0.5, max = 5, value = 3, step = 0.5),
         uiOutput("a_q2")),
    card(card_header(txt("h_q3")), uiOutput("a_q3"))
  )
)

trial_panel <- nav_panel(
  txt("tab_trial"), value = "trial",
  tags$p(class = "text-muted", txt("trial_help")),
  tags$div(class = "alert alert-warning py-2", txt("trial_warn")),
  layout_columns(
    fill = FALSE,
    value_box(txt("vb_rho"), uiOutput("vb_rho"), theme = "secondary"),
    value_box(txt("vb_n0"), uiOutput("vb_n0"), theme = "secondary"),
    value_box(txt("vb_n1"), uiOutput("vb_n1"), theme = "primary"),
    value_box(txt("vb_saved"), uiOutput("vb_saved"), theme = "success")),
  layout_columns(
    col_widths = c(4, 8),
    card(card_header(txt("h_design")),
         radioButtons("ep", txt("endpoint"), endpoint_choices("en")),
         conditionalPanel("input.ep == 'bio'",
                          selectInput("ep_bio", txt("q_bio"), bio_choices("en"))),
         sliderInput("ep_x", txt("ep_year"), min = min(trial_years),
                     max = max(trial_years), value = 2, step = 1),
         conditionalPanel("input.ep == 'bio'",
                          sliderInput("eff_d", txt("effect_d"), min = 0.1, max = 1,
                                      value = 0.3, step = 0.05)),
         conditionalPanel("input.ep == 'event'",
                          sliderInput("eff_rrr", txt("effect_rrr"), min = 10, max = 60,
                                      value = 30, step = 5)),
         layout_columns(
           selectInput("alpha", txt("alpha"), c(0.05, 0.025, 0.01)),
           selectInput("power", txt("power"), c("80%" = 0.8, "90%" = 0.9)))),
    card(card_header(txt("h_trial_tab")), uiOutput("trial_table"),
         card_footer(tags$small(class = "text-muted", uiOutput("trial_note"))))
  ),
  layout_columns(
    card(card_header(txt("h_scatter")), plotOutput("trial_scatter", height = 330),
         card_footer(tags$small(class = "text-muted", txt("cap_scatter")))),
    card(card_header(txt("h_ncurve")), plotOutput("trial_ncurve", height = 330),
         card_footer(tags$small(class = "text-muted", txt("cap_ncurve"))))
  )
)

ui <- page_sidebar(
  title = txt("title"),
  window_title = "BJM Digital Twin",
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
    # the trial design tab uses all of pbc3, not one patient
    conditionalPanel(
      "input.tabs != 'trial'",
      radioButtons("source", txt("source"), source_choices("en")),
      conditionalPanel(
        "input.source == 'demo'",
        selectInput("pid", txt("pid"), levels(pbc3$id), selected = "2")),
      conditionalPanel(
        "input.source == 'upload'",
        fileInput("file", txt("file"), accept = ".csv", buttonLabel = txt("browse")),
        downloadLink("template", txt("template"))),
      sliderInput("landmark", txt("landmark"), min = 0, max = 10, value = 3, step = 0.25,
                  animate = animationOptions(interval = 1500)),
      sliderInput("H", txt("H"), min = 1, max = 8, value = 5, step = 1),
      selectInput("n_sim", txt("n_sim"), c(200, 500, 1000), selected = 500),
      selectInput("bio", txt("bio"), bio_choices("en"))),
    conditionalPanel(
      "input.tabs == 'trial'",
      tags$p(class = "text-muted", txt("sidebar_trial"))),
    tags$small(class = "text-muted", txt("disclaimer"))
  ),
  navset_card_underline(id = "tabs", replay_panel, ask_panel, trial_panel)
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

  # choice labels can't be outputs, so relabel them in place
  observeEvent(lang(), ignoreInit = TRUE, {
    updateRadioButtons(session, "source", choices = source_choices(lang()),
                       selected = input$source)
    updateSelectInput(session, "bio", choices = bio_choices(lang()), selected = input$bio)
    updateSelectInput(session, "q_bio", choices = bio_choices(lang()), selected = input$q_bio)
    updateRadioButtons(session, "q_dir", choices = dir_choices(lang()),
                       selected = input$q_dir, inline = TRUE)
    updateRadioButtons(session, "ep", choices = endpoint_choices(lang()), selected = input$ep)
    updateSelectInput(session, "ep_bio", choices = bio_choices(lang()), selected = input$ep_bio)
  })

  patient <- reactive({
    if (input$source == "demo") {
      prepare_patient(pbc3[pbc3$id == input$pid, ])
    } else {
      req(input$file)
      tryCatch(prepare_patient(read.csv(input$file$datapath), isolate(lang())),
               error = function(e) {
                 showNotification(conditionMessage(e), type = "error"); NULL })
    }
  })

  # the real outcome, known for the example patients only
  actual <- reactive({
    if (input$source != "demo") return(NULL)
    p <- pbc3[pbc3$id == input$pid, ][1, ]
    list(time = p$years, event = p$status3 == 1, type = p$status4)
  })

  # keep the landmark slider within the patient's visits
  observeEvent(patient(), {
    yr <- patient()$year
    if (is.null(yr) || !length(yr)) return()
    last <- floor(max(yr, na.rm = TRUE) * 4) / 4
    updateSliderInput(session, "landmark", max = max(last, 0.25),
                      value = min(isolate(input$landmark), last))
  })

  # default threshold: the patient's last value of that biomarker that the
  # twin has seen (not one measured after the year it is built at)
  observeEvent(list(patient(), input$q_bio), {
    p <- patient()
    v <- p[[input$q_bio]][p$year <= isolate(input$landmark) + 1e-8]
    v <- v[!is.na(v)]
    if (length(v)) updateNumericInput(session, "q_thr", value = round(v[length(v)], 1))
  })

  # the "within" sliders cannot reach past the forecast horizon
  observeEvent(twin(), {
    span <- twin()$end - twin()$landmark
    for (id in c("q1_x", "q2_x"))
      updateSliderInput(session, id, max = span, value = min(input[[id]], span))
  })

  twin <- reactive({
    p <- req(patient())
    validate(need(any(p$year <= input$landmark + 1e-8), tr("err_nodata", lang())))
    withProgress(message = tr("progress", isolate(lang())), {
      build_twin(p, input$landmark, input$H, as.integer(input$n_sim))
    })
  })

  ## --- Replay ----------------------------------------------------

  window_probs <- reactive(event_prob(twin()$futures, twin()$end))
  output$vb_any <- renderUI(pct(window_probs()[["any"]]))
  output$vb_death <- renderUI(pct(window_probs()[["death"]]))
  output$vb_tx <- renderUI(pct(window_probs()[["tx"]]))
  output$vb_actual <- renderUI({
    a <- actual(); lg <- lang()
    if (is.null(a)) return(tags$span(style = "font-size: 1rem", tr("actual_na", lg)))
    key <- if (!a$event) "actual_cens" else if (a$type == 1) "actual_tx" else "actual_death"
    tags$span(style = "font-size: 1.1rem", sprintf(tr(key, lg), a$time))
  })

  output$replay_bio <- renderPlot({
    tw <- twin(); lg <- lang(); b <- input$bio
    vals <- tw$sims[!is.na(tw$sims[[b]]), ]
    bands <- do.call(rbind, lapply(split(vals[[b]], vals$year), function(v)
      data.frame(lower = quantile(v, 0.05), median = median(v), upper = quantile(v, 0.95))))
    bands$year <- as.numeric(rownames(bands))
    measured <- patient()[!is.na(patient()[[b]]), c("year", b)]
    names(measured)[2] <- "value"
    measured$used <- ifelse(measured$year <= tw$landmark + 1e-8,
                            tr("used", lg), tr("later", lg))
    shapes <- setNames(c(16, 1), c(tr("used", lg), tr("later", lg)))
    ggplot(bands, aes(x = year)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2980b9", alpha = 0.2) +
      geom_line(aes(y = median), colour = "#2980b9", linewidth = 1) +
      geom_point(data = measured, aes(y = value, shape = used), size = 2.5) +
      scale_shape_manual(values = shapes, name = NULL) +
      geom_vline(xintercept = tw$landmark, linetype = 2, colour = "grey40") +
      labs(x = tr("x_time", lg), y = bio_labels[[lg]][match(b, bio_vars)]) +
      theme_bw(base_size = 14) + theme(legend.position = "bottom")
  })

  output$replay_cif <- renderPlot({
    tw <- twin(); lg <- lang()
    probs <- t(vapply(tw$times, function(t) event_prob(tw$futures, t), numeric(3)))
    curves <- data.frame(
      year = rep(tw$times, 3),
      prob = c(probs[, "any"], probs[, "death"], probs[, "tx"]),
      event = factor(rep(c(tr("ev_any", lg), tr("ev_death", lg), tr("ev_tx", lg)),
                         each = length(tw$times)),
                     levels = c(tr("ev_any", lg), tr("ev_death", lg), tr("ev_tx", lg))))
    cols <- setNames(c("#2c3e50", "#c0392b", "#8e44ad"), levels(curves$event))
    g <- ggplot(curves, aes(year, prob, colour = event)) +
      geom_step(linewidth = 1) +
      scale_colour_manual(values = cols, name = NULL) +
      guides(colour = guide_legend(nrow = 2)) +
      scale_y_continuous(tr("y_prob", lg), labels = scales::percent,
                         limits = c(0, NA)) +
      labs(x = tr("x_time", lg)) +
      theme_bw(base_size = 14) + theme(legend.position = "bottom")
    a <- actual()
    if (!is.null(a) && a$time >= tw$landmark && a$time <= tw$end) {
      key <- if (!a$event) "actual_cens" else if (a$type == 1) "actual_tx" else "actual_death"
      g <- g + geom_vline(xintercept = a$time, linetype = 3, linewidth = 0.8) +
        annotate("label", x = a$time, y = Inf, vjust = 1.2, size = 3.5,
                 label = sprintf(tr(key, lg), a$time))
    }
    g
  })

  ## --- Ask the twin ----------------------------------------------

  output$spaghetti <- renderPlot({
    tw <- twin(); lg <- lang(); b <- input$bio
    shown <- tw$sims[tw$sims$sim <= 40 & !is.na(tw$sims[[b]]), ]
    ends <- shown[!duplicated(shown$sim, fromLast = TRUE), ]
    ends$outcome <- with(ends, ifelse(event == 1 & event_time <= tw$end + 1e-8,
                                      ifelse(status4 == 1, tr("ev_tx", lg), tr("ev_death", lg)),
                                      tr("ev_free", lg)))
    vals <- tw$sims[!is.na(tw$sims[[b]]), ]
    med <- aggregate(vals[[b]], list(year = vals$year), median)
    hist <- patient()[patient()$year <= tw$landmark + 1e-8 & !is.na(patient()[[b]]), ]
    ggplot() +
      geom_line(data = shown, aes(year, .data[[b]], group = sim),
                colour = "grey60", alpha = 0.6) +
      geom_line(data = med, aes(year, x), colour = "#2980b9", linewidth = 1.2) +
      geom_point(data = ends, aes(year, .data[[b]], colour = outcome), size = 2.2) +
      scale_colour_manual(values = event_colours(lg), name = tr("outcome", lg)) +
      geom_point(data = hist, aes(year, .data[[b]]), size = 2.5) +
      geom_line(data = hist, aes(year, .data[[b]]), colour = "grey30") +
      geom_vline(xintercept = tw$landmark, linetype = 2, colour = "grey40") +
      labs(x = tr("x_time", lg), y = bio_labels[[lg]][match(b, bio_vars)]) +
      theme_bw(base_size = 14) + theme(legend.position = "bottom")
  })

  answer_row <- function(label, value)
    tags$tr(tags$td(label), tags$td(class = "text-end", tags$strong(value)))

  output$a_q1 <- renderUI({
    tw <- twin(); lg <- lang()
    by <- tw$landmark + input$q1_x
    p <- event_prob(tw$futures, by)
    tagList(
      tags$p(class = "text-muted mb-1",
             sprintf(tr("a_q1", lg), fmt_year(input$q1_x), fmt_year(by))),
      tags$table(class = "table table-sm mb-0", tags$tbody(
        answer_row(tr("ev_any", lg), pct(p[["any"]])),
        answer_row(tr("ev_death", lg), pct(p[["death"]])),
        answer_row(tr("ev_tx", lg), pct(p[["tx"]])))))
  })

  output$a_q2 <- renderUI({
    tw <- twin(); lg <- lang(); b <- input$q_bio
    thr <- req(input$q_thr)
    by <- tw$landmark + input$q2_x
    crosses <- function(v) if (input$q_dir == "above") v > thr else v < thr
    future <- tw$sims[tw$sims$year > tw$landmark + 1e-8 & tw$sims$year <= by + 1e-8, ]
    p_cross <- mean(vapply(split(future[[b]], future$sim), function(v)
      any(crosses(v[!is.na(v)])), logical(1)))
    # the simulated time closest to `by`
    at <- tw$sims[abs(tw$sims$year - tw$times[which.min(abs(tw$times - by))]) < 1e-8, ]
    p_end <- mean(!is.na(at[[b]]) & crosses(at[[b]]))
    dir <- tr(paste0("dir_", input$q_dir), lg)
    tags$table(class = "table table-sm mb-0", tags$tbody(
      answer_row(sprintf(tr("a_cross", lg), dir, thr, fmt_year(input$q2_x)), pct(p_cross)),
      answer_row(sprintf(tr("a_end", lg), fmt_year(by), dir, thr), pct(p_end))))
  })

  output$a_q3 <- renderUI({
    tw <- twin(); lg <- lang()
    f <- tw$futures
    ev <- f$event_time[f$event == 1]
    tags$table(class = "table table-sm mb-0", tags$tbody(
      answer_row(sprintf(tr("a_median", lg), fmt_year(max_follow_up)),
                 if (length(ev)) sprintf(tr("yr_fmt", lg), median(ev)) else tr("none", lg)),
      answer_row(sprintf(tr("a_free", lg), fmt_year(tw$end)),
                 pct(1 - event_prob(f, tw$end)[["any"]]))))
  })

  ## --- Trial design ----------------------------------------------

  trial <- reactive({
    d <- trial_data(input$ep, input$ep_x, input$ep_bio)
    alpha <- as.numeric(input$alpha); power <- as.numeric(input$power)
    p0 <- mean(d$outcome)
    if (input$ep == "bio") {
      n_for <- function(rho, eff) n_continuous(eff, alpha, power, rho)
      eff <- input$eff_d
    } else {
      n_for <- function(rho, eff) n_binary(p0, eff / 100, alpha, power, rho)
      eff <- input$eff_rrr
    }
    rho <- c(none = 0, base = cor(d$base, d$outcome), twin = cor(d$twin, d$outcome))
    list(d = d, p0 = p0, rho = rho, n = vapply(rho, n_for, numeric(1), eff = eff),
         n_for = n_for)
  })

  output$vb_rho <- renderUI(sprintf("%.2f", trial()$rho[["twin"]]))
  output$vb_n0 <- renderUI(format(trial()$n[["none"]], big.mark = ","))
  output$vb_n1 <- renderUI(format(trial()$n[["twin"]], big.mark = ","))
  output$vb_saved <- renderUI(pct(1 - trial()$n[["twin"]] / trial()$n[["none"]]))

  output$trial_table <- renderUI({
    tl <- trial(); lg <- lang()
    base_label <- if (input$ep == "bio") tr("an_base", lg) else tr("an_base_bili", lg)
    row <- function(label, k, class = NULL)
      tags$tr(class = class, tags$td(label),
              tags$td(class = "text-end", if (k == "none") "-" else sprintf("%.2f", tl$rho[[k]])),
              tags$td(class = "text-end", format(tl$n[[k]], big.mark = ",")),
              tags$td(class = "text-end",
                      if (k == "none") "" else sprintf("-%s", pct(1 - tl$n[[k]] / tl$n[["none"]]))))
    tagList(
      tags$table(class = "table table-sm mb-2",
        tags$thead(tags$tr(tags$th(tr("col_analysis", lg)),
                           tags$th(class = "text-end", tr("col_rho", lg)),
                           tags$th(class = "text-end", tr("col_n", lg)),
                           tags$th())),
        tags$tbody(row(tr("an_none", lg), "none"), row(base_label, "base"),
                   row(tr("an_twin", lg), "twin", class = "table-primary"))),
      tags$small(class = "text-muted", tr("trial_base_note", lg)))
  })

  output$trial_note <- renderUI({
    tl <- trial(); lg <- lang()
    if (input$ep == "bio") sprintf(tr("trial_n_note_bio", lg), nrow(tl$d))
    else sprintf(tr("trial_n_note", lg), nrow(tl$d), pct(tl$p0))
  })

  output$trial_scatter <- renderPlot({
    tl <- trial(); lg <- lang(); x <- input$ep_x
    if (input$ep == "bio") {
      lab <- bio_labels[[lg]][match(input$ep_bio, bio_vars)]
      ggplot(tl$d, aes(twin, outcome)) +
        geom_abline(linetype = 2, colour = "grey50") +
        geom_point(colour = "#2980b9", alpha = 0.6, size = 2) +
        labs(x = paste0(tr("x_twin_pred", lg), ": ", lab),
             y = paste0(sprintf(tr("y_observed", lg), x), ": ", lab)) +
        theme_bw(base_size = 14)
    } else {
      levs <- c(sprintf(tr("obs_free", lg), x), sprintf(tr("obs_event", lg), x))
      tl$d$observed <- factor(levs[tl$d$outcome + 1], levels = levs)
      ggplot(tl$d, aes(twin, observed)) +
        geom_boxplot(outlier.shape = NA, width = 0.5) +
        geom_jitter(aes(colour = observed), height = 0.15, width = 0, alpha = 0.6,
                    size = 2, show.legend = FALSE) +
        scale_colour_manual(values = c("#7f8c8d", "#c0392b")) +
        scale_x_continuous(tr("x_twin_risk", lg), labels = scales::percent) +
        labs(y = NULL) +
        theme_bw(base_size = 14)
    }
  })

  output$trial_ncurve <- renderPlot({
    tl <- trial(); lg <- lang()
    if (input$ep == "bio") {
      grid <- seq(0.1, 1, by = 0.05); chosen <- input$eff_d; xlab <- tr("effect_d", lg)
    } else {
      grid <- seq(10, 60, by = 1); chosen <- input$eff_rrr; xlab <- tr("effect_rrr", lg)
    }
    base_label <- if (input$ep == "bio") tr("an_base", lg) else tr("an_base_bili", lg)
    labels <- c(tr("an_none", lg), base_label, tr("an_twin", lg))
    curves <- data.frame(
      effect = rep(grid, 3),
      n = unlist(lapply(tl$rho, function(r) vapply(grid, function(e) tl$n_for(r, e), numeric(1)))),
      analysis = factor(rep(labels, each = length(grid)), levels = labels))
    ggplot(curves, aes(effect, n, colour = analysis)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = chosen, linetype = 2, colour = "grey40") +
      scale_colour_manual(values = setNames(c("#7f8c8d", "#e67e22", "#2980b9"), labels),
                          name = NULL) +
      scale_y_log10(tr("y_n", lg), labels = scales::comma) +
      labs(x = xlab) +
      guides(colour = guide_legend(nrow = 2)) +
      theme_bw(base_size = 14) + theme(legend.position = "bottom")
  })

  output$template <- downloadHandler(
    filename = "bjm_twin_patient_template.csv",
    content = function(file) write.csv(template, file, row.names = FALSE))
}

shinyApp(ui, server)
