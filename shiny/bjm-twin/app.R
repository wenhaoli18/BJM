# BJM digital twin demo (Shiny)
#
# Run from R:  shiny::runApp("shiny/bjm-twin")
#
# A patient's twin is a set of futures drawn by simulateTrajectory(),
# conditional on the measurements up to the year the twin is built at.
# "Replay" rebuilds the twin as that year moves through follow-up;
# "Ask the twin" answers questions from the drawn futures.
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
  err_nodata   = c(en = "No visits up to this year", zh = "截至该年没有随访记录")
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

replay_panel <- nav_panel(
  txt("tab_replay"),
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
  txt("tab_ask"),
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
    selectInput("bio", txt("bio"), bio_choices("en")),
    tags$small(class = "text-muted", txt("disclaimer"))
  ),
  navset_card_underline(replay_panel, ask_panel)
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

  output$template <- downloadHandler(
    filename = "bjm_twin_patient_template.csv",
    content = function(file) write.csv(template, file, row.names = FALSE))
}

shinyApp(ui, server)
