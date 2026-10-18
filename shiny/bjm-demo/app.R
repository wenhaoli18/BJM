# BJM dynamic prediction demo (Shiny)
#
# Run from R:  shiny::runApp("shiny/bjm-demo")
#
# The models are fit once on pbc3 at start-up; each prediction then only
# calls predictRisk()/predictLongitudinal() for one patient, which is fast.
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
i18n <- list(
  title       = c(en = "BJM Dynamic Prediction Demo", zh = "BJM 动态预测演示"),
  source      = c(en = "Data source", zh = "数据来源"),
  src_demo    = c(en = "Example patient (pbc3)", zh = "示例患者 (pbc3)"),
  src_upload  = c(en = "Upload CSV", zh = "上传 CSV"),
  pid         = c(en = "Patient ID", zh = "患者 ID"),
  file        = c(en = "Follow-up data for one patient (.csv)", zh = "单个患者的随访数据 (.csv)"),
  browse      = c(en = "Browse...", zh = "选择文件..."),
  template    = c(en = "Download CSV template", zh = "下载 CSV 模板"),
  t           = c(en = "Prediction time (years)", zh = "预测时点 (年)"),
  H           = c(en = "Prediction window (years)", zh = "预测窗口 (年)"),
  go         = c(en = "Predict", zh = "开始预测"),
  disclaimer  = c(en = "Models are fit on the PBC data set. Research demo only; not for clinical decisions.",
                  zh = "模型基于 PBC 数据集拟合，仅用于研究演示，不可用于临床决策。"),
  risk_end    = c(en = "Event risk at end of window", zh = "窗口末的事件风险"),
  n_visits    = c(en = "Visits used", zh = "已用随访次数"),
  h_risk      = c(en = "Cumulative event risk", zh = "累积事件风险"),
  h_bio       = c(en = "Biomarkers: observed history + predicted future",
                  zh = "生物标志物：历史观测 + 未来预测"),
  h_pred      = c(en = "Predicted biomarker values at end of window",
                  zh = "窗口末的生物标志物预测值"),
  col_bio     = c(en = "Biomarker", zh = "生物标志物"),
  col_last    = c(en = "Last observed", zh = "最近观测值"),
  col_pred    = c(en = "Predicted at year ", zh = "预测值：第 "),
  col_pred2   = c(en = "", zh = " 年"),
  col_ci      = c(en = "95% interval", zh = "95% 区间"),
  h_data      = c(en = "Follow-up data used", zh = "使用的随访数据"),
  x_time      = c(en = "Follow-up time (years)", zh = "随访时间 (年)"),
  y_risk      = c(en = "Cumulative event risk", zh = "事件累积风险"),
  y_bio       = c(en = "Biomarker value", zh = "指标值"),
  caption     = c(en = "Points = observed; blue line/band = predicted mean and 95% interval. Predictions are conditional on the patient still being event-free at that future time.",
                  zh = "点 = 历史观测；蓝线/带 = 预测均值与 95% 区间。预测以患者在该未来时点仍未发生事件为条件。"),
  progress    = c(en = "Running BJM...", zh = "BJM 计算中…"),
  err_missing = c(en = "Missing columns: ", zh = "缺少列: "),
  err_nodata  = c(en = "No visits before the prediction time", zh = "预测时点之前没有随访记录")
)
tr <- function(key, lang) unname(i18n[[key]][lang])

## ---------------------------------------------------------------
## Model fitting (once, at start-up)
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
                        Surv(years, status3) ~ age + sex, NULL)
long_fit <- longitudinalSub(pbc3[pbc3$status3 == 1, ], long_fixed, long_random)
trans    <- survivalTrans(c(1, 3, 5, 7))

required_cols <- c("year", "age", "sex", bio_vars)

## ---------------------------------------------------------------
## Prediction helpers
## ---------------------------------------------------------------
# risk of the event in (t, t + h], for a grid of horizons h
risk_curve <- function(dat, t, horizons) {
  data.frame(horizon = horizons, time = t + horizons,
             risk = vapply(horizons, function(h)
               unname(predictRisk(dat, long_fit, surv_fit,
                                  prediction_time = t, horizon = h,
                                  time_variable = "year",
                                  trans$survival_variable_all,
                                  trans$survival_trans_function,
                                  bandcount1 = 10, bandcount2 = 10)$risk_prob_1),
               numeric(1)))
}

# predicted mean and 95% interval of one biomarker at future times
bio_curve <- function(dat, t, horizons, bio_i) {
  do.call(rbind, lapply(horizons, function(h) {
    p <- predictLongitudinal(bio_i = bio_i, dat, long_fit, surv_fit,
                             prediction_time = t, horizon = h,
                             time_variable = "year",
                             trans$survival_variable_all,
                             trans$survival_trans_function,
                             bandcount2 = 40, bandcount3 = 400)
    cdf <- cumsum(p$Y_density[, 1]); cdf <- cdf / cdf[length(cdf)]
    q <- function(a) p$Y_all[which(cdf >= a)[1]]
    data.frame(time = t + h, mean = unname(p$Y_predict),
               lo = q(0.025), hi = q(0.975))
  }))
}

# make an uploaded / selected patient look like BJM's prediction input
prepare_patient <- function(dat, lang = "en") {
  missing <- setdiff(required_cols, names(dat))
  if (length(missing))
    stop(tr("err_missing", lang), paste(missing, collapse = ", "))
  dat <- dat[order(dat$year), required_cols]
  dat$id <- "patient"
  dat$years <- NA_real_  # event time is unknown; BJM integrates over it
  dat
}

template <- pbc3[pbc3$id == 2 & pbc3$year <= 3, required_cols]

## ---------------------------------------------------------------
## UI
## ---------------------------------------------------------------
# every static piece of text is an output, so it follows input$lang
# (renderUI rather than renderText: renderText goes through cat(), which
# escapes Chinese to <U+XXXX> in a non-UTF-8 locale)
txt <- function(key) uiOutput(paste0("txt_", key), inline = TRUE)

source_choices <- function(lang)
  setNames(c("demo", "upload"), c(tr("src_demo", lang), tr("src_upload", lang)))

ui <- page_sidebar(
  title = txt("title"),
  window_title = "BJM Dynamic Prediction",
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
    sliderInput("t", txt("t"), min = 0, max = 10, value = 3, step = 0.25),
    sliderInput("H", txt("H"), min = 1, max = 5, value = 3, step = 0.5),
    actionButton("go", txt("go"), class = "btn-primary"),
    tags$small(class = "text-muted", txt("disclaimer"))
  ),
  layout_columns(
    fill = FALSE,
    value_box(txt("risk_end"), textOutput("risk_end"), theme = "danger"),
    value_box(txt("n_visits"), textOutput("n_visits"), theme = "primary")
  ),
  layout_columns(
    col_widths = c(5, 7),
    card(card_header(txt("h_risk")), plotOutput("risk_plot", height = 300)),
    card(card_header(txt("h_pred")), uiOutput("pred_table"))
  ),
  card(card_header(txt("h_bio")), plotOutput("bio_plot", height = 320),
       card_footer(tags$small(class = "text-muted", txt("caption")))),
  card(card_header(txt("h_data")), tableOutput("data_table"))
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
  })

  patient_raw <- reactive({
    if (input$source == "demo") {
      pbc3[pbc3$id == input$pid, ]
    } else {
      req(input$file)
      read.csv(input$file$datapath)
    }
  })

  # keep the prediction-time slider within the patient's follow-up
  observeEvent(patient_raw(), {
    yr <- patient_raw()$year
    if (is.null(yr) || !length(yr)) return()
    last <- max(yr, na.rm = TRUE)
    updateSliderInput(session, "t", max = max(ceiling(last), 1),
                      value = min(input$t, round(last * 4) / 4))
  })

  result <- eventReactive(input$go, {
    lg <- isolate(lang())
    dat <- tryCatch(prepare_patient(patient_raw(), lg),
                    error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    req(dat)
    dat <- dat[dat$year <= input$t, ]
    validate(need(nrow(dat) > 0, tr("err_nodata", lg)))

    horizons <- seq(0.25, input$H, by = 0.25)
    bio_horizons <- unique(c(horizons[horizons %% 0.5 == 0], input$H))
    withProgress(message = tr("progress", lg), value = 0, {
      risk <- risk_curve(dat, input$t, horizons)
      incProgress(0.25)
      bio <- do.call(rbind, lapply(seq_along(bio_vars), function(i) {
        out <- cbind(bio_i = i, bio_curve(dat, input$t, bio_horizons, i))
        incProgress(0.75 / length(bio_vars))
        out
      }))
    })
    list(dat = dat, t = input$t, H = input$H, risk = risk, bio = bio)
  })

  output$risk_end <- renderText({
    r <- result()$risk
    sprintf("%.1f%%", 100 * r$risk[nrow(r)])
  })
  output$n_visits <- renderText(nrow(result()$dat))

  output$risk_plot <- renderPlot({
    r <- result()
    ggplot(rbind(data.frame(horizon = 0, time = r$t, risk = 0), r$risk),
           aes(time, risk)) +
      geom_step(linewidth = 1, colour = "#c0392b") +
      geom_vline(xintercept = r$t, linetype = 2) +
      scale_y_continuous(tr("y_risk", lang()), labels = scales::percent) +
      labs(x = tr("x_time", lang())) +
      theme_bw(base_size = 14)
  })

  # predicted mean / 95% interval of every biomarker at t + H
  # (built as HTML rather than renderTable, which escapes Chinese headers)
  output$pred_table <- renderUI({
    r <- result(); lg <- lang()
    end <- r$bio[abs(r$bio$time - (r$t + r$H)) < 1e-8, ]
    last_obs <- vapply(bio_vars, function(b) {
      v <- r$dat[[b]][!is.na(r$dat[[b]])]
      if (length(v)) v[length(v)] else NA_real_
    }, numeric(1))
    head <- c(tr("col_bio", lg), tr("col_last", lg),
              paste0(tr("col_pred", lg), format(r$t + r$H), tr("col_pred2", lg)),
              tr("col_ci", lg))
    rows <- lapply(seq_along(bio_vars), function(i) {
      e <- end[end$bio_i == i, ]
      tags$tr(tags$td(bio_labels[[lg]][i]),
              tags$td(sprintf("%.2f", last_obs[i])),
              tags$td(tags$strong(sprintf("%.2f", e$mean))),
              tags$td(sprintf("(%.2f, %.2f)", e$lo, e$hi)))
    })
    tags$table(class = "table table-sm mb-0",
               tags$thead(tags$tr(lapply(head, tags$th))),
               tags$tbody(rows))
  })

  output$bio_plot <- renderPlot({
    r <- result(); lg <- lang()
    lab <- function(i) factor(bio_labels[[lg]][i], levels = bio_labels[[lg]])
    pred <- transform(r$bio, label = lab(bio_i))
    obs <- do.call(rbind, lapply(seq_along(bio_vars), function(i)
      data.frame(label = lab(i), time = r$dat$year, value = r$dat[[bio_vars[i]]])))
    obs <- obs[!is.na(obs$value), ]
    ggplot() +
      geom_ribbon(data = pred, aes(time, ymin = lo, ymax = hi),
                  fill = "#2980b9", alpha = 0.2) +
      geom_line(data = pred, aes(time, mean), colour = "#2980b9", linewidth = 1) +
      geom_point(data = obs, aes(time, value), size = 2.5) +
      geom_line(data = obs, aes(time, value), colour = "grey40") +
      geom_vline(xintercept = r$t, linetype = 2) +
      facet_wrap(~ label, nrow = 1, scales = "free_y") +
      labs(x = tr("x_time", lg), y = tr("y_bio", lg)) +
      theme_bw(base_size = 14)
  })

  output$data_table <- renderTable(result()$dat[, required_cols], digits = 2)

  output$template <- downloadHandler(
    filename = "bjm_patient_template.csv",
    content = function(file) write.csv(template, file, row.names = FALSE))
}

shinyApp(ui, server)
