library(tradepolicy)
library(fixest)
library(dplyr)
library(knitr)
library(broom)

trade <- tradepolicy::agtpa_applications %>% 
  rename(fta = rta) %>% 
  mutate_if(is.character, tolower)

# step 1 ----

## stage 1 ----

trade2 <- trade %>%
  select(exporter, importer, pair_id, year, trade, dist, cntg, lang, clny, fta) %>%
  filter(year %in% seq(1990, 2006, 4)) %>%
  mutate(
    log_dist = log(dist),
    intl = ifelse(exporter != importer, 1, 0),
    exporter = ifelse(exporter == "deu", "0-deu", exporter),
    importer = ifelse(importer == "deu", "0-deu", importer)
  ) %>%
  
  # Create Yit
  group_by(exporter, year) %>%
  mutate(y = sum(trade)) %>%
  
  # Create Eit
  group_by(importer, year) %>%
  mutate(e = sum(trade)) %>%
  
  # Create Er
  group_by(year) %>%
  mutate(e_r = max(ifelse(importer == "0-deu", e, NA), na.rm = T)) %>% 
  arrange(importer)

trade2 <- trade2 %>%
  mutate(
    exp_year = paste0(exporter, year),
    imp_year = paste0(importer, year),
    pair_id_2 = ifelse(exporter == importer, "0-intra", as.character(pair_id))
  )

trade2 <- trade2 %>%
  group_by(pair_id) %>%
  mutate(sum_trade = sum(trade, na.rm = T)) %>%
  ungroup()

fit_baseline <- fepois(
  trade ~ fta | exp_year + imp_year + pair_id_2,
  data = filter(trade2, sum_trade > 0),
  glm.iter = 500
)

fit_baseline

options(scipen = 999)

# kable(tidy(fit_baseline) %>% 
#         mutate_if(is.numeric, function(x) round(x,5)),
#       format = "latex", booktabs = T, longtable = T,
#       caption = "Stage I model summary, part I")
# 
# kable(glance(fit_baseline) %>% 
#         select(pseudo.r.squared, nobs:logLik) %>% 
#         mutate_if(is.numeric, function(x) round(x,5)),
#       format = "latex", booktabs = T, longtable = T,
#       caption = "Stage I model summary, part II")

(exp(fit_baseline$coefficients) - 1) * 100

trade2 <- trade2 %>%
  mutate(
    fe_exporter_bln = fixef(fit_baseline)$exp_year[exp_year],
    fe_importer_bln = fixef(fit_baseline)$imp_year[imp_year],
    fe_pair_id_2_bln = fixef(fit_baseline)$pair_id_2[pair_id_2]
  )

## stage 2 ----

trade2 <- trade2 %>%
  mutate(
    tij_bar = exp(fe_pair_id_2_bln),
    tij_bln = exp(fe_pair_id_2_bln + fit_baseline$coefficients["fta"] * fta)
  )

trade2 %>% 
  mutate(tij_bar_na = ifelse(is.na(tij_bar), 1L, 0L)) %>% 
  group_by(tij_bar_na) %>% 
  count()

quantile(trade2$tij_bar, na.rm = T)

trade2_2006 <- trade2 %>%
  filter(year %in% c(2006), exporter != importer)

fit_costs <- fepois(
  tij_bar ~ log_dist + cntg + lang + clny | exporter + importer,
  data = trade2_2006,
  glm.iter = 500
)

fit_costs

etable(fit_costs, tex = T)

trade2_2006 <- trade2_2006 %>%
  mutate(tij_no_fta = predict(fit_costs, trade2_2006)) %>%
  select(exporter, importer, tij_no_fta)

trade2 <- trade2 %>%
  filter(year %in% c(2006)) %>%
  left_join(trade2_2006, by = c("exporter", "importer")) %>%
  mutate(
    tij_bar = ifelse(is.na(tij_bar), tij_no_fta, tij_bar),
    tij_bln = ifelse(is.na(tij_bln), tij_bar * exp(fit_baseline$coefficients["fta"] * fta), tij_bln)
  ) %>%
  select(-tij_no_fta)

fit_constrained <- fepois(
  trade ~ 0 | exporter + importer,
  data = trade2,
  offset = ~log(tij_bln),
  glm.iter = 500
)

trade2 <- trade2 %>%
  mutate(tradehat_bln = predict(fit_constrained, trade2)) %>%
  group_by(exporter) %>%
  mutate(xi_bln = sum(tradehat_bln * (exporter != importer))) %>%
  ungroup()

trade2 <- trade2 %>%
  mutate(
    fe_exporter_cns = fixef(fit_constrained)$exporter[exporter],
    fe_importer_cns = fixef(fit_constrained)$importer[importer]
  )

trade2 <- trade2 %>%
  mutate(
    omr_bln = y * e_r/ exp(fe_exporter_cns),
    imr_bln = e / (exp(fe_importer_cns) * e_r)
  )

# kable(tidy(fit_constrained) %>% 
#         mutate_if(is.numeric, function(x) round(x,5)),
#       format = "latex", booktabs = T, longtable = T,
#       caption = "Stage I model summary, part I")
# 
# kable(glance(fit_baseline) %>% 
#         select(pseudo.r.squared, nobs:logLik) %>% 
#         mutate_if(is.numeric, function(x) round(x,5)),
#       format = "latex", booktabs = T, longtable = T,
#       caption = "Stage I model summary, part II")

etable(fit_baseline, tex = T)

# step 2 ----

ftas <- c("chl", "chn", "usa")

trade2 <- trade2 %>%
  mutate(
    fta_no_ccu = ifelse(exporter %in% ftas & importer %in% ftas, 0L, fta),
    tij_cfl = tij_bar * exp(fit_baseline$coefficients["fta"] * fta_no_ccu)
  )

# step 3 ----

fit_counterfactual <- fepois(
  trade ~ 0 | exporter + importer,
  data = trade2,
  offset = ~log(tij_cfl),
  glm.iter = 500
)

trade2 <- trade2 %>%
  mutate(
    fe_exporter_cfl = fixef(fit_counterfactual)$exporter[exporter],
    fe_importer_cfl = fixef(fit_counterfactual)$importer[importer]
  )

trade2 <- trade2 %>%
  mutate(
    omr_cfl = y * e_r / exp(fe_exporter_cfl),
    imr_cfl = e / (exp(fe_importer_cfl) * e_r)
  )

trade2 <- trade2 %>%
  mutate(tradehat_cfl = predict(fit_counterfactual, trade2)) %>%
  group_by(exporter) %>%
  mutate(xi_cfl = sum(tradehat_cfl * (exporter != importer))) %>%
  ungroup()

sigma <- 7

trade2 <- trade2 %>%
  mutate(
    change_tij = tij_cfl / tij_bln,
    phi = ifelse(importer == exporter, e / y, 0)
  ) %>%
  group_by(exporter) %>%
  mutate(phi = max(phi)) %>%
  ungroup()

trade2 <- trade2 %>%
  group_by(exporter) %>%
  mutate(change_p_i = ((exp(fe_exporter_cfl) / e_r) / (exp(fe_exporter_cns) / e_r))^(1 /(1 - sigma))) %>%
  ungroup() %>%
  
  group_by(importer) %>%
  mutate(
    change_p_j = ifelse(importer == exporter, change_p_i, 0),
    change_p_j = max(change_p_j)
  ) %>%
  ungroup()

trade2 <- trade2 %>%
  mutate(trade_cfl = tradehat_cfl * change_p_i * change_p_j)

trade2 <- trade2 %>%
  mutate(
    omr_cfl_0 = omr_cfl,
    imr_cfl_0 = imr_cfl,
    change_imr_full_0 = 1,
    change_omr_full_0 = 1,
    change_p_i_0 = change_p_i,
    change_p_j_0 = change_p_j,
    fe_exporter_cfl_0 = fe_exporter_cfl,
    fe_importer_cfl_0 = fe_importer_cfl,
    tradehat_0 = tradehat_cfl,
    e_r_cfl_0 = e_r
  )

# set parameters
max_dif <- 1
sd_dif <- 1
change_price_i_old <- 0

i2 <- 1
while(sd_dif > 1e-3 | max_dif > 1e-3) {
  print(i2)
  trade2 <- trade2 %>%
    # mutate(trade_1 = tradehat_0 * change_p_i_0 * change_p_j_0 / (change_omr_full_0 * change_imr_full_0))
    mutate(trade_1 = tradehat_0 * change_p_i_0 * change_p_j_0 / (change_omr_full_0 * change_imr_full_0))
  
  # repeat the counterfactual model
  fit_counterfactual_app2_2 <- fepois(
    trade_1 ~ 0 | exporter + importer,
    data = trade2,
    offset = ~log(tij_cfl),
    glm.iter = 500
  )
  
  trade2 <- trade2 %>%
    mutate(
      fe_exporter_cfl = fixef(fit_counterfactual_app2_2)$exporter[exporter],
      fe_importer_cfl = fixef(fit_counterfactual_app2_2)$importer[importer]
    )
  
  # compute the conditional general equilibrium effects of trade
  trade2 <- trade2 %>%
    mutate(tradehat_1 = predict(fit_counterfactual_app2_2, trade2)) %>%
    group_by(exporter) %>%
    mutate(y_cfl_1 = sum(tradehat_1)) %>%
    ungroup() %>%
    
    mutate(e_cfl_1 = ifelse(importer == exporter, phi * y_cfl_1, 0)) %>%
    group_by(importer) %>%
    mutate(e_cfl_1 = max(e_cfl_1)) %>%
    ungroup() %>%
    
    mutate(
      e_r_cfl_1 = ifelse(importer == "0-deu", e_cfl_1, 0),
      e_r_cfl_1 = max(e_r_cfl_1)
    )
  
  # compute the change in prices for exporters and importers
  trade2 <- trade2 %>%
    mutate(change_p_i_1 = ((exp(fe_exporter_cfl) / e_r_cfl_1) /
                             (exp(fe_exporter_cfl_0) / e_r_cfl_0))^(1 / (1 - sigma)))
  
  # compute the change in prices for exporters and importers
  trade2 <- trade2 %>%
    group_by(importer) %>%
    mutate(
      change_p_j_1 = ifelse(importer == exporter, change_p_i_1, 0),
      change_p_j_1 = max(change_p_j_1)
    ) %>%
    ungroup()
  
  # compute both outward and inward multilateral resistance
  trade2 <- trade2 %>%
    mutate(
      omr_cfl_1 = (y_cfl_1 * e_r_cfl_1) / exp(fe_exporter_cfl),
      imr_cfl_1 = e_cfl_1 / (exp(fe_importer_cfl) * e_r_cfl_1)
    )
  
  # update the differences
  max_dif <- abs(max(trade2$change_p_i_0 - change_price_i_old))
  sd_dif <- sd(trade2$change_p_i_0 - change_price_i_old)
  change_price_i_old <- trade2$change_p_i_0
  
  # compute changes in outward and inward multilateral resistance
  trade2 <- trade2 %>%
    mutate(
      change_omr_full_1 = omr_cfl_1 / omr_cfl_0,
      change_imr_full_1 = imr_cfl_1 / imr_cfl_0,
      omr_cfl_0 = omr_cfl_1,
      imr_cfl_0 = imr_cfl_1,
      change_omr_full_0 = change_omr_full_1,
      change_imr_full_0 = change_imr_full_1,
      change_p_i_0 = change_p_i_1,
      change_p_j_0 = change_p_j_1,
      fe_exporter_cfl_0 = fe_exporter_cfl,
      fe_importer_cfl_0 = fe_importer_cfl,
      tradehat_0 = tradehat_1,
      e_r_cfl_0 = e_r_cfl_1
    ) %>%
    select(-fe_exporter_cfl, -fe_importer_cfl)
  
  i2 <- i2 + 1
}

trade2 <- trade2 %>%
  mutate(
    change_p_i_full = ((exp(fe_exporter_cfl_0) / e_r_cfl_0) /
                         (exp(fe_exporter_cns) / e_r))^(1 / (1 - sigma)),
    change_p_j_full = change_p_i_full * (exporter == importer)
  ) %>%
  group_by(importer) %>%
  mutate(change_p_j_full = max(change_p_j_full)) %>%
  ungroup() %>%
  mutate(y_full = change_p_i_full * y)

trade2 <- trade2 %>%
  mutate(e_full = change_p_j_full * e * (exporter == importer)) %>%
  group_by(importer) %>%
  mutate(e_full = max(e_full, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    e_full_r = e_full * (importer == "0-DEU"),
    e_full_r = max(e_full_r)
  )

trade2 <- trade2 %>%
  mutate(
    omr_full = y_full * e_r_cfl_0 / exp(fe_exporter_cfl_0),
    imr_full = e_cfl_1 / (exp(fe_importer_cfl_0) * e_r_cfl_0)
  )

trade2 <- trade2 %>%
  mutate(x_full = (y_full * e_full * tij_cfl) / (imr_full * omr_full)) %>%
  group_by(exporter) %>%
  mutate(xi_full = sum(x_full * (importer != exporter))) %>%
  ungroup()

exporter_indexes <- trade2 %>%
  select(
    exporter, starts_with("omr_"), change_p_i_full,
    starts_with("xi_"), y, y_full
  ) %>%
  distinct() %>%
  mutate(exporter = ifelse(exporter == "0-DEU", "DEU", exporter)) %>%
  arrange(exporter) %>%
  mutate(
    change_p_i_full = (1 - change_p_i_full) * 100,
    change_omr_cfl = ((omr_bln / omr_cfl)^(1 / (1-sigma)) - 1) * 100,
    change_omr_full = ((omr_bln / omr_full)^(1 / (1-sigma)) - 1) * 100,
    change_xi_cfl = (xi_bln / xi_cfl - 1) * 100,
    change_xi_full = (xi_bln / xi_full - 1) * 100
  ) %>%
  select(exporter, starts_with("change"), starts_with("y"))

importer_indexes <- trade2 %>%
  select(importer, imr_bln, imr_cfl, imr_full) %>%
  distinct() %>%
  mutate(importer = ifelse(importer == "0-DEU", "DEU", importer)) %>%
  arrange(importer) %>%
  mutate(
    change_imr_cfl = ((imr_bln / imr_cfl)^(1 / (1 - sigma)) - 1) * 100,
    change_imr_full = ((imr_bln / imr_full)^(1 / (1 - sigma)) - 1) * 100
  )

indexes_final <- exporter_indexes %>%
  left_join(importer_indexes, by = c("exporter" = "importer")) %>%
  mutate(
    rgdp_bln = y / (imr_bln^(1 / (1 - sigma))),
    rgdp_full = y_full / (imr_full^(1 / (1 - sigma))),
    change_rgdp_full = (rgdp_bln / rgdp_full - 1) * 100
  ) %>%
  select(exporter, change_xi_cfl, change_xi_full,
         change_rgdp_full, change_imr_full, change_omr_full, change_p_i_full)

indexes_final <- indexes_final %>%
  mutate_if(is.numeric, function(x) round(x, 2))

indexes_final

kable(indexes_final, format = "latex",
      caption = "Conditional GE and Full Endownment GE simulation results",
      booktabs = T,
      longtable = T)
