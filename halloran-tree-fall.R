# Mangrove surface elevation loss after tree fall during extreme weather 
# Vicki Bennion, Jack W Hill, Catherine E Lovelock

##### 0. setup workspace #####

library(betareg)
library(ggforce)
library(gratia)
library(lubridate)
library(mgcv)
library(ragg)
library(readxl)
library(tidyverse) 

# custom ggplot theme
theme_tree_fall <- function() {
  theme_classic(base_size = 12, base_family = "sans") %+replace%
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(size = rel(0.75)),
      axis.text.x = element_text(margin = margin(2, 0, 3, 0)),
      axis.text.y = element_text(margin = margin(0, 2, 0, 1)),
      axis.ticks = element_line(colour = "black"),
      axis.line = element_line(colour = "black"),
      legend.position = "none"
    )
}

#
##### 1. import data #####

halloran_raw <- read_excel("data/halloran-tidied.xlsx", 
                           sheet = "cumulative")

halloran_all <- halloran_raw %>% 
  pivot_longer(cols = -c(habitat, site),
               names_to = "date",
               values_to = "elev") %>% 
  mutate(across(c(habitat, site),
                tolower),
         # convert Excel dates to calendar dates
         date = as.Date(as.numeric(date), 
                        origin = "1899-12-30")) %>% 
  filter(!is.na(elev))

# drop pin sites to keep only RSET data
halloran_rsets <- halloran_all %>% 
  filter(str_detect(site, "pin", negate = TRUE))

storm_markers <- tibble(date = as.Date(c("20/03/2021", "01/02/2022"), 
                                       format = "%d/%m/%Y"),
                        y = c(27, 22),
                        yend = c(22, 17))

root_bags_raw <- read_excel("data/halloran-tidied.xlsx", 
                            sheet = "root-bags") 

root_bags <- root_bags_raw  %>% 
  group_by(date, site) %>% 
  summarise(mean_root_prop = mean(root_prop),
            se_root_prop = sd(root_prop)/sqrt(n())) %>% 
  ungroup() %>% 
  mutate(site = fct_relevel(as_factor(site),
                            "TD",
                            "M3",
                            "M2",
                            "M1"))

halloran_pins <- halloran_all %>% 
  filter(str_detect(site, "pin"))

# td = 'tree downed', i.e. the pins near fallen tree
td_pins <- halloran_pins %>%
  filter(site == "tdpins") %>% 
  mutate(trunk_side = str_split_i(habitat, "_", 1),
         depth = str_split_i(habitat, "_", 2))

interval_dat_raw <- read_excel("data/halloran-tidied.xlsx", 
                               sheet = "interval")

interval_dat <- interval_dat_raw %>% 
  pivot_longer(cols = -c(type, site),
               names_to = "date",
               values_to = "elev_interval") %>% 
  mutate(across(c(type, site),
                tolower),
         date = as.Date(as.numeric(date), 
                        origin = "1899-12-30"))  

elevation_dat <- interval_dat %>%
  filter(!is.na(elev_interval),
         type != "root",
         between(date, as.Date("2020-01-03"), 
                 as.Date("2022-12-10"))) %>% 
  group_by(site, type) %>% 
  # should already be arranged by date
  arrange(date, .by_group = TRUE) %>% 
  mutate(elev = cumsum(elev_interval)) %>% 
  ungroup() %>% 
  mutate(mang_or_salt = case_when(str_detect(site, "mang") ~ "mang",
                                  str_detect(site, "salt") ~ "salt",
                                  TRUE ~ "downed"))

#
##### 2. analyse -- RSETs #####

# create relative decimal dates -- although downstream code should produce the 
# same output with decimal dates created directly from calendar dates
halloran_gam_dat <- halloran_rsets %>% 
  mutate(year = year(date),
         dec_year = case_when(year == 2013 ~ 0,
                              year == 2014 ~ 1,
                              year == 2015 ~ 2,
                              year == 2016 ~ 3,
                              year == 2017 ~ 4,
                              year == 2018 ~ 5,
                              year == 2019 ~ 6,
                              year == 2020 ~ 7,
                              year == 2021 ~ 8,
                              year == 2022 ~ 9,
                              year == 2023 ~ 10),
         dec_date_real_year = decimal_date(date),
         dec_date = dec_year + (dec_date_real_year - year),
         site = as_factor(site))

halloran_gam <- gam(elev ~ s(dec_date, by = site, k = 9) + site,
                    method = "REML", data = halloran_gam_dat)
summary(halloran_gam)
# result: GAM edfs

gam_dat_raw <- smooth_estimates(halloran_gam) %>% 
  add_confint() # default is a 95% CI

gam_dat <- gam_dat_raw %>% 
  # convert decimal date back to calendar date
  mutate(decimal_part = dec_date %% 1,
         whole_part = dec_date - decimal_part,
         year = case_when(whole_part == 0 ~ 2013,
                          whole_part == 1 ~ 2014,
                          whole_part == 2 ~ 2015,
                          whole_part == 3 ~ 2016,
                          whole_part == 4 ~ 2017,
                          whole_part == 5 ~ 2018,
                          whole_part == 6 ~ 2019,
                          whole_part == 7 ~ 2020,
                          whole_part == 8 ~ 2021,
                          whole_part == 9 ~ 2022,
                          whole_part == 10 ~ 2023),
         month = month(date_decimal(decimal_part)),
         day = day(date_decimal(decimal_part)),
         text_date = paste0(year, "/",
                            month, "/",
                            day),
         date = ymd(text_date),
         # correct y intercept of each RSET site GAM fit using model coefficients
         across(c(est, lower_ci, upper_ci),
                ~ case_when(
                  str_detect(smooth, "mang1") ~ .x + coef(halloran_gam)[[1]],
                  str_detect(smooth, "mang2") ~ .x + coef(halloran_gam)[[1]] + coef(halloran_gam)[[2]],
                  str_detect(smooth, "mang3") ~ .x + coef(halloran_gam)[[1]] + coef(halloran_gam)[[3]],
                  str_detect(smooth, "salt") ~ .x + coef(halloran_gam)[[1]] + coef(halloran_gam)[[4]],
                  .default = .x)))

# GAM derivatives after Feher et al. (2023)
gam_deriv_raw <- derivatives(halloran_gam,
                             n = 500)

gam_deriv_raw %>% 
  group_by(site) %>% 
  summarise(mean(derivative),
            sd(derivative)/sqrt(n()))
# result: surface elevation change over full monitoring period 

# decimal date used to split GAM derivatives into before and after tree fall
dec_date_after_dec2020_measurement <- decimal_date(as_date("2020/12/05")) - 2020 + 7 
# +7 because it's the seventh year in dataset, as per relative decimal dates;
# see decimal date code and comments above

gam_deriv_raw %>% 
  mutate(event = case_when(data < dec_date_after_dec2020_measurement ~
                             "before",
                           data > dec_date_after_dec2020_measurement ~
                             "after")) %>% 
  group_by(site, event) %>% 
  summarise(mean(derivative),
            sd(derivative)/sqrt(n()),
            n()) %>% 
  arrange(event)
# result: surface elevation change before and after tree fall

#
##### 3. plot -- RSETs #####

plot_rset_gams <- ggplot() +
  geom_line(data = halloran_rsets, 
            aes(x = date,
                y = elev,
                colour = site),
            size = 0.4,
            linetype = "dashed") +
  geom_ribbon(data = gam_dat,
              aes(x = date,
                  ymin = lower_ci,
                  ymax = upper_ci,
                  fill = site),
              alpha = 0.2,
              colour = NA) +
  geom_line(data = gam_dat,
            aes(x = date,
                y = est,
                colour = site),
            size = 1) +
  geom_segment(data = storm_markers,
               mapping = aes(x = date, xend = date,
                             y = y, yend = yend),
               colour = "black",
               size = rel(1.1),
               lineend = "butt", linejoin = "mitre",
               arrow = arrow(angle = 30,
                             length = unit(0.03,
                                           units = "npc"),
                             type = "closed")) +
  scale_fill_manual(values = c(mang1 = "mediumorchid2",
                               mang2 = "olivedrab3",
                               mang3 = "dodgerblue2",
                               salt = "coral2")) +
  scale_color_manual(values = c(mang1 = "mediumorchid2",
                                mang2 = "olivedrab3",
                                mang3 = "dodgerblue2",
                                salt = "coral2")) +
  scale_y_continuous(name = "Cumulative surface elevation change (mm)") +
  scale_x_date(name = "Year",
               date_labels = "%Y",
               breaks = as.Date(c("01/01/2014", "01/01/2016", "01/01/2018",
                                  "01/01/2020", "01/01/2021", "01/01/2022"), 
                                format = "%d/%m/%Y"),
               limits = as.Date(c('30/11/2013', NA), format = "%d/%m/%Y"),
               expand = expansion()) +
  facet_zoom(xlim = c(as.Date('01/03/2020', format = "%d/%m/%Y"),
                      as.Date('30/09/2022', format = "%d/%m/%Y")),
             zoom.size = 1) +
  theme_tree_fall() %+replace%
  theme(panel.border = element_rect(fill = NA, colour = "black"), 
        strip.background = element_rect(fill = "grey85", 
                                        colour = "black",
                                        linewidth = rel(1.25)))

agg_png("fig-output/fig-2_left-panel.png",
        res = 1080, scaling = 1.25, units = "cm",
        width = 15, height = 18)
plot_rset_gams
dev.off()

#
##### 4. analyse -- root bags #####

summary(betareg(root_prop ~ site * date, 
                data = root_bags_raw %>% filter(root_prop != 0)))
# result: root growth over time and across sites

# test for elevation change in each root bag deployment period using pin data.
# the 'raw data' are the four measurements at each time point -- 
# seaward (deep and shallow) and landward (deep and shallow) -- which are each
# the average of four pin measurements

# interval change in period 25/03/22 - 15/07/22
pins_to_july <- c(-4.40, -4.25, -7.50, -12.75) 

# interval change in period 16/07/22 - 10/12/22
pins_to_dec <- c(7.20, 10.58, 5.50, 10.50)

t.test(pins_to_july, pins_to_dec)
# result: root zone soil thickness change across time

#
##### 5. plot -- root bags #####

plot_root_bags <- ggplot(data = root_bags,
                         aes(x = date,
                             y = mean_root_prop,
                             fill = site)) +
  geom_col(colour = "black",
           size = 0.2,
           position = position_dodge(0.9),
           width = 0.75) +
  geom_errorbar(aes(ymin = mean_root_prop - se_root_prop,
                    ymax = mean_root_prop + se_root_prop), 
                width = 0, 
                position = position_dodge(0.9)) +
  scale_fill_manual(values = c(M1 = "mediumorchid2",
                               M2 = "olivedrab3",
                               M3 = "dodgerblue2",
                               TD = "grey70")) +
  scale_y_continuous(name = "Root growth (g of dry roots per g of soil)",
                     labels = scales::label_number(),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(name = "Root bag deployment period",
                   labels = c("Dec 2021 – Jun 2022",
                              "Jun 2022 – Dec 2022")) +
  theme_tree_fall() %+replace%
  theme(axis.ticks.x = element_blank())

agg_png("fig-output/fig-3.png",
        res = 1080, scaling = 0.75, units = "cm",
        width = 10, height = 8)
plot_root_bags
dev.off()

#
##### 6. analyse -- fallen tree sediment pins #####

td_pins_gam_dat <- td_pins %>% 
  mutate(year = year(date),
         dec_year = case_when(year == 2021 ~ 1,
                              year == 2022 ~ 2),
         dec_date_real_year = decimal_date(date),
         dec_date = dec_year + (dec_date_real_year - year),
         across(c(trunk_side, depth),
                ~ as_factor(.)))

td_pin_gam <- gam(elev ~ s(dec_date, by = trunk_side, k = 7) + trunk_side + depth,
                  method = "REML",
                  data = td_pins_gam_dat)
summary(td_pin_gam)
# result: landward and seaward root zone soil thickness change around fallen tree

#
##### 7. plot -- fallen tree sediment pins #####

plot_td_pins <- ggplot(td_pins,
                       aes(x = date, 
                           y = elev,
                           colour = trunk_side,
                           group = interaction(trunk_side, depth))) +
  geom_line(linetype = "dashed",
            size = 0.4) +
  geom_smooth(aes(group = trunk_side),
              method = "gam",
              formula = y ~ s(x, k = 7),
              se = FALSE) +
  scale_colour_manual(values = c("#ED7D31",
                                 "turquoise3")) +
  annotate("text", label = "landward",
           x = as.Date("07/12/2021",
                       format = "%d/%m/%Y"),
           y = 9.1,
           size = 3, lineheight = 0.8,
           hjust = 1, vjust = 0,
           colour = "#ED7D31", fontface = "bold") +
  annotate("text", label = "seaward",
           x = as.Date("28/01/2022",
                       format = "%d/%m/%Y"),
           y = -1.8,
           size = 3, lineheight = 0.75,
           hjust = 1, vjust = 0,
           colour = "turquoise3", fontface = "bold") +
  annotate("rect", 
           xmin = as.Date("24/12/2021",
                          format = "%d/%m/%Y"),
           xmax = as.Date("13/07/2022",
                          format = "%d/%m/%Y"),
           ymin = -4.9, ymax = -5.7,
           colour = NA, fill = "grey90") +
  annotate("text", label = "deployment one",
           x = as.Date("02/04/2022",
                       format = "%d/%m/%Y"),
           y = -5.45,
           size = 2.3,
           hjust = 0.5, vjust = 0,
           colour = "grey20", fontface = "bold") +
  annotate("rect", 
           xmin = as.Date("16/07/2022",
                          format = "%d/%m/%Y"),
           xmax = as.Date("10/12/2022",
                          format = "%d/%m/%Y"),
           ymin = -4.9, ymax = -5.7,
           colour = NA, fill = "grey90") +
  annotate("text", label = "deployment two",
           x = as.Date("28/09/2022",
                       format = "%d/%m/%Y"),
           y = -5.45,
           size = 2.3,
           hjust = 0.5, vjust = 0,
           colour = "grey20", fontface = "bold") +
  scale_y_continuous(name = "Change in root zone soil thickness (cumulative mm)",
                     breaks = c(-4, 0, 4, 8, 12),
                     expand = expansion(add = c(0, 0.5))) +
  scale_x_date(name = "Date",
               date_labels = "%b %Y",
               breaks = as.Date(c("01/05/2021", "01/11/2021", 
                                  "01/05/2022", "01/11/2022"),
                                format = "%d/%m/%Y"),
               expand = expansion()) +
  labs(tag = "b") +
  theme_tree_fall() %+replace%
  theme(plot.tag = element_text(face = "bold",
                                colour = "grey70"),
        plot.tag.position = c(0.14, 0.99))

#
##### 8. plot -- fallen tree pin and root bag correlation #####

root_corro_dat <- root_bags %>% 
  filter(site == "TD") %>% 
  select(-site) %>% 
  rename(period = date) %>% 
  mutate(group = c("time1", "time2"))

pin_corro_dat <- tibble(period = rep(c("24/12/22-15/07/22", "16/07/22-10/12/22"),
                                     each = 4),
                        group = rep(c("time1", "time2"),
                                    each = 4),
                        pin_change = c(pins_to_july, pins_to_dec)) %>% 
  group_by(group, period) %>% 
  summarise(mean_pin_change = mean(pin_change),
            se_pin_change = sd(pin_change)/sqrt(n())) %>% 
  ungroup()

full_corro_dat <- left_join(root_corro_dat, pin_corro_dat,
                            by = c("group")) %>% 
  select(-contains("period"))

plot_corro <- ggplot(data = full_corro_dat,
                     aes(x = mean_root_prop,
                         y = mean_pin_change,
                         ymin = mean_pin_change - se_pin_change, 
                         ymax = mean_pin_change + se_pin_change,
                         xmin = mean_root_prop - se_root_prop,
                         xmax = mean_root_prop + se_root_prop)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             colour = "grey70") +
  geom_errorbar(width = 0) +
  geom_errorbarh(height = 0) +
  geom_point() +
  annotate("text", label = "deployment\none",
           x = 0.00001,
           y = -6.63,
           size = 3, lineheight = 0.8,
           hjust = 0, vjust = 1,
           colour = "grey40", fontface = "bold") +
  annotate("text", label = "deployment\ntwo",
           x = 0.000156,
           y = 7.75,
           size = 3, lineheight = 0.8,
           hjust = 0, vjust = 1, 
           colour = "grey40", fontface = "bold") +
  scale_x_continuous(name = "Root growth (g of dry roots per g of soil)",
                     labels = scales::label_number(),
                     limits = c(0, NA), 
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(name = "Change in root zone soil thickness (mm)",
                     limits = c(-10, 10),
                     # breaks = c(-4, 0, 4, 8, 12),
                     expand = expansion(add = c(0.5, 0.5))) +
  labs(tag = "a") +
  theme_tree_fall() %+replace%
  theme(plot.tag = element_text(face = "bold",
                                colour = "grey70"),
        plot.tag.position = c(0.14, 0.99))

agg_png("fig-output/fig-4.png",
        res = 1080, scaling = 0.75, units = "cm",
        width = 15, height = 8)
plot_corro + plot_td_pins
dev.off()

#
##### 9. plot -- RSET sediment pins #####

plot_pins <- ggplot(elevation_dat %>% 
                      filter(type == "pinsg",
                             mang_or_salt == "mang"),
                    aes(x = date, 
                        y = elev,
                        colour = site)) +
  geom_line(size = 0.4) +
  annotate("text", label = "furthest\nRSET    ",
           x = as.Date("30/08/2021",
                       format = "%d/%m/%Y"),
           y = 4,
           size = 3, lineheight = 0.8,
           hjust = 1, vjust = 0,
           colour = "mediumorchid2", fontface = "bold") +
  annotate("text", label = "nearest RSET",
           x = as.Date("20/10/2021",
                       format = "%d/%m/%Y"),
           y = -2.4, angle = 3,
           size = 3, lineheight = 0.75,
           hjust = 1, vjust = 0,
           colour = "dodgerblue2", fontface = "bold") +
  annotate("text", label = "middle RSET",
           x = as.Date("20/10/2021",
                       format = "%d/%m/%Y"),
           y = -4.6, angle = -3,
           size = 3, lineheight = 0.75,
           hjust = 1, vjust = 0, 
           colour = "olivedrab3", fontface = "bold") +
  scale_color_manual(values = c(mang1 = "mediumorchid2",
                                mang2 = "olivedrab3",
                                mang3 = "dodgerblue2")) +
  scale_y_continuous(name = "Change in root zone soil thickness (cumulative mm)",
                     breaks = c(-4, 0, 4, 8, 12),
                     expand = expansion(add = c(0.5, 0.5))) +
  scale_x_date(name = "Date",
               date_labels = "%b %Y",
               limits = c(as.Date("07/04/2021", format = "%d/%m/%Y"),
                          NA),
               breaks = as.Date(c("01/05/2021", "01/11/2021", 
                                  "01/05/2022", "01/11/2022"),
                                format = "%d/%m/%Y"),
               expand = expansion()) +
  theme_tree_fall() 

agg_png("fig-output/fig-s3.png",
        res = 1080, scaling = 0.75, units = "cm",
        width = 10, height = 8)
plot_pins
dev.off()

#