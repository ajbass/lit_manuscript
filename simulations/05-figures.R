###
# 05-figures.R creates figures for simulation results
###

# plotting function from gap::qqunif package (minor changes)
QQ_setup <- function(p) {
  n <- length(p)
  c <- abs(qnorm(0.05/2))
  m <- (1:n)/(n + 1)
  v <- (1:n) * (n - (1:n) + 1)/(n + 1)^2/(n + 2)
  s <- sqrt(v)
  lcl <- m - c * s
  ucl <- m + c * s
  lid <- (lcl > 0)
  uid <- (ucl <= 1)
  id <- lid & uid
  p <- p[id]
  lcl <- lcl[id]
  ucl <- ucl[id]
  df <- data.frame(observed = -log10(sort(p)),
                   theoretical = -log10(1:length(p) / (length(p) + 1)),
                   lower = -log10(lcl),
                   upper = -log10(ucl))
  df
  }

# load plotting functions
library(tidyverse)
library(cowplot)

########### 01-lit-null.R plots ###########
df = NULL
fnames <- list.files("./data/null/")
for (i in fnames) {
  tmp <- readRDS(paste0("./data/null/", i))
  tmp <- tmp %>%
    group_by(method, error, correlation, N, num_traits, seed) %>%
    mutate(error = as.character(error)) %>%
    summarise(threshold1 = mean(p.value < 1e-3))
  df <- dplyr::bind_rows(tmp, df)
}
df$error = factor(df$error, labels =  c("Chi-squared", "Gaussian","t"))

# Plotting functions
cbbPalette <- c( "#009E73", "#0072B2", "#D55E00", "#CC79A7")
df$Method <- factor(df$method,
                      labels = c("aLIT", "wLIT", "uLIT"))

expSup <- function(w, digits=1) {
  sprintf(paste0('"%.1f"',"*x*10^%d"), w/10^floor(log10(abs(w))), floor(log10(abs(w))))
}
p1 <- df %>%
  ggplot(aes(x = Method, y = threshold1, fill = as.factor(correlation))) +
  geom_boxplot(size = .53, alpha = 0.6,outlier.colour = NA) +
  theme_bw(base_size = 11) +
  xlab("") +
  scale_fill_manual(values = cbbPalette, name = "Correlation")  +
  facet_grid(num_traits~error) +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_abline(intercept = 0.001, slope = 0, linetype = "dashed", color = "firebrick3",  size = .7) +
  ylab("False positive rate") +
  scale_y_continuous(breaks=c(0.000,0.0005, 0.001, 0.0015, 0.002), labels=c("0", parse(text=c(expSup(c(0.0005 , 0.001, 0.0015, 0.002)))))) +
  theme(plot.title = element_text(size=12, vjust = -1.7), axis.text.x = element_text(color = "black", size = 11.5))

ggsave(p1,
       width = 7,
       height = 3.5,
       filename = "./figures/simulations-type-i-error.png")

# QQ-plot of all p-values
df = NULL
fnames <- list.files("./data/null/")
for (i in fnames) {
  tmp <- readRDS(paste0("./data/null/", i))
  tmp <- tmp %>%
    group_by(method, error, correlation, N, num_traits, seed) %>%
    mutate(error = as.character(error))
  df <- dplyr::bind_rows(tmp, df)
}

cbbPalette <- c( "#009E73", "#0072B2", "#D55E00", "#CC79A7")
df$Method <- factor(df$method,
                    labels = c("aLIT", "wLIT", "uLIT"))

df_plot <- df %>%
  group_by(N,
           Method,
           error,
           correlation,
           num_traits) %>%
  do(QQ_setup(.$p.value))

dat2 <- df_plot %>%
  ungroup() %>%
  select(theoretical, upper, lower) %>%
  distinct()

expSup <- function(w, digits=1) {
  sprintf(paste0('"%.1f"',"*x*10^%d"), w/10^floor(log10(abs(w))), floor(log10(abs(w))))
}

df_plot$error = factor(df_plot$error, labels =  c("Chi-squared", "Gaussian","t"))
df_plot$correlation = factor(df_plot$correlation, labels =  c("0.25", "0.50", "0.75"))
p1 <- df_plot %>%
  ggplot(aes(x = theoretical, y = observed, shape = Method, color = correlation)) +
  geom_point( size = .5, alpha = 0.75) +
  geom_abline(slope = 1)  +
  geom_line(aes(x = theoretical, y = upper),  color = "black", linetype = "dashed") +
  geom_line(aes(x = theoretical, y = lower),color = "black", linetype = "dashed") +
  scale_color_manual(values = cbbPalette, name = "Correlation") +
  theme_bw(base_size = 12) +
  coord_fixed() +
  facet_grid(num_traits ~ error) +
  xlab(expression(-log[10](theoretical)))+
  ylab(expression(-log[10](observed)))+
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())+
  theme(plot.title = element_text(size=12, vjust = -1.7), axis.text.x = element_text(color = "black", size = 11.5))

ggsave(p1,
       filename = "./figures/simulations-type-1-error-qq.png")
###

########### 02-lit-alt.R plots ###########
# positive pleiotropy
out <- readRDS("./data/02-lit-alt-5-traits.rds")
out2 <- readRDS("./data/02-lit-alt-10-traits.rds")
out <- rbind(out, out2)

out %>% filter(positive == TRUE,
               counteract == FALSE,
               method %in% c("aLIT", "Marginal_SQ", "Marginal")) %>%
  filter(prop_null == 0, heritability == 0.1) %>%
  group_by(num_traits, method) %>%
  summarise(mean(p.value < 5e-8))

## discovery comparison
out %>% filter(positive == TRUE,
               counteract == FALSE,
               method %in% c("aLIT", "Marginal_SQ", "Marginal")) %>%
  filter(prop_null == .5, num_traits == 10, heritability %in% c(0.35, .6)) %>%
  group_by(num_traits, method, heritability) %>%
  summarise(mean(p.value < 5e-8))

out_plot <- out %>%
  filter(positive == TRUE,
         counteract == FALSE,
         method %in% c("aLIT", "uLIT", "wLIT"))
out_plot$method <- factor(out_plot$method,
                          labels = c("aLIT", "uLIT", "wLIT"))

cbbPalette <- c( "black", "#0072B2", "#009E73")

out_plot$baseline = out_plot$heritability + out_plot$env
p_positive <- out_plot %>%
  group_by(N, baseline, counteract, prop_null, positive, num_traits, method) %>%
  summarise(power = mean(p.value < 5e-8)) %>%
  ggplot(aes(x = 1 - prop_null, y = power, col = method, group = rev(method))) +
  geom_point(alpha = 0.5, size = 1.5) + geom_line(size = 0.65) +
  facet_grid(num_traits ~ baseline, labeller = label_parsed) +
  xlab("Proportion of traits with shared interaction effects") +
  ylab("Power") +
  scale_colour_manual(values=cbbPalette) +
  labs(col = "Method") +
  scale_x_continuous(sec.axis = sec_axis(~ ., name = "Trait correlation", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Number of traits", breaks = NULL, labels = NULL)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())

# positive/negative pleiotropy
out_plot <- out %>%
  filter(positive == FALSE,
         counteract == FALSE,
         method %in% c("aLIT", "uLIT", "wLIT"))
out_plot$method <- factor(out_plot$method,
                          labels = c("aLIT", "uLIT", "wLIT"))

cbbPalette <- c( "black", "#0072B2","#009E73")
out_plot$baseline = out_plot$heritability + out_plot$env
p_pn <- out_plot %>%
  group_by(N, baseline, counteract, prop_null, positive, num_traits, method) %>%
  summarise(power = mean(p.value < 5e-8)) %>%
  ggplot(aes(x = 1 - prop_null, y = power, col = method, group = rev(method))) +
  geom_point(alpha = 0.5, size = 1.5) + geom_line(size = 0.65) +
  facet_grid(num_traits ~ baseline, labeller = label_parsed) +
  xlab("Proportion of traits with shared interaction effects") +
  ylab("Power") +
  scale_colour_manual(values = cbbPalette) +
  labs(col = "Method") +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Trait correlation", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of traits", breaks = NULL, labels = NULL)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())

# combine plots
prow <- plot_grid(
  p_positive + theme(legend.position="none"),
  p_pn + theme(legend.position="none"),
  labels = c("A", "B"),
  nrow = 2,
  label_size = 12
)

legend <- get_legend(
  p_positive + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p <- plot_grid(prow, legend, rel_widths = c(3.5, .42))

ggsave(p,  width = 7.0, height = 7.35, filename = "./figures/alt-simulations-lit.png")

# positive pleiotropy w/ counteracting effect size direction
out_plot <- out %>%
  filter(positive == TRUE,
         counteract == TRUE,
         method %in% c("aLIT", "uLIT", "wLIT"))
out_plot$method <- factor(out_plot$method,
                          labels = c("aLIT", "uLIT", "wLIT"))

cbbPalette <- c( "black", "#0072B2","#009E73")
out_plot$baseline = out_plot$heritability + out_plot$env
p_positive <- out_plot %>%
  group_by(N,  baseline, counteract, prop_null, positive, num_traits, method) %>%
  summarise(power = mean(p.value < 5e-8)) %>%
  ggplot(aes(x = 1 - prop_null, y = power, col = method, group = rev(method))) +
  geom_point(alpha = 0.5, size = 1.5) + geom_line(size = 0.65) +
  facet_grid(num_traits ~ baseline, labeller = label_parsed) +
  xlab("Proportion of traits with shared interaction effects") +
  ylab("Power") +
  scale_colour_manual(values=cbbPalette) +
  labs(col = "Method") +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Trait correlation", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of traits", breaks = NULL, labels = NULL)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())

# positive/negative pleiotropy w/ counteracting effect size
out_plot <- out %>%
  filter(positive == FALSE,
         counteract == TRUE,
         method %in% c("aLIT", "uLIT", "wLIT"))
out_plot$method <- factor(out_plot$method,
                          labels = c("aLIT", "uLIT", "wLIT"))

cbbPalette <- c( "black", "#0072B2","#009E73")
out_plot$baseline = out_plot$heritability + out_plot$env
p_pn <- out_plot %>%
  group_by(N,  baseline, counteract, prop_null, positive, num_traits, method) %>%
  summarise(power = mean(p.value < 5e-8)) %>%
  ggplot(aes(x = 1 - prop_null, y = power, col = method, group = rev(method))) +
  geom_point(alpha = 0.5, size = 1.5) + geom_line(size = 0.65) +
  facet_grid(num_traits ~ baseline, labeller = label_parsed) +
  xlab("Proportion of traits with shared interaction effects") +
  ylab("Power") +
  scale_colour_manual(values=cbbPalette) +
  labs(col = "Method") +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Trait correlation", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of traits", breaks = NULL, labels = NULL)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())

# combine plots
prow <- plot_grid(
  p_positive + theme(legend.position="none"),
  p_pn + theme(legend.position="none"),
  labels = c("A", "B"),
  nrow = 2,
  label_size = 12
)

legend <- get_legend(
  p_positive + theme(legend.box.margin = margin(0, 0, 0, 12))
)

p <- plot_grid(prow, legend, rel_widths = c(3.5, .42))
ggsave(p,  width = 7.0, height = 7.35, filename = "./figures/alt-simulations-lit-counteract.png")

# positive pleiotropy for comparison with marginal procedures
out_plot <- out %>%
  filter(positive == TRUE,
         counteract == FALSE,
         method %in% c("aLIT", "Marginal", "Marginal_SQ"))
out_plot$method <- factor(out_plot$method,
                          labels = c("aLIT", "Marginal (SQ/CP)", "Marginal (SQ)"))

cbbPalette <- c("black", "#8c8c8c", "#cccccc")
out_plot$baseline = out_plot$heritability + out_plot$env
p_positive <- out_plot %>%
  group_by(N, baseline, counteract, prop_null, positive, num_traits, method) %>%
  summarise(power = mean(p.value < 5e-8)) %>%
  ggplot(aes(x = 1 - prop_null, y = power, col = method, group = rev(method))) +
  geom_point(alpha = 0.5, size = 1.5) + geom_line(size = 0.65) +
  facet_grid(num_traits ~ baseline, labeller = label_parsed) +
  xlab("Proportion of traits with shared interaction effects") +
  ylab("Power") +
  scale_colour_manual(values=cbbPalette) +
  labs(col = "Method") +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Trait correlation", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of traits", breaks = NULL, labels = NULL)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())

# positive/negative pleiotropy for comparison with marginal procedures
out_plot <- out %>%
  filter(positive == FALSE,
         counteract == FALSE,
         method %in% c("aLIT", "Marginal", "Marginal_SQ"))
out_plot$method <- factor(out_plot$method,
                          labels = c("aLIT", "Marginal (SQ/CP)", "Marginal (SQ)"))

cbbPalette <- c("black", "#8c8c8c", "#cccccc")#c( "black", "#0072B2","#009E73")
out_plot$baseline = out_plot$heritability + out_plot$env
p_pn <- out_plot %>%
  group_by(N, baseline, counteract, prop_null, positive, num_traits, method) %>%
  summarise(power = mean(p.value < 5e-8)) %>%
  ggplot(aes(x = 1 - prop_null, y = power, col = method, group = rev(method))) +
  geom_point(alpha = 0.5, size = 1.5) + geom_line(size = 0.65) +
  facet_grid(num_traits ~ baseline, labeller = label_parsed) +
  xlab("Proportion of traits with shared interaction effects") +
  ylab("Power") +
  scale_colour_manual(values=cbbPalette) +
  labs(col = "Method") +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Trait correlation", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of traits", breaks = NULL, labels = NULL)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())

# combine plots
prow <- plot_grid(
  p_positive + theme(legend.position="none"),
  p_pn + theme(legend.position="none"),
  labels = c("A", "B"),
  nrow = 2,
  label_size = 12
)

legend <- get_legend(
  p_positive + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p <- plot_grid(prow, legend, rel_widths = c(3, .82))

ggsave(p,  width = 7.0, height = 7.35, filename = "./figures/alt-simulations-lit-vs-marginal.png")

# positive pleiotropy with counteracting effect sizes for comparison with marginal procedures
out_plot <- out %>%
  filter(positive == TRUE,
         counteract == TRUE,
         method %in% c("aLIT", "Marginal", "Marginal_SQ"))
out_plot$method <- factor(out_plot$method,
                          labels = c("aLIT", "Marginal (SQ/CP)", "Marginal (SQ)"))

cbbPalette <- c("black", "grey49", "grey70")#c( "black", "#0072B2","#009E73")
out_plot$heritability = out_plot$heritability + out_plot$env
p_positive <- out_plot %>%
  group_by(N,  heritability, counteract, prop_null, positive, num_traits, method) %>%
  summarise(power = mean(p.value < 5e-8)) %>%
  ggplot(aes(x = 1 - prop_null, y = power, col = method, group = rev(method))) +
  geom_point(alpha = 0.5, size = 1.5) + geom_line(size = 0.65) +
  facet_grid(num_traits ~heritability, labeller = label_parsed) +
  xlab("Proportion of traits with shared interaction effects") +
  ylab("Power") +
  scale_colour_manual(values=cbbPalette) +
  labs(col = "Method") +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Trait correlation", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of traits", breaks = NULL, labels = NULL)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())

# positive/negative pleiotropy with counteracting effect sizes for comparison with marginal procedures
out_plot <- out %>%
  filter(positive == FALSE,
         counteract == TRUE,
         method %in% c("aLIT", "Marginal", "Marginal_SQ"))
out_plot$method <- factor(out_plot$method,
                          labels = c("aLIT", "Marginal (SQ/CP)", "Marginal (SQ)"))

out_plot$baseline = out_plot$heritability + out_plot$env
p_pn <- out_plot %>%
  group_by(N, baseline, counteract, prop_null, positive, num_traits, method) %>%
  summarise(power = mean(p.value < 5e-8)) %>%
  ggplot(aes(x = 1 - prop_null, y = power, col = method, group = rev(method))) +
  geom_point(alpha = 0.5, size = 1.5) + geom_line(size = 0.65) +
  facet_grid(num_traits ~ baseline, labeller = label_parsed) +
  xlab("Proportion of traits with shared interaction effects") +
  ylab("Power") +
  scale_colour_manual(values=cbbPalette) +
  labs(col = "Method") +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Trait correlation", breaks = NULL, labels = NULL)) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of traits", breaks = NULL, labels = NULL)) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = NA, color = "black"), axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank())

# combine plots
prow <- plot_grid(
  p_positive + theme(legend.position="none"),
  p_pn + theme(legend.position="none"),
  labels = c("A", "B"),
  nrow = 2,
  label_size = 12
)

legend <- get_legend(
  p_positive + theme(legend.box.margin = margin(0, 0, 0, 12))
)

p <- plot_grid(prow, legend, rel_widths = c(3, .82))
ggsave(p,  width = 7.0, height = 7.35, filename = "./figures/alt-simulations-counteract-lit-vs-marginal.png")
###

########### 03-lit-time.R plots ###########
out <- readRDS(file = "./data/03-lit-time.rds")
cbbPalette <- c("#009E73", "#0072B2", "#D55E00", "#CC79A7", "green")

p <- out %>% group_by(N, num_traits) %>%
  summarise(time = mean(time)) %>%
  ggplot(aes(x = N / 1000 , y = time, color = as.factor(num_traits))) +
  geom_point() +
  geom_line(size = 0.5)+
  xlab("Sample size / 10") +
  ylab("Time (s)") +
  theme_bw(base_size = 12) +
  xlab("Sample size (x1000)") +
  scale_colour_manual(values=cbbPalette) +
  labs(color = "Traits")

ggsave(p, width = 4, height = 2.5 ,filename = "./figures/computational_time.png")

########### 04-PC.R plots ###########
out <- readRDS(file = "./data/04-PC.rds")
p <- out %>%
  group_by(PC, heritability, prop_null) %>%
  filter(as.factor(prop_null) %in% as.factor(seq(0.0, 1,0.2))) %>%
  summarise(power = mean(p < 5e-8)) %>%
  mutate(prop_null = 1-prop_null, ID = ifelse(PC <= 1, "PC 1",
                                              ifelse(PC <= 5, "PC 2-5",
                                                     ifelse(PC <=12, "PC 6-12",
                                                            ifelse(PC <= 25, "PC 13-25",
                                                                   ifelse(PC <= 44, "PC 26-44"))))))%>%
  mutate(ID = factor(ID, levels =c('PC 1', "PC 2-5", "PC 6-12", "PC 13-25", "PC 26-44") )) %>%
  ggplot(aes(x = heritability + 0.15, y = power, group = PC)) +
  geom_line(size = .8, alpha = 0.4) +
  geom_point(alpha = 0.5, size = 1.0)+
  facet_grid(ID~prop_null) +
  theme_bw() +
  xlab("Trait correlation") +
  ylab("Empirical power") +
  ggtitle("Proportion of traits with shared interaction effects")+theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank()) +
  scale_color_brewer(palette = "Set1")

ggsave(p, width = 7, height = 7 ,filename = "./figures/PC.png")
