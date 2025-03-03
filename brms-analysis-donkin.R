
library("tidyverse")
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position="bottom", 
                  panel.grid.major.x = element_blank()))
library("brms")
options(mc.cores = parallel::detectCores())
library("patchwork")
library("bayestestR")
library("emmeans")

## data from: https://www2.psy.unsw.edu.au/Users/CDonkin/data.htm
## Donkin, Tran & Nosofsky (2014): Experiment 2
df <- list.files("DTN_2014_APP_Exp2/", full.names = TRUE)
dl <- map_dfr(df, read_table) %>% 
  mutate(change_prob = factor(Prop))

dl

dl %>% 
  summarise(n_distinct(Subj))
## should be 20

## three set sizes: 2, 5, 8
p1 <- dl %>% 
  group_by(Subj, SetSize) %>% 
  summarise(acc = mean(ACC)) %>% 
  ggplot(aes(x = factor(SetSize), y = acc)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.3) +
  stat_summary(colour = "red")

## three change probabilities: .15,.3,.5,.7,.85
p2 <- dl %>% 
  group_by(Subj, change_prob) %>% 
  summarise(acc = mean(ACC)) %>% 
  ggplot(aes(x = change_prob, y = acc)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.3) +
  stat_summary(colour = "red")

p1/p2


p3 <- dl %>% 
  group_by(Subj, CorrResp, SetSize) %>% 
  summarise(resp = mean(keyIdx)) %>% 
  ggplot(aes(x = factor(SetSize), y = resp)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.3) +
  stat_summary(colour = "red") +
  facet_wrap(vars(CorrResp))

p4 <- dl %>% 
  group_by(Subj, CorrResp, change_prob) %>% 
  summarise(resp = mean(keyIdx)) %>% 
  ggplot(aes(x = change_prob, y = resp)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.3) +
  stat_summary(colour = "red") +
  facet_wrap(vars(CorrResp))
p3/p4

##---------------------------------------------------------------
##                          Modelling                           -
##---------------------------------------------------------------


##---------------
##  Slots Model  
##---------------

source("slots-model-brms.R")

options(contrasts=c('contr.equalprior_deviations', 'contr.poly'))
fit_slot <- brm(
  brmsformula(
    keyIdx | vint(CorrResp, SetSize) ~ 1 + (1|p|Subj), ## formula for first model parameter (m)
    ## because set size is part of m, only RI necessary
    g ~ change_prob + (change_prob|p|Subj) ## for g, we need effect of change probability
    ## "p" in between "|" indicates that we estimate correlations for random terms across terms 
  ), data = dl,
  family = slots, stanvars = stanvars_slots
)
save(fit_slot, file = "fit-slot.rda", compress = "xz")
load("fit-slot.rda")
fit_slot

emmeans(fit_slot, "change_prob", dpar = "g", type = "response")
emmeans(fit_slot, "1", dpar = "mu", type = "response")

pred_slots <- posterior_epred(fit_slot)
str(pred_slots)

dl$pred_slots <- apply(pred_slots, 2, mean)

dl

dl %>% 
  group_by(Subj, CorrResp, change_prob, SetSize) %>% 
  summarise(resp = mean(keyIdx),
            pred = mean(pred_slots)) %>% 
  ggplot(aes(x = resp, y = pred)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.3) +
  coord_fixed(xlim = c(0, 1), ylim = c(0,1)) +
  facet_grid(rows = vars(SetSize, CorrResp), cols = vars(change_prob))


##------------------------------
##  Slots Model with Attention  
##------------------------------

source("slots2-model-brms.R")

fit_slot2_1 <- brm(
  brmsformula(
    keyIdx | vint(CorrResp, SetSize) ~ 1 + (1|p|Subj),
    g ~ change_prob + (change_prob|p|Subj),
    a ~ 1 + (1|p|Subj)
  ), data = dl,
  family = slots2, stanvars = stanvars_slots2
)
fit_slot2_1
save(fit_slot2_1, file = "fit-slot2-1.rda", compress = "xz")
load("fit-slot2-1.rda")

emmeans(fit_slot2_1, "1", dpar = "mu", type = "response")
emmeans(fit_slot2_1, "1", dpar = "a", type = "response")

emmeans(fit_slot2_1, "change_prob", dpar = "g", type = "response")

pred_slots2 <- posterior_epred(fit_slot2_1)
dl$pred_slots2 <- apply(pred_slots2, 2, mean)

psize <- 3
p3 <- dl %>% 
  group_by(Subj, CorrResp, SetSize) %>% 
  summarise(resp = mean(keyIdx),
            pred = mean(pred_slots),
            pred2 = mean(pred_slots2)) %>% 
  ggplot(aes(x = factor(SetSize), y = resp)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.3) +
  stat_summary(colour = "green") + 
  stat_summary(aes(y = pred), colour = "blue", shape = 4, 
               geom = "point", size = psize) +
  stat_summary(aes(y = pred2), colour = "orange", shape = 5
               , geom = "point", size = psize) +
  facet_wrap(vars(CorrResp))

p4 <- dl %>% 
  group_by(Subj, CorrResp, change_prob) %>% 
  summarise(resp = mean(keyIdx),
            pred = mean(pred_slots),
            pred2 = mean(pred_slots2)) %>% 
  ggplot(aes(x = change_prob, y = resp)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.3) +
  stat_summary(colour = "green") +
  stat_summary(aes(y = pred), colour = "blue", shape = 4
               , geom = "point", size = psize) +
  stat_summary(aes(y = pred2), colour = "orange", shape = 5
               , geom = "point", size = psize) +
  facet_wrap(vars(CorrResp))
p3/p4

pc1 <- dl %>% 
  group_by(Subj, CorrResp, change_prob, SetSize) %>% 
  summarise(resp = mean(keyIdx),
            pred = mean(pred_slots)) %>% 
  ggplot(aes(x = resp, y = pred)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.3) +
  coord_fixed(xlim = c(0, 1), ylim = c(0,1)) +
  facet_grid(rows = vars(SetSize, CorrResp), cols = vars(change_prob)) 
pc2 <- dl %>% 
  group_by(Subj, CorrResp, change_prob, SetSize) %>% 
  summarise(resp = mean(keyIdx),
            pred2 = mean(pred_slots2)) %>% 
  ggplot(aes(x = resp, y = pred2)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point(alpha = 0.3) +
  coord_fixed(xlim = c(0, 1), ylim = c(0,1)) +
  facet_grid(rows = vars(SetSize, CorrResp), cols = vars(change_prob)) 
pc1+pc2


#######
dl

dl %>% 
  count(Blk, Trl) %>% 
  print(n = Inf)

dl %>% 
  count(Blk, change_prob) %>% 
  print(n = Inf)

dl %>% 
  group_by(Blk) %>% 
  summarise(acc = mean(ACC)) %>% 
  ggplot(aes(x = Blk, y = acc)) +
  geom_line(aes(group = 1)) +
  geom_smooth() 

dl %>% 
  group_by(Blk, Trl) %>% 
  summarise(acc = mean(ACC)) %>% 
  ggplot(aes(x = Trl, y = acc)) +
  geom_line(aes(group = 1)) +
  geom_smooth() +
  facet_wrap(vars(Blk))

fit_slot2_2 <- brm(
  brmsformula(
    keyIdx | vint(CorrResp, SetSize) ~ Blk + (Blk|p|Subj),
    g ~ change_prob + Blk + (change_prob+Blk|p|Subj),
    a ~ Blk + (Blk|p|Subj)
  ), data = dl,
  family = slots2, stanvars = stanvars_slots2
)
fit_slot2_2
save(fit_slot2_2, file = "fit-slot2-2.rda", compress = "xz")
load("fit-slot2-2.rda")

