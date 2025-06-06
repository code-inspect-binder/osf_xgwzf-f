rm(list = ls())

# change to project directory
setwd("/Volumes/GoogleDrive/My Drive/grad-school/scaffold_behav/manuscript/psychsci_submission/psychsci_scaffold")

# setup
library(Rmisc)
library(tidyverse)
library(ggbeeswarm)
library(effsize)
library(lme4)
source("scripts/utils.R")

# set plot theme
plot_theme <- theme_light() + theme(panel.grid = element_blank(), legend.position = 'none')
theme_set(plot_theme)
cond_cols <- c('salmon2','dodgerblue2')

data_dir <- "clean_data"
plot_dir <- "plots"

# ----------------
# PRETRAINING
# ----------------

## 1. pretraining study RTs for pred v. rand sequences

# load data
pretrain_study <- read.csv(paste0(data_dir, "/exp1_pretrain_study_clean.csv"), stringsAsFactors = F)

# remove outlier RTs
all_rts <- filter(pretrain_study, !is.na(response))$rt
sd_cutoff <- mean(all_rts) + (3 * sd(all_rts))
pretrain_study <- filter(pretrain_study, rt <= sd_cutoff)

# mean RTs by condition
pretrain_study_rt <- pretrain_study %>%
  filter(block == 3) %>%
  group_by(id, condition) %>%
  summarise(mean_rt = mean(rt)) %>%
  ungroup()

# summarize across participants
pretrain_study_rt_group <- pretrain_study_rt %>%
  group_by(condition) %>%
  summarise(group_mean_rt = mean(mean_rt)) %>%
  ungroup()
pretrain_study_rt_group$group_sem_rt <- summarySEwithin(data = pretrain_study, measurevar = "rt",
                                                        withinvars = "condition", idvar = "id")$se

# plot - fig. 2b
ggplot(pretrain_study_rt, aes(x = condition, y = mean_rt, fill = condition)) +
  geom_bar(data = pretrain_study_rt_group, aes(x = condition, y = group_mean_rt),
           stat = 'identity', width = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.2, color = 'white') +
  geom_errorbar(data = pretrain_study_rt_group, width = 0.1,
                aes(x = condition, y = group_mean_rt,
                    ymin = group_mean_rt - group_sem_rt, ymax = group_mean_rt + group_sem_rt)) +
  scale_fill_manual(values = cond_cols) +
  labs(x = 'condition', y = 'RT (ms)')
ggsave(paste0(plot_dir, "/exp1_pretrainRT.pdf"), width = 4, height = 3.5)

# stats
wilcox.test(mean_rt ~ condition, data = pretrain_study_rt, paired = T)
effsize::cliff.delta(mean_rt ~ condition, data = pretrain_study_rt, paired = T)

## 2. pretraining test accuracy for pred sequences

pretrain_test <- read.csv(paste0(data_dir, "/exp1_pretrain_test_clean.csv"), stringsAsFactors = F)

# mean acc by condition and test repetition
pretrain_test_acc <- pretrain_test %>%
  group_by(id, test_rep) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()

# summarize across participants
pretrain_test_acc_group <- pretrain_test_acc %>%
  group_by(test_rep) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
pretrain_test_acc_group$group_sem_acc <- summarySEwithin(data = pretrain_test_acc, measurevar = "mean_acc",
                                                         withinvars = "test_rep", idvar = "id")$se

# plot - fig. 2b
ggplot(pretrain_test_acc, aes(x = test_rep, y = mean_acc)) +
  geom_jitter(width = 0.2, height = 0.02, fill = 'dodgerblue1', color = 'white',
              alpha = 0.4, shape = 21, size = 2.3) +
  geom_line(data = pretrain_test_acc_group, aes(x = test_rep, y = group_mean_acc),
            size = 1.5, color = 'dodgerblue4') +
  geom_errorbar(data = pretrain_test_acc_group,
                aes(x = test_rep, y = group_mean_acc,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc),
                size = 1.5, width = 0, color = 'dodgerblue4') +
  scale_x_continuous(breaks = 1.6) +
  labs(x = 'sequence test repetition', y = 'accuracy')
ggsave(paste0(plot_dir, "/exp1_pretrain_testAcc.pdf"), width = 5, height = 3.5)

# show get descriptive stats (how many/what proportion of participants reach each accuracy level)
pretrain_test_acc %>%
  group_by(test_rep) %>%
  summarise(n_80 = sum(mean_acc > 0.8),
            n_100 = sum(mean_acc == 1),
            prop_100 = sum(mean_acc == 1) / length(mean_acc),
            prop_80 = sum(mean_acc > 0.8) / length(mean_acc))
pretrain_test_acc %>%
  filter(test_rep == max(pretrain_test_acc$test_rep)) %>%
  summarise(group_mean_acc = mean(mean_acc), group_sd_acc = sd(mean_acc))


# ----------------
# ENCODING
# ----------------

## 1. encoding accuracy for pred v. rand events

enc <- read.csv(paste0(data_dir, "/exp1_enc_clean.csv"), stringsAsFactors = F)

# mean acc by condition
enc_resp_acc <- enc %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()

# summarize across participants
enc_resp_acc_group <- enc_resp_acc %>%
  group_by(condition) %>%
  summarise(group_mean_acc = mean(mean_acc), group_sd_acc = sd(mean_acc)) %>%
  ungroup()
enc_resp_acc_group

# stats
t.test(mean_acc ~ condition, data = enc_resp_acc, paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = enc_resp_acc, paired = T)

# mean acc by condition and block/list number
enc_resp_acc <- enc %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id, condition, list) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()

# stats
ttest_multiple_fdr(enc_resp_acc, equation = 'mean_acc ~ condition', grouping = 'list')

# 2. encoding RTs for pred v. rand events

# remove outlier RTs
all_rts <- enc$rt
sd_cutoff <- mean(all_rts) + (3 * sd(all_rts))
enc <- filter(enc, rt <= sd_cutoff)

# mean rt by condition
enc_resp_rt <- enc %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id, condition) %>%
  summarise(mean_rt = mean(rt)) %>%
  ungroup()

# stats
wilcox.test(mean_rt ~ condition, data = enc_resp_rt, paired = T)
effsize::cliff.delta(mean_rt ~ condition, data = enc_resp_rt, paired = T)

# ----------------
# ORDER RECONSTRUCTION
# ----------------

recon <- read.csv(paste0(data_dir, "/exp1_recon_clean.csv"), stringsAsFactors = F)

# 1. order recon accuracy for pred v. rand events — ordinal accuracy

# order memory by condition
recon_acc <- recon %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
recon_acc_group <- recon_acc %>%
  group_by(condition) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
recon_acc_group$group_sem_acc <- summarySEwithin(data = recon_acc, measurevar = "mean_acc",
                                                 withinvars = "condition", idvar = "id")$se

# plot - fig. 3a
ggplot(recon_acc, aes(x = condition, y = mean_acc, fill = condition)) +
  geom_bar(data = recon_acc_group, aes(x = condition, y = group_mean_acc),
           stat = 'identity', width = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.5, color = 'white') +
  geom_errorbar(data = recon_acc_group, width = 0.1,
                aes(x = condition, y = group_mean_acc,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc)) +
  geom_hline(yintercept = 1/6, color = 'grey20', linetype = 'dashed') +
  scale_fill_manual(values = cond_cols) +
  labs(x = 'condition', y = 'ordinal accuracy')
ggsave(paste0(plot_dir, "/exp1_reconAcc.pdf"), width = 4.5, height = 3.5)

t.test(mean_acc ~ condition, data = recon_acc, paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = recon_acc, paired = T)

# 2. order recon accuracy for pred v. rand events — levenshtein distance

# get dataframe with this measure 
recon_leven_data <- suppressMessages(get_leven_dist(recon))
# implemented in utils.R script for readability
# without message suppression, this will print out a bunch of warnings, because participants
#     don't always select all 6 of the items during each recon trial

# levenshtein distance by condition
recon_leven <- recon_leven_data %>%
  group_by(id, condition) %>%
  summarise(mean_leven_dist = mean(leven_dist)) %>%
  ungroup()
recon_leven_group <- recon_leven %>%
  group_by(condition) %>%
  summarise(group_mean_leven_dist = mean(mean_leven_dist)) %>%
  ungroup()
recon_leven_group$group_sem_leven_dist <- summarySEwithin(data = recon_leven,
                                                          measurevar = "mean_leven_dist",
                                                          withinvars = "condition", idvar = "id")$se

# stats
t.test(mean_leven_dist ~ condition, data = recon_leven, paired = T)
effsize::cohen.d(mean_leven_dist ~ condition | Subject(id), data = recon_leven, paired = T)

# 2. proportion of items remembered in correct order within pred v. rand events

recon_event <- recon %>%
  group_by(id, list, event) %>%
  summarise(n_correct = sum(accuracy), condition = get_mode(condition)) %>%
  group_by(id, condition) %>%
  summarise(prop_6 = sum(n_correct == 6) / length(n_correct),
            prop_5 = sum(n_correct == 5) / length(n_correct),
            prop_4 = sum(n_correct == 4) / length(n_correct),
            prop_3 = sum(n_correct == 3) / length(n_correct),
            prop_2 = sum(n_correct == 2) / length(n_correct),
            prop_1 = sum(n_correct == 1) / length(n_correct),
            prop_0 = sum(n_correct == 0) / length(n_correct)) %>%
  pivot_longer(cols = starts_with('prop_'), names_to = 'num_correct', values_to = 'proportion') %>%
  mutate(num_correct = as.numeric(gsub('prop_', '', num_correct)),
         condition = as.factor(condition)) %>%
  ungroup()
recon_event_group <- recon_event %>%
  group_by(condition, num_correct) %>%
  summarise(group_mean_prop = mean(proportion)) %>%
  ungroup()
recon_event_group$group_sem_prop <- summarySEwithin(data = recon_event, measurevar = "proportion",
                                                    withinvars = c("condition", "num_correct"),
                                                    idvar = "id")$se

# plot - fig. 3b
ggplot(recon_event_group, aes(x = num_correct, y = group_mean_prop, fill = condition)) +
  geom_bar(stat = 'identity', position = position_dodge(0.9)) +
  geom_point(data = recon_event, aes(x = num_correct, y = proportion, color = condition),
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0.01, dodge.width = 0.9),
             pch = 21, color = 'white', size = 2) +
  geom_errorbar(aes(ymin = group_mean_prop - group_sem_prop,
                    ymax = group_mean_prop + group_sem_prop),
                position = position_dodge(0.9), width = 0) +
  labs(x = 'number of items selected in correct order (per event)', y = 'proportion of events') +
  scale_x_continuous(breaks = seq(0, 6)) +
  scale_fill_manual(values = cond_cols) + scale_color_manual(values = cond_cols)
ggsave(paste0(plot_dir, "/exp1_reconAcc_event.pdf"), width = 4.5, height = 3.5)

# stats
ttest_multiple_fdr(recon_event, equation = "proportion ~ condition", grouping = "num_correct")

# 3. order recon for pred v. rand events as a function of sequence position

# order memory by condition and seq pos
recon_acc <- recon %>%
  group_by(id, condition, seq_pos) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
recon_acc_group <- recon_acc %>%
  group_by(condition, seq_pos) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
recon_acc_group$group_sem_acc <- summarySEwithin(data = recon_acc, measurevar = "mean_acc",
                                                 withinvars = c("condition", "seq_pos"), idvar = "id")$se

# plot - fig. 3c
ggplot(recon_acc_group, aes(x = seq_pos, y = group_mean_acc, color = condition)) +
  geom_hline(yintercept = 1/6, color = 'grey20', linetype = 'dashed') +
  geom_jitter(data = recon_acc, aes(x = seq_pos, y = mean_acc, fill = condition),
              width = 0.1, height = 0.02, color = 'white', alpha = 0.4, shape = 21, size = 2) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(x = seq_pos,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc),
                width = 0, size = 1) +
  scale_x_continuous(breaks = seq(1, 6)) +
  scale_color_manual(values = cond_cols) + scale_fill_manual(values = cond_cols) +
  labs(x = 'sequence position', y = 'ordinal accuracy')
ggsave(paste0(plot_dir, "/exp1_reconAcc_seqpos.pdf"), width = 4.5, height = 3.5)

# stats - ttest
ttest_multiple_fdr(recon_acc, equation = "mean_acc ~ condition", grouping = "seq_pos")

# stats - logistic regression model
recon_model_data <- mutate(recon,
                           condition_e = ifelse(condition == 'pred', 0.5, -0.5), # effect coded
                           seq_pos_z = scale(seq_pos, center = T, scale = F)) # center pos = 0
recon_model <- glmer(accuracy ~ condition_e * seq_pos_z + (condition_e + seq_pos_z || id),
                     data = recon_model_data, family = 'binomial')
summary(recon_model)
confint(recon_model) # takes awhile to run


# ----------------
# SPATIAL MEMORY
# ----------------

spatial <- read.csv(paste0(data_dir, "/exp1_spatial_clean.csv"), stringsAsFactors = F)

# 1. spatial memory for pred v. rand events

# mean acc by condition
spatial_acc <- spatial %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
spatial_acc_group <- spatial_acc %>%
  group_by(condition) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
spatial_acc_group$group_sem_acc <- summarySEwithin(data = spatial_acc, measurevar = "mean_acc",
                                                   withinvars = "condition", idvar = "id")$se

# plot - fig. 4a
ggplot(spatial_acc, aes(x = condition, y = mean_acc, fill = condition)) +
  geom_bar(data = spatial_acc_group, aes(x = condition, y = group_mean_acc),
           stat = 'identity', width = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.5, color = 'white') +
  geom_errorbar(width = 0.1, data = spatial_acc_group,
                aes(x = condition, y = group_mean_acc,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc)) +
  geom_hline(yintercept = 1/4, color = 'grey20', linetype = 'dashed') +
  scale_fill_manual(values = cond_cols) +
  labs(x = 'condition', y = 'spatial memory accuracy')
ggsave(paste0(plot_dir, "/exp1_spatialAcc.pdf"), width = 4.5, height = 3.5)

# stats
t.test(mean_acc ~ condition, data = spatial_acc, paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = spatial_acc, paired = T)

# 2. spatial memory for pred v. rand events as a function of sequence position

# order memory by condition and seq pos
spatial_acc <- spatial %>%
  group_by(id, condition, seq_pos) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
spatial_acc_group <- spatial_acc %>%
  group_by(condition, seq_pos) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
spatial_acc_group$group_sem_acc <- summarySEwithin(data = spatial_acc, measurevar = "mean_acc",
                                                   withinvars = c("condition", "seq_pos"), idvar = "id")$se

# plot - fig. 4b
ggplot(spatial_acc_group, aes(x = seq_pos, y = group_mean_acc, color = condition)) +
  geom_hline(yintercept = 1/4, color = 'grey20', linetype = 'dashed') +
  geom_jitter(data = spatial_acc, aes(x = seq_pos, y = mean_acc, fill = condition),
              width = 0.1, height = 0.02, color = 'white', alpha = 0.4, shape = 21, size = 3) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(x = seq_pos,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc),
                width = 0, size = 1) +
  scale_color_manual(values = cond_cols) + scale_fill_manual(values = cond_cols) +
  scale_x_continuous(breaks = seq(1,6)) +
  labs(x = 'sequence position', y = 'spatial memory accuracy')
ggsave(paste0(plot_dir, "/exp1_spatialAcc_seqpos.pdf"), width = 5, height = 3.5)

# stats - ttest
ttest_multiple_fdr(spatial_acc, equation = "mean_acc ~ condition", grouping = "seq_pos")
