#########################################################################################
####                 Tools_analysis_simulations                                      ####
#### This script provides functions that are helpful for the analysis of simulations ####
#### and examples on their usage.                                                    ####
#########################################################################################

# Load libraries

library(ggplot2)
library(tidyverse)
library(magrittr)
library(cowplot)
library(grid)
library(gridExtra)
library(Cairo)

#### Define functions ####
# A function that will receive the folder with the six tables of fixed and proposed mutations to summarize them
# The user has to provide the other parameters used
# Returns a list with:
# - the energy values for the trajectory plots
# - the ddGs for the scatterplots and supplementary figures on mutational effects 
condensed_energies <- function(folder, num_reps, num_subs, scenario, beta, N){
  setwd(folder)
  
  outlist <- list()
  
  # Load the fixed mutations
  het_fixed <- read.table('all_deltaGs_heterodimer_all_reps.tab', h = T, sep = '\t') 
  het_fixed %<>%
    mutate(Replicate = rep(1:num_reps, each = (nrow(het_fixed)/num_reps)),
           Substitution = rep(0:num_subs, nrow(het_fixed)/(num_subs+1)),
           ComplexType = 'Heterodimer AB',
           Scenario = scenario,
           Beta = beta,
           N = N)
  
  homo_AA_fixed <- read.table('all_deltaGs_homodimer_A_all_reps.tab', h = T, sep = '\t') 
  homo_AA_fixed %<>%
    mutate(Replicate = rep(1:num_reps, each = (nrow(homo_AA_fixed)/num_reps)),
           Substitution = rep(0:num_subs, nrow(homo_AA_fixed)/(num_subs+1)),
           ComplexType = 'Homodimer AA',
           Scenario = scenario,
           Beta = beta,
           N = N)
  
  homo_BB_fixed <- read.table('all_deltaGs_homodimer_B_all_reps.tab', h = T, sep = '\t')
  homo_BB_fixed %<>%
    mutate(Replicate = rep(1:num_reps, each = (nrow(homo_BB_fixed)/num_reps)),
           Substitution = rep(0:num_subs, nrow(homo_BB_fixed)/(num_subs+1)),
           ComplexType = 'Homodimer BB',
           Scenario = scenario,
           Beta = beta,
           N = N)
  
  outlist$all_energies <- rbind(het_fixed, homo_AA_fixed, homo_BB_fixed)
  
  previous_HET <- rep(NA, nrow(het_fixed))
  previous_HM_AA <- rep(NA, nrow(het_fixed))
  previous_HM_BB <- rep(NA, nrow(het_fixed))
  Verdict = rep("REJECTED", nrow(het_fixed))
  
  for(i in 1:nrow(het_fixed)){
    previous_HET[i] <- ifelse(het_fixed$Substitution[i] == 0, NA, het_fixed$Binding_energy[i-1])
    previous_HM_AA[i] <- ifelse(homo_AA_fixed$Substitution[i] == 0, NA, homo_AA_fixed$Binding_energy[i-1])
    previous_HM_BB[i] <- ifelse(homo_BB_fixed$Substitution[i] == 0, NA, homo_BB_fixed$Binding_energy[i-1])
    
    # Recall the selection scenario
    if(scenario == "Selection on heterodimer AB"){
      Verdict <- het_fixed$Verdict
    }else if(scenario == "Selection on both homodimers"){
      if(is.na(homo_AA_fixed$Verdict[i]) || is.na(homo_BB_fixed$Verdict[i])){
        Verdict[i] = NA
      }
      else if(and(homo_AA_fixed$Verdict[i] == 'ACCEPTED', homo_BB_fixed$Verdict[i] == 'ACCEPTED')){
        Verdict[i] = 'ACCEPTED'
      }
    }else if(scenario == "Selection on homodimer AA"){
      Verdict <- homo_AA_fixed$Verdict
    }else if(scenario == "Selection on homodimer BB"){
      Verdict <- homo_BB_fixed$Verdict
    }
  }
  # Load the proposed mutations
  het_proposed <- read.table('all_deltaGs_proposed_heterodimer_all_reps.tab', h = T, sep = '\t')
  het_proposed %<>%
    mutate(Replicate = rep(1:num_reps, each = (nrow(het_proposed)/num_reps)),
           Substitution = rep(0:num_subs, nrow(het_proposed)/(num_subs+1)),
           ComplexType = 'Heterodimer AB',
           Scenario = scenario,
           Beta = beta,
           N = N,
           prev_binding_energy = previous_HET,
           Verdict = Verdict) %>%
    filter(!(is.na(prev_binding_energy))) %>%
    mutate(ddG = Binding_energy - prev_binding_energy)

  
  homo_AA_proposed <- read.table('all_deltaGs_proposed_homodimer_A_all_reps.tab', h = T, sep = '\t')
  homo_AA_proposed %<>%
    mutate(Replicate = rep(1:num_reps, each = (nrow(homo_AA_proposed)/num_reps)),
           Substitution = rep(0:num_subs, nrow(homo_AA_proposed)/(num_subs+1)),
           ComplexType = 'Homodimer AA',
           Scenario = scenario,
           Beta = beta,
           N = N, 
           prev_binding_energy = previous_HM_AA)  %>%
    filter(!(is.na(prev_binding_energy))) %>%
    mutate(ddG = Binding_energy - prev_binding_energy)
  
  homo_BB_proposed <- read.table('all_deltaGs_proposed_homodimer_B_all_reps.tab', h = T, sep = '\t')
  homo_BB_proposed %<>%
    mutate(Replicate = rep(1:num_reps, each = (nrow(homo_BB_proposed)/num_reps)),
           Substitution = rep(0:num_subs, nrow(homo_BB_proposed)/(num_subs+1)),
           ComplexType = 'Homodimer BB',
           Scenario = scenario,
           Beta = beta,
           N = N,
           prev_binding_energy = previous_HM_BB)  %>%
    filter(!(is.na(prev_binding_energy))) %>%
    mutate(ddG = Binding_energy - prev_binding_energy)
  
  outlist$all_ddgs <- het_proposed %>% select(Chain_A_subs_initial, Chain_A_subs_final, Chain_A_subs_position,
                                              Chain_B_subs_initial, Chain_B_subs_final, Chain_B_subs_position,
                                              Replicate, Substitution, Scenario, Beta, N, ddG, Verdict) %>%
    mutate(ddG_HM_AA = homo_AA_proposed$ddG,
           ddG_HM_BB = homo_BB_proposed$ddG)
  
  return(outlist)
  
}

# This function will extract the single mutants from the dataset of ddGs for the scatterplots (figure S13)
extract_single_mutants <- function(df){
  
  # Check the selection scenario
  df_single <- df %>% filter(xor(Chain_A_subs_initial == Chain_A_subs_final, Chain_B_subs_initial == Chain_B_subs_final))
  
  df_single %<>% mutate(het_ddg_be = ddG,
                        homo_ddg_be = ifelse(Chain_A_subs_initial == Chain_A_subs_final, ddG_HM_BB, ddG_HM_AA))
  
}

# This function will prepare the double mutants from the dataset of ddGs for figure S15
extract_double_mutants <- function(df, scenario){
  df_double <- df %>% filter(and(Chain_A_subs_initial != Chain_A_subs_final, Chain_B_subs_initial != Chain_B_subs_final))
  
  df_double %<>% mutate(quadrant = ifelse(and(ddG_HM_AA > 0, ddG_HM_BB > 0), 'Destabilizing both HMs',
                                                                   ifelse(and(ddG_HM_AA <= -0, ddG_HM_BB > 0), 'Opposite direction',
                                                                          ifelse(and(ddG_HM_AA <= -0, ddG_HM_BB <= -0), 'Stabilizing both HMs',
                                                                                 ifelse(and(ddG_HM_AA > 0, ddG_HM_BB <= -0), 'Opposite direction', 'Neutral')))))
  
  accepted_double_sel_hm <- df_double %>% filter(Verdict == 'ACCEPTED')
  
  summary_sel_hm_double_mutants <- df_double %>% group_by(quadrant) %>%
    summarise(total_count = n()) %>% mutate(total_percentage = round(100*total_count / sum(total_count), 2))
  
  summary_accepted_sel_hm_double_mutants <- accepted_double_sel_hm %>% group_by(quadrant) %>%
    summarise(total_count = n()) %>% mutate(total_percentage = round(100*total_count / sum(total_count), 2))
  
  summary_sel_hm_double_mutants$accepted_count <- summary_accepted_sel_hm_double_mutants$total_count
  summary_sel_hm_double_mutants$accepted_percentage <- summary_accepted_sel_hm_double_mutants$total_percentage
  summary_sel_hm_double_mutants$selection <- rep(scenario, nrow(summary_sel_hm_double_mutants))
  
  summary_final <- summary_sel_hm_double_mutants %>% mutate(accepted_percentage_of_total = 100*accepted_count / total_count)
  
  conf_intervals <- c()
  for(entry_num in 1:nrow(summary_final)){
    accepted <- summary_final$accepted_count[entry_num]
    total <- summary_final$total_count[entry_num]
    ratio <- summary_final$accepted_percentage_of_total[entry_num] / 100
    
    a <- binom.test(accepted, total, p = ratio)
    conf_intervals <- rbind(conf_intervals, a$conf.int)
  }
  
  conf_intervals <- as.data.frame(conf_intervals)
  colnames(conf_intervals) <- c('min', 'max')
  summary_final$min <- conf_intervals$min * 100
  summary_final$max <- conf_intervals$max * 100
  
  return(summary_final)
}

#### Use the functions to produce an example of the plots ####

# Load the data
sel_AA_BB <- condensed_energies(<folder_sel_AA_BB>,
                                         50, 200, 'Selection on both homodimers', 10, 1000)
sel_AB <- condensed_energies(<folder_sel_AB>,
                                          50, 200, 'Selection on heterodimer AB', 10, 1000)


# Put them together and plot the trajectories
full_dataset <- rbind(sel_AA, sel_AB)

final_dataset <- full_dataset %>% rowwise() %>%
  mutate(params = paste('Beta = ', toString(Beta), ', N = ', toString(N), sep = ''))

p_all <- final_dataset %>% ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(Scenario ~ params, labeller = label_wrap_gen(width = 30)) +
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("purple","blue", "violet"),
                      name = "")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=20), axis.text.y= element_text(size=20), 
        axis.title.x = element_text(size=25, face = 'bold'), axis.title.y = element_text(size=25, face = 'bold'))+
  theme(
    legend.key = element_rect(fill = "white"),
    legend.text=element_text(), legend.position = "none",
    legend.justification = "center",
    strip.text = element_text(size = 25))+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-55,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=4),
        strip.text.x = element_text(face = 'bold'),
        strip.text.y = element_text(face = 'bold'),
        strip.background = element_rect(fill = 'white')
  ) 


legend <- get_legend(p_all + theme(legend.position = 'top', 
                                   legend.text=element_text(size = 25),
                                   legend.justification = "center"
) +
  guides(color = guide_legend(override.aes = list(size=10))))

fig_s9 <- plot_grid(legend, p_all,
                    rel_heights = c(0.05, 3), ncol = 1)

# fig_s9 <- plot_grid(header_pdb, p_all,
#                     rel_heights = c(0.5, 3), ncol = 1)

# The one for all the parameters
ggsave(filename = <path>,
       width = 28, height = 28, dpi = 500, device = cairo_pdf, plot = fig_s9)

#### Show an example of the scatterplots ####

single_mutants_sel_AA_BB <- extract_single_mutants(sel_AA_BB$all_ddgs)
single_mutants_sel_AB <- extract_single_mutants(sel_AB$all_ddgs)

# Put them together for the plot
full_dataset <- rbind(single_mutants_sel_AA_BB, single_mutants_sel_AB) %>%
  mutate(PDB = '1M38')

# Define limits for points out of bounds

y_min <- -5
y_max <- 10
x_min <- -5
x_max <- 15

p <- full_dataset %>% 
  mutate(out_of_bounds = or(het_ddg_be < y_min, or(het_ddg_be > y_max, or(homo_ddg_be < x_min, homo_ddg_be > x_max))),
         het_ddg_be = ifelse(het_ddg_be < y_min, y_min,
                             ifelse(het_ddg_be > y_max, y_max, het_ddg_be)),
         homo_ddg_be = ifelse(homo_ddg_be < x_min, x_min,
                              ifelse(homo_ddg_be > x_max, x_max, homo_ddg_be)),
         Verdict = gsub(x = Verdict, pattern = 'ACCEPTED', replacement = 'Fixed'),
         Verdict = gsub(x = Verdict, pattern = 'REJECTED', replacement = 'Lost')) %>%
  ggplot(aes(y = het_ddg_be, x = homo_ddg_be, color = Verdict, shape = out_of_bounds)) + 
  facet_grid(Scenario ~ PDB, labeller = label_wrap_gen(width = 15)) +
  scale_color_manual(values = c('red', 'black')) +
  scale_shape_manual(values = c(16, 17), guide = 'none') +
  geom_point(aes(alpha = Verdict)) +
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.line = element_line(),
        axis.text.x= element_text(face = 'bold'), axis.text.y= element_text(face = 'bold'), 
        axis.title.x = element_text(), axis.title.y = element_text(),
        strip.text = element_text(face = 'bold'), 
        strip.background = element_blank())+
  theme(plot.title = element_text(face = "bold",family="Arial", hjust = 0.5),
        legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(), legend.position = "none",
        legend.justification = "center")+
  guides(color = guide_legend(override.aes = list())) +
  ylab("ΔΔG HETs (Kcal/mol)") + xlab("ΔΔG HMs (Kcal/mol)") +
  stat_cor(method = 'pearson', label.x.npc = 0.35,
           label.y.npc = 0.1, show.legend = FALSE,
           inherit.aes = FALSE,
           aes(x = het_ddg_be, y = homo_ddg_be, label = gsub(x = ..label.., replacement = 'r', pattern = 'R')),
           size = 7) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  xlim(x_min, x_max) + ylim(y_min, y_max) +
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=4))
p

legend <- get_legend(p + theme(legend.position = 'top', 
                                            legend.text=element_text()) +
                       guides(color = guide_legend(override.aes = list())))

fig_s13 <- plot_grid(legend, p, rel_heights = c(0.2, 3), ncol = 1)

ggsave(filename = <path>,
       width = 7, height = 14, dpi = 500, device = cairo_pdf, plot = fig_s13)

#### Getting the double mutants ####

# Load the data on sel on AA_BB
double_mutants_sel_AA_BB <- extract_double_mutants(sel_AA_BB$all_ddgs, 'Both homodimers')
double_mutants_sel_AB <- extract_double_mutants(sel_AB$all_ddgs, 'Heterodimer')

summary_final <- rbind(double_mutants_sel_AA_BB, double_mutants_sel_AB)

fig_S15 <- summary_final %>%
  ggplot(aes(x = quadrant, y = accepted_percentage_of_total, fill = selection)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_errorbar(aes(ymin = min, ymax = max), 
                width = 0.5, size = 0.3,
                position = position_dodge(0.9)) +
  scale_fill_manual(values = c('#737373', '#d9d9d9')) + 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold'), axis.title.y = element_text(face = 'bold'),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.5)) +
  xlab("Effect on HMs") + ylab("Fixed mutations (%)") + labs(fill = 'Selection on')
p15

