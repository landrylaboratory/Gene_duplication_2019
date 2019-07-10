#########################################################################################
####                 Tools_analysis_simulations                                      ####
#### This script provides functions that are helpful for the analysis of simulations ####
#### and uses them to produce the following figures from the paper:                  ####
#### - Figure 4                                                                      ####
#### - Figure 4- figure supplements 2 to 4                                           ####
#### - Figure 5                                                                      ####
#### - Figure 5- figure supplements 1 to 3                                           ####
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
condensed_energies <- function(folder, num_reps, num_subs, scenario, beta, N, PDB_structure){
  original_wd <- getwd()
  
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
           N = N, 
           PDB = PDB_structure)
  
  homo_AA_fixed <- read.table('all_deltaGs_homodimer_A_all_reps.tab', h = T, sep = '\t') 
  homo_AA_fixed %<>%
    mutate(Replicate = rep(1:num_reps, each = (nrow(homo_AA_fixed)/num_reps)),
           Substitution = rep(0:num_subs, nrow(homo_AA_fixed)/(num_subs+1)),
           ComplexType = 'Homodimer AA',
           Scenario = scenario,
           Beta = beta,
           N = N,
           PDB = PDB_structure)
  
  homo_BB_fixed <- read.table('all_deltaGs_homodimer_B_all_reps.tab', h = T, sep = '\t')
  homo_BB_fixed %<>%
    mutate(Replicate = rep(1:num_reps, each = (nrow(homo_BB_fixed)/num_reps)),
           Substitution = rep(0:num_subs, nrow(homo_BB_fixed)/(num_subs+1)),
           ComplexType = 'Homodimer BB',
           Scenario = scenario,
           Beta = beta,
           N = N,
           PDB = PDB_structure)
  
  outlist$all_energies <- rbind(het_fixed, homo_AA_fixed, homo_BB_fixed)
  
  # Load the proposed mutations
  het_proposed <- read.table('all_deltaGs_proposed_heterodimer_all_reps.tab', h = T, sep = '\t')
  homo_AA_proposed <- read.table('all_deltaGs_proposed_homodimer_A_all_reps.tab', h = T, sep = '\t')
  homo_BB_proposed <- read.table('all_deltaGs_proposed_homodimer_B_all_reps.tab', h = T, sep = '\t')
  
  previous_HET <- rep(NA, nrow(het_fixed))
  previous_HM_AA <- rep(NA, nrow(het_fixed))
  previous_HM_BB <- rep(NA, nrow(het_fixed))
  verdict_final = rep("REJECTED", nrow(het_fixed))
  
  for(i in 1:nrow(het_fixed)){

    previous_HET[i] <- ifelse(het_fixed$Substitution[i] == 0, NA, het_fixed$Binding_energy[i-1])
    previous_HM_AA[i] <- ifelse(homo_AA_fixed$Substitution[i] == 0, NA, homo_AA_fixed$Binding_energy[i-1])
    previous_HM_BB[i] <- ifelse(homo_BB_fixed$Substitution[i] == 0, NA, homo_BB_fixed$Binding_energy[i-1])
    
    # Recall the selection scenario
    if(scenario == "Selection on HET AB"){
      verdict_final <- het_proposed$Verdict
    }else if(scenario == "Selection on both HMs"){
      if(is.na(homo_AA_proposed$Verdict[i]) || is.na(homo_BB_proposed$Verdict[i])){
        verdict_final[i] = 'REJECTED'
      }
      else if(homo_AA_proposed$Verdict[i] == 'ACCEPTED' && homo_BB_proposed$Verdict[i] == 'ACCEPTED'){
        verdict_final[i] = 'ACCEPTED'
      }
    }else if(scenario == "Selection on HM AA"){
      verdict_final <- homo_AA_proposed$Verdict
    }else if(scenario == "Selection on HM BB"){
      verdict_final <- homo_BB_proposed$Verdict
    }else if(scenario == "Neutral evolution"){
      verdict_final <- rep('ACCEPTED', nrow(het_fixed))
    }
  }
  
  het_proposed %<>%
    mutate(Replicate = rep(1:num_reps, each = (nrow(het_proposed)/num_reps)),
           Substitution = rep(0:num_subs, nrow(het_proposed)/(num_subs+1)),
           ComplexType = 'Heterodimer AB',
           Scenario = scenario,
           Beta = beta,
           N = N,
           prev_binding_energy = previous_HET,
           Verdict = verdict_final,
           PDB = PDB_structure) %>%
    filter(!(is.na(prev_binding_energy))) %>%
    mutate(ddG = Binding_energy - prev_binding_energy)
  
  homo_AA_proposed %<>%
    mutate(Replicate = rep(1:num_reps, each = (nrow(homo_AA_proposed)/num_reps)),
           Substitution = rep(0:num_subs, nrow(homo_AA_proposed)/(num_subs+1)),
           ComplexType = 'Homodimer AA',
           Scenario = scenario,
           Beta = beta,
           N = N, 
           prev_binding_energy = previous_HM_AA,
           Verdict = verdict_final,
           PDB = PDB_structure)  %>%
    filter(!(is.na(prev_binding_energy))) %>%
    mutate(ddG = Binding_energy - prev_binding_energy)
  
  homo_BB_proposed %<>%
    mutate(Replicate = rep(1:num_reps, each = (nrow(homo_BB_proposed)/num_reps)),
           Substitution = rep(0:num_subs, nrow(homo_BB_proposed)/(num_subs+1)),
           ComplexType = 'Homodimer BB',
           Scenario = scenario,
           Beta = beta,
           N = N,
           prev_binding_energy = previous_HM_BB,
           Verdict = verdict_final,
           PDB = PDB_structure)  %>%
    filter(!(is.na(prev_binding_energy))) %>%
    mutate(ddG = Binding_energy - prev_binding_energy)
  
  outlist$all_ddgs <- het_proposed %>% select(Chain_A_subs_initial, Chain_A_subs_final, Chain_A_subs_position,
                                              Chain_B_subs_initial, Chain_B_subs_final, Chain_B_subs_position,
                                              Replicate, Substitution, Scenario, Beta, N, ddG, Verdict, PDB) %>%
    mutate(ddG_HM_AA = homo_AA_proposed$ddG,
           ddG_HM_BB = homo_BB_proposed$ddG,
           Verdict = gsub(pattern = 'ACCEPTED', replacement = 'Fixed', x = Verdict),
           Verdict = gsub(pattern = 'REJECTED', replacement = 'Lost', x = Verdict))
  
  setwd(original_wd)
  
  return(outlist)
  
}

# This function will extract the single mutants used for the figure on effect sizes
extract_single_mutants <- function(df){
  
  # Check the selection scenario
  df_single <- df %>% filter(xor(Chain_A_subs_initial == Chain_A_subs_final, Chain_B_subs_initial == Chain_B_subs_final))
  
  df_single %<>% mutate(het_ddg_be = ddG,
                        homo_ddg_be = ifelse(Chain_A_subs_initial == Chain_A_subs_final, ddG_HM_BB, ddG_HM_AA))
  
}

# A function that prepares the single mutants for the figure on effect sizes
assign_quadrants_single_mutants <- function(df){
  quadrant_labels <- c('Destabilizing both',
                       'Destabilizing HM',
                       'Stabilizing both',
                       'Destabilizing HET')
  
  # Without the neutral label
  df$quadrant <- ifelse(df$het_ddg_be > 0 & df$homo_ddg_be > 0, quadrant_labels[1],
                        ifelse(df$het_ddg_be <= -0 & df$homo_ddg_be > 0, quadrant_labels[2],
                               ifelse(df$het_ddg_be <= -0 & df$homo_ddg_be <= -0, quadrant_labels[3],
                                      ifelse(df$het_ddg_be > 0 & df$homo_ddg_be <= -0, quadrant_labels[4], 'Neutral')
                               )
                               
                        )
                        
  )
  
  df_plot <- df %>%
    select(het_ddg_be, homo_ddg_be, quadrant) %>%
    gather(key = Complex, value = ddg, -quadrant) %>%
    mutate(Complex = gsub(x = Complex, pattern = "het_ddg_be", replacement = "HET"),
           Complex = gsub(x = Complex, pattern = "homo_ddg_be", replacement = "HM"))
  
  df_plot$Complex <- factor(df_plot$Complex, levels = c('HM', 'HET'))
  
  return(df_plot)
}

# This function will extract the double mutants from the dataset of ddGs for the scatterplots
extract_double_mutants <- function(df){
  df_double <- df %>% filter(and(Chain_A_subs_initial != Chain_A_subs_final, Chain_B_subs_initial != Chain_B_subs_final))
  
  return(df_double)
}

# This function works with double mutants to prepare them for the plot on fixation rates
get_fixation_rates <- function(df_double, scenario){
  df_double %<>% mutate(quadrant = ifelse(and(ddG_HM_AA > 0, ddG_HM_BB > 0), 'Destabilizing both HMs',
                                          ifelse(and(ddG_HM_AA <= -0, ddG_HM_BB > 0), 'Opposite direction',
                                                 ifelse(and(ddG_HM_AA <= -0, ddG_HM_BB <= -0), 'Stabilizing both HMs',
                                                        ifelse(and(ddG_HM_AA > 0, ddG_HM_BB <= -0), 'Opposite direction', 'Neutral')))))
  
  accepted_double_sel_hm <- df_double %>% filter(Verdict == 'Fixed')
  
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

# Define the function to get the linear regressions
lm_eqn <- function(df){
  y <- df$ddG
  x <- df$exp_ddg_het
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

# Define a formula to get a p-value based on the z-score of the differences of two binomial distributions
# Discussed here: https://stats.stackexchange.com/questions/113602/test-if-two-binomial-distributions-are-statistically-different-from-each-other
get_pvals_binomial <- function(p1, p2, n1, n2){
  phat <- (n1*p1 + n2*p2)/(n1 + n2)
  z_scores <- (p1 - p2) / (phat*(1 - phat)*(1/n1 + 1/n2))
  pvals <- ifelse(pnorm(-abs(z_scores)) < 2.2e-16, 'p < 2.2e-16', 
                  toString(format(pnorm(-abs(z_scores)), exponential = TRUE)))
  return(pvals)
}

# Define limits for points out of bounds

y_min <- -5
y_max <- 20
x_min <- -5
x_max <- 20

#### Use the functions to produce an example of the plots ####

setwd(<path_to_scripts_for_simulations>)

#### Load the data ####

# Main set of parameters #

# 1A82 #
nosel_1a82 <- condensed_energies('Data/1a82/nosel', 50, 200, 'Neutral evolution', 10, 1000, '1A82')
sel_AA_BB_1a82 <- condensed_energies('Data/1a82/sel_AA_BB/', 50, 200, 'Selection on both HMs', 10, 1000, '1A82')
sel_AB_1a82 <- condensed_energies('Data/1a82/sel_AB', 50, 200, 'Selection on HET AB', 10, 1000, '1A82')
sel_AA_1a82 <- condensed_energies('Data/1a82/sel_AA', 50, 200, 'Selection on HM AA', 10, 1000, '1A82')
sel_BB_1a82 <- condensed_energies('Data/1a82/sel_BB', 50, 200, 'Selection on HM BB', 10, 1000, '1A82')

# 1M38 #
nosel_1m38 <- condensed_energies('Data/1m38/nosel', 50, 200, 'Neutral evolution', 10, 1000, '1M38')
sel_AA_BB_1m38 <- condensed_energies('Data/1m38/sel_AA_BB/', 50, 200, 'Selection on both HMs', 10, 1000, '1M38')
sel_AB_1m38 <- condensed_energies('Data/1m38/sel_AB', 50, 200, 'Selection on HET AB', 10, 1000, '1M38')
sel_AA_1m38 <- condensed_energies('Data/1m38/sel_AA', 50, 200, 'Selection on HM AA', 10, 1000, '1M38')
sel_BB_1m38 <- condensed_energies('Data/1m38/sel_BB', 50, 200, 'Selection on HM BB', 10, 1000, '1M38')

# 2JKY #
nosel_2jky <- condensed_energies('Data/2jky/nosel', 50, 200, 'Neutral evolution', 10, 1000, '2JKY')
sel_AA_BB_2jky <- condensed_energies('Data/2jky/sel_AA_BB/', 50, 200, 'Selection on both HMs', 10, 1000, '2JKY')
sel_AB_2jky <- condensed_energies('Data/2jky/sel_AB', 50, 200, 'Selection on HET AB', 10, 1000, '2JKY')
sel_AA_2jky <- condensed_energies('Data/2jky/sel_AA', 50, 200, 'Selection on HM AA', 10, 1000, '2JKY')
sel_BB_2jky <- condensed_energies('Data/2jky/sel_BB', 50, 200, 'Selection on HM BB', 10, 1000, '2JKY')

# 2O1V #
nosel_2o1v <- condensed_energies('Data/2o1v/nosel', 50, 200, 'Neutral evolution', 10, 1000, '2O1V')
sel_AA_BB_2o1v <- condensed_energies('Data/2o1v/sel_AA_BB/', 50, 200, 'Selection on both HMs', 10, 1000, '2O1V')
sel_AB_2o1v <- condensed_energies('Data/2o1v/sel_AB', 50, 200, 'Selection on HET AB', 10, 1000, '2O1V')
sel_AA_2o1v <- condensed_energies('Data/2o1v/sel_AA', 50, 200, 'Selection on HM AA', 10, 1000, '2O1V')
sel_BB_2o1v <- condensed_energies('Data/2o1v/sel_BB', 50, 200, 'Selection on HM BB', 10, 1000, '2O1V')

# 3D8X #
nosel_3d8x <- condensed_energies('Data/3d8x/nosel', 50, 200, 'Neutral evolution', 10, 1000, '3D8X')
sel_AA_BB_3d8x <- condensed_energies('Data/3d8x/sel_AA_BB/', 50, 200, 'Selection on both HMs', 10, 1000, '3D8X')
sel_AB_3d8x <- condensed_energies('Data/3d8x/sel_AB', 50, 200, 'Selection on HET AB', 10, 1000, '3D8X')
sel_AA_3d8x <- condensed_energies('Data/3d8x/sel_AA', 50, 200, 'Selection on HM AA', 10, 1000, '3D8X')
sel_BB_3d8x <- condensed_energies('Data/3d8x/sel_BB', 50, 200, 'Selection on HM BB', 10, 1000, '3D8X')

# 4FGW #
nosel_4fgw <- condensed_energies('Data/4fgw/nosel', 50, 200, 'Neutral evolution', 10, 1000, '4FGW')
sel_AA_BB_4fgw <- condensed_energies('Data/4fgw/sel_AA_BB/', 50, 200, 'Selection on both HMs', 10, 1000, '4FGW')
sel_AB_4fgw <- condensed_energies('Data/4fgw/sel_AB', 50, 200, 'Selection on HET AB', 10, 1000, '4FGW')
sel_AA_4fgw <- condensed_energies('Data/4fgw/sel_AA', 50, 200, 'Selection on HM AA', 10, 1000, '4FGW')
sel_BB_4fgw <- condensed_energies('Data/4fgw/sel_BB', 50, 200, 'Selection on HM BB', 10, 1000, '4FGW')

#### Figure 4 ####

# Put them together to plot the trajectories
full_dataset <- rbind(nosel_1a82$all_energies, sel_AA_BB_1a82$all_energies, sel_AB_1a82$all_energies, sel_AA_1a82$all_energies, sel_BB_1a82$all_energies,
                      nosel_1m38$all_energies, sel_AA_BB_1m38$all_energies, sel_AB_1m38$all_energies, sel_AA_1m38$all_energies, sel_BB_1m38$all_energies,
                      nosel_2jky$all_energies, sel_AA_BB_2jky$all_energies, sel_AB_2jky$all_energies, sel_AA_2jky$all_energies, sel_BB_2jky$all_energies,
                      nosel_2o1v$all_energies, sel_AA_BB_2o1v$all_energies, sel_AB_2o1v$all_energies, sel_AA_2o1v$all_energies, sel_BB_2o1v$all_energies,
                      nosel_3d8x$all_energies, sel_AA_BB_3d8x$all_energies, sel_AB_3d8x$all_energies, sel_AA_3d8x$all_energies, sel_BB_3d8x$all_energies,
                      nosel_4fgw$all_energies, sel_AA_BB_4fgw$all_energies, sel_AB_4fgw$all_energies, sel_AA_4fgw$all_energies, sel_BB_4fgw$all_energies
) %>% filter(Substitution != 0) %>%
  mutate(ComplexType = gsub(pattern = 'Heterodimer AB', replacement = 'HET AB', x = ComplexType),
         ComplexType = gsub(pattern = 'Homodimer AA', replacement = 'HM AA', x = ComplexType),
         ComplexType = gsub(pattern = 'Homodimer BB', replacement = 'HM BB', x = ComplexType))

# Load the diagram for the simulations
fig_4A <- ggdraw() + draw_image('Data/Figure4-A.png')

full_dataset$ComplexType = factor(full_dataset$ComplexType, levels = c('HM AA', 'HET AB', 'HM BB'))

# Neutral evolution #

p_1M38_nosel <- full_dataset %>% 
  filter(PDB == '1M38', Scenario == 'Neutral evolution') %>%
  ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(~Scenario) + 
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "pink"),
                      name = "Complex")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = 'bold'), axis.title.y = element_text(size=15, face = 'bold'))+
  theme(legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "top",
        legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-35,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  )

# Selection on both HMs #

p_1M38_sel_AA_BB <- full_dataset %>% 
  filter(PDB == '1M38', Scenario == 'Selection on both HMs') %>%
  ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(~Scenario) +
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "pink"),
                      name = "Complex")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = 'bold'), axis.title.y = element_text(size=15, face = 'bold'))+
  theme(legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "top",
        legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-35,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  )

# Selection on HET #

p_1M38_sel_AB <- full_dataset %>% 
  filter(PDB == '1M38', Scenario == 'Selection on HET AB') %>%
  ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(~Scenario) + 
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "pink"),
                      name = "Complex")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = 'bold'), axis.title.y = element_text(size=15, face = 'bold'))+
  theme(legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "top",
        legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-35,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  )

# Selection on HM AA #

p_1M38_sel_AA <- full_dataset %>% 
  filter(PDB == '1M38', Scenario == 'Selection on HM AA') %>%
  ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(~Scenario) +
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "pink"),
                      name = "Complex")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = 'bold'), axis.title.y = element_text(size=15, face = 'bold'))+
  theme(legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "top",
        legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-35,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  )

# Selection on HM BB #

p_1M38_sel_BB <- full_dataset %>% 
  filter(PDB == '1M38', Scenario == 'Selection on HM BB') %>%
  ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(~Scenario) +
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "pink"),
                      name = "Complex")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = 'bold'), axis.title.y = element_text(size=15, face = 'bold'))+
  theme(legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "top",
        legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-35,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  )

# Save the panels for figure 4 #
final_fig4 <- plot_grid(fig_4A, p_1M38_nosel, p_1M38_sel_AA_BB,
                        p_1M38_sel_AA, p_1M38_sel_BB,p_1M38_sel_AB, 
                        nrow = 3, 
                        labels = c('A', 'B', 'C', 'D', 'E', 'F'))

ggsave(filename = 'Figures/Figure4.pdf',
       width = 14, height = 21, dpi = 500, plot = final_fig4, device = cairo_pdf)
ggsave(filename = 'Figures/Figure4.png',
       width = 14, height = 21, dpi = 300, plot = final_fig4, device = 'png')

#### Figure 5 ####

#### Getting the double mutants ####

double_mutants_sel_AA_BB <- extract_double_mutants(sel_AA_BB_1m38$all_ddgs)
double_mutants_sel_AB <- extract_double_mutants(sel_AB_1m38$all_ddgs) 

#### Plot deviation of observations from expectations ####

double_mutants_sel_AA_BB %<>% 
  mutate(exp_ddg_het = (ddG_HM_AA + ddG_HM_BB)/2)

double_mutants_sel_AB %<>%
  mutate(exp_ddg_het = (ddG_HM_AA + ddG_HM_BB)/2)

# Relevel to show red on top
double_mutants_sel_AA_BB$Verdict <- factor(double_mutants_sel_AA_BB$Verdict, levels = c('Lost', 'Fixed'))
double_mutants_sel_AB$Verdict <- factor(double_mutants_sel_AB$Verdict, levels = c('Lost', 'Fixed'))

fig_5A <- double_mutants_sel_AA_BB %>%
  mutate(out_of_bounds = ifelse(or(exp_ddg_het > x_max,
                                   or(exp_ddg_het < x_min,
                                      or(ddG > y_max, ddG < y_min))),
                                1, 0)) %>%
  mutate(exp_ddg_het = ifelse(exp_ddg_het > x_max, x_max,
                              ifelse(exp_ddg_het < x_min, x_min, exp_ddg_het)),
         ddG = ifelse(ddG > y_max, y_max, 
                      ifelse(ddG < y_min, y_min, ddG)),
         out_of_bounds = as.factor(out_of_bounds)
  ) %>%
  ggplot(aes(x = exp_ddg_het, y = ddG)) +
  facet_grid(~Scenario) +
  geom_point(aes(colour = Verdict, alpha = Verdict, shape = out_of_bounds)) +
  scale_shape_manual(values = c(16, 17), guide = 'none') +
  scale_colour_manual(values = c("black", "red")) +
  scale_alpha_manual(values = c(0.2, 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5),
        legend.title = element_blank(),
        legend.position = "top", 
        legend.justification = "center",
        legend.text=element_text(size=15),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = "bold"), axis.title.y = element_text(size=15, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  xlab('Expected HET ΔΔG (kcal/mol)') + ylab('Observed HET ΔΔG (kcal/mol)') + labs(fill = '') +
  xlim(x_min, x_max) + ylim(y_min, y_max) +
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_smooth(method = 'lm', aes(group = Verdict, colour = Verdict),  show.legend = FALSE) +
  annotate(geom = 'text', x = x_max - 5, y = y_min + 2, label = lm_eqn(double_mutants_sel_AA_BB %>% filter(Verdict == 'Fixed')),
           parse = TRUE, colour = 'red') +
  annotate(geom = 'text', x = x_max - 5, y = y_min + 3, label = lm_eqn(double_mutants_sel_AA_BB %>% filter(Verdict == 'Lost')),
           parse = TRUE, colour = 'black')
fig_5A

fig_5B <- double_mutants_sel_AB %>%
  mutate(out_of_bounds = ifelse(or(exp_ddg_het > x_max,
                                   or(exp_ddg_het < x_min,
                                      or(ddG > y_max, ddG < y_min))),
                                1, 0)) %>%
  mutate(exp_ddg_het = ifelse(exp_ddg_het > x_max, x_max,
                              ifelse(exp_ddg_het < x_min, x_min, exp_ddg_het)),
         ddG = ifelse(ddG > y_max, y_max, 
                      ifelse(ddG < y_min, y_min, ddG)),
         out_of_bounds = as.factor(out_of_bounds)
  ) %>%
  ggplot(aes(x = exp_ddg_het, y = ddG)) +
  facet_grid(~Scenario) +
  geom_point(aes(colour = Verdict, alpha = Verdict, shape = out_of_bounds)) +
  scale_shape_manual(values = c(16, 17), guide = 'none') +
  scale_colour_manual(values = c("black", "red")) +
  scale_alpha_manual(values = c(0.2, 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5),
        legend.title = element_blank(),
        legend.position = "top", 
        legend.justification = "center",
        legend.text=element_text(size=15),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = "bold"), axis.title.y = element_text(size=15, face = "bold")) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  xlab('Expected HET ΔΔG (kcal/mol)') + ylab('Observed HET ΔΔG (kcal/mol)') + labs(fill = '') +
  xlim(x_min, x_max) + ylim(y_min, y_max) +
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=1),
        strip.text.x = element_text(size = 15, face = 'bold'),
        strip.text.y = element_text(size = 15, face = 'bold'),
        strip.background = element_rect(fill = 'white')) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_smooth(method = 'lm', aes(group = Verdict, colour = Verdict),  show.legend = FALSE) +
  annotate(geom = 'text', x = x_max - 5, y = y_min + 2, label = lm_eqn(double_mutants_sel_AB %>% filter(Verdict == 'Fixed')),
           parse = TRUE, colour = 'red') +
  annotate(geom = 'text', x = x_max - 5, y = y_min + 3, label = lm_eqn(double_mutants_sel_AB %>% filter(Verdict == 'Lost')),
           parse = TRUE, colour = 'black')
fig_5B

fig_5 <- plot_grid(fig_5A, fig_5B, ncol = 2, labels = c('A', 'B'))

ggsave(filename = 'Figures/Figure5.pdf',
       width = 14, height = 7, dpi = 500, device = cairo_pdf, plot = fig_5)
ggsave(filename = 'Figures/Figure5.png',
       width = 14, height = 7, dpi = 300, device = 'png', plot = fig_5)

# Test for the difference in the magnitude of deviations from expectations between 
# mutations fixed under selection for the HET and under selction for the two HMs
double_mutants_sel_AA_BB %<>% mutate(diff_exp_obs = ddG - exp_ddg_het)
double_mutants_sel_AB %<>% mutate(diff_exp_obs = ddG - exp_ddg_het)

data_fig_5C <- rbind(double_mutants_sel_AA_BB, double_mutants_sel_AB)

comparison_list <- list(c('Selection on both HMs', 'Selection on HET AB'))

# Violin plot with all comparisons
fig_5C <- data_fig_5C %>%
  filter(Verdict == 'Fixed') %>%
  ggplot(aes(x = Scenario, 
             y = diff_exp_obs, fill = Scenario)) + 
  facet_wrap(~Verdict) +
  scale_fill_manual(values = c('gray', 'gray')) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
  geom_violin(position = position_dodge(), alpha = 0.5) + 
  stat_compare_means(method = 't.test', comparisons = comparison_list, paired = FALSE,
                     position = 'top', 
                     label.y = 12.5) +
  theme(legend.position = 'none',
        axis.title = element_text(face = 'bold')) +
  ylab('Observed ΔΔG - Expected ΔΔG') + xlab('Selection on') +
  scale_x_discrete(labels = c('Both HMs', 'HET AB')) +
  ylim(-20, 14)
fig_5C


#### Figure 4 - figure supplement 2 ####

# Load the figures of the PDB structures
fig_1a82 <- ggdraw() + draw_image('Data/1A82_new.png')
fig_2o1v <- ggdraw() + draw_image('Data/2O1V_new.png')
fig_1m38 <- ggdraw() + draw_image('Data/1M38_new.png')
fig_4fgw <- ggdraw() + draw_image('Data/4FGW_new.png')
fig_3d8x <- ggdraw() + draw_image('Data/3D8X_new.png')
fig_2jky <- ggdraw() + draw_image('Data/4FGW_new.png')

# Put them together and plot the trajectories (data loaded for figure 4)
p_all <- full_dataset %>% ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(Scenario ~ PDB, labeller = label_wrap_gen(width = 30)) +
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "violet"),
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

header_pdb <- plot_grid(NULL, fig_1a82, fig_1m38, fig_2jky, fig_2o1v, fig_3d8x, fig_4fgw, NULL, nrow = 1,
                        rel_widths = c(0.2, 1, 1, 1, 1, 1, 1, 0.1))

fig4_suppl2 <- plot_grid(header_pdb, legend, p_all,
                     rel_heights = c(0.5, 0.05, 3), ncol = 1)

ggsave(filename = 'Figures/Figure 4-figure_supplement_2.pdf',
       width = 28, height = 28, dpi = 500, device = cairo_pdf, plot = fig4_suppl2)
ggsave(filename = 'Figures/Figure 4-figure_supplement_2.png',
       width = 28, height = 28, dpi = 300, device = 'png', plot = fig4_suppl2)


#### Figure 4 - figure supplement 3 ####

# Load results from simulations with new parameters

# 1M38, beta = 1, N = 1000 #
sel_AA_BB_1m38_b1_n1000 <- condensed_energies('Data/1m38_other_params/1m38_b1_n1000/sel_AA_BB/', 50, 200, 'Selection on both HMs', 1, 1000, '1M38')
sel_AB_1m38_b1_n1000 <- condensed_energies('Data/1m38_other_params/1m38_b1_n1000/sel_AB', 50, 200, 'Selection on HET AB', 1, 1000, '1M38')
sel_AA_1m38_b1_n1000 <- condensed_energies('Data/1m38_other_params/1m38_b1_n1000/sel_AA', 50, 200, 'Selection on HM AA', 1, 1000, '1M38')
sel_BB_1m38_b1_n1000 <- condensed_energies('Data/1m38_other_params/1m38_b1_n1000/sel_BB', 50, 200, 'Selection on HM BB', 1, 1000, '1M38')

# 1M38, beta = 10, N = 100 #
sel_AA_BB_1m38_b10_n100 <- condensed_energies('Data/1m38_other_params/1m38_b10_n100/sel_AA_BB/', 50, 200, 'Selection on both HMs', 10, 100, '1M38')
sel_AB_1m38_b10_n100 <- condensed_energies('Data/1m38_other_params/1m38_b10_n100/sel_AB', 50, 200, 'Selection on HET AB', 10, 100, '1M38')
sel_AA_1m38_b10_n100 <- condensed_energies('Data/1m38_other_params/1m38_b10_n100/sel_AA', 50, 200, 'Selection on HM AA', 10, 100, '1M38')
sel_BB_1m38_b10_n100 <- condensed_energies('Data/1m38_other_params/1m38_b10_n100/sel_BB', 50, 200, 'Selection on HM BB', 10, 100, '1M38')

# 1M38, beta = 10, N = 1000, 500 subs #
sel_AA_BB_1m38_500_subs <- condensed_energies('Data/1m38_other_params/1m38_b10_n1000_500_subs/sel_AA_BB/', 50, 500, 'Selection on both HMs', 10.1, 1000, '1M38')
sel_AB_1m38_500_subs <- condensed_energies('Data/1m38_other_params/1m38_b10_n1000_500_subs/sel_AB', 50, 500, 'Selection on HET AB', 10.1, 1000, '1M38')
sel_AA_1m38_500_subs <- condensed_energies('Data/1m38_other_params/1m38_b10_n1000_500_subs/sel_AA', 50, 500, 'Selection on HM AA', 10.1, 1000, '1M38')
sel_BB_1m38_500_subs <- condensed_energies('Data/1m38_other_params/1m38_b10_n1000_500_subs/sel_BB', 50, 500, 'Selection on HM BB', 10.1, 1000, '1M38')

# 1M38, beta = 10, N = 10000 #
sel_AA_BB_1m38_b10_n10000 <- condensed_energies('Data/1m38_other_params/1m38_b10_n10000/sel_AA_BB/', 50, 200, 'Selection on both HMs', 10, 10000, '1M38')
sel_AB_1m38_b10_n10000 <- condensed_energies('Data/1m38_other_params/1m38_b10_n10000/sel_AB', 50, 200, 'Selection on HET AB', 10, 10000, '1M38')
sel_AA_1m38_b10_n10000 <- condensed_energies('Data/1m38_other_params/1m38_b10_n10000/sel_AA', 50, 200, 'Selection on HM AA', 10, 10000, '1M38')
sel_BB_1m38_b10_n10000 <- condensed_energies('Data/1m38_other_params/1m38_b10_n10000/sel_BB', 50, 200, 'Selection on HM BB', 10, 10000, '1M38')

# 1M38, beta = 20, N = 1000 #
sel_AA_BB_1m38_b20_n1000 <- condensed_energies('Data/1m38_other_params/1m38_b20_n1000/sel_AA_BB/', 50, 200, 'Selection on both HMs', 20, 1000, '1M38')
sel_AB_1m38_b20_n1000 <- condensed_energies('Data/1m38_other_params/1m38_b20_n1000/sel_AB', 50, 200, 'Selection on HET AB', 20, 1000, '1M38')
sel_AA_1m38_b20_n1000 <- condensed_energies('Data/1m38_other_params/1m38_b20_n1000/sel_AA', 50, 200, 'Selection on HM AA', 20, 1000, '1M38')
sel_BB_1m38_b20_n1000 <- condensed_energies('Data/1m38_other_params/1m38_b20_n1000/sel_BB', 50, 200, 'Selection on HM BB', 20, 1000, '1M38')

final_dataset <- rbind(
  sel_AA_BB_1m38_b1_n1000$all_energies, sel_AB_1m38_b1_n1000$all_energies, sel_AA_1m38_b1_n1000$all_energies, sel_BB_1m38_b1_n1000$all_energies,
  sel_AA_BB_1m38_b10_n100$all_energies, sel_AB_1m38_b10_n100$all_energies, sel_AA_1m38_b10_n100$all_energies, sel_BB_1m38_b10_n100$all_energies,
  sel_AA_BB_1m38_500_subs$all_energies, sel_AB_1m38_500_subs$all_energies, sel_AA_1m38_500_subs$all_energies, sel_BB_1m38_500_subs$all_energies,
  sel_AA_BB_1m38_b10_n10000$all_energies, sel_AB_1m38_b10_n10000$all_energies, sel_AA_1m38_b10_n10000$all_energies, sel_BB_1m38_b10_n10000$all_energies,
  sel_AA_BB_1m38_b20_n1000$all_energies, sel_AB_1m38_b20_n1000$all_energies, sel_AA_1m38_b20_n1000$all_energies, sel_BB_1m38_b20_n1000$all_energies) %>% 
  filter(Substitution != 0)
  
final_dataset %<>% rowwise() %>%
  mutate(params = ifelse(Beta == 10.1, paste('Beta = 10', ', N = ', toString(N), ", 500 substitutions", sep = ''),
                         paste('Beta = ', toString(Beta), ', N = ', toString(N), sep = ''))) %>%
  mutate(Beta = ifelse(Beta == 10.1, 10, Beta)) %>%
  mutate(ComplexType = gsub(pattern = 'Heterodimer AB', replacement = 'HET AB', x = ComplexType),
         ComplexType = gsub(pattern = 'Homodimer AA', replacement = 'HM AA', x = ComplexType),
         ComplexType = gsub(pattern = 'Homodimer BB', replacement = 'HM BB', x = ComplexType))

final_dataset$ComplexType = factor(final_dataset$ComplexType, levels = c('HM AA', 'HET AB', 'HM BB'))

p_all <- final_dataset %>% 
  ggplot(aes(x=Substitution,y=Binding_energy, group=interaction(Replicate, ComplexType), color=ComplexType)) + 
  facet_grid(Scenario ~ params, labeller = label_wrap_gen(width = 30), scales = "free") +
  geom_line(size=1, alpha=0.1) +
  scale_colour_manual(values = c("blue", "purple", "violet"),
                      name = "")+
  scale_alpha_manual(values = c(1, 0.3)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), 
        axis.title.x = element_text(size=20, face = 'bold'), axis.title.y = element_text(size=20, face = 'bold'))+
  theme(
    legend.key = element_rect(fill = "white"),
    legend.text=element_text(size=20), legend.position = "none",
    legend.justification = "center")+
  ylab("Binding energy (Kcal/mol)") + xlab("Time (substitutions)") +
  stat_summary(aes(group=ComplexType), fun.y=mean, geom="line", size=2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  ylim(c(-55,20))+
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size=4),
        strip.text.x = element_text(size = 20, face = 'bold'),
        strip.text.y = element_text(size = 20, face = 'bold'),
        strip.background = element_rect(fill = 'white')
  ) 


legend <- get_legend(p_all + theme(legend.position = 'top', 
                                   legend.text=element_text(size = 25),
                                   legend.justification = "center"
) +
  guides(color = guide_legend(override.aes = list(size=10))))

fig4_suppl3 <- plot_grid(legend, p_all,
                    rel_heights = c(0.05, 3), ncol = 1)

ggsave(filename = 'Figures/Figure 4-figure_supplement_3.pdf',
       width = 21, height = 21, dpi = 500, device = cairo_pdf, plot = fig4_suppl3)
ggsave(filename = 'Figures/Figure 4-figure_supplement_3.png',
       width = 21, height = 21, dpi = 300, device = 'png', plot = fig4_suppl3)

#### Figure 4 - figure supplement 4 ####

# Put the data together (loaded for figure 4)
full_dataset <- rbind(sel_AA_BB_1a82$all_ddgs, sel_AB_1a82$all_ddgs,
                      sel_AA_BB_1m38$all_ddgs, sel_AB_1m38$all_ddgs,
                      sel_AA_BB_2jky$all_ddgs, sel_AB_2jky$all_ddgs,
                      sel_AA_BB_2o1v$all_ddgs, sel_AB_2o1v$all_ddgs,
                      sel_AA_BB_3d8x$all_ddgs, sel_AB_3d8x$all_ddgs,
                      sel_AA_BB_4fgw$all_ddgs, sel_AB_4fgw$all_ddgs)

single_mutants_all <- extract_single_mutants(full_dataset)

single_mutants_all$Verdict <- factor(single_mutants_all$Verdict, levels = c('Lost', 'Fixed'))

p_all_facet <- single_mutants_all %>% 
  mutate(out_of_bounds = ifelse(or(homo_ddg_be > x_max,
                                   or(homo_ddg_be < x_min,
                                      or(het_ddg_be > y_max, het_ddg_be < y_min))),
                                1, 0)) %>%
  mutate(homo_ddg_be = ifelse(homo_ddg_be > x_max, x_max,
                              ifelse(homo_ddg_be < x_min, x_min, homo_ddg_be)),
         ddG = ifelse(het_ddg_be > y_max, y_max, 
                      ifelse(het_ddg_be < y_min, y_min, het_ddg_be)),
         out_of_bounds = as.factor(out_of_bounds)
  ) %>%
  ggplot(aes(y = het_ddg_be, x = homo_ddg_be, color = Verdict, shape = out_of_bounds)) + 
  facet_grid(Scenario ~ PDB, labeller = label_wrap_gen(width = 15)) +
  scale_color_manual(values = c('black', 'red')) +
  scale_shape_manual(values = c(16, 17), guide = 'none') +
  geom_point(aes(alpha = Verdict), size = 2) +
  scale_alpha_manual(values = c(0.2, 1)) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), 
        axis.title.x = element_text(size=20, face = "bold"), axis.title.y = element_text(size=20, face = "bold"))+
  theme(plot.title = element_text(size = 20, face = "bold",family="Arial", hjust = 0.5),
        legend.key = element_rect(fill = "white"), legend.title = element_blank(),
        legend.text=element_text(size=12), legend.position = "none",
        legend.justification = "center")+
  guides(color = guide_legend(override.aes = list(size=3))) +
  ylab("ΔΔG HET (Kcal/mol)") + xlab("ΔΔG HM (Kcal/mol)") + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  xlim(x_min, x_max) + ylim(y_min, y_max) +
  stat_cor(method = 'pearson', label.x.npc = 0.35,
           label.y.npc = 0.1, show.legend = FALSE,
           inherit.aes = FALSE,
           aes(x = het_ddg_be, y = homo_ddg_be, label = gsub(x = ..label.., replacement = 'r', pattern = 'R')),
           size = 6
  ) +
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size = 1),
        strip.text.x = element_text(size = 20, face = 'bold'),
        strip.text.y = element_text(size = 20, face = 'bold'),
        strip.background = element_rect(fill = 'white')) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed')

header_pdb <- plot_grid(NULL, fig_1a82, fig_1m38, fig_2jky, fig_2o1v, fig_3d8x, fig_4fgw, NULL, nrow = 1,
                        rel_widths = c(0.2, 1, 1, 1, 1, 1, 1, 0.1))

legend <- get_legend(p_all_facet + theme(legend.position = 'top', 
                                         legend.text=element_text(size=20)) +
                       guides(color = guide_legend(override.aes = list(size=4))))

fig4_suppl4 <- plot_grid(header_pdb, legend, p_all_facet, rel_heights = c(1, 0.2, 4), ncol = 1)

ggsave(filename = 'Figures/Figure 4-figure_supplement_4.pdf',
       width = 28, height = 14, dpi = 500, device = cairo_pdf, plot = fig4_suppl4)
ggsave(filename = 'Figures/Figure 4-figure_supplement_4.png',
       width = 28, height = 14, dpi = 300, device = 'png', plot = fig4_suppl4)

#### Figure 5 - figure supplement 1 #### 

# Put the data together (loaded for figure 4)
full_dataset <- rbind(sel_AA_BB_1a82$all_ddgs, sel_AB_1a82$all_ddgs,
                      sel_AA_BB_1m38$all_ddgs, sel_AB_1m38$all_ddgs,
                      sel_AA_BB_2jky$all_ddgs, sel_AB_2jky$all_ddgs,
                      sel_AA_BB_2o1v$all_ddgs, sel_AB_2o1v$all_ddgs,
                      sel_AA_BB_3d8x$all_ddgs, sel_AB_3d8x$all_ddgs,
                      sel_AA_BB_4fgw$all_ddgs, sel_AB_4fgw$all_ddgs)

double_mutants_all <- extract_double_mutants(full_dataset) 

double_mutants_all %<>% 
  mutate(exp_ddg_het = (ddG_HM_AA + ddG_HM_BB)/2)

double_mutants_all %<>%
  mutate(exp_ddg_het = (ddG_HM_AA + ddG_HM_BB)/2)

# Prepare the labels for the equation and the R squared

labels_fixed_hm <- as.data.frame(cbind(as.character(unique(double_mutants_all$PDB)),
                                       rep('Fixed', 6), rep('Selection on both HMs', 6)))
labels_lost_hm <- as.data.frame(cbind(as.character(unique(double_mutants_all$PDB)),
                                      rep('Lost', 6), rep('Selection on both HMs', 6)))
labels_list_fixed <- rep('A', 6)
labels_list_lost <- rep('A', 6)

colnames(labels_fixed_hm) <- c('PDB', 'Type_mutation', 'Scenario')
colnames(labels_lost_hm) <- c('PDB', 'Type_mutation', 'Scenario')

for(i in 1:nrow(labels_fixed_hm)){
  labels_list_fixed[i] <- lm_eqn(double_mutants_all %>% filter(Verdict == 'Fixed', Scenario == 'Selection on both HMs', PDB == labels_fixed_hm$PDB[i]))
  labels_list_lost[i] <- lm_eqn(double_mutants_all %>% filter(Verdict == 'Lost', Scenario == 'Selection on both HMs', PDB == labels_lost_hm$PDB[i]))
}

labels_fixed_hm$Label <- labels_list_fixed
labels_lost_hm$Label <- labels_list_lost

labels_fixed_het <- as.data.frame(cbind(as.character(unique(double_mutants_all$PDB)),
                                        rep('Fixed', 6), rep('Selection on HET AB', 6)))
labels_lost_het <- as.data.frame(cbind(as.character(unique(double_mutants_all$PDB)),
                                       rep('Lost', 6), rep('Selection on HET AB', 6)))
labels_list_fixed <- rep('A', 6)
labels_list_lost <- rep('A', 6)

colnames(labels_fixed_het) <- c('PDB', 'Type_mutation', 'Scenario')
colnames(labels_lost_het) <- c('PDB', 'Type_mutation', 'Scenario')

for(i in 1:nrow(labels_fixed_het)){
  labels_list_fixed[i] <- lm_eqn(double_mutants_all %>% filter(Verdict == 'Fixed', Scenario == 'Selection on HET AB', PDB == labels_fixed_het$PDB[i]))
  labels_list_lost[i] <- lm_eqn(double_mutants_all %>% filter(Verdict == 'Lost', Scenario == 'Selection on HET AB', PDB == labels_lost_het$PDB[i]))
  
}

labels_fixed_het$Label <- labels_list_fixed
labels_lost_het$Label <- labels_list_lost

labels_fixed <- rbind(labels_fixed_het, labels_fixed_hm)
labels_lost <- rbind(labels_lost_het, labels_lost_hm)

double_mutants_all$Verdict <- factor(double_mutants_all$Verdict, levels = c('Lost', 'Fixed'))

p_facet_double_mutants <- double_mutants_all %>%
  mutate(out_of_bounds = ifelse(or(exp_ddg_het > x_max,
                                   or(exp_ddg_het < x_min,
                                      or(ddG > y_max, ddG < y_min))),
                                1, 0)) %>%
  mutate(exp_ddg_het = ifelse(exp_ddg_het > x_max, x_max,
                              ifelse(exp_ddg_het < x_min, x_min, exp_ddg_het)),
         ddG = ifelse(ddG > y_max, y_max, 
                         ifelse(ddG < y_min, y_min, ddG)),
         out_of_bounds = as.factor(out_of_bounds)
  ) %>%
  ggplot(aes(x = exp_ddg_het, y = ddG)) +
  facet_grid(Scenario ~ PDB) +
  geom_point(aes(colour = Verdict, alpha = Verdict, shape = out_of_bounds)) +
  scale_shape_manual(values = c(16, 17), guide = 'none') +
  scale_colour_manual(values = c("black", "red")) +
  scale_alpha_manual(values = c(0.2, 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        legend.text=element_text(size=15),
        axis.text.x= element_text(size=12), axis.text.y= element_text(size=12), 
        axis.title.x = element_text(size=15, face = "bold"), axis.title.y = element_text(size=15, face = "bold"),
        axis.line = element_line(size = 0.5),
        legend.title = element_blank(),
        legend.position = "none", 
        legend.justification = "center") +
  xlab('Expected HET ΔΔG (Kcal/mol)') + ylab('Observed HET ΔΔG (Kcal/mol)') + labs(fill = '') +
  xlim(x_min, x_max) + ylim(y_min, y_max) +
  theme(panel.border = element_rect(linetype = "solid", colour = "gray50", size = 1),
        strip.text.x = element_text(size = 20, face = 'bold'),
        strip.text.y = element_text(size = 20, face = 'bold'),
        strip.background = element_rect(fill = 'white')) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50"),
        axis.text.x= element_text(size=15), axis.text.y= element_text(size=15), 
        axis.title.x = element_text(size=20, face = "bold"), axis.title.y = element_text(size=20, face = "bold")) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_smooth(method = 'lm', aes(group = Verdict, colour = Verdict),  show.legend = FALSE) +
  geom_text(data = labels_fixed, x = x_max - 8, y = y_min + 2, aes(label = Label), parse = TRUE, colour = 'red') +
  geom_text(data = labels_lost, x = x_max - 8, y = y_min + 3, aes(label = Label), parse = TRUE, colour = 'black')

legend <- get_legend(p_facet_double_mutants + theme(legend.position = 'top', 
                                                    legend.text=element_text(size=20)) +
                       guides(color = guide_legend(override.aes = list(size=4))))

header_pdb <- plot_grid(NULL, fig_1a82, fig_1m38, fig_2jky, fig_2o1v, fig_3d8x, fig_4fgw, NULL, nrow = 1,
                        rel_widths = c(0.2, 1, 1, 1, 1, 1, 1, 0.1))

fig5_suppl1 <- plot_grid(header_pdb, legend, 
                     p_facet_double_mutants, rel_heights = c(0.8, 0.2, 3), ncol = 1)

ggsave(filename = 'Figures/Figure 5-figure_supplement_1.pdf',
       width = 28, height = 14, dpi = 500, device = cairo_pdf, plot = fig5_suppl1)
ggsave(filename = 'Figures/Figure 5-figure_supplement_1.png',
       width = 28, height = 14, dpi = 300, device = 'png', plot = fig5_suppl1)

#### Figure 5 - figure supplement 2 ####

data_fig5_suppl2A <- assign_quadrants_single_mutants(df = single_mutants_sel_AA_BB)

fig5_suppl2A <- data_fig5_suppl2A %>% 
  ggplot(aes(x = quadrant, y = ddg, fill = Complex)) +
  geom_point(aes(colour = Complex),
             position = position_jitterdodge(jitter.width = 0.25)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +
  scale_fill_manual(values = c('#737373', '#d9d9d9')) +
  scale_colour_manual(values = c('#737373', '#d9d9d9'), guide = 'none') +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold', size = 15),
        axis.title.y = element_text(face = 'bold', size = 15),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 12)) +
  xlab('Class of mutation') + ylab('ΔΔG (kcal/mol)') + labs(fill = '') +
  ylim(-5,10)
fig5_suppl2A

data_fig5_suppl2B <- assign_quadrants_single_mutants(df = single_mutants_sel_AB)

fig5_suppl2B <- data_fig5_suppl2B %>% 
  ggplot(aes(x = quadrant, y = ddg, fill = Complex)) +
  geom_point(aes(colour = Complex),
             position = position_jitterdodge(jitter.width = 0.25)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +
  scale_fill_manual(values = c('#737373', '#d9d9d9')) +
  scale_colour_manual(values = c('#737373', '#d9d9d9'), guide = 'none') +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'), 
        panel.background = element_rect(fill = "white"),
        axis.title.x = element_text(face = 'bold', size = 15), axis.title.y = element_text(face = 'bold', size = 15),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5), 
        axis.text.y = element_text(size = 12),
        axis.line = element_line(size = 0.5)) +
  xlab('Class of mutation') + ylab('ΔΔG (kcal/mol)') + labs(fill = '') +
  ylim(-5,10)
fig5_suppl2B

p_effect_sizes <- plot_grid(fig5_suppl2A, NULL, fig5_suppl2B, labels = c('A', NULL, 'B'), nrow = 3, rel_heights = c(1, 0.1, 1))

ggsave(filename = 'Figures/Figure 5-figure_supplement_2.pdf',
       device = cairo_pdf, width = 7, height = 14, plot = p_effect_sizes, dpi = 500)
ggsave(filename = 'Figures/Figure 5-figure_supplement_2.png',
       device = 'png', width = 7, height = 14, plot = p_effect_sizes, dpi = 300)

#### Figure 5 - figure supplement 3 ####

double_mutants_sel_AA_BB <- extract_double_mutants(sel_AA_BB_1m38$all_ddgs)
double_mutants_sel_AB <- extract_double_mutants(sel_AB_1m38$all_ddgs)

double_mutants_sel_AA_BB <- get_fixation_rates(double_mutants_sel_AA_BB, 'Both homodimers')
double_mutants_sel_AB <- get_fixation_rates(double_mutants_sel_AB, 'Heterodimer')

summary_final <- rbind(double_mutants_sel_AA_BB, double_mutants_sel_AB) %>%
  mutate(rejected_count = total_count - accepted_count)

# Get the z-scores for the p-values
z_score_table <- summary_final %>%
  select(quadrant, accepted_percentage_of_total, total_count, selection) %>%
  unite(col = values, accepted_percentage_of_total, total_count, sep = ':') %>%
  spread(key = selection, value = values) %>%
  separate(col = `Both homodimers`, into = c('Acceptance_prob_HMs', 'Total_attempts_HMs'), sep = ':') %>%
  separate(col = Heterodimer, into = c('Acceptance_prob_HETs', 'Total_attempts_HETs'), sep = ':') %>%
  mutate(Acceptance_prob_HMs = as.numeric(Acceptance_prob_HMs)/ 100,
         Acceptance_prob_HETs = as.numeric(Acceptance_prob_HETs)/ 100,
         Total_attempts_HMs = as.numeric(Total_attempts_HMs),
         Total_attempts_HETs = as.numeric(Total_attempts_HETs)
         )

# Get the binomial p-values
p_vals <- get_pvals_binomial(p1 = z_score_table$Acceptance_prob_HMs, p2 = z_score_table$Acceptance_prob_HETs, 
                             n1 = z_score_table$Total_attempts_HMs, n2 = z_score_table$Total_attempts_HETs)


fig5_suppl3 <- summary_final %>%
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
  xlab("Effect on HMs") + ylab("Fixed mutations (%)") + labs(fill = 'Selection on') +
  geom_segment(aes(x = 0.75, xend = 1.25, y = 7, yend = 7), size = 0.5) +
  geom_segment(aes(x = 1.75, xend = 2.25, y = 12, yend = 12), size = 0.5) +
  geom_segment(aes(x = 2.75, xend = 3.25, y = 22.5, yend = 22.5), size = 0.5) +
  annotate(geom = 'text', label = p_vals[1], x = 1, y = 8) +
  annotate(geom = 'text', label = p_vals[2], x = 2, y = 13) +
  annotate(geom = 'text', label = p_vals[3], x = 3, y = 23.5)
fig5_suppl3

ggsave(filename = 'Figures/Figure 5-figure_supplement_3.pdf',
       device = cairo_pdf, width = 10, height = 7, plot = fig5_suppl3, dpi = 500)
ggsave(filename = 'Figures/Figure 5-figure_supplement_3.png',
       device = 'png', width = 10, height = 7, plot = fig5_suppl3, dpi = 300)
