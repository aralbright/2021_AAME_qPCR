#### Variability of knockdown in single Drosophila embryos using 
#### the maternal-Gal4 shRNA system and qRT-PCR
#### 
#### Author: Ashley Albright

## Load Packages  

library(readxl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(gridExtra)

## Load excel file into R as dataframe 
qPCR <- read_excel("/Users/aralbright/Box Sync/Eisen Things/qPCR/112217-qPCR/111217_qpcr.xlsx")
View(qPCR)

## Wrangle your data (have two separate Cp columns, one for each gene)

lid_df <- qPCR %>%
  dplyr::filter(Gene=='Lid') %>%
  dplyr::rename('lid_Cp' = 'Cp') %>%
  dplyr::select(-Category)

act_df <- qPCR %>% 
  dplyr::filter(Gene=='Act5c') %>%
  dplyr::rename('act_Cp' = 'Cp')

df <- dplyr::full_join(lid_df, act_df, by = c("Name","Replicate"))

# Calculate average control lid and actin Cp values, then dCp (used for all plots)

lid_control_mean <- mean(df[df$Category == 'nanos',]$lid_Cp)
act_control_mean <- mean(df[df$Category == 'nanos',]$act_Cp)

dCp_control <- (lid_control_mean-act_control_mean)

################## Figure 1A ################## 

# Calculate average Cp value per gene per condition

lid_line1_mean <- mean(df[df$Category == 'nanos 35706',]$lid_Cp)
act_line1_mean <- mean(df[df$Category == 'nanos 35706',]$act_Cp)
dCp_line1 <- (lid_line1_mean-act_line1_mean)


lid_line2_mean <- mean(df[df$Category == 'nanos 36652',]$lid_Cp)
act_line2_mean <- mean(df[df$Category == 'nanos 36652',]$act_Cp)
dCp_line2 <- (lid_line2_mean-act_line2_mean)

# Put those values in new dataframe

cats <- c('nanos', 'nanos 35706', 'nanos 36652')
dCp <- c(dCp_control, dCp_line1, dCp_line2)

p1_df <- data.frame(cats, dCp)

# Calculate ddCp, then 2^-ddCp

p1_df$ddCp <- (p1_df$dCp-dCp_control)

p1_df$'2^-ddCp' <- (2^-(p1_df$ddCp))

## Generate figure 

p1 <- ggplot(data=p1_df, aes(x=cats, y=`2^-ddCp`)) +
  geom_bar(stat="identity", fill="#88CCEE", width = 0.5)+
  ylab('Relative Gene Expression') + #changes y axis label
  xlab('') + #changes x axis label +
  ggtitle('A. pseudo-bulk before 2^-ddCp') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), #these 3 lines get rid of background
        axis.line = element_line(color="black", size = 0.5),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18)
        )

p1

################## Figure 1B ################## 

df$dCp <- df$lid_Cp-df$act_Cp

df$ddCp <- df$dCp-dCp_control

df$`2^-ddCp` <- (2^-(df$ddCp))

# Mean and Standard Error Per Condition

means_se <- df %>% 
  group_by(Category) %>% # Group the data  
  summarize(means=mean(`2^-ddCp`), # Create variable with mean  cty per group
            sds=sd(`2^-ddCp`), # Create variable with sd of cty per group
            N_=n(), # Create new variable N of cty per group
            se=sds/sqrt(N_), # Create variable with se of cty per group
            upper_limit=means+se, # Upper limit
            lower_limit=means-se # Lower limit
  )

means_se

# Significance 

stat.test <- compare_means(
  `2^-ddCp` ~ Category, data = df,
  method = "t.test")

stat.test <- stat.test[-c(3),] # remove unwanted row

stat.test <- stat.test %>%
  mutate(y.position = c(1.3, 1.6)) # determines where plotted in figure 

# Generate Figure 

p2 <- ggplot(data=means_se, aes(x=Category, y=means)) +
  geom_bar(stat="identity", fill="#88CCEE", width = 0.5)+
  geom_errorbar(aes(ymin=lower_limit, ymax=upper_limit), width = 0.1) +
  ylab('') + #changes y axis label
  xlab('Driver and shRNA Line') + #changes x axis label
  ggtitle('B. pseudo-bulk after 2^-ddCp') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
  stat_pvalue_manual(stat.test, label = "p.signif") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), #these 3 lines get rid of background
        axis.line = element_line(color="black", size = 0.5), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 16), 
        plot.title = element_text(size = 18)
        )

p2

################## Figure 1C ################## 


## Generate figure for single embryos 

p3 <- ggplot(df, aes(x = as.factor(`Category`), y = `2^-ddCp`)) +
  geom_jitter(color = "#88CCEE", width = 0.1) +
  guides(fill = FALSE) + #gets rid of legend
  ylab('') + #changes y axis label
  xlab('') +
  ggtitle('C. single embryo') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="pointrange", color="#117733") +
  stat_pvalue_manual(stat.test, label = "p.signif") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), #these 3 lines get rid of background
        axis.line = element_line(color="black", size = 0.5), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 18)
        )

p3


################## Arrange and Save Figure 1 ################## 

#test
grid.arrange(p1, p2, p3, nrow = 1)

#save
pdf("plot1.pdf", width = 13, height = 4)
print(grid.arrange(p1, p2, p3, nrow = 1))
dev.off()


################## Figure 2 ################## 

muted24 <- c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#999933", "#44AA99", "#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#999933", "#44AA99", "#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#999933", "#44AA99")

p4 <- ggplot(df, aes(x = as.factor(`Category`), y = `2^-ddCp`)) +
  geom_dotplot(aes(fill = as.factor(Name)), binaxis = "y", stackdir = "center") +
  scale_fill_manual(values=muted24) +
  guides(fill = FALSE) + #gets rid of legend
  ylab('Relative Gene Expression') + #changes y axis label
  xlab('Driver and shRNA Line') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), #these 3 lines get rid of background
        axis.line = element_line(color="black", size = 0.5), 
        plot.title = element_text(size = 18)
        )

p4

#save
pdf("plot2.pdf", width = 4, height = 3)
print(p4)
dev.off()
