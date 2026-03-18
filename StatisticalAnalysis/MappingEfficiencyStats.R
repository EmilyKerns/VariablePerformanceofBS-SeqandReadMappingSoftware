### Analyze variation in mapping efficiency for each alignment method

## RRBS
# Threespine stickleback, Gasterosteus aculeatus (this study)
# Cichlid, Astatotilapia calliptera (Vernaz et al., 2021)
# Great Tit, Parus major (Viitaniemi et al., 2019)
# Sea Urchin, Strongylocentrotus purpuratus (Bogan et al., 2023)

## WGBS
# Threespine stickleback, Gasterosteus aculeatus (this study)
# Cichlid, Astatotilapia calliptera (Vernaz et al., 2021)
# Coral, Acropora nana (Guerrero & Bay, 2024)
# Stick Insects, Timema cristinae (de Carvalho et al., 2023)




# RRBS --------------------------------------------------------------




## All --------------------------------------------------------------



setwd("/MethylMethods")

dat <- read_excel("QCStats_1.xlsx", sheet = "BamtoolsStats")
{
dat <- filter(dat, !AlignmentMethod == "bwa mem")

dat$PercentFailedQC[is.na(dat$PercentFailedQC)] <- 0

dat$PercMappedErr <- dat$PercentMapped - dat$PercentFailedQC

dat$AlignmentMethod <- factor(dat$AlignmentMethod, levels = c("Bismark","BismarkLocal","BisulfiteBolt", "Biscuit", "BWA meth"))
}

summary <- dat %>%
  group_by(Organism, AlignmentMethod) %>%
  summarise(
    n = n(),
    mean = mean(PercMappedErr, na.rm = TRUE),
    median = median(PercMappedErr, na.rm = TRUE),
    sd = sd(PercMappedErr, na.rm = TRUE),
    se = (sd(PercMappedErr, na.rm = TRUE))/(sqrt(n())),
    min = min(PercMappedErr, na.rm = TRUE),
    max = max(PercMappedErr, na.rm = TRUE)
  )
summary
# Organism    AlignmentMethod     n  mean median      sd       se   min   max
# 1 Cichlid     Bismark            35 0.842  0.844 0.0156  0.00264  0.801 0.864
# 2 Cichlid     BismarkLocal       35 0.852  0.855 0.0156  0.00264  0.812 0.874
# 3 Cichlid     BisulfiteBolt      35 0.882  0.885 0.0131  0.00222  0.847 0.899
# 4 Cichlid     Biscuit            35 0.995  0.996 0.00103 0.000174 0.992 0.997
# 5 Cichlid     BWA meth           35 0.987  0.989 0.00564 0.000953 0.967 0.994
# 6 Great Tit   Bismark           112 0.367  0.483 0.191   0.0181   0.047 0.593
# 7 Great Tit   BismarkLocal      112 0.813  0.812 0.0238  0.00225  0.745 0.894
# 8 Great Tit   BisulfiteBolt     112 0.930  0.941 0.0258  0.00244  0.865 0.960
# 9 Great Tit   Biscuit           112 0.716  0.733 0.0912  0.00862  0.454 0.838
# 10 Great Tit   BWA meth          112 0.636  0.653 0.0737  0.00697  0.454 0.763
# 11 Sea Urchin  Bismark            12 0.374  0.374 0.00619 0.00179  0.363 0.383
# 12 Sea Urchin  BismarkLocal       12 0.654  0.655 0.00252 0.000726 0.648 0.658
# 13 Sea Urchin  BisulfiteBolt      12 0.837  0.837 0.00308 0.000890 0.832 0.842
# 14 Sea Urchin  Biscuit            12 0.977  0.977 0.00361 0.00104  0.969 0.981
# 15 Sea Urchin  BWA meth           12 0.931  0.931 0.00490 0.00142  0.921 0.938
# 16 Stickleback Bismark            34 0.544  0.560 0.0435  0.00746  0.419 0.611
# 17 Stickleback BismarkLocal       34 0.652  0.668 0.0555  0.00951  0.518 0.74 
# 18 Stickleback BisulfiteBolt      34 0.798  0.810 0.0438  0.00752  0.708 0.861
# 19 Stickleback Biscuit            34 0.980  0.981 0.00415 0.000711 0.969 0.987
# 20 Stickleback BWA meth           34 0.853  0.844 0.0323  0.00553  0.802 0.909

# ANOVA
ANOVA <- aov(PercentMapped ~ AlignmentMethod,
             data = dat
)
ANOVA

# Call:
#   aov(formula = PercentMapped ~ AlignmentMethod, data = dat)
# 
# Terms:
#   AlignmentMethod Residuals
# Sum of Squares         16.23077  20.35532
# Deg. of Freedom               4       887
# 
# Residual standard error: 0.1514876
# Estimated effects may be unbalanced
# 73 observations deleted due to missingness

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# Shapiro-Wilk normality test
# 
# data:  ANOVA$residuals
# W = 0.9715, p-value = 3.338e-12

#### Does not meet assumption of normality, moving to Friedman test

dat_wide <- dat %>%
  select(SampleID, AlignmentMethod, PercentMapped) %>%
  pivot_wider(names_from = AlignmentMethod, 
              values_from = PercentMapped)

dat_complete <- dat_wide %>%
  filter(complete.cases(.)) %>%  # Remove rows with any NAs
  pivot_longer(cols = -SampleID,
               names_to = "AlignmentMethod",
               values_to = "PercentMapped") %>%
  mutate(AlignmentMethod = factor(AlignmentMethod, 
                                  levels = c("Bismark", "BismarkLocal", 
                                             "BisulfiteBolt", "Biscuit", "BWA meth")))

# Check how many complete samples you have
length(unique(dat_complete$SampleID))

## Friedman test
res.fried <- friedman_test(dat_complete, PercentMapped ~ AlignmentMethod |SampleID)
res.fried
#     .y.               n   statistic    df        p      method       
#   1 PercentMapped   139      286.      4    1.30e-60    Friedman test

dat_complete %>% friedman_effsize(PercentMapped ~ AlignmentMethod |SampleID)
#      .y.               n  effsize     method      magnitude
#   1 PercentMapped   139   0.514     Kendall W     large

# post hoc test
pwc <- dat_complete %>%
  wilcox_test(PercentMapped ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
# .y.           group1        group2           n1    n2 statistic        p    p.adj p.adj.signif
# PercentMapped Bismark       BismarkLocal    139   139         0 1.38e-24 1.38e-23 ****        
# PercentMapped Bismark       BisulfiteBolt   139   139         0 1.49e-24 1.49e-23 ****        
# PercentMapped Bismark       Biscuit         139   139        48 4.21e-24 4.21e-23 ****        
# PercentMapped Bismark       BWA meth        139   139        48 4.21e-24 4.21e-23 ****        
# PercentMapped BismarkLocal  BisulfiteBolt   139   139         1 1.53e-24 1.53e-23 ****        
# PercentMapped BismarkLocal  Biscuit         139   139      2507 7.17e- 7 7.17e- 6 ****        
# PercentMapped BismarkLocal  BWA meth        139   139      3285 8.97e- 4 9   e- 3 **          
# PercentMapped BisulfiteBolt Biscuit         139   139      6099 1   e- 2 9.5 e- 2 ns          
# PercentMapped BisulfiteBolt BWA meth        139   139      6497 6.03e- 4 6   e- 3 **          
# PercentMapped Biscuit       BWA meth        139   139      5230 1.24e- 7 1.24e- 6 **** 

dat_complete %>% wilcox_effsize(PercentMapped ~ AlignmentMethod, paired = TRUE)
# .y.           group1        group2        effsize    n1    n2 magnitude
# PercentMapped Bismark       BismarkLocal    0.868   139   139 large    
# PercentMapped Bismark       BisulfiteBolt   0.868   139   139 large    
# PercentMapped Bismark       Biscuit         0.859   139   139 large    
# PercentMapped Bismark       BWA meth        0.859   139   139 large    
# PercentMapped BismarkLocal  BisulfiteBolt   0.867   139   139 large    
# PercentMapped BismarkLocal  Biscuit         0.421   139   139 moderate 
# PercentMapped BismarkLocal  BWA meth        0.282   139   139 small    
# PercentMapped BisulfiteBolt Biscuit         0.220   139   139 small    
# PercentMapped BisulfiteBolt BWA meth        0.291   139   139 small    
# PercentMapped Biscuit       BWA meth        0.444   139   139 moderate



## Threesping stickleback --------------------------------------------------


stickle <- dat %>% 
  filter(Organism == "Stickleback")


summary <- stickle %>%
  group_by(Batch, AlignmentMethod) %>%
  summarise(
    n = n(),
    mean = mean(PercMappedErr, na.rm = TRUE),
    median = median(PercMappedErr, na.rm = TRUE),
    sd = sd(PercMappedErr, na.rm = TRUE),
    se = (sd(PercMappedErr, na.rm = TRUE))/(sqrt(n())),
    min = min(PercMappedErr, na.rm = TRUE),
    max = max(PercMappedErr, na.rm = TRUE)
  )
summary
# Batch AlignmentMethod     n  mean median      sd       se   min   max
# 1     Bismark            22 0.545  0.560 0.0447  0.00952  0.419 0.611
# 1     BismarkLocal       22 0.668  0.683 0.0523  0.0111   0.518 0.74 
# 1     BisulfiteBolt      22 0.825  0.831 0.0253  0.00539  0.765 0.861
# 1     Biscuit            22 0.982  0.983 0.00268 0.000571 0.977 0.987
# 1     BWA meth           22 0.832  0.833 0.0156  0.00333  0.802 0.863
# 2     Bismark            12 0.542  0.556 0.0431  0.0124   0.466 0.589
# 2     BismarkLocal       12 0.618  0.625 0.0478  0.0138   0.546 0.676
# 2     BisulfiteBolt      12 0.749  0.751 0.0211  0.00608  0.708 0.787
# 2     Biscuit            12 0.976  0.977 0.00343 0.000991 0.969 0.980
# 2     BWA meth           12 0.892  0.893 0.0122  0.00351  0.870 0.909

# ANOVA
ANOVA <- aov(PercMappedErr ~ AlignmentMethod,
             data = stickle
)
ANOVA

# AlignmentMethod Residuals
# Sum of Squares         3.966472  0.259136
# Deg. of Freedom               4       164
# 
# Residual standard error: 0.03703051
# Estimated effects may be unbalanced
# 1 observation deleted due to missingness

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# Shapiro-Wilk normality test
# 
# data:  ANOVA$residuals
# W = 0.97302, p-value = 0.002232

#### Does not meet assumption of normality, moving to Friedman test

stickle_wide <- stickle %>%
  select(SampleID, AlignmentMethod, PercMappedErr) %>%
  pivot_wider(names_from = AlignmentMethod, 
              values_from = PercMappedErr)

stickle_complete <- stickle_wide %>%
  filter(complete.cases(.)) %>%  # Remove rows with any NAs
  pivot_longer(cols = -SampleID,
               names_to = "AlignmentMethod",
               values_to = "PercMappedErr") %>%
  mutate(AlignmentMethod = factor(AlignmentMethod, 
                                  levels = c("Bismark", "BismarkLocal", 
                                             "BisulfiteBolt", "Biscuit", "BWA meth")))


## Friedman test
res.fried <- friedman_test(stickle_complete, PercMappedErr ~ AlignmentMethod |SampleID)
res.fried
#     .y.               n statistic    df          p  method     
#   PercMappedErr    33      127.     4 1.92e-26 Friedman test

stickle_complete %>% friedman_effsize(PercMappedErr ~ AlignmentMethod |SampleID)
#     .y.               n   effsize   method      magnitude
#   PercMappedErr    33   0.960 Kendall W large    

# post hoc test
pwc <- stickle_complete %>%
  wilcox_test(PercMappedErr ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
# .y.           group1        group2           n1    n2 statistic        p        p.adj p.adj.signif
#   1 PercMappedErr Bismark       BismarkLocal     33    33         0 5.63e- 7      5.63e-6 ****        
#   2 PercMappedErr Bismark       BisulfiteBolt    33    33         0 2.33e-10      2.33e-9 ****        
#   3 PercMappedErr Bismark       Biscuit          33    33         0 2.33e-10      2.33e-9 ****        
#   4 PercMappedErr Bismark       BWA meth         33    33         0 2.33e-10      2.33e-9 ****        
#   5 PercMappedErr BismarkLocal  BisulfiteBolt    33    33         0 2.33e-10      2.33e-9 ****        
#   6 PercMappedErr BismarkLocal  Biscuit          33    33         0 2.33e-10      2.33e-9 ****        
#   7 PercMappedErr BismarkLocal  BWA meth         33    33         0 2.33e-10      2.33e-9 ****        
#   8 PercMappedErr BisulfiteBolt Biscuit          33    33         0 2.33e-10      2.33e-9 ****        
#   9 PercMappedErr BisulfiteBolt BWA meth         33    33       106 1   e- 3      1.3 e-2 *           
#   10 PercMappedErr Biscuit       BWA meth         33    33       561 2.33e-10      2.33e-9 ****  

stickle_complete %>% wilcox_effsize(PercMappedErr ~ AlignmentMethod, paired = TRUE)
# .y.           group1        group2        effsize    n1    n2 magnitude
# * <chr>         <chr>         <chr>           <dbl> <int> <int> <ord>    
#   1 PercMappedErr Bismark       BismarkLocal    0.873    33    33 large    
# 2 PercMappedErr Bismark       BisulfiteBolt   0.872    33    33 large    
# 3 PercMappedErr Bismark       Biscuit         0.872    33    33 large    
# 4 PercMappedErr Bismark       BWA meth        0.872    33    33 large    
# 5 PercMappedErr BismarkLocal  BisulfiteBolt   0.872    33    33 large    
# 6 PercMappedErr BismarkLocal  Biscuit         0.872    33    33 large    
# 7 PercMappedErr BismarkLocal  BWA meth        0.872    33    33 large    
# 8 PercMappedErr BisulfiteBolt Biscuit         0.872    33    33 large    
# 9 PercMappedErr BisulfiteBolt BWA meth        0.543    33    33 large    
# 10 PercMappedErr Biscuit       BWA meth        0.872    33    33 large  

kruskal.test(PercMappedErr ~ Batch, data = stickle)
# Kruskal-Wallis chi-squared = 0.46156, df = 1, p-value = 0.4969


## Cichlid --------------------------------------------------


cich <- dat %>% 
  filter(Organism == "Cichlid") %>% 
  select("AlignmentMethod", "PercentMapped", "Batch", "Population", "Environment", "SampleID")


# ANOVA
ANOVA <- aov(PercentMapped ~ AlignmentMethod,
             data = cich
)
ANOVA

# Call:
#   aov(formula = PercentMapped ~ AlignmentMethod, data = dat)
# 
#                 AlignmentMethod Residuals
# Sum of Squares        0.7558351 0.0235296
# Deg. of Freedom               4       168
# 
# Residual standard error: 0.01183457
# Estimated effects may be unbalanced
# 2 observations deleted due to missingness

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# Shapiro-Wilk normality test
# 
# data:  ANOVA$residuals
#W = 0.92967, p-value = 1.873e-07

#### Does not meet assumption of normality, moving to Friedman test

cich_wide <- cich %>%
  select(SampleID, AlignmentMethod, PercentMapped) %>%
  pivot_wider(names_from = AlignmentMethod, 
              values_from = PercentMapped)

cich_complete <- cich_wide %>%
  filter(complete.cases(.)) %>%  # Remove rows with any NAs
  pivot_longer(cols = -SampleID,
               names_to = "AlignmentMethod",
               values_to = "PercentMapped") %>%
  mutate(AlignmentMethod = factor(AlignmentMethod, 
                                  levels = c("Bismark", "BismarkLocal", 
                                             "BisulfiteBolt", "Biscuit", "BWA meth")))

## Friedman test
res.fried <- friedman_test(cich_complete, PercentMapped ~ AlignmentMethod |SampleID)
res.fried
#         .y.               n    statistic     df         p   method       
#   1 PercentMapped        33         131.     4    2.13e-27  Friedman test

cich_complete %>% friedman_effsize(PercentMapped ~ AlignmentMethod |SampleID)
#     .y.               n   effsize   method      magnitude
#   1 PercentMapped    33     0.994   Kendall W   large

# post hoc test
pwc <- cich_complete %>%
  wilcox_test(PercentMapped ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
# .y.           group1        group2           n1    n2 statistic        p         p.adj p.adj.signif
#   1 PercentMapped Bismark       BismarkLocal     33    33         0 1.24e- 7 0.00000124    ****        
#   2 PercentMapped Bismark       BisulfiteBolt    33    33         0 2.33e-10 0.00000000233 ****        
#   3 PercentMapped Bismark       Biscuit          33    33         0 2.33e-10 0.00000000233 ****        
#   4 PercentMapped Bismark       BWA meth         33    33         0 2.33e-10 0.00000000233 ****        
#   5 PercentMapped BismarkLocal  BisulfiteBolt    33    33         1 4.66e-10 0.00000000466 ****        
#   6 PercentMapped BismarkLocal  Biscuit          33    33         0 2.33e-10 0.00000000233 ****        
#   7 PercentMapped BismarkLocal  BWA meth         33    33         0 2.33e-10 0.00000000233 ****        
#   8 PercentMapped BisulfiteBolt Biscuit          33    33         0 2.33e-10 0.00000000233 ****        
#   9 PercentMapped BisulfiteBolt BWA meth         33    33         0 2.33e-10 0.00000000233 ****        
#   10 PercentMapped Biscuit       BWA meth         33    33       561 2.33e-10 0.00000000233 ****

cich_complete %>% wilcox_effsize(PercentMapped ~ AlignmentMethod, paired = TRUE)
# .y.           group1        group2        effsize    n1    n2 magnitude
# 1 PercentMapped Bismark       BismarkLocal    0.922    33    33 large    
# 2 PercentMapped Bismark       BisulfiteBolt   0.872    33    33 large    
# 3 PercentMapped Bismark       Biscuit         0.872    33    33 large    
# 4 PercentMapped Bismark       BWA meth        0.872    33    33 large    
# 5 PercentMapped BismarkLocal  BisulfiteBolt   0.869    33    33 large    
# 6 PercentMapped BismarkLocal  Biscuit         0.872    33    33 large    
# 7 PercentMapped BismarkLocal  BWA meth        0.872    33    33 large    
# 8 PercentMapped BisulfiteBolt Biscuit         0.872    33    33 large    
# 9 PercentMapped BisulfiteBolt BWA meth        0.872    33    33 large    
# 10 PercentMapped Biscuit       BWA meth        0.872    33    33 large 




## Great Tit --------------------------------------------------


GT <- dat %>% 
  filter(Organism == "Great Tit") %>% 
  select("AlignmentMethod", "PercentMapped", "Batch", "Population", "Environment", "SampleID")


# ANOVA
ANOVA <- aov(PercentMapped ~ AlignmentMethod,
             data = GT
)
ANOVA

# Call:
#   aov(formula = PercentMapped ~ AlignmentMethod, data = dat)
# 
#                 AlignmentMethod Residuals
# Sum of Squares        17.340011  4.572371
# Deg. of Freedom               4       485
# 
# Residual standard error: 0.09709567
# Estimated effects may be unbalanced
# 70 observations deleted due to missingness

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# Shapiro-Wilk normality test
# 
# data:  ANOVA$residuals
# W = 0.88831, p-value < 2.2e-16

#### Does not meet assumption of normality, moving to Friedman test

GT_wide <- GT %>%
  select(SampleID, AlignmentMethod, PercentMapped) %>%
  pivot_wider(names_from = AlignmentMethod, 
              values_from = PercentMapped)

GT_complete <- GT_wide %>%
  filter(complete.cases(.)) %>%  # Remove rows with any NAs
  pivot_longer(cols = -SampleID,
               names_to = "AlignmentMethod",
               values_to = "PercentMapped") %>%
  mutate(AlignmentMethod = factor(AlignmentMethod, 
                                  levels = c("Bismark", "BismarkLocal", 
                                             "BisulfiteBolt", "Biscuit", "BWA meth")))

## Friedman test
res.fried <- friedman_test(GT_complete, PercentMapped ~ AlignmentMethod |SampleID)
res.fried
#     .y.               n statistic    df         p     method       
#   1 PercentMapped    61      226.     4   1.03e-47    Friedman test

GT_complete %>% friedman_effsize(PercentMapped ~ AlignmentMethod |SampleID)
# .y.                   n   effsize   method      magnitude
#   1 PercentMapped    61     0.926   Kendall W   large   

# post hoc test
pwc <- GT_complete %>%
  wilcox_test(PercentMapped ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
#     .y.           group1        group2           n1    n2 statistic        p    p.adj p.adj.signif
#   1 PercentMapped Bismark       BismarkLocal     61    61         0 1.14e-11 1.14e-10 ****        
#   2 PercentMapped Bismark       BisulfiteBolt    61    61         0 1.14e-11 1.14e-10 ****        
#   3 PercentMapped Bismark       Biscuit          61    61        48 1.17e-10 1.17e- 9 ****        
#   4 PercentMapped Bismark       BWA meth         61    61        48 1.17e-10 1.17e- 9 ****        
#   5 PercentMapped BismarkLocal  BisulfiteBolt    61    61         0 1.14e-11 1.14e-10 ****        
#   6 PercentMapped BismarkLocal  Biscuit          61    61      1891 1.14e-11 1.14e-10 ****        
#   7 PercentMapped BismarkLocal  BWA meth         61    61      1891 1.14e-11 1.14e-10 ****        
#   8 PercentMapped BisulfiteBolt Biscuit          61    61      1891 1.14e-11 1.14e-10 ****        
#   9 PercentMapped BisulfiteBolt BWA meth         61    61      1891 1.14e-11 1.14e-10 ****        
#   10 PercentMapped Biscuit       BWA meth         61    61       703 1.19e- 7 1.19e- 6 ****

GT_complete %>% wilcox_effsize(PercentMapped ~ AlignmentMethod, paired = TRUE)
# .y.           group1        group2        effsize    n1    n2 magnitude
# 1 PercentMapped Bismark       BismarkLocal    0.870    61    61 large    
# 2 PercentMapped Bismark       BisulfiteBolt   0.870    61    61 large    
# 3 PercentMapped Bismark       Biscuit         0.825    61    61 large    
# 4 PercentMapped Bismark       BWA meth        0.825    61    61 large    
# 5 PercentMapped BismarkLocal  BisulfiteBolt   0.870    61    61 large    
# 6 PercentMapped BismarkLocal  Biscuit         0.870    61    61 large    
# 7 PercentMapped BismarkLocal  BWA meth        0.870    61    61 large    
# 8 PercentMapped BisulfiteBolt Biscuit         0.870    61    61 large    
# 9 PercentMapped BisulfiteBolt BWA meth        0.870    61    61 large    
# 10 PercentMapped Biscuit       BWA meth        0.756    61    61 large



## Urchin --------------------------------------------------


Urch <- dat %>% 
  filter(Organism == "Sea Urchin") %>% 
  select("AlignmentMethod", "PercentMapped", "Batch", "Population", "Environment", "SampleID")


# ANOVA
ANOVA <- aov(PercentMapped ~ AlignmentMethod,
             data = Urch
)
ANOVA

# Call:
#   aov(formula = PercentMapped ~ AlignmentMethod, data = dat)
# 
#                 AlignmentMethod Residuals
# Sum of Squares        3.0047492 0.0009889
# Deg. of Freedom               4        55
# 
# Residual standard error: 0.004240184
# Estimated effects may be unbalanced

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# Shapiro-Wilk normality test
# 
# data:  ANOVA$residuals
# W = 0.98225, p-value = 0.5306 

#### normally distributed, but using Friedman test for consistency

Urch_wide <- Urch %>%
  select(SampleID, AlignmentMethod, PercentMapped) %>%
  pivot_wider(names_from = AlignmentMethod, 
              values_from = PercentMapped)

Urch_complete <- Urch_wide %>%
  filter(complete.cases(.)) %>%  # Remove rows with any NAs
  pivot_longer(cols = -SampleID,
               names_to = "AlignmentMethod",
               values_to = "PercentMapped") %>%
  mutate(AlignmentMethod = factor(AlignmentMethod, 
                                  levels = c("Bismark", "BismarkLocal", 
                                             "BisulfiteBolt", "Biscuit", "BWA meth")))

## Friedman test
res.fried <- friedman_test(Urch_complete, PercentMapped ~ AlignmentMethod |SampleID)
res.fried
#     .y.               n statistic    df        p method       
#   1 PercentMapped    12        48     4 9.44e-10 Friedman test

Urch_complete %>% friedman_effsize(PercentMapped ~ AlignmentMethod |SampleID)
#     .y.               n effsize method    magnitude
#   1 PercentMapped    12       1 Kendall W large

# post hoc test
pwc <- Urch_complete %>%
  wilcox_test(PercentMapped ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
# .y.           group1        group2           n1    n2 statistic        p p.adj p.adj.signif
# * <chr>         <chr>         <chr>         <int> <int>     <dbl>    <dbl> <dbl> <chr>       
#   1 PercentMapped Bismark       BismarkLocal     12    12         0 0.002    0.025 *           
#   2 PercentMapped Bismark       BisulfiteBolt    12    12         0 0.000488 0.005 **          
#   3 PercentMapped Bismark       Biscuit          12    12         0 0.000488 0.005 **          
#   4 PercentMapped Bismark       BWA meth         12    12         0 0.000488 0.005 **          
#   5 PercentMapped BismarkLocal  BisulfiteBolt    12    12         0 0.000488 0.005 **          
#   6 PercentMapped BismarkLocal  Biscuit          12    12         0 0.000488 0.005 **          
#   7 PercentMapped BismarkLocal  BWA meth         12    12         0 0.000488 0.005 **          
#   8 PercentMapped BisulfiteBolt Biscuit          12    12         0 0.000488 0.005 **          
#   9 PercentMapped BisulfiteBolt BWA meth         12    12         0 0.000488 0.005 **          
#   10 PercentMapped Biscuit       BWA meth         12    12        78 0.000488 0.005 **

Urch_complete %>% wilcox_effsize(PercentMapped ~ AlignmentMethod, paired = TRUE)
# .y.           group1        group2        effsize    n1    n2 magnitude
# 1 PercentMapped Bismark       BismarkLocal    0.884    12    12 large    
# 2 PercentMapped Bismark       BisulfiteBolt   0.883    12    12 large    
# 3 PercentMapped Bismark       Biscuit         0.883    12    12 large    
# 4 PercentMapped Bismark       BWA meth        0.883    12    12 large    
# 5 PercentMapped BismarkLocal  BisulfiteBolt   0.883    12    12 large    
# 6 PercentMapped BismarkLocal  Biscuit         0.883    12    12 large    
# 7 PercentMapped BismarkLocal  BWA meth        0.883    12    12 large    
# 8 PercentMapped BisulfiteBolt Biscuit         0.883    12    12 large    
# 9 PercentMapped BisulfiteBolt BWA meth        0.883    12    12 large    
# 10 PercentMapped Biscuit       BWA meth        0.883    12    12 large   




# WGBS ------------------------------------------------------------


## All ---------------------------------------------------------------------


dat <- read_excel("QCStats_1.xlsx", sheet = "WBGSBamtoolsStats")
{
  dat <- filter(dat, !AlignmentMethod == "bwa mem")
  
  dat$PercentFailedQC[is.na(dat$PercentFailedQC)] <- 0
  
  dat$PercMappedErr <- dat$PercentMapped - dat$PercentFailedQC
  
  dat$AlignmentMethod <- factor(dat$AlignmentMethod, levels = c("Bismark","BismarkLocal","BisulfiteBolt", "Biscuit", "bwa meth"))
}

summary <- dat %>%
  group_by(Organism, AlignmentMethod) %>%
  summarise(
    n = n(),
    mean = mean(PercMappedErr, na.rm = TRUE),
    median = median(PercMappedErr, na.rm = TRUE),
    sd = sd(PercMappedErr, na.rm = TRUE),
    se = (sd(PercMappedErr, na.rm = TRUE))/(sqrt(n())),
    min = min(PercMappedErr, na.rm = TRUE),
    max = max(PercMappedErr, na.rm = TRUE)
  )
summary
# Organism    AlignmentMethod     n  mean median       sd       se   min   max
# 1 Cichlid     Bismark             8 0.527  0.518 0.0210   0.00741  0.502 0.565
# 2 Cichlid     BismarkLocal        8 0.770  0.767 0.0128   0.00452  0.755 0.797
# 3 Cichlid     BisulfiteBolt       8 0.881  0.880 0.00447  0.00158  0.875 0.888
# 4 Cichlid     Biscuit             8 0.997  0.997 0.000623 0.000220 0.996 0.998
# 5 Cichlid     bwa meth            8 0.825  0.828 0.00972  0.00344  0.808 0.837
# 6 Coral       Bismark            48 0.207  0.210 0.0253   0.00365  0.04  0.223
# 7 Coral       BismarkLocal       48 0.628  0.639 0.0689   0.00994  0.181 0.674
# 8 Coral       BisulfiteBolt      48 0.940  0.940 0.0105   0.00151  0.917 0.987
# 9 Coral       Biscuit            48 0.963  0.967 0.0159   0.00230  0.884 0.981
# 10 Coral       bwa meth           48 0.598  0.596 0.0271   0.00391  0.509 0.687
# 11 StickInsect Bismark            24 0.415  0.417 0.0170   0.00347  0.372 0.451
# 12 StickInsect BismarkLocal       24 0.769  0.774 0.0175   0.00357  0.72  0.785
# 13 StickInsect BisulfiteBolt      24 0.934  0.931 0.00778  0.00159  0.926 0.948
# 14 StickInsect Biscuit            24 0.984  0.987 0.00855  0.00174  0.954 0.990
# 15 StickInsect bwa meth           24 0.864  0.870 0.0197   0.00403  0.804 0.886
# 16 Stickleback Bismark             4 0.387  0.392 0.0124   0.00618  0.369 0.396
# 17 Stickleback BismarkLocal        4 0.888  0.893 0.0153   0.00763  0.865 0.899
# 18 Stickleback BisulfiteBolt       4 0.941  0.941 0.00117  0.000585 0.940 0.943
# 19 Stickleback Biscuit             4 0.991  0.991 0.00125  0.000624 0.990 0.992
# 20 Stickleback bwa meth            4 0.941  0.943 0.00463  0.00231  0.934 0.945

# ANOVA
ANOVA <- aov(PercentMapped ~ AlignmentMethod,
             data = dat
)
ANOVA

# Call:
#   aov(formula = PercentMapped ~ AlignmentMethod, data = dat)
# 
# Terms:
#   AlignmentMethod Residuals
# Sum of Squares        26.811252  1.964853
# Deg. of Freedom               4       365
# 
# Residual standard error: 0.07337001
# Estimated effects may be unbalanced
# 50 observations deleted due to missingness

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# Shapiro-Wilk normality test
# 
# data:  ANOVA$residuals
# W = 0.8794, p-value < 2.2e-16

#### Does not meet assumption of normality, moving to Friedman test

wgbs_wide <- dat %>%
  select(SampleID, AlignmentMethod, PercentMapped) %>%
  pivot_wider(names_from = AlignmentMethod, 
              values_from = PercentMapped)

wgbs_complete <- wgbs_wide %>%
  filter(complete.cases(.)) %>%  
  pivot_longer(cols = -SampleID,
               names_to = "AlignmentMethod",
               values_to = "PercentMapped") %>%
  mutate(AlignmentMethod = factor(AlignmentMethod, 
                                  levels = c("Bismark", "BismarkLocal", 
                                             "BisulfiteBolt", "Biscuit", "bwa meth"))) %>% 
  filter(!is.na(AlignmentMethod)) 
 

## Friedman test
res.fried <- friedman_test(wgbs_complete, PercentMapped ~ AlignmentMethod |SampleID)
res.fried
#     .y.               n statistic    df        p method       
#   1 PercentMapped    47      186.     4 3.04e-39 Friedman test

wgbs_complete %>% friedman_effsize(PercentMapped ~ AlignmentMethod |SampleID)
#     .y.               n effsize method    magnitude
#   1 PercentMapped    47   0.992 Kendall W large

# post hoc test
pwc <- wgbs_complete %>%
  wilcox_test(PercentMapped ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
#     .y.           group1        group2           n1    n2 statistic        p    p.adj p.adj.signif
#   1 PercentMapped Bismark       BismarkLocal     47    47         0 2.47e- 9 2.47e- 8 ****        
#   2 PercentMapped Bismark       BisulfiteBolt    47    47         0 1.42e-14 1.42e-13 ****        
#   3 PercentMapped Bismark       Biscuit          47    47         0 1.42e-14 1.42e-13 ****        
#   4 PercentMapped Bismark       bwa meth         47    47         0 1.42e-14 1.42e-13 ****        
#   5 PercentMapped BismarkLocal  BisulfiteBolt    47    47         0 1.42e-14 1.42e-13 ****        
#   6 PercentMapped BismarkLocal  Biscuit          47    47         0 1.42e-14 1.42e-13 ****        
#   7 PercentMapped BismarkLocal  bwa meth         47    47         0 1.42e-14 1.42e-13 ****        
#   8 PercentMapped BisulfiteBolt Biscuit          47    47        43 1.85e-10 1.85e- 9 ****        
#   9 PercentMapped BisulfiteBolt bwa meth         47    47         0 1.42e-14 1.42e-13 ****        
#   10 PercentMapped Biscuit       bwa meth         47    47         0 1.42e-14 1.42e-13 **** 

wgbs_complete %>% wilcox_effsize(PercentMapped ~ AlignmentMethod, paired = TRUE)
#   .y.           group1        group2        effsize    n1    n2 magnitude
# 1 PercentMapped Bismark       BismarkLocal    0.871    47    47 large    
# 2 PercentMapped Bismark       BisulfiteBolt   0.871    47    47 large    
# 3 PercentMapped Bismark       Biscuit         0.871    47    47 large    
# 4 PercentMapped Bismark       bwa meth        0.871    47    47 large    
# 5 PercentMapped BismarkLocal  BisulfiteBolt   0.871    47    47 large    
# 6 PercentMapped BismarkLocal  Biscuit         0.871    47    47 large    
# 7 PercentMapped BismarkLocal  bwa meth        0.871    47    47 large    
# 8 PercentMapped BisulfiteBolt Biscuit         0.804    47    47 large    
# 9 PercentMapped BisulfiteBolt bwa meth        0.871    47    47 large    
# 10 PercentMapped Biscuit       bwa meth        0.871    47    47 large 



## Threespine Stickleback --------------------------------------------------------


stickle <- dat %>% 
  filter(Organism == "Stickleback")



# ANOVA
ANOVA <- aov(PercentMapped ~ AlignmentMethod,
             data = stickle
)
ANOVA

# Call:
#   aov(formula = PercentMapped ~ AlignmentMethod, data = dat)
# 
# Terms:
#                   AlignmentMethod Residuals
# Sum of Squares        1.0588639 0.0011732
# Deg. of Freedom               4        15
# 
# Residual standard error: 0.008843719
# Estimated effects may be unbalanced

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# Shapiro-Wilk normality test
# 
# data:  ANOVA$residuals
# W = 0.78819, p-value = 0.000579

#### Does not meet assumption of normality, moving to Friedman test

stickle_wide <- stickle %>%
  select(SampleID, AlignmentMethod, PercentMapped) %>%
  pivot_wider(names_from = AlignmentMethod, 
              values_from = PercentMapped)

stickle_complete <- stickle_wide %>%
  filter(complete.cases(.)) %>%  
  pivot_longer(cols = -SampleID,
               names_to = "AlignmentMethod",
               values_to = "PercentMapped") %>%
  mutate(AlignmentMethod = factor(AlignmentMethod, 
                                  levels = c("Bismark", "BismarkLocal", 
                                             "BisulfiteBolt", "Biscuit", "bwa meth"))) %>% 
  filter(!is.na(AlignmentMethod)) 

## Friedman test
res.fried <- friedman_test(stickle_complete, PercentMapped ~ AlignmentMethod |SampleID)
res.fried
#     .y.               n statistic    df       p method       
#   1 PercentMapped     4        16     4 0.00302 Friedman test

stickle_complete %>% friedman_effsize(PercentMapped ~ AlignmentMethod |SampleID)
#     .y.               n effsize method    magnitude
#   1 PercentMapped     4       1 Kendall W large

# post hoc test
pwc <- stickle_complete %>%
  wilcox_test(PercentMapped ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
#   .y.           group1        group2           n1    n2 statistic     p p.adj p.adj.signif
# 1 PercentMapped Bismark       BismarkLocal      4     4         0 0.125     1 ns          
# 2 PercentMapped Bismark       BisulfiteBolt     4     4         0 0.125     1 ns          
# 3 PercentMapped Bismark       Biscuit           4     4         0 0.125     1 ns          
# 4 PercentMapped Bismark       bwa meth          4     4         0 0.125     1 ns          
# 5 PercentMapped BismarkLocal  BisulfiteBolt     4     4         0 0.125     1 ns          
# 6 PercentMapped BismarkLocal  Biscuit           4     4         0 0.125     1 ns          
# 7 PercentMapped BismarkLocal  bwa meth          4     4         0 0.125     1 ns          
# 8 PercentMapped BisulfiteBolt Biscuit           4     4         0 0.125     1 ns          
# 9 PercentMapped BisulfiteBolt bwa meth          4     4         0 0.125     1 ns          
# 10 PercentMapped Biscuit       bwa meth          4     4         0 0.125     1 ns

stickle_complete %>% wilcox_effsize(PercentMapped ~ AlignmentMethod, paired = TRUE)
#   .y.           group1        group2        effsize    n1    n2 magnitude
# 1 PercentMapped Bismark       BismarkLocal    0.913     4     4 large    
# 2 PercentMapped Bismark       BisulfiteBolt   0.913     4     4 large    
# 3 PercentMapped Bismark       Biscuit         0.913     4     4 large    
# 4 PercentMapped Bismark       bwa meth        0.913     4     4 large    
# 5 PercentMapped BismarkLocal  BisulfiteBolt   0.913     4     4 large    
# 6 PercentMapped BismarkLocal  Biscuit         0.913     4     4 large    
# 7 PercentMapped BismarkLocal  bwa meth        0.913     4     4 large    
# 8 PercentMapped BisulfiteBolt Biscuit         0.913     4     4 large    
# 9 PercentMapped BisulfiteBolt bwa meth        0.913     4     4 large    
# 10 PercentMapped Biscuit       bwa meth        0.913     4     4 large




## Cichlid --------------------------------------------------


cich <- dat %>% 
  filter(Organism == "Cichlid")

# ANOVA
ANOVA <- aov(PercentMapped ~ AlignmentMethod,
             data = cich
)
ANOVA

# Call:
#   aov(formula = PercentMapped ~ AlignmentMethod, data = dat)
# 
#                 AlignmentMethod Residuals
# Sum of Squares        1.0837403 0.0039029
# Deg. of Freedom               4        31
# 
# Residual standard error: 0.01122057
# Estimated effects may be unbalanced
# 4 observations deleted due to missingness

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# Shapiro-Wilk normality test
# 
# data:  ANOVA$residuals
# W = 0.83058, p-value = 7.011e-05

#### Does not meet assumption of normality, moving to Friedman test

cich_wide <- cich %>%
  select(SampleID, AlignmentMethod, PercentMapped) %>%
  pivot_wider(names_from = AlignmentMethod, 
              values_from = PercentMapped)

cich_complete <- cich_wide %>%
  filter(complete.cases(.)) %>%  
  pivot_longer(cols = -SampleID,
               names_to = "AlignmentMethod",
               values_to = "PercentMapped") %>%
  mutate(AlignmentMethod = factor(AlignmentMethod, 
                                  levels = c("Bismark", "BismarkLocal", 
                                             "BisulfiteBolt", "Biscuit", "bwa meth"))) %>% 
  filter(!is.na(AlignmentMethod)) 


## Friedman test
res.fried <- friedman_test(cich_complete, PercentMapped ~ AlignmentMethod |SampleID)
res.fried
#     .y.               n statistic    df         p method       
#   1 PercentMapped     6        24     4 0.0000799 Friedman test

cich_complete %>% friedman_effsize(PercentMapped ~ AlignmentMethod |SampleID)
#     .y.               n effsize method    magnitude
#   1 PercentMapped     6       1 Kendall W large

# post hoc test
pwc <- cich_complete %>%
  wilcox_test(PercentMapped ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
# .y.           group1        group2           n1    n2 statistic     p p.adj p.adj.signif
# 1 PercentMapped Bismark       BismarkLocal      6     6         0 0.035 0.355 ns          
# 2 PercentMapped Bismark       BisulfiteBolt     6     6         0 0.031 0.313 ns          
# 3 PercentMapped Bismark       Biscuit           6     6         0 0.031 0.313 ns          
# 4 PercentMapped Bismark       bwa meth          6     6         0 0.031 0.313 ns          
# 5 PercentMapped BismarkLocal  BisulfiteBolt     6     6         0 0.031 0.313 ns          
# 6 PercentMapped BismarkLocal  Biscuit           6     6         0 0.031 0.313 ns          
# 7 PercentMapped BismarkLocal  bwa meth          6     6         0 0.031 0.313 ns          
# 8 PercentMapped BisulfiteBolt Biscuit           6     6         0 0.031 0.313 ns          
# 9 PercentMapped BisulfiteBolt bwa meth          6     6         0 0.031 0.313 ns          
# 10 PercentMapped Biscuit       bwa meth          6     6         0 0.031 0.313 ns

cich_complete %>% wilcox_effsize(PercentMapped ~ AlignmentMethod, paired = TRUE)
# .y.           group1        group2        effsize    n1    n2 magnitude
# 1 PercentMapped Bismark       BismarkLocal    0.901     6     6 large    
# 2 PercentMapped Bismark       BisulfiteBolt   0.899     6     6 large    
# 3 PercentMapped Bismark       Biscuit         0.899     6     6 large    
# 4 PercentMapped Bismark       bwa meth        0.899     6     6 large    
# 5 PercentMapped BismarkLocal  BisulfiteBolt   0.899     6     6 large    
# 6 PercentMapped BismarkLocal  Biscuit         0.899     6     6 large    
# 7 PercentMapped BismarkLocal  bwa meth        0.899     6     6 large    
# 8 PercentMapped BisulfiteBolt Biscuit         0.899     6     6 large    
# 9 PercentMapped BisulfiteBolt bwa meth        0.899     6     6 large    
# 10 PercentMapped Biscuit       bwa meth        0.899     6     6 large



## Coral --------------------------------------------------


Coral <- dat %>% 
  filter(Organism == "Coral")


# ANOVA
ANOVA <- aov(PercentMapped ~ AlignmentMethod,
             data = Coral
)
ANOVA

# Call:
#   aov(formula = PercentMapped ~ AlignmentMethod, data = dat)
# 
#                 AlignmentMethod Residuals
# Sum of Squares        20.756208  0.262638
# Deg. of Freedom               4       214
# 
# Residual standard error: 0.03503255
# Estimated effects may be unbalanced
# 21 observations deleted due to missingness

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# Shapiro-Wilk normality test
# 
# data:  ANOVA$residuals
# W = 0.27825, p-value < 2.2e-16

#### Does not meet assumption of normality, moving to Friedman test

coral_wide <- Coral %>%
  select(SampleID, AlignmentMethod, PercentMapped) %>%
  pivot_wider(names_from = AlignmentMethod, 
              values_from = PercentMapped)

coral_complete <- coral_wide %>%
  filter(complete.cases(.)) %>%  
  pivot_longer(cols = -SampleID,
               names_to = "AlignmentMethod",
               values_to = "PercentMapped") %>%
  mutate(AlignmentMethod = factor(AlignmentMethod, 
                                  levels = c("Bismark", "BismarkLocal", 
                                             "BisulfiteBolt", "Biscuit", "bwa meth"))) %>% 
  filter(!is.na(AlignmentMethod)) 


## Friedman test
res.fried <- friedman_test(coral_complete, PercentMapped ~ AlignmentMethod |SampleID)
res.fried
#     .y.               n statistic    df        p method       
#   1 PercentMapped    30      119.     4 1.11e-24 Friedman test

coral_complete %>% friedman_effsize(PercentMapped ~ AlignmentMethod |SampleID)
#     .y.               n effsize method    magnitude
#   1 PercentMapped    30   0.988 Kendall W large

# post hoc test
pwc <- coral_complete %>%
  wilcox_test(PercentMapped ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
#     .y.           group1        group2           n1    n2 statistic             p        p.adj p.adj.signif
#   1 PercentMapped Bismark       BismarkLocal     30    30         0 0.00000181    0.0000181    ****        
#   2 PercentMapped Bismark       BisulfiteBolt    30    30         0 0.00000000186 0.0000000186 ****        
#   3 PercentMapped Bismark       Biscuit          30    30         0 0.00000000186 0.0000000186 ****        
#   4 PercentMapped Bismark       bwa meth         30    30         0 0.00000000186 0.0000000186 ****        
#   5 PercentMapped BismarkLocal  BisulfiteBolt    30    30         0 0.00000000186 0.0000000186 ****        
#   6 PercentMapped BismarkLocal  Biscuit          30    30         0 0.00000000186 0.0000000186 ****        
#   7 PercentMapped BismarkLocal  bwa meth         30    30         0 0.00000000186 0.0000000186 ****        
#   8 PercentMapped BisulfiteBolt Biscuit          30    30        32 0.00000514    0.0000514    ****        
#   9 PercentMapped BisulfiteBolt bwa meth         30    30         0 0.00000000186 0.0000000186 ****        
#   10 PercentMapped Biscuit       bwa meth         30    30         0 0.00000000186 0.0000000186 ****

coral_complete %>% wilcox_effsize(PercentMapped ~ AlignmentMethod, paired = TRUE)
#   .y.           group1        group2        effsize    n1    n2 magnitude
# 1 PercentMapped Bismark       BismarkLocal    0.873    30    30 large    
# 2 PercentMapped Bismark       BisulfiteBolt   0.873    30    30 large    
# 3 PercentMapped Bismark       Biscuit         0.873    30    30 large    
# 4 PercentMapped Bismark       bwa meth        0.873    30    30 large    
# 5 PercentMapped BismarkLocal  BisulfiteBolt   0.873    30    30 large    
# 6 PercentMapped BismarkLocal  Biscuit         0.873    30    30 large    
# 7 PercentMapped BismarkLocal  bwa meth        0.873    30    30 large    
# 8 PercentMapped BisulfiteBolt Biscuit         0.753    30    30 large    
# 9 PercentMapped BisulfiteBolt bwa meth        0.873    30    30 large    
# 10 PercentMapped Biscuit       bwa meth        0.873    30    30 large



## Stickbug  --------------------------------------------------


SB <- dat %>% 
  filter(Organism == "StickInsect")


# ANOVA
ANOVA <- aov(PercentMapped ~ AlignmentMethod,
             data = SB
)
ANOVA

# Call:
#   aov(formula = PercentMapped ~ AlignmentMethod, data = dat)
# 
#                 AlignmentMethod Residuals
# Sum of Squares         4.975236  0.015663
# Deg. of Freedom               4        90
# 
# Residual standard error: 0.01319234
# Estimated effects may be unbalanced
# 25 observations deleted due to missingness

# QQ plot
qqPlot(ANOVA$residuals,
       id = FALSE # id = FALSE to remove point identification
)

# histogram
hist(ANOVA$residuals)

# Shapiro-Wilk normality test
shapiro.test(ANOVA$residuals)
# Shapiro-Wilk normality test
# 
# data:  ANOVA$residuals
# W = 0.77721, p-value = 1.095e-10

#### Does not meet assumption of normality, moving to Friedman test

sb_wide <- SB %>%
  select(SampleID, AlignmentMethod, PercentMapped) %>%
  pivot_wider(names_from = AlignmentMethod, 
              values_from = PercentMapped)

sb_complete <- sb_wide %>%
  filter(complete.cases(.)) %>%  
  pivot_longer(cols = -SampleID,
               names_to = "AlignmentMethod",
               values_to = "PercentMapped") %>%
  mutate(AlignmentMethod = factor(AlignmentMethod, 
                                  levels = c("Bismark", "BismarkLocal", 
                                             "BisulfiteBolt", "Biscuit", "bwa meth"))) %>% 
  filter(!is.na(AlignmentMethod)) 

## Friedman test
res.fried <- friedman_test(sb_complete, PercentMapped ~ AlignmentMethod |SampleID)
res.fried
#     .y.               n statistic    df         p method       
#   1 PercentMapped     7        28     4 0.0000125 Friedman test

sb_complete %>% friedman_effsize(PercentMapped ~ AlignmentMethod |SampleID)
#     .y.               n effsize method    magnitude
#   1 PercentMapped     7       1 Kendall W large

# post hoc test
pwc <- sb_complete %>%
  wilcox_test(PercentMapped ~ AlignmentMethod, paired = TRUE, p.adjust.method = "bonferroni")
pwc
#   .y.           group1        group2           n1    n2 statistic     p p.adj p.adj.signif
# 1 PercentMapped Bismark       BismarkLocal      7     7         0 0.016 0.156 ns          
# 2 PercentMapped Bismark       BisulfiteBolt     7     7         0 0.016 0.156 ns          
# 3 PercentMapped Bismark       Biscuit           7     7         0 0.016 0.156 ns          
# 4 PercentMapped Bismark       bwa meth          7     7         0 0.016 0.156 ns          
# 5 PercentMapped BismarkLocal  BisulfiteBolt     7     7         0 0.016 0.156 ns          
# 6 PercentMapped BismarkLocal  Biscuit           7     7         0 0.016 0.156 ns          
# 7 PercentMapped BismarkLocal  bwa meth          7     7         0 0.016 0.156 ns          
# 8 PercentMapped BisulfiteBolt Biscuit           7     7         0 0.016 0.156 ns          
# 9 PercentMapped BisulfiteBolt bwa meth          7     7         0 0.016 0.156 ns          
# 10 PercentMapped Biscuit       bwa meth          7     7         0 0.016 0.156 ns

sb_complete %>% wilcox_effsize(PercentMapped ~ AlignmentMethod, paired = TRUE)
# .y.           group1        group2        effsize    n1    n2 magnitude
# 1 PercentMapped Bismark       BismarkLocal    0.894     7     7 large    
# 2 PercentMapped Bismark       BisulfiteBolt   0.894     7     7 large    
# 3 PercentMapped Bismark       Biscuit         0.894     7     7 large    
# 4 PercentMapped Bismark       bwa meth        0.894     7     7 large    
# 5 PercentMapped BismarkLocal  BisulfiteBolt   0.894     7     7 large    
# 6 PercentMapped BismarkLocal  Biscuit         0.894     7     7 large    
# 7 PercentMapped BismarkLocal  bwa meth        0.894     7     7 large    
# 8 PercentMapped BisulfiteBolt Biscuit         0.894     7     7 large    
# 9 PercentMapped BisulfiteBolt bwa meth        0.894     7     7 large    
# 10 PercentMapped Biscuit       bwa meth        0.894     7     7 large

