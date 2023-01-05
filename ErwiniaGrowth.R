library(growthcurver)
library(tidyverse)
library(dplyr)
library(tidyr)


# Read in CSVs
setwd("~/Erwinia2022/GrowthAssay")
temp <- list.files(pattern="*.csv")

for (i in 1:length(temp)){
  assign(temp[i], read.csv(temp[i]))
} 

temp

# Convert the "time" column from days to hours, omit NAs, 
# and trim the first 10 hours of abs data. It is wacky on Synergy machine, and still in lag phase. 
d1 <-`ErwiniaGrowthProtocol1.1 (Synergy).csv` %>% mutate(time = time*24) %>% filter(time > 10) %>% drop_na()
d2 <- `ErwiniaGrowthProtocol1.2 (Cytation).csv` %>% mutate(time = time*24) %>% filter(time > 10) %>% drop_na()
d3 <-`ErwiniaGrowthProtocol2 (Cytation).csv` %>% mutate(time = time*24) %>% filter(time > 10) %>% drop_na()
d4 <-`ErwiniaGrowthProtocol3 (Synergy).csv` %>% mutate(time = time*24) %>% filter(time > 10) %>% drop_na()


# Make sure that you have a column called "time"  .
# Now, we'll use Growthcurver to summarize the growth curve data for the entire plate 

gc_out.1 <- SummarizeGrowthByPlate(d1)
gc_out.2 <- SummarizeGrowthByPlate(d2)
gc_out.3 <- SummarizeGrowthByPlate(d3)
gc_out.4 <- SummarizeGrowthByPlate(d4)


# You can also summarize growth by plate while simultaneously making a pdf of your plots. 

gc_out1 <- d1 %>% SummarizeGrowthByPlate(plot_fit = TRUE, 
                                 plot_file = "gc_plot1.pdf")
gc_out2 <- d2 %>% SummarizeGrowthByPlate(plot_fit = TRUE, 
                                  plot_file = "gc_plot2.pdf")
gc_out3 <- d3 %>% SummarizeGrowthByPlate(plot_fit = TRUE, 
                                  plot_file = "gc_plot3.pdf")
gc_out4 <- d4 %>% SummarizeGrowthByPlate(plot_fit = TRUE, 
                                  plot_file = "gc_plot4.pdf")


## Testing differences in growth profiles between each machine
# ran identical plates to test for differences in abs between machines 
t.test(gc_out.1$r, gc_out.2$r, paired = TRUE)
t.test(gc_out.1$k, gc_out.2$k, paired = TRUE)

# Paired t-test
# data:  gc_out.1$r and gc_out.2$r
# t = 0.29425, df = 95, p-value = 0.7692
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#  -0.03930705  0.05298679
# sample estimates:
#  mean difference 
# 0.006839866 

# Paired t-test
# data:  gc_out.1$k and gc_out.2$k
# t = -0.1373, df = 95, p-value = 0.8911
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#   -477.5657  415.7825
# sample estimates:
# mean difference 
# -30.8916 

# Looks good for both r and K, no significant difference



# Let's create an output data frame to store the results in. 
# We'll create it so that it is the right size (it's faster this way!), 
# but leave it empty.
num_analyses <- length(names(d1)) - 1
d_gc <- data.frame(sample = character(num_analyses),
                   k = numeric(num_analyses),
                   n0  = numeric(num_analyses),
                   r = numeric(num_analyses),
                   t_mid = numeric(num_analyses),
                   t_gen = numeric(num_analyses),
                   auc_l = numeric(num_analyses),
                   auc_e = numeric(num_analyses),
                   sigma = numeric(num_analyses),
                   stringsAsFactors = FALSE)





# Check if Growthcurver provided any notes in a plate of growthcurves returned 
# from SummarizeGrowthByPlate
gc_out.1 %>% filter(note != "") 
gc_out.2 %>% filter(note != "") 
gc_out.3 %>% filter(note != "") 
gc_out.4 %>% filter(note != "") 



# Load dplyr and the sample output data
library(dplyr)

is.data.frame(gc_out.1) #TRUE
is.data.frame(gc_out.2) #TRUE
is.data.frame(gc_out.3) #TRUE
is.data.frame(gc_out.4) #TRUE




# Load ggplot2, and the sample data
library(ggplot2)

pca_gc_out.1 <- gc_out.1 
pca_gc_out.2 <- gc_out.2 
pca_gc_out.3 <- gc_out.3 
pca_gc_out.4 <- gc_out.4 


# Prepare the gc_out data for the PCA
rownames(pca_gc_out.1) <- pca_gc_out.1$sample
rownames(pca_gc_out.2) <- pca_gc_out.2$sample
rownames(pca_gc_out.3) <- pca_gc_out.3$sample
rownames(pca_gc_out.4) <- pca_gc_out.4$sample




# Perform the PCA
pca.res.1 <- prcomp(pca_gc_out.1 %>% select(k:sigma), center=TRUE, scale=TRUE)
pca.res.2 <- prcomp(pca_gc_out.2 %>% select(k:sigma), center=TRUE, scale=TRUE)
pca.res.3 <- prcomp(pca_gc_out.3 %>% select(k:sigma), center=TRUE, scale=TRUE)
pca.res.4 <- prcomp(pca_gc_out.4 %>% select(k:sigma), center=TRUE, scale=TRUE)


# Plot the results
as_data_frame(list(PC1=pca.res.1$x[,1],
                   PC2=pca.res.1$x[,2],
                   samples = pca_gc_out.1$sample)) %>% 
  ggplot(aes(x=PC1,y=PC2, label=samples)) + 
  geom_text(size = 3)


as_data_frame(list(PC1=pca.res.2$x[,1],
                   PC2=pca.res.2$x[,2],
                   samples = pca_gc_out.2$sample)) %>% 
  ggplot(aes(x=PC1,y=PC2, label=samples)) + 
  geom_text(size = 3)


as_data_frame(list(PC1=pca.res.3$x[,1],
                   PC2=pca.res.3$x[,2],
                   samples = pca_gc_out.3$sample)) %>% 
  ggplot(aes(x=PC1,y=PC2, label=samples)) + 
  geom_text(size = 3)



as_data_frame(list(PC1=pca.res.4$x[,1],
                   PC2=pca.res.4$x[,2],
                   samples = pca_gc_out.4$sample)) %>% 
  ggplot(aes(x=PC1,y=PC2, label=samples)) + 
  geom_text(size = 3)





# Get r, and K  values for each antagonist 
ant_table.1 <- gc_out.1 %>% 
  mutate(ant_name = case_when(sample== c("A1","B1","C1")~"Ant1_7.5",
                                                        sample== c("A2","B2","C2")~"Ant1_15",
                                                        sample== c("A3","B3","C3")~"Ant3_7.5",
                                                        sample== c("A4","B4","C4")~"Ant3_15",
                                                        sample== c("A5","B5","C5")~"Ant5_7.5",
                                                        sample== c("A6","B6","C6")~"Ant5_15",
                                                        sample== c("A7","B7","C7")~"Ant7_7.5",
                                                        sample== c("A8","B8","C8")~"Ant7_15",
                                                        sample== c("A9","B9","C9")~"Ant9_7.5",
                                                        sample== c("A10","B10","C10")~"Ant9_15",
                                                        sample== c("A11","B11","C11")~"Ant11_7.5",
                                                        sample== c("A12","B12","C12")~"Ant11_15",
                                                        sample== c("D1","E1","F1")~"Ant2_7.5",
                                                        sample== c("D2","E2","F2")~"Ant2_15",
                                                        sample== c("D3","E3","F3")~"Ant4_7.5",
                                                        sample== c("D4","E4","F4")~"Ant4_15",
                                                        sample== c("D5","E5","F5")~"Ant6_7.5",
                                                        sample== c("D6","E6","F6")~"Ant6_15",
                                                        sample== c("D7","E7","F7")~"Ant8_7.5",
                                                        sample== c("D8","E8","F8")~"Ant8_15",
                                                        sample== c("D9","E9","F9")~"Ant10_7.5",
                                                        sample== c("D10","E10","F10")~"Ant10_15",
                                                        sample== c("D11","E11","F11")~"Ant12_7.5",
                                                        sample== c("D12","E12","F12")~"Ant12_15",
                                                        sample== c("G1","G2","G3")~"Ant13_7.5",
                                                        sample== c("H1","H2","H3")~"Ant13_15",
                                                        sample== c("G4","G5","G6")~"Ant14_7.5",
                                                        sample== c("H4","H5","H6")~"Ant14_15",
                                                        sample== c("G7","G8","G9")~"Ant15_7.5",
                                                        sample== c("H7","H8","H9")~"Ant15_15", sample== c("G10","G11","G12","H10","H11", "H12") ~ "Control")) %>%
  filter(note != "questionable fit (k < n0)", note != "questionable fit") %>%
  group_by(ant_name) %>% summarize(avg_growth_rate = mean(r), avg_carry_cap = mean(k))
 


                                                       
ant_table.2 <- gc_out.2 %>% 
  mutate(ant_name = case_when(sample== c("A1","B1","C1")~"Ant1_7.5",
                                                        sample== c("A2","B2","C2")~"Ant1_15",
                                                        sample== c("A3","B3","C3")~"Ant3_7.5",
                                                        sample== c("A4","B4","C4")~"Ant3_15",
                                                        sample== c("A5","B5","C5")~"Ant5_7.5",
                                                        sample== c("A6","B6","C6")~"Ant5_15",
                                                        sample== c("A7","B7","C7")~"Ant7_7.5",
                                                        sample== c("A8","B8","C8")~"Ant7_15",
                                                        sample== c("A9","B9","C9")~"Ant9_7.5",
                                                        sample== c("A10","B10","C10")~"Ant9_15",
                                                        sample== c("A11","B11","C11")~"Ant11_7.5",
                                                        sample== c("A12","B12","C12")~"Ant11_15",
                                                        sample== c("D1","E1","F1")~"Ant2_7.5",
                                                        sample== c("D2","E2","F2")~"Ant2_15",
                                                        sample== c("D3","E3","F3")~"Ant4_7.5",
                                                        sample== c("D4","E4","F4")~"Ant4_15",
                                                        sample== c("D5","E5","F5")~"Ant6_7.5",
                                                        sample== c("D6","E6","F6")~"Ant6_15",
                                                        sample== c("D7","E7","F7")~"Ant8_7.5",
                                                        sample== c("D8","E8","F8")~"Ant8_15",
                                                        sample== c("D9","E9","F9")~"Ant10_7.5",
                                                        sample== c("D10","E10","F10")~"Ant10_15",
                                                        sample== c("D11","E11","F11")~"Ant12_7.5",
                                                        sample== c("D12","E12","F12")~"Ant12_15",
                                                        sample== c("G1","G2","G3")~"Ant13_7.5",
                                                        sample== c("H1","H2","H3")~"Ant13_15",
                                                        sample== c("G4","G5","G6")~"Ant14_7.5",
                                                        sample== c("H4","H5","H6")~"Ant14_15",
                                                        sample== c("G7","G8","G9")~"Ant15_7.5",
                                                        sample== c("H7","H8","H9")~"Ant15_15", sample== c("G10","G11","G12","H10","H11", "H12") ~ "Control"))%>%
filter(note != "questionable fit (k < n0)", note != "questionable fit") %>%
  group_by(ant_name) %>% summarize(avg_growth_rate = mean(r), avg_carry_cap = mean(k))




ant_table.3 <- gc_out.3 %>% 
  mutate(ant_name = case_when(sample== c("A1","B1","C1")~"Ant1_7.5",
                                                        sample== c("A2","B2","C2")~"Ant1_15",
                                                        sample== c("A3","B3","C3")~"Ant3_7.5",
                                                        sample== c("A4","B4","C4")~"Ant3_15",
                                                        sample== c("A5","B5","C5")~"Ant5_7.5",
                                                        sample== c("A6","B6","C6")~"Ant5_15",
                                                        sample== c("A7","B7","C7")~"Ant7_7.5",
                                                        sample== c("A8","B8","C8")~"Ant7_15",
                                                        sample== c("A9","B9","C9")~"Ant9_7.5",
                                                        sample== c("A10","B10","C10")~"Ant9_15",
                                                        sample== c("A11","B11","C11")~"Ant11_7.5",
                                                        sample== c("A12","B12","C12")~"Ant11_15",
                                                        sample== c("D1","E1","F1")~"Ant2_7.5",
                                                        sample== c("D2","E2","F2")~"Ant2_15",
                                                        sample== c("D3","E3","F3")~"Ant4_7.5",
                                                        sample== c("D4","E4","F4")~"Ant4_15",
                                                        sample== c("D5","E5","F5")~"Ant6_7.5",
                                                        sample== c("D6","E6","F6")~"Ant6_15",
                                                        sample== c("D7","E7","F7")~"Ant8_7.5",
                                                        sample== c("D8","E8","F8")~"Ant8_15",
                                                        sample== c("D9","E9","F9")~"Ant10_7.5",
                                                        sample== c("D10","E10","F10")~"Ant10_15",
                                                        sample== c("D11","E11","F11")~"Ant12_7.5",
                                                        sample== c("D12","E12","F12")~"Ant12_15",
                                                        sample== c("G1","G2","G3")~"Ant13_7.5",
                                                        sample== c("H1","H2","H3")~"Ant13_15",
                                                        sample== c("G4","G5","G6")~"Ant14_7.5",
                                                        sample== c("H4","H5","H6")~"Ant14_15",
                                                        sample== c("G7","G8","G9")~"Ant15_7.5",
                                                        sample== c("H7","H8","H9")~"Ant15_15", sample== c("G10","G11","G12","H10","H11", "H12") ~ "Control"))%>%
  filter(note != "questionable fit (k < n0)", note != "questionable fit") %>%
  group_by(ant_name) %>% summarize(avg_growth_rate = mean(r), avg_carry_cap = mean(k))


ant_table.4 <- gc_out.4 %>% 
  mutate(ant_name = case_when(sample== c("A1","B1","C1")~"Ant1_7.5",
                                                        sample== c("A2","B2","C2")~"Ant1_15",
                                                        sample== c("A3","B3","C3")~"Ant3_7.5",
                                                        sample== c("A4","B4","C4")~"Ant3_15",
                                                        sample== c("A5","B5","C5")~"Ant5_7.5",
                                                        sample== c("A6","B6","C6")~"Ant5_15",
                                                        sample== c("A7","B7","C7")~"Ant7_7.5",
                                                        sample== c("A8","B8","C8")~"Ant7_15",
                                                        sample== c("A9","B9","C9")~"Ant9_7.5",
                                                        sample== c("A10","B10","C10")~"Ant9_15",
                                                        sample== c("A11","B11","C11")~"Ant11_7.5",
                                                        sample== c("A12","B12","C12")~"Ant11_15",
                                                        sample== c("D1","E1","F1")~"Ant2_7.5",
                                                        sample== c("D2","E2","F2")~"Ant2_15",
                                                        sample== c("D3","E3","F3")~"Ant4_7.5",
                                                        sample== c("D4","E4","F4")~"Ant4_15",
                                                        sample== c("D5","E5","F5")~"Ant6_7.5",
                                                        sample== c("D6","E6","F6")~"Ant6_15",
                                                        sample== c("D7","E7","F7")~"Ant8_7.5",
                                                        sample== c("D8","E8","F8")~"Ant8_15",
                                                        sample== c("D9","E9","F9")~"Ant10_7.5",
                                                        sample== c("D10","E10","F10")~"Ant10_15",
                                                        sample== c("D11","E11","F11")~"Ant12_7.5",
                                                        sample== c("D12","E12","F12")~"Ant12_15",
                                                        sample== c("G1","G2","G3")~"Ant13_7.5",
                                                        sample== c("H1","H2","H3")~"Ant13_15",
                                                        sample== c("G4","G5","G6")~"Ant14_7.5",
                                                        sample== c("H4","H5","H6")~"Ant14_15",
                                                        sample== c("G7","G8","G9")~"Ant15_7.5",
                                                        sample== c("H7","H8","H9")~"Ant15_15", sample== c("G10","G11","G12","H10","H11", "H12") ~ "Control"))%>%
filter(note != "questionable fit (k < n0)", note != "questionable fit")%>% 
  group_by(ant_name) %>% summarize(avg_growth_rate = mean(r), avg_carry_cap = mean(k))



