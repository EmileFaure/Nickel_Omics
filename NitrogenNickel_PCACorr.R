require(PermCor)
require(vegan)

#### Investigating the relationship between nitrogen-related variables and Ni measures in the ACE dataset ####

### Load data ###

meta <- read.table(file = "/Users/emifaure/Documents/ACE/MetaData/ACE_metadata/Metadata_For_Submission/meta_CTD_ForSubmission.csv", 
                   header=TRUE, sep = ";", dec = ".", stringsAsFactors = FALSE)
meta <- meta %>% 
  mutate_if(.predicate = is.character, .funs = as.factor)

# Load list of metagenomes
metagenomes <- read.table("/Users/emifaure/Documents/ACE/Metagenomics/ACEsamples_With_Ace_seq_name.tsv", header = TRUE)

# Load Nickel dataset
Nickel=read.table("/Users/emifaure/Documents/ACE/NolwennCollab/NickelData.txt", sep="\t", header=T, dec=",")
meta_nickel = meta %>% inner_join(Nickel, by=c("TM_station_number","Depth_m")) %>%
  filter(!is.na(Ni_60_58_D_DELTA_BOTTLE))
# 85 ACE_seq_Name have a nickel values, how many in out sequenced metagenomes ?

meta_nickel = meta_nickel[which(meta_nickel$ACE_seq_name %in% metagenomes$ace_seq_name),]
# we have 48 metagenomes that can be associated with an isotopic value

### Pre-process data ###

meta_nickel_rda = select(meta_nickel,
                            c("ACE_seq_name",
                              "NOx_µmol.L", "Nitrite_µmol.L", "Ammonium_µmol.L", "Nitrate_µmol.L", "PON.µmol", "d15N", 
                              "Ni_60_58_D_DELTA_BOTTLE", "Ni_D_CONC_BOTTLE"))

row.names(meta_nickel_rda)=meta_nickel_rda$ACE_seq_name
meta_nickel_rda=select(meta_nickel_rda,-c("ACE_seq_name"))

# Water samples for which multiple filters were taken are duplicated lines in this table, we need to remove them
meta_nickel_rda <- distinct(meta_nickel_rda)
meta_clean <- na.omit(meta_nickel_rda)
meta_clean=data.frame(scale(meta_clean))

### PCA analysis ###

PCA_nickel <- rda(meta_clean, na.action = "na.omit")

samples = PCA_nickel$CA$u
samples = as.data.frame(samples)
samples = merge(samples,meta_nickel, by.x="row.names", by.y="ACE_seq_name")
row.names(samples)=samples[,1]
samples=samples[,-1]

# Retrieve envi scores
enviscore=PCA_nickel$CA$v

eig_vals <- PCA_nickel$CA$eig
var_exp <- round(100 * eig_vals / sum(eig_vals), 1)

#Plot the triplot
ggplot() +   geom_hline(yintercept = 0, linetype='dotted') +
  geom_vline(xintercept = 0, linetype='dotted') +
  labs(x = paste0("PCA1 (",
                  var_exp[1], "%)"),
       y = paste0("PCA2 (",
                  var_exp[2], "%)")) +
  theme(plot.title=element_text(hjust=0.5)) +
  geom_point(aes(samples[,1], samples[,2], col=samples$Ni_60_58_D_DELTA_BOTTLE), size = 5, alpha = 0.9) +
  geom_text(aes(samples[,1], samples[,2], label=samples$TM_station_number), col="white", size=3) +
  geom_segment(data=as.data.frame(enviscore),aes(xend = PC1, yend = PC2),x=0,y=0,linewidth = 0.5, linetype="F1",color = 'orange',arrow = arrow(length = unit(0.2,"cm"))) +
  geom_text(aes(enviscore[,1]*1.05, enviscore[,2]*1.05, label=rownames(enviscore)), color="orange") +
  labs(col = expression(delta^60 * Ni ~ "(‰)")) +
  theme_bw() 

write.table(meta_nickel_rda, "Documents/ACE/NolwennCollab/Post_Review/Explore_NitrogenCycle/Metadata_Nickel_Nitrogen.tsv", sep="\t", quote=F, row.names = T)

### Regression/correlation analysis ###

ggplot(data=meta_nickel_rda, aes(x=PON.µmol, y=Ni_60_58_D_DELTA_BOTTLE)) +
  geom_point() +
  geom_smooth(method="lm", col="black") +
  labs(x="Particulate Organic Nitrogen (µmol.L-1)", y="∂60Ni (‰)") +
  theme_bw() 
summary(lm(Ni_60_58_D_DELTA_BOTTLE~PON.µmol, data=meta_nickel_rda))
cor(meta_clean$Ni_60_58_D_DELTA_BOTTLE,meta_clean$PON.µmol, method = "spearman")
perm_test(meta_clean$Ni_60_58_D_DELTA_BOTTLE,meta_clean$PON.µmol,method = "Spearman", alternative = "two.sided", B=1000)$p.value

ggplot(data=meta_nickel_rda, aes(x=Nitrate_µmol.L, y=Ni_60_58_D_DELTA_BOTTLE)) +
  geom_point() +
  geom_smooth(method="lm", col="black") +
  labs(x="Nitrate (µmol.L-1)", y="∂60Ni (‰)") +
  theme_bw() 
summary(lm(Ni_60_58_D_DELTA_BOTTLE~Nitrate_µmol.L, data=meta_nickel_rda))
cor(meta_clean$Ni_60_58_D_DELTA_BOTTLE,meta_clean$Nitrate_µmol.L, method = "spearman")
perm_test(meta_clean$Ni_60_58_D_DELTA_BOTTLE,meta_clean$Nitrate_µmol.L,method = "Spearman", alternative = "two.sided", B=1000)$p.value
