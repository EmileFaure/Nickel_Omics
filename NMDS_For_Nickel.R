library(vegan)
library(tidyverse)
library(viridis)
library(patchwork)

#Attached
NMDSData = readRDS("/Users/emifaure/Library/Mobile Documents/com~apple~CloudDocs/Documents/Figure_NMDS_ACE/AGN_NMDS_ATT_WOOutlier_T60MAX_NZV20.rds")

Nickel=read.table("/Users/emifaure/Documents/ACE/NolwennCollab/NickelData.txt", sep="\t", header=T, dec=",")

meta <- read.table(file = "/Users/emifaure/Documents/ACE/MetaData/ACE_metadata/Metadata_For_Submission/meta_CTD_ForSubmission.csv", 
                   header=TRUE, sep = ";", dec = ".", stringsAsFactors = FALSE)
meta <- meta %>% 
  mutate_if(.predicate = is.character, .funs = as.factor)

meta_nickel = meta %>% inner_join(Nickel, by=c("TM_station_number","Depth_m")) %>%
  filter(!is.na(Ni_60_58_D_DELTA_BOTTLE))

df.NMDS <- scores(NMDSData, display = "sites") %>% as.data.frame()

df.NMDS.nickel <- merge(df.NMDS,meta_nickel,by.x="row.names",by.y="ACE_seq_name")
df.NMDS.nickel$Size_fraction <- factor(df.NMDS.nickel$Size_fraction, levels = c("0.2-3 µm", "0.2-40 µm"))
df.NMDS.nickel$TM_station_number = as.factor(df.NMDS.nickel$TM_station_number)

df.NMDS <- merge(df.NMDS,meta,by.x="row.names",by.y="ACE_seq_name")
df.NMDS$Size_fraction <- factor(df.NMDS$Size_fraction, levels = c("0.2-3 µm", "0.2-40 µm", "> 3µm"))


plot.NMDS.att <- ggplot() +
  geom_point(data = df.NMDS,
             aes(x = NMDS1, y = NMDS2, shape = MertzGlacier),
             size=6, alpha=0.2, col="grey15") +
  geom_point(data = df.NMDS.nickel,
             aes(x = NMDS1, y = NMDS2, shape = MertzGlacier, 
                col=TM_station_number), size=6) +
  labs(col = "Station", shape = "Mertz Glacier Sample") +
  scale_color_manual(values=c("orange", "springgreen1", "springgreen3",
                              "brown4", "darkblue", "deepskyblue2",
                              "lightblue","red")) +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.box.background= element_rect(colour="black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  guides(shape = guide_legend(ncol=2))
plot.NMDS.att

#FreeLiving
NMDSDatafl = readRDS("/Users/emifaure/Library/Mobile Documents/com~apple~CloudDocs/Documents/Figure_NMDS_ACE/AGN_NMDS_FL_T60MAX_NZV20.rds")

Nickel=read.table("/Users/emifaure/Documents/ACE/NolwennCollab/NickelData.txt", sep="\t", header=T, dec=",")

meta <- read.table(file = "/Users/emifaure/Documents/ACE/MetaData/ACE_metadata/Metadata_For_Submission/meta_CTD_ForSubmission.csv", 
                   header=TRUE, sep = ";", dec = ".", stringsAsFactors = FALSE)
meta <- meta %>% 
  mutate_if(.predicate = is.character, .funs = as.factor)

meta_nickel = meta %>% inner_join(Nickel, by=c("TM_station_number","Depth_m")) %>%
  filter(!is.na(Ni_60_58_D_DELTA_BOTTLE))

df.NMDS <- scores(NMDSDatafl, display = "sites") %>% as.data.frame()

df.NMDS.nickel <- merge(df.NMDS,meta_nickel,by.x="row.names",by.y="ACE_seq_name")
df.NMDS.nickel$Size_fraction <- factor(df.NMDS.nickel$Size_fraction, levels = c("0.2-3 µm", "0.2-40 µm"))
df.NMDS.nickel$TM_station_number = as.factor(df.NMDS.nickel$TM_station_number)

df.NMDS <- merge(df.NMDS,meta,by.x="row.names",by.y="ACE_seq_name")
df.NMDS$Size_fraction <- factor(df.NMDS$Size_fraction, levels = c("0.2-3 µm", "0.2-40 µm"))


plot.NMDS.fl <- ggplot() +
  geom_point(data = df.NMDS,
             aes(x = NMDS1, y = NMDS2, shape = MertzGlacier, size=Size_fraction),
             alpha=0.2, col="grey15") +
  geom_point(data = df.NMDS.nickel,
             aes(x = NMDS1, y = NMDS2, shape = MertzGlacier, size=Size_fraction, col=TM_station_number)) +
  labs(col = "Station", shape = "Mertz Glacier Sample") +
  scale_color_manual(values=c("orange", "springgreen1", "springgreen3",
                              "brown4", "darkblue", "deepskyblue2",
                              "lightblue","red")) +
  scale_size_manual(values=c(3,5)) +
  theme_bw() +
  theme(legend.background = element_blank(),
        legend.box.background= element_rect(colour="black"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) +
  guides(shape = guide_legend(ncol=2))
plot.NMDS.fl

plot.NMDS.fl | plot.NMDS.att
