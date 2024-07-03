# Code by Lucas Ustick to analyze nifH data

# import tree data ####
library(ggtree)

u_tree <- read.tree("/nifH_microdiversity/phylogenetic_analysis/tree/itol_ucyna_newick.txt")
t_tree <- read.tree("/nifH_microdiversity/phylogenetic_analysis/tree/itol_tricho_newick.txt")

u_labels <- u_tree$tip.label
t_labels <- t_tree$tip.label

u_metadata <- data.frame(u_labels)
t_metadata <- data.frame(t_labels)

# identify each clade
u_metadata$clade <- "0"
t_metadata$clade <- "0"

# first get index of start and stop points
u1_s <- which(u_metadata$u_labels == "Cornejo-Castillo_OTU04")
u1_e <- which(u_metadata$u_labels == "Cornejo-Castillo_OTU00")

u2_s <- which(u_metadata$u_labels == "Cornejo-Castillo_OTU06")
u2_e <- which(u_metadata$u_labels == "UCYNA-nifH-ASV0040")

u3_s <- which(u_metadata$u_labels == "UCYNA-nifH-ASV0033")
u3_e <- which(u_metadata$u_labels == "UCYNA-nifH-ASV0024")

u4_s <- which(u_metadata$u_labels == "Cornejo-Castillo_OTU03")
u4_e <- which(u_metadata$u_labels == "UCYNA-nifH-ASV0058")

u5_s <- which(u_metadata$u_labels == "UCYNA-nifH-ASV0045")
u5_e <- which(u_metadata$u_labels == "UCYNA-nifH-ASV0047")

u6_s <- which(u_metadata$u_labels == "Turk-Kubo_oligo24_A6")
u6_e <- which(u_metadata$u_labels == "UCYNA-nifH-ASV0072")

# first get index of start and stop points
t1_s <- which(t_metadata$t_labels == "Tricho-nifH-ASV0086")
t1_e <- which(t_metadata$t_labels == "Tricho-nifH-ASV0025")

t2_s <- which(t_metadata$t_labels == "Tricho-nifH-ASV0026")
t2_e <- which(t_metadata$t_labels == "Tricho-nifH-ASV0010")

t3_s <- which(t_metadata$t_labels == "Tricho-nifH-ASV0018")
t3_e <- which(t_metadata$t_labels == "Tricho-nifH-ASV0052")

t4_s <- which(t_metadata$t_labels == "Tricho-nifH-ASV0013")
t4_e <- which(t_metadata$t_labels == "Tricho-nifH-ASV0017")

# assign the clade to each
u_metadata$clade[u1_s:u1_e] <- "1"
u_metadata$clade[u2_s:u2_e] <- "2"
u_metadata$clade[u3_s:u3_e] <- "3"
u_metadata$clade[u4_s:u4_e] <- "4"
u_metadata$clade[u5_s:u5_e] <- "5"
u_metadata$clade[u6_s:u6_e] <- "6"
u_metadata$clade <- factor(u_metadata$clade)

# assign the clade to each
t_metadata$clade[t1_s:t1_e] <- "I"
t_metadata$clade[t2_s:t2_e] <- "UIV"
t_metadata$clade[t3_s:t3_e] <- "III"
t_metadata$clade[t4_s:t4_e] <- "UII"
t_metadata$clade <- factor(t_metadata$clade)

# combine tree data, metadata, & observations ####

u_asv <- read.csv("/nifH_microdiversity/data_files/sdata2_UCYNA_ASV_24.01.11.csv")
t_asv <- read.csv("/nifH_microdiversity/data_files/sdata1_Tricho_ASV_24.01.11.csv")

u_full_OTU <- merge(u_asv,u_metadata,by.x="ASV_ID",by.y="u_labels",all=FALSE)
u_full_OTU <- u_full_OTU[,c(1,2,dim(u_full_OTU)[2],3:(dim(u_full_OTU)[2]-1))]

t_full_OTU <- merge(t_asv,t_metadata,by.x="ASV_ID",by.y="t_labels",all=FALSE)
t_full_OTU <- t_full_OTU[,c(1,2,dim(u_full_OTU)[2],3:(dim(u_full_OTU)[2]-1))]

# create relative abundance for each clade
t_clade <- t(t_full_OTU[1:5,5:dim(u_full_OTU)[2]])
colnames(t_clade)<- c("Tricho_I","Tricho_III","Tricho_UIV","Tricho_UII","Tricho_Total")
t_clade[] <- 0

u_clade <- t(u_full_OTU[1:7,5:dim(u_full_OTU)[2]])
colnames(u_clade)<- c("UCYN-A_1","UCYN-A_2","UCYN-A_3","UCYN-A_4","UCYN-A_5","UCYN-A_6","UCYN-A_Total")
u_clade[] <- 0

for (i in 1:dim(t_clade)[1]) {
  t_clade[i,1] <- sum(t_full_OTU[t_full_OTU$clade=="I",i+4])
  t_clade[i,2] <- sum(t_full_OTU[t_full_OTU$clade=="III",i+4])
  t_clade[i,3] <- sum(t_full_OTU[t_full_OTU$clade=="UIV",i+4])
  t_clade[i,4] <- sum(t_full_OTU[t_full_OTU$clade=="UII",i+4])
  t_clade[i,5] <- sum(t_full_OTU[,i+4])
  
  u_clade[i,1] <- sum(u_full_OTU[u_full_OTU$clade=="1",i+4])
  u_clade[i,2] <- sum(u_full_OTU[u_full_OTU$clade=="2",i+4])
  u_clade[i,3] <- sum(u_full_OTU[u_full_OTU$clade=="3",i+4])
  u_clade[i,4] <- sum(u_full_OTU[u_full_OTU$clade=="4",i+4])
  u_clade[i,5] <- sum(u_full_OTU[u_full_OTU$clade=="5",i+4])
  u_clade[i,6] <- sum(u_full_OTU[u_full_OTU$clade=="6",i+4])
  u_clade[i,7] <- sum(u_full_OTU[,i+4])
  
  
}

clade_master <- merge(t_clade,u_clade,by="row.names")
clade_master$nifH_total <- rowSums(clade_master[, c("Tricho_Total", "UCYN-A_Total")])

# add in lat lon data
sample_metadata <- read.csv("/nifH_microdiversity/data_files/sdata3_sample_metadata_24.01.11.csv")

clade_master <- merge(clade_master,sample_metadata[,1:3],by.x="Row.names",by.y="Sample_ID")

# calculate relative abundances
clade_relative_abundance <- clade_master[, 2:5] / clade_master[, 6] #tricho
clade_relative_abundance[,5:10] <- clade_master[, 7:12] / clade_master[, 13] #UCYNA
clade_relative_abundance[,11:12] <- clade_master[,c(6,13)] / clade_master[, 14] #TvA
clade_relative_abundance[,13:15] <- clade_master[,c(1,15:16)]

# Figure 1 sample map ####

library(ggplot2)
library(rnaturalearthdata)

theme1 <-theme(panel.background = element_blank() #remove background
               ,panel.border=element_rect(fill=NA,size=1) #make box around axis
               ,panel.grid.major = element_blank()
               ,panel.grid.minor = element_blank()
               ,strip.background=element_blank()
               ,axis.text.x=element_text(colour="black") #make all text black
               ,axis.text.y=element_text(colour="black") #make all text black
               ,axis.ticks=element_line(colour="black") #make all text black
               ,plot.title = element_text(hjust = 0.5) #center the title
)

size_obj <- 3

col_v <- c("#ffff99","#fdbf6f","#ff7f00","#b15928","#cab2d6","#6a3d9a","#a6cee3","#1f78b4","#b2df8a","#33a02c","#175014")

ggplot() +
  geom_path(data = rnaturalearthdata::coastline50, aes(x = long, y = lat, group = group)) +
  geom_point(data = sample_metadata[!is.na(clade_relative_abundance$Tricho_Total),], aes(x = lon_E, y = lat_N, color = Cruise), size = size_obj) + 
  xlim(-180,180) + ylim(-90,90) +
  coord_cartesian(xlim =c(-165,165), ylim = c(-77,77)) +
  #scale_colour_continuous(type = "viridis") +
  #scale_color_brewer(palette = "Paired") +
  scale_color_manual(values=col_v) +
  theme1

#ggsave("/nifH_microdiversity/images/figure_1/figure_1_24.01.11.pdf",width=12,height=6)

# Figure 2C&D Top ASV barplot ####
library(ggplot2)

# percent of total reads
# Tricho top 100
100*( sum(t_asv$ASV_count[1:100])/sum(t_asv$ASV_count[1:dim(t_asv)[1]]) )
# UCYN-A top 100
100*( sum(u_asv$ASV_count[1:100])/sum(u_asv$ASV_count[1:dim(t_asv)[1]]) )
# Tricho most abundant ASV's
100*( t_asv$ASV_count[1]/sum(t_asv$ASV_count[1:dim(t_asv)[1]]) )
100*( t_asv$ASV_count[2]/sum(t_asv$ASV_count[1:dim(t_asv)[1]]) )
# Tricho 1 & 2
100*( sum(t_full_OTU$ASV_count[t_full_OTU$clade == "I" | t_full_OTU$clade == "UII"]) / sum(t_asv$ASV_count[1:dim(t_asv)[1]]) )
# UCYN-A most abundant ASV
100*( u_asv$ASV_count[1]/sum(u_asv$ASV_count[1:dim(t_asv)[1]]) )
100*( u_asv$ASV_count[2]/sum(u_asv$ASV_count[1:dim(t_asv)[1]]) )
100*( u_asv$ASV_count[3]/sum(u_asv$ASV_count[1:dim(t_asv)[1]]) )
100*( u_asv$ASV_count[4]/sum(u_asv$ASV_count[1:dim(t_asv)[1]]) )
# UCYN-A top asv out of total A1
100*( u_asv$ASV_count[1] / sum(u_full_OTU$ASV_count[u_full_OTU$clade == 1]) )
# UCYN-A clades 4/5/6 abundance
100*( sum(u_full_OTU$ASV_count[u_full_OTU$clade == "4" | u_full_OTU$clade == "5" | u_full_OTU$clade == "6" | u_full_OTU$clade == "0"]) / sum(u_asv$ASV_count[1:dim(u_asv)[1]]) )

#make custom colors for each
col_v1 <- c("#cab2d6","#6a3d9a","#fdbf6f","#ff7f00")
col_v2 <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c")

t_full_OTU$clade <- factor(t_full_OTU$clade, levels=c("I","UII","III","UIV"))

p1 <- ggplot(data=t_full_OTU[1:25,], aes(x=(1:25),y=ASV_count,fill=clade)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=col_v1) +
  labs(title="Top Trichodesmium ASV's",x ="ASV", y = "count") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme_classic() +
  theme(legend.position = c(0.9, 0.5))

p2 <- ggplot(data=u_full_OTU[1:25,], aes(x=(1:25),y=ASV_count,fill=clade)) +
  geom_bar(stat="identity") + 
  scale_fill_manual(values=col_v2) +
  labs(title="Top UCYN-A ASV's",x ="ASV", y = "count") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme_classic() +
  theme(legend.position = c(0.9, 0.5))

library(ggpubr)
ggarrange(
  p1,p2, labels = c("C", "D")
)

#ggsave("/nifH_microdiversity/images/figure_2/figure_2B_24.01.10.pdf",width=12,height=3)

# Figure 2 Diversity indexes ####

library(vegan)

t_asv_meta_df <- merge(sample_metadata , t(t_asv[1:100,4:dim(t_asv)[2]]) , by.x="Sample_ID",by.y="row.names")
u_asv_meta_df <- merge(sample_metadata , t(u_asv[1:100,4:dim(u_asv)[2]]) , by.x="Sample_ID",by.y="row.names")

complete_otu_t <- !apply(t_asv_meta_df[, 15:dim(t_asv_meta_df)[2]], 1, function(x) all(x == 0))
complete_otu_u <- !apply(u_asv_meta_df[, 15:dim(u_asv_meta_df)[2]], 1, function(x) all(x == 0))

t_shannon <- diversity(t_asv_meta_df[ complete_otu_t , 15:dim(t_asv_meta_df)[2] ])
u_shannon <- diversity(u_asv_meta_df[ complete_otu_u , 15:dim(u_asv_meta_df)[2] ])

shannon_df <- data.frame(species=character(),
                      shannon=numeric())

shannon_df[1:length(t_shannon),2] <- t_shannon
shannon_df[1:length(t_shannon),1] <- "Trichodesmium"
shannon_df[(length(t_shannon)+1):(length(t_shannon)+length(u_shannon)),2] <- u_shannon
shannon_df[(length(t_shannon)+1):(length(t_shannon)+length(u_shannon)),1] <- "UCYN-A"


library(ggplot2)

theme1 <-theme(panel.background = element_blank() #remove background
               ,panel.border=element_rect(fill=NA,size=1) #make box around axis
               ,panel.grid.major = element_blank()
               ,panel.grid.minor = element_blank()
               ,strip.background=element_blank()
               ,axis.text.x=element_text(colour="black") #make all text black
               ,axis.text.y=element_text(colour="black") #make all text black
               ,axis.ticks=element_line(colour="black") #make all text black
               ,plot.title = element_text(hjust = 0.5) #center the title
)

library(ggpubr)

p3 <- ggplot(shannon_df, aes(x=species, y=shannon)) + 
  geom_violin() +
  labs(x ="species",y ="Shannon index")+
  stat_summary(fun.y=mean, geom="point", size=2, color="darkred")+
  stat_compare_means(method = "t.test", label.y = 2, paired = FALSE) +
  theme1


ggarrange(
  p1,p2,p3, labels = c("C", "D", "E"),
  ncol=3,nrow=1, widths = c(1, 1,0.5)
)


#ggsave("/nifH_microdiversity/images/figure_2/figure_2CDE_24.01.19.pdf",width=8.5,height=2)

# Figure 3 biogeography ####

library(ggplot2)
library(rnaturalearthdata)

theme1 <-theme(panel.background = element_blank() #remove background
               ,panel.border=element_rect(fill=NA,size=1) #make box around axis
               ,panel.grid.major = element_blank()
               ,panel.grid.minor = element_blank()
               ,strip.background=element_blank()
               ,axis.text.x=element_text(colour="black") #make all text black
               ,axis.text.y=element_text(colour="black") #make all text black
               ,axis.ticks=element_line(colour="black") #make all text black
               ,plot.title = element_text(hjust = 0.5) #center the title
)

size_obj <- 2

t0 <- ggplot() +
  geom_path(data = rnaturalearthdata::coastline50, aes(x = long, y = lat, group = group)) +
  geom_point(data = clade_relative_abundance[!is.na(clade_relative_abundance$Tricho_Total),], aes(x = lon_E, y = lat_N, color = Tricho_Total), size = size_obj) + 
  xlim(-180,180) + ylim(-90,90) +
  coord_cartesian(xlim =c(-165,165), ylim = c(-77,77)) +
  labs(title="Trichodesmium vs. UCYN-A relative abundance",x ="longitude",y ="latitude")+
  #scale_colour_continuous(type = "viridis") +
  scale_color_distiller(palette = "RdBu") +
  theme1

library(ggpubr)
ggarrange(
  t0, labels = c("A"))

#ggsave("/nifH_microdiversity/images/figure_3/figure_3A_24.01.19.pdf",width=8.5,height=4)

theme1 <-theme(panel.background = element_blank() #remove background
               ,panel.border=element_rect(fill=NA,size=1) #make box around axis
               ,panel.grid.major = element_blank()
               ,panel.grid.minor = element_blank()
               ,strip.background=element_blank()
               ,plot.title = element_text(hjust = 0.5)  #center the title
               ,axis.text.x=element_blank() #remove x axis labels
               ,axis.ticks.x=element_blank() #remove x axis ticks
               ,axis.text.y=element_blank()  #remove y axis labels
               ,axis.ticks.y=element_blank()  #remove y axis ticks
               ,axis.title.x = element_blank()
               ,axis.title.y = element_blank()
              
)


size_obj <- 0.5

t1 <- ggplot() +
  geom_path(data = rnaturalearthdata::coastline50, aes(x = long, y = lat, group = group)) +
  geom_point(data = clade_relative_abundance[!is.na(clade_relative_abundance$Tricho_I),], aes(x = lon_E, y = lat_N, color = Tricho_I), size = size_obj) + 
  xlim(-180,180) + ylim(-90,90) +
  coord_cartesian(xlim =c(-165,165), ylim = c(-77,77)) +
  labs(title="Tricho I",x ="longitude",y ="latitude")+
  scale_color_distiller(palette = "Reds", direction = 1) +
  guides(color = FALSE)+
  theme1

t3 <- ggplot() +
  geom_path(data = rnaturalearthdata::coastline50, aes(x = long, y = lat, group = group)) +
  geom_point(data = clade_relative_abundance[!is.na(clade_relative_abundance$Tricho_I),], aes(x = lon_E, y = lat_N, color = Tricho_III), size = size_obj) + 
  xlim(-180,180) + ylim(-90,90) +
  coord_cartesian(xlim =c(-165,165), ylim = c(-77,77)) +
  labs(title="Tricho III",x ="longitude",y ="latitude")+
  scale_color_distiller(palette = "Reds", direction = 1) +
  guides(color = FALSE)+
  theme1

tu4 <- ggplot() +
  geom_path(data = rnaturalearthdata::coastline50, aes(x = long, y = lat, group = group)) +
  geom_point(data = clade_relative_abundance[!is.na(clade_relative_abundance$Tricho_I),], aes(x = lon_E, y = lat_N, color = Tricho_UIV), size = size_obj) + 
  xlim(-180,180) + ylim(-90,90) +
  coord_cartesian(xlim =c(-165,165), ylim = c(-77,77)) +
  labs(title="Tricho UIV",x ="longitude",y ="latitude")+
  scale_color_distiller(palette = "Reds", direction = 1) +
  guides(color = FALSE)+
  theme1

tu2 <- ggplot() +
  geom_path(data = rnaturalearthdata::coastline50, aes(x = long, y = lat, group = group)) +
  geom_point(data = clade_relative_abundance[!is.na(clade_relative_abundance$Tricho_I),], aes(x = lon_E, y = lat_N, color = Tricho_UII), size = size_obj) + 
  xlim(-180,180) + ylim(-90,90) +
  coord_cartesian(xlim =c(-165,165), ylim = c(-77,77)) +
  labs(title="Tricho UII",x ="longitude",y ="latitude")+
  scale_color_distiller(palette = "Reds", direction = 1) +
  guides(color = FALSE)+
  theme1

u1 <- ggplot() +
  geom_path(data = rnaturalearthdata::coastline50, aes(x = long, y = lat, group = group)) +
  geom_point(data = clade_relative_abundance[!is.na(clade_relative_abundance$`UCYN-A_1`),], aes(x = lon_E, y = lat_N, color = `UCYN-A_1`), size = size_obj) + 
  xlim(-180,180) + ylim(-90,90) +
  coord_cartesian(xlim =c(-165,165), ylim = c(-77,77)) +
  labs(title="UCYN-A1",x ="longitude",y ="latitude")+
  scale_color_distiller(palette = "Blues", direction = 1) +
  guides(color = FALSE)+
  theme1

u2 <- ggplot() +
  geom_path(data = rnaturalearthdata::coastline50, aes(x = long, y = lat, group = group)) +
  geom_point(data = clade_relative_abundance[!is.na(clade_relative_abundance$`UCYN-A_1`),], aes(x = lon_E, y = lat_N, color = `UCYN-A_2`), size = size_obj) + 
  xlim(-180,180) + ylim(-90,90) +
  coord_cartesian(xlim =c(-165,165), ylim = c(-77,77)) +
  labs(title="UCYN-A2",x ="longitude",y ="latitude")+
  scale_color_distiller(palette = "Blues", direction = 1) +
  guides(color = FALSE)+
  theme1

u3 <- ggplot() +
  geom_path(data = rnaturalearthdata::coastline50, aes(x = long, y = lat, group = group)) +
  geom_point(data = clade_relative_abundance[!is.na(clade_relative_abundance$`UCYN-A_1`),], aes(x = lon_E, y = lat_N, color = `UCYN-A_3`), size = size_obj) + 
  xlim(-180,180) + ylim(-90,90) +
  coord_cartesian(xlim =c(-165,165), ylim = c(-77,77)) +
  labs(title="UCYN-A3",x ="longitude",y ="latitude")+
  scale_color_distiller(palette = "Blues", direction = 1) +
  guides(color = FALSE)+
  theme1

u4 <- ggplot() +
  geom_path(data = rnaturalearthdata::coastline50, aes(x = long, y = lat, group = group)) +
  geom_point(data = clade_relative_abundance[!is.na(clade_relative_abundance$`UCYN-A_1`),], aes(x = lon_E, y = lat_N, color = `UCYN-A_4`), size = size_obj) + 
  xlim(-180,180) + ylim(-90,90) +
  coord_cartesian(xlim =c(-165,165), ylim = c(-77,77)) +
  labs(title="UCYN-A4",x ="longitude",y ="latitude")+
  scale_color_distiller(palette = "Blues", direction = 1) +
  guides(color = FALSE)+
  theme1


library(ggpubr)
ggarrange(
  t1,tu2,t3,tu4,u1,u2,u3, labels = c("B","C","D","E","F","G","H","I"),
  ncol = 4,
  nrow = 2
)

#ggsave("/nifH_microdiversity/images/figure_3/figure_3B-I_24.07.02.pdf",width=7.65,height=3)


# Figure 4 PCOA ####
set.seed(1994)
library(vegan)

t_asv_meta_df <- merge(sample_metadata , t(t_asv[1:100,4:dim(t_asv)[2]]) , by.x="Sample_ID",by.y="row.names")
u_asv_meta_df <- merge(sample_metadata , t(u_asv[1:100,4:dim(u_asv)[2]]) , by.x="Sample_ID",by.y="row.names")

complete_meta <- !(is.na(t_asv_meta_df$Nutricline_1uM_Interp) |
                     is.na(t_asv_meta_df$Temp_SST) |
                     is.na(t_asv_meta_df$Cruise) |
                     is.na(t_asv_meta_df$Omega_P) |
                     is.na(t_asv_meta_df$Omega_N) |
                     is.na(t_asv_meta_df$Omega_Fe))

complete_otu_t <- !apply(t_asv_meta_df[, 15:dim(t_asv_meta_df)[2]], 1, function(x) all(x == 0))
complete_otu_u <- !apply(u_asv_meta_df[, 15:dim(u_asv_meta_df)[2]], 1, function(x) all(x == 0))

complete_t <- (complete_meta & complete_otu_t)
complete_u <- (complete_meta & complete_otu_u)

dist_t <- vegdist(t_asv_meta_df[ complete_t , 15:dim(t_asv_meta_df)[2] ], method = "bray")
dist_u <- vegdist(u_asv_meta_df[ complete_u , 15:dim(u_asv_meta_df)[2] ], method = "bray")

PCOA_t <- wcmdscale(dist_t,eig = TRUE, add = "lingoes")
PCOA_u <- wcmdscale(dist_u,eig = TRUE, add = "lingoes")

PCOA_dat_t <- as.data.frame(scores(PCOA_t,choices=1:2))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
PCOA_dat_u <- as.data.frame(scores(PCOA_u,choices=1:2))  #Using the scores function from vegan to extract the site scores and convert to a data.frame

PCOA_dat_t[,3:16] <- t_asv_meta_df[complete_t,1:14]
PCOA_dat_u[,3:16] <- u_asv_meta_df[complete_u,1:14]

PCOA_dat_t <- merge(PCOA_dat_t,clade_relative_abundance[,c(1:4,11:13)],by.x="Sample_ID",by.y="Row.names")
PCOA_dat_u <- merge(PCOA_dat_u,clade_relative_abundance[,c(5:8,11:13)],by.x="Sample_ID",by.y="Row.names")

env_t <- envfit(PCOA_t ~ Temp_SST+Nutricline_1uM_Interp+Omega_P+Omega_Fe+Omega_N+Tricho_I+Tricho_III+Tricho_UIV+Tricho_UII,
                   data = PCOA_dat_t,
                   choices = 1:2,
                   scaling = "symmetric",
                   permutations = 1000)
env_u <- envfit(PCOA_u ~ Temp_SST+Nutricline_1uM_Interp+Omega_P+Omega_Fe+Omega_N+`UCYN-A_1`+`UCYN-A_2`+`UCYN-A_3`,
                data = PCOA_dat_u,
                choices = 1:2,
                scaling = "symmetric",
                permutations = 1000)

env_dft <- as.data.frame(scores(env_t, display = "vectors"))
env_dfu <- as.data.frame(scores(env_u, display = "vectors"))

env_dft$data <- c("SST","Nut","P","Fe","N","TI","TIII","TUIV","TUII")
env_dfu$data <- c("SST","Nut","P","Fe","N","A1","A2","A3")


library(ggplot2)

col_v1 <- c("#fdbf6f","#ff7f00","#b15928","#cab2d6","#6a3d9a","#a6cee3","#1f78b4","#b2df8a","#33a02c")
col_v2 <- c("#ffff99","#fdbf6f","#ff7f00","#b15928","#cab2d6","#6a3d9a","#a6cee3","#1f78b4","#b2df8a","#33a02c")
theme1 <-theme(panel.background = element_blank() #remove background
               ,panel.border=element_rect(fill=NA,size=1) #make box around axis
               ,panel.grid.major = element_blank()
               ,panel.grid.minor = element_blank()
               ,strip.background=element_blank()
               ,axis.text.x=element_text(colour="black") #make all text black
               ,axis.text.y=element_text(colour="black") #make all text black
               ,axis.ticks=element_line(colour="black") #make all text black
               ,plot.title = element_text(hjust = 0.5) #center the title
)

p1 <- ggplot(PCOA_dat_t, aes(x=Dim1, y=Dim2)) + 
  geom_point(aes(col=Cruise)) +
  scale_color_manual(values=col_v1) +
  geom_segment(data = env_dft[6:9,],#add in the lines for envfit
               aes(x = 0, xend = Dim1*0.65, y = 0, yend = Dim2*0.65),
               arrow = arrow(length = unit(0.1, "cm")), colour = "darkgrey") +
  geom_text(data = env_dft[6:9,], aes(x = Dim1*0.75, y = Dim2*0.75, label = data), size = 5, colour="black") +
  geom_segment(data = env_dft[1:5,],#add in the lines for envfit
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_text(data = env_dft[1:5,], aes(x = Dim1*1.2, y = Dim2*1.2, label = data), size = 5, fontface = "bold") +
  guides(col = FALSE)+
  labs(title="Trichodesmium",x ="PCoA 1",y ="PCoA 2")+
  theme1
  
p2 <- ggplot(PCOA_dat_u, aes(x=Dim1, y=Dim2)) + 
  geom_point(aes(col=Cruise)) +
  scale_color_manual(values=col_v2) +
  geom_segment(data = env_dfu[6:9,],#add in the lines for envfit
               aes(x = 0, xend = Dim1*0.6, y = 0, yend = Dim2*0.6),
               arrow = arrow(length = unit(0.1, "cm")), colour = "darkgrey") +
  geom_text(data = env_dfu[6:9,], aes(x = Dim1*0.7, y = Dim2*0.7, label = data), size = 5, colour="black") +
  geom_segment(data = env_dfu[1:5,],#add in the lines for envfit
               aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "black") +
  geom_text(data = env_dfu[1:5,], aes(x = Dim1*1.2, y = Dim2*1.2, label = data), size = 5, fontface = "bold") +
  guides(col = FALSE)+
  labs(title="UCYN-A",x ="PCoA 1",y ="PCoA 2")+
  theme1

library(ggpubr)
ggarrange(
  p1,p2, labels = c("A", "B")
)
#ggsave("/nifH_microdiversity/images/figure_4/figure_4AB_24.07.02.pdf",
#       width=7.5,height=3.75)

# Figure 4 Spearman correlation ####

t_asv_meta_df <- merge(sample_metadata , t(t_asv[1:100,4:dim(t_asv)[2]]) , by.x="Sample_ID",by.y="row.names")
u_asv_meta_df <- merge(sample_metadata , t(u_asv[1:100,4:dim(u_asv)[2]]) , by.x="Sample_ID",by.y="row.names")


complete_meta <- !(is.na(t_asv_meta_df$Nutricline_1uM_Interp) |
                     is.na(t_asv_meta_df$Temp_SST) |
                     is.na(t_asv_meta_df$Cruise) |
                     is.na(t_asv_meta_df$Omega_P) |
                     is.na(t_asv_meta_df$Omega_N) |
                     is.na(t_asv_meta_df$Omega_Fe))

complete_otu_t <- !apply(t_asv_meta_df[, 15:dim(t_asv_meta_df)[2]], 1, function(x) all(x == 0))
complete_otu_u <- !apply(u_asv_meta_df[, 15:dim(u_asv_meta_df)[2]], 1, function(x) all(x == 0))

complete_otu <- (complete_otu_t | complete_otu_u)

complete_t <- (complete_meta & complete_otu_t)
complete_u <- (complete_meta & complete_otu_u)
complete <- (complete_meta & complete_otu)

t_rel_meta <- merge(t_asv_meta_df[complete_t,1:14],clade_relative_abundance,by.x="Sample_ID",by.y="Row.names")
u_rel_meta <- merge(u_asv_meta_df[complete_u,1:14],clade_relative_abundance,by.x="Sample_ID",by.y="Row.names")
rel_meta <- merge(t_asv_meta_df[complete,1:14],clade_relative_abundance,by.x="Sample_ID",by.y="Row.names")

spearman_corr <- data.frame(corr=character(),
                            meta=character(),
                            rho=numeric(),
                            pval=numeric())
# tricho vs UCYN-A relative abundance spearman correlations
i=0
spearman_corr[(1+i):(5+i),1] <- "Tricho vs. UCYN-A"
spearman_corr[(1+i):(5+i),2] <- c("SST","Nutricline","Omega P","Omega Fe","Omega N")
spearman_corr[1+i,3:4] <- cor.test(rel_meta$Tricho_Total,rel_meta$Temp_SST,method="spearman",exact=F)[4:3]
spearman_corr[2+i,3:4] <- cor.test(rel_meta$Tricho_Total,rel_meta$Nutricline_1uM_Interp ,method="spearman",exact=F)[4:3]
spearman_corr[3+i,3:4] <- cor.test(rel_meta$Tricho_Total,rel_meta$Omega_P,method="spearman",exact=F)[4:3]
spearman_corr[4+i,3:4] <- cor.test(rel_meta$Tricho_Total,rel_meta$Omega_Fe,method="spearman",exact=F)[4:3]
spearman_corr[5+i,3:4] <- cor.test(rel_meta$Tricho_Total,rel_meta$Omega_N,method="spearman",exact=F)[4:3]

# tricho clades spearman correlations
i=5
spearman_corr[(1+i):(5+i),1] <- "Trichodesmium I"
spearman_corr[(1+i):(5+i),2] <- c("SST","Nutricline","Omega P","Omega Fe","Omega N")
spearman_corr[1+i,3:4] <- cor.test(t_rel_meta$Tricho_I,t_rel_meta$Temp_SST,method="spearman",exact=F)[4:3]
spearman_corr[2+i,3:4] <- cor.test(t_rel_meta$Tricho_I,t_rel_meta$Nutricline_1uM_Interp ,method="spearman",exact=F)[4:3]
spearman_corr[3+i,3:4] <- cor.test(t_rel_meta$Tricho_I,t_rel_meta$Omega_P,method="spearman",exact=F)[4:3]
spearman_corr[4+i,3:4] <- cor.test(t_rel_meta$Tricho_I,t_rel_meta$Omega_Fe,method="spearman",exact=F)[4:3]
spearman_corr[5+i,3:4] <- cor.test(t_rel_meta$Tricho_I,t_rel_meta$Omega_N,method="spearman",exact=F)[4:3]

i=10
spearman_corr[(1+i):(5+i),1] <- "Trichodesmium II"
spearman_corr[(1+i):(5+i),2] <- c("SST","Nutricline","Omega P","Omega Fe","Omega N")
spearman_corr[1+i,3:4] <- cor.test(t_rel_meta$Tricho_UII,t_rel_meta$Temp_SST,method="spearman",exact=F)[4:3]
spearman_corr[2+i,3:4] <- cor.test(t_rel_meta$Tricho_UII,t_rel_meta$Nutricline_1uM_Interp ,method="spearman",exact=F)[4:3]
spearman_corr[3+i,3:4] <- cor.test(t_rel_meta$Tricho_UII,t_rel_meta$Omega_P,method="spearman",exact=F)[4:3]
spearman_corr[4+i,3:4] <- cor.test(t_rel_meta$Tricho_UII,t_rel_meta$Omega_Fe,method="spearman",exact=F)[4:3]
spearman_corr[5+i,3:4] <- cor.test(t_rel_meta$Tricho_UII,t_rel_meta$Omega_N,method="spearman",exact=F)[4:3]

i=15
spearman_corr[(1+i):(5+i),1] <- "Trichodesmium III"
spearman_corr[(1+i):(5+i),2] <- c("SST","Nutricline","Omega P","Omega Fe","Omega N")
spearman_corr[1+i,3:4] <- cor.test(t_rel_meta$Tricho_III,t_rel_meta$Temp_SST,method="spearman",exact=F)[4:3]
spearman_corr[2+i,3:4] <- cor.test(t_rel_meta$Tricho_III,t_rel_meta$Nutricline_1uM_Interp ,method="spearman",exact=F)[4:3]
spearman_corr[3+i,3:4] <- cor.test(t_rel_meta$Tricho_III,t_rel_meta$Omega_P,method="spearman",exact=F)[4:3]
spearman_corr[4+i,3:4] <- cor.test(t_rel_meta$Tricho_III,t_rel_meta$Omega_Fe,method="spearman",exact=F)[4:3]
spearman_corr[5+i,3:4] <- cor.test(t_rel_meta$Tricho_III,t_rel_meta$Omega_N,method="spearman",exact=F)[4:3]

i=20
spearman_corr[(1+i):(5+i),1] <- "Trichodesmium IV"
spearman_corr[(1+i):(5+i),2] <- c("SST","Nutricline","Omega P","Omega Fe","Omega N")
spearman_corr[1+i,3:4] <- cor.test(t_rel_meta$Tricho_UIV,t_rel_meta$Temp_SST,method="spearman",exact=F)[4:3]
spearman_corr[2+i,3:4] <- cor.test(t_rel_meta$Tricho_UIV,t_rel_meta$Nutricline_1uM_Interp ,method="spearman",exact=F)[4:3]
spearman_corr[3+i,3:4] <- cor.test(t_rel_meta$Tricho_UIV,t_rel_meta$Omega_P,method="spearman",exact=F)[4:3]
spearman_corr[4+i,3:4] <- cor.test(t_rel_meta$Tricho_UIV,t_rel_meta$Omega_Fe,method="spearman",exact=F)[4:3]
spearman_corr[5+i,3:4] <- cor.test(t_rel_meta$Tricho_UIV,t_rel_meta$Omega_N,method="spearman",exact=F)[4:3]

# UCYN-A clades spearman correlations
i=25
spearman_corr[(1+i):(5+i),1] <- "UCYN-A1"
spearman_corr[(1+i):(5+i),2] <- c("SST","Nutricline","Omega P","Omega Fe","Omega N")
spearman_corr[1+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_1`,u_rel_meta$Temp_SST,method="spearman",exact=F)[4:3]
spearman_corr[2+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_1`,u_rel_meta$Nutricline_1uM_Interp ,method="spearman",exact=F)[4:3]
spearman_corr[3+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_1`,u_rel_meta$Omega_P,method="spearman",exact=F)[4:3]
spearman_corr[4+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_1`,u_rel_meta$Omega_Fe,method="spearman",exact=F)[4:3]
spearman_corr[5+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_1`,u_rel_meta$Omega_N,method="spearman",exact=F)[4:3]

i=30
spearman_corr[(1+i):(5+i),1] <- "UCYN-A2"
spearman_corr[(1+i):(5+i),2] <- c("SST","Nutricline","Omega P","Omega Fe","Omega N")
spearman_corr[1+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_2`,u_rel_meta$Temp_SST,method="spearman",exact=F)[4:3]
spearman_corr[2+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_2`,u_rel_meta$Nutricline_1uM_Interp ,method="spearman",exact=F)[4:3]
spearman_corr[3+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_2`,u_rel_meta$Omega_P,method="spearman",exact=F)[4:3]
spearman_corr[4+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_2`,u_rel_meta$Omega_Fe,method="spearman",exact=F)[4:3]
spearman_corr[5+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_2`,u_rel_meta$Omega_N,method="spearman",exact=F)[4:3]

i=35
spearman_corr[(1+i):(5+i),1] <- "UCYN-A3"
spearman_corr[(1+i):(5+i),2] <- c("SST","Nutricline","Omega P","Omega Fe","Omega N")
spearman_corr[1+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_3`,u_rel_meta$Temp_SST,method="spearman",exact=F)[4:3]
spearman_corr[2+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_3`,u_rel_meta$Nutricline_1uM_Interp ,method="spearman",exact=F)[4:3]
spearman_corr[3+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_3`,u_rel_meta$Omega_P,method="spearman",exact=F)[4:3]
spearman_corr[4+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_3`,u_rel_meta$Omega_Fe,method="spearman",exact=F)[4:3]
spearman_corr[5+i,3:4] <- cor.test(u_rel_meta$`UCYN-A_3`,u_rel_meta$Omega_N,method="spearman",exact=F)[4:3]

spearman_corr$meta <- factor(spearman_corr$meta, levels=c("SST","Nutricline","Omega P","Omega Fe","Omega N"))

spearman_corr$sig <- ""
spearman_corr$sig[spearman_corr$pval < 0.05] <- "*"
spearman_corr$sig[spearman_corr$pval < 0.01] <- "**"
spearman_corr$sig[spearman_corr$pval < 0.001] <- "***"


# make figures
library(ggplot2)

theme2 <-theme(panel.background = element_blank() #remove background
               ,panel.border=element_rect(fill=NA,size=1) #make box around axis
               ,panel.grid.major = element_line(colour="gray90") #grey grid
               ,panel.grid.minor = element_line(colour="gray90") #grey grid
               ,axis.text.x=element_text(colour="black") #make all text black
               ,axis.text.y=element_text(colour="black") #make all text black
               ,axis.ticks=element_line(colour="black") #make all text black
               ,plot.title = element_text(hjust = 0.5) #center the title
)
p0 <- ggplot(data=spearman_corr[1:5,], aes(x=corr,y=rho,fill=meta))+
  geom_bar(stat="identity", position=position_dodge(),color="black")+
  labs(title="Trichodesmium vs. UCYN-A",y ="Spearman rho")+
  guides(fill = FALSE)+
  ylim(-0.7,0.63)+
  theme2+
  scale_fill_manual(values=c("#785EF0","#D55E00","#648FFF","#DC267F","#FFB000"))

p1 <- ggplot(data=spearman_corr[6:25,], aes(x=corr,y=rho,fill=meta))+
  geom_bar(stat="identity", position=position_dodge(),width=0.7,color="black")+
  labs(title="Trichodesmium",x ="clade",y ="Spearman rho")+
  scale_x_discrete(labels=c("I","UII","III","UIV"))+
  guides(fill = FALSE)+
  ylim(-0.7,0.63)+
  theme2+
  scale_fill_manual(values=c("#785EF0","#D55E00","#648FFF","#DC267F","#FFB000"))

p2 <- ggplot(data=spearman_corr[26:40,], aes(x=corr,y=rho,fill=meta))+
  geom_bar(stat="identity", position=position_dodge(),width=0.7,color="black")+
  labs(title="UCYN-A",x ="clade",y ="Spearman rho")+
  scale_x_discrete(labels=c("1","2","3"))+
  guides(fill = FALSE)+
  ylim(-0.7,0.63)+
  theme2+
  scale_fill_manual(values=c("#785EF0","#D55E00","#648FFF","#DC267F","#FFB000"))

library(ggpubr)
ggarrange(
  p1,p2, labels = c("D", "E"),
  ncol = 2,
  nrow = 1
)
#ggsave("/nifH_microdiversity/images/figure_4/figure_4DE_24.07.02.pdf",
#       width=6,height=3)

ggarrange(
  p0,labels = c("C")
)
 #ggsave("/nifH_microdiversity/images/figure_4/figure_4C_24.01.18.pdf",
#       width=1.5,height=3)
