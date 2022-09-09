rm(list = ls())

trout <- read.csv("~/work/Trout/BAMs/Atha-Buffalo_Prairie.s.OhnologSummary.csv")
jared <- read.csv("~/work/Trout/BAMs/Atha-Buffalo_Prairie.s.JaredClean.OhnologSummary.csv")

tetrasomic_chroms <- c("OmyA_01","OmyA_02","OmyA_03","OmyA_06","OmyA_07","OmyA_10","OmyA_12","OmyA_13","OmyA_15","OmyA_17","OmyA_18","OmyA_19","OmyA_21","OmyA_23","OmyA_26")
jared$tetrasome <- jared$CHROM%in%tetrasomic_chroms
trout$tetrasome <- trout$CHROM%in%tetrasomic_chroms

library(ggplot2)

plot( jared$total_mappings ~ trout$total_mappings )
jared$data_removed <- 1- jared$total_mappings/trout$total_mappings

##############################

ggplot(data = jared,
       aes(x = START/1e6,
           y = (total_mappings*70)/(END-START)))+
  geom_point(aes(fill  = tetrasome), 
             alpha = 0.5, pch = 21, size= 1.2)+
  ylab("Coverage")+
  xlab("Position in genome (Mbp)")+
  scale_y_continuous(limits = c(0,30))+
  geom_smooth(aes(col = "Filtered"), span = 0.01, method = "loess", se = FALSE)+
  geom_smooth(data = trout, aes(col = "Unfiltered"),  span = 0.1, method = "loess")+
  facet_wrap(~CHROM,
             scales = "free_x")+
  scale_fill_manual("Predicted tetrasomic",
                     values = c("black","#D55E00"))+
  scale_color_manual("Source",
                    values = c( "#56B4E9","#009E73"))+
#  scale_y_continuous(limits = c(0,1))+
  theme_bw()

################################

ggplot(data = jared,
       aes(x = START/1e6,
           y =data_removed))+
  geom_point(aes(fill  = tetrasome), 
             alpha = 0.5, pch = 21, size= 1.2)+
  ggtitle("The effect of filtering")+
  ylab("Proportion of data removed by filtering")+
  xlab("Position in genome (Mbp)")+
  scale_y_continuous(limits = c(0,1))+
  geom_smooth(aes(col = "Filtered"), span = 0.1, method = "loess", se = FALSE)+
 # geom_smooth(data = trout, aes(col = "Unfiltered"),  span = 0.1, method = "loess")+
  facet_wrap(~CHROM,
             scales = "free_x")+
  scale_fill_manual("Predicted tetrasomic",
                    values = c("black","#D55E00"))+
  scale_color_manual("Source",
                     values = c( "#56B4E9","#009E73"))+
  #  scale_y_continuous(limits = c(0,1))+
  theme_bw()

##############################

ggplot(data = trout,
       aes(x = START/1e6,
           y = (primary_alignment_wrong)/(total_mappings),
           col  = tetrasome))+
  geom_point(alpha = 0.5)+
  ggtitle("The frequency of incorrect alignments (before filtering)")+
  ylab("Proportion of all reads that are incorrectly mapped")+
  xlab("Position in genome (Mbp)")+
#  geom_hline(aes(yintercept = 12), col = "red")+
  geom_smooth(col = "#56B4E9", span = 0.1, method = "loess")+
  facet_wrap(~CHROM,
             scales = "free_x")+
  scale_color_manual("Predicted tetrasomic",
                     values = c("black","#D55E00"))+
  theme_bw()

ggplot(data = jared,
       aes(x = START/1e6,
           y = (primary_alignment_wrong)/(total_mappings),
           col  = tetrasome))+
  geom_point(alpha = 0.5)+
  ggtitle("The frequency of incorrect alignments (after filtering)")+
  ylab("Proportion of all reads that are incorrectly mapped")+
  xlab("Position in genome (Mbp)")+
  #  geom_hline(aes(yintercept = 12), col = "red")+
  geom_smooth(col = "#56B4E9", span = 0.1, method = "loess")+
  facet_wrap(~CHROM,
             scales = "free_x")+
  scale_color_manual("Predicted tetrasomic",
                     values = c("black","#D55E00"))+
  theme_bw()



##############################


trout$single_ohno_maps_removed <- trout$single_ohno_maps - jared$single_ohno_maps 

trout$ohno_maps <- (trout$single_ohno_maps + trout$double_ohno_maps + trout$multi_ohno_maps) - (jared$single_ohno_maps + jared$double_ohno_maps + jared$multi_ohno_maps) 



ggplot(data = trout,
       aes(x = START/1e6,
           y = ohno_maps/(total_mappings)))+
  geom_point()+
  facet_wrap(~CHROM, scales = "free_x")+
#  geom_smooth(col = "#56B4E9", span = 0.1, method = "loess")+
  theme_bw()


ggplot(data = trout,
       aes(x = START,
           y = single_maps))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  xlab("Number of single Ohnolog mappings removed")+
  ylab("Total number of Ohnolog mappings removed")+
  facet_wrap(~CHROM)+
  theme_bw()


ggplot(data = jared,
       aes(x = single_maps,
           y = single_ohno_maps))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point()+
  facet_wrap(~CHROM)+
  theme_bw()
