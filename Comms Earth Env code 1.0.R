rm(list = ls())


# Jez note - check to see if all packages are relevant

#BiocManager::install("variancePartition")
library(variancePartition) # missing
library(openxlsx)
library(tibble)
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(emmeans)
library(sjPlot)
library(MuMIn)
library(car)
library(tidyxl)
library(DHARMa)
library(visdat) # missing
library(partR2) # missing


### Davis et al (2021) Communications Earth & Environment Volume 2, Article number: 105 (https://www.nature.com/articles/s43247-021-00168)
    # Data availability
    # The authors declare that the data supporting the findings of this study and its source data are available within the paper and 
    # its Supplementary Information Files. The metadata is also available on the SEANOE database at https://doi.org/10.17882/80022.
    # Code availability
    # The authors declare that the R code supporting the findings of this study are available within the paper.

### Wolfe & Roff (2021) "Global predictions of coral reef dissolution in the Anthropocene" data analysis 

# 1) load data files
#  a) data file from Communications Earth & Environment (Supplementary Data 1) ------ https://static-content.springer.com/esm/art%3A10.1038%2Fs43247-021-00168-w/MediaObjects/43247_2021_168_MOESM3_ESM.xlsx
#  b) data file from SEANOE (DOI: 10.17882/80022) ----------------------------------- https://www.seanoe.org/data/00688/80022/data/83009.xlsx

data <- read.xlsx("https://www.seanoe.org/data/00688/80022/data/83009.xlsx") # replace throughout
data <- read.xlsx("~/Dropbox/Data - Carbonate R Package/caRb/43247_2021_168_MOESM3_ESM.xlsx") 
formats <- xlsx_formats( "~/Dropbox/Data - Carbonate R Package/caRb/43247_2021_168_MOESM3_ESM.xlsx" )
cells <- xlsx_cells( "~/Dropbox/Data - Carbonate R Package/caRb/43247_2021_168_MOESM3_ESM.xlsx" )

# 2) extract colours for key from the first row of the xlsx file
figure3cols <- cells[cells$local_format_id %in% which(formats$local$fill$patternFill$fgColor$rgb == "FFFFC000"), "row"] %>% distinct %>%as.data.frame()
figure4cols <- cells[cells$local_format_id %in% which(formats$local$fill$patternFill$fgColor$rgb == "FF4472C4"), "row"] %>% distinct %>% as.data.frame()
figure3and4cols <- cells[cells$local_format_id %in% which( formats$local$fill$patternFill$fgColor$rgb == "FF70AD47"), "row" ] %>% distinct %>% as.data.frame()
figure3cols$key <- 'fig3'
figure4cols$key <- 'fig4'
figure3and4cols$key <- 'fig3and4'
rowcolours <- rbind(figure3cols, figure4cols, figure3and4cols) %>% arrange(row) %>% filter(row>1)

# 3) clean up data and column names from xlsx format
colnames(data) <- data %>% dplyr::slice(1)
data$row <- seq(2,(nrow(data)+1))
data <- data %>% select("Short Reference", "Location", "Year Studied", "row", "Heat Type","Lat", "Wave Action", "Duration (days)", 
                         "Net dissolution at Night? (study average)" , "Night measured?","avg Ωarag", "Temp avg",	"Temp Stdev", 
                        "Depth",	"% coral cover ", "% Coralline",	"% calcifiers", "Health", 
                        "Gnet (mmol m-2 d-1)", "Gnet error",  "Pnet (mmol m-2 d-1)", "Pnet error",  "Method") %>% 
                 dplyr::rename("Short_ref" = "Short Reference", "Coral" = "% coral cover ", "Coralline" =  "% Coralline",  "Calcifiers"= 
                        "% calcifiers", "Gnet" =  "Gnet (mmol m-2 d-1)", "Pnet" = "Pnet (mmol m-2 d-1)",  "Pnet error",
                        "Night_d" = "Net dissolution at Night? (study average)", "Night" = "Night measured?", "Aragonite" = "avg Ωarag", "Duration" = "Duration (days)", 
                        "Season" = "Heat Type", "Year" = "Year Studied", "Wave" = "Wave Action",
                        "Temp" = "Temp avg", "TempStd" = "Temp Stdev",  "Gnet.error" = "Gnet error", "Pnet.error" = "Pnet error",) %>%
                 slice(3:118) %>% 
                #filter(!Season == "A") %>%
                 replace_na(list(Health = "H")) %>% 
                 mutate(Health = dplyr::recode(Health,'Dead' = "D")) %>%
                 mutate(Night_d = dplyr::recode(Night_d,'y' = "Y")) %>%
                 replace_na(list(Night_d = "Y")) %>%
                 mutate(Night = dplyr::recode(Night,'yes' = "Yes")) %>%
                 type_convert() %>% mutate_if(is.character,as.factor) %>%
                 mutate(Lat = abs(Lat)) %>%
                 left_join(., rowcolours, by="row")

# 4) Amend original data file with missing data points from original studies:
data$Depth[data$Short_ref == 'Barnes, 1988. Proc 6th Intl Coral Reef Symp' & is.na(data$Depth)] <- '1.11'
data$Depth[data$Short_ref == 'Kayanne et al., 2005. Glob Biogeochem Cycles' & data$Year == "1994" & is.na(data$Depth)] <- '2.0'
data$Depth[data$Short_ref == 'Kayanne et al., 2005. Glob Biogeochem Cycles' & data$Year == "2000" & is.na(data$Depth)] <- '1.55'
data$Depth[data$Short_ref == 'Pisapia et al., 2019. Front Mar Sci' & is.na(data$Depth)] <- '2.3'
data$Depth[data$Short_ref == 'Silverman et al., 2014. GCA' & is.na(data$Depth)] <- '1.5'

# Duration missing in data but 5-10 days in Figure 4 Davis et al 2021 
data$Duration[data$Short_ref == 'Kinsey & Davies, 1979. Stud. Env. Sci'] <- 5


# 5) subset the data from from Davis et al (2021) to include only observations which occurred using the same site and
#     seasonal bin over different years. Recreate Figure 4 (Long-term changes in coral reef ecosystem calcification for 
#     well-studied reefs) and linear mixed effects model to estimate decline in Gnet ("If future change continues at the 
#     current rate of decline, we can expect average global net-zero calcification around 2054")

fig4subset <- data %>% filter(key=="fig3and4" | key=="fig4") %>% 
                       mutate(Depth = as.numeric(Depth)) %>%
                       mutate(across(Location, factor, levels=c("Kanehoe Bay Hawaii","Shiraho Reef, Japan", "Palau", "Davies Reef GBR", "Lizard Island, GBR", "One Tree Island GBR")))

# 6) Add Corrections from COMMSENV-21-0437 Davis et al Response to Matters Arising (replace with new DOI at later date

      # For the One Tree Island values presented in Matters Arising, Table 1 & Figure 1f: All One Tree Island
      # datapoints were collected in the Spring (September-November). In the Davis et al metadata, the
      # incorrect Kinsey 1977 datapoint is highlighted as green, the actual value used in our analyses was the
      # spring (September 1975) datapoint. The Kwiatkowsi datapoint in the metadata table states ‘Spring
      # (Sept-Oct 2014)’, though the ‘Season’ column mistakenly states ‘Autumn’. 

      # change Kinsey 1977 datapoint to the spring (September 1975) datapoint:
      fig4subset$Season[fig4subset$Short_ref == 'Kinsey, 1977. Proc 3rd Intl Coral Reef Symp'] <- "C" # season
      fig4subset$Gnet[fig4subset$Short_ref == 'Kinsey, 1977. Proc 3rd Intl Coral Reef Symp'] <- 125 # gnet for September 1975
      fig4subset$Gnet.error[fig4subset$Short_ref == 'Kinsey, 1977. Proc 3rd Intl Coral Reef Symp'] <- 12.5 # gnet for September 1975
      
      # change Kwiatkowsi datapoint from 'autumn' to 'spring' ("Hot" to "Cold")
      fig4subset$Season[fig4subset$Short_ref == 'Kwiatkowski et al., 2016. GRL'] <- "C"

# 7) reanalyse Figure 4 ("Long-term changes in coral reef ecosystem calcification (Gnet) for well-studied reefs") from Davis et al (2021)

fm1.subset <- lmer(Gnet ~ Year + (1|Location), data=fig4subset)
Anova(fm1.subset)
r.squaredGLMM(fm1.subset)
plot(simulationOutput <- simulateResiduals(fittedModel = fm1.subset))

newdata <- with(fig4subset, list(Year = seq(1975, 2200)))
preds  <-  emmeans(fm1.subset, ~ Year, data=fig4subset, at=newdata, type = 'response') %>% as.data.frame

ggplot() + theme_bw() + # linear predictions with 95% CI
  geom_point(data=fig4subset, aes(x=Year, y=Gnet, size=Duration)) + scale_size_binned(breaks=c(1,5,10,20,30)) +
  geom_errorbar(data=fig4subset, aes(x=Year, ymin=Gnet-Gnet.error, ymax=Gnet+Gnet.error)) + 
  geom_line(data=preds, aes(x=Year, y=emmean)) +
  geom_ribbon(data=preds, aes(x=Year, ymin=lower.CL, max=upper.CL), alpha=0.2) +
  geom_vline(xintercept = 2054) + coord_cartesian(ylim = c(0, 400), xlim = c(1970, 2170))
  
# 8) Plot global changes in coral reef ecosystem calcification (Gnet) over time from subset data (Figure 1 in Wolfe & Roff (2021))

ggplot() + theme_bw() + facet_wrap(~Location) + xlim(1970,2020) +
  theme(legend.position="right") +
  geom_smooth(data=fig4subset %>% filter(Season=="A"), aes(x=Year, y=Gnet), color="green", method = lm, se=FALSE, size=0.5) +
  geom_smooth(data=fig4subset %>% filter(Season=="C"), aes(x=Year, y=Gnet), color="lightblue", method = lm, se=FALSE, size=0.5) +
  geom_smooth(data=fig4subset %>% filter(Season=="H") %>% filter(Wave=="P"), aes(x=Year, y=Gnet), color="red", method = lm, se=FALSE, size=0.5) +
  geom_smooth(data=fig4subset %>% filter(Season=="H") %>% filter(Wave=="M"), aes(x=Year, y=Gnet), color="red", method = lm, se=FALSE, size=0.5) +
  geom_smooth(data=fig4subset %>% filter(Season=="H") %>% filter(Wave=="E"), aes(x=Year, y=Gnet), color="red", method = lm, se=FALSE, size=0.5) +
  geom_errorbar(data=fig4subset, aes(x=Year, ymin=Gnet-Gnet.error, ymax=Gnet+Gnet.error)) +
  geom_point(data=fig4subset %>% filter(Season=="A"), aes(x=Year, y=Gnet), fill="green", shape=21, size=3) +
  geom_point(data=fig4subset %>% filter(Season=="C") %>% filter(Wave=="P"), aes(x=Year, y=Gnet), fill="lightblue", shape=21, size=3, show.legend=TRUE) + 
  geom_point(data=fig4subset %>% filter(Season=="C") %>% filter(Wave=="M"), aes(x=Year, y=Gnet), fill="lightblue", shape=24, size=3, show.legend=TRUE) + 
  geom_point(data=fig4subset %>% filter(Season=="C") %>% filter(Wave=="E"), aes(x=Year, y=Gnet), fill="lightblue", shape=22, size=3, show.legend=TRUE) + 
  geom_point(data=fig4subset %>% filter(Season=="H") %>% filter(Wave=="P"), aes(x=Year, y=Gnet), fill="red", shape=21, size=3, show.legend=TRUE) +
  geom_point(data=fig4subset %>% filter(Season=="H") %>% filter(Wave=="M"), aes(x=Year, y=Gnet), fill="red", shape=24, size=3, show.legend=TRUE) +
  geom_point(data=fig4subset %>% filter(Season=="H") %>% filter(Wave=="E"), aes(x=Year, y=Gnet), fill="red", shape=22, size=3, show.legend=TRUE)
  


### Kenny - see below

fm1.lizard.full <- lm(Gnet ~ Year, data=fig4subset %>% filter(Location=="Lizard Island, GBR"))
Anova(fm1.lizard.full)
r.squaredGLMM(fm1.lizard.full)

fm1.lizard <- lm(Gnet ~ Year, data=fig4subset %>% filter(Location=="Lizard Island, GBR") %>% filter(Year > 2000))
Anova(fm1.lizard)
r.squaredGLMM(fm1.lizard)

fm1.onetree <- lm(Gnet ~ Year, data=fig4subset %>% filter(Location=="One Tree Island GBR"))
Anova(fm1.lizard)
r.squaredGLMM(fm1.onetree)
