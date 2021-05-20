library(ggplot2)
library(rio)
library(lubridate)
library(seasonal)
library(reshape2)
library(tidyr)
library(scales)
library(tseries)
library(stringr)

Sys.setlocale("LC_TIME", "C")
dat <- import("/Users/rutra/ВШЭ/Магистратура/Thesis/Data/Aggregated data/Data Quarterly.xlsx", sheet=3)
dat <- dat[NROW(dat):1,]
dat <- na.omit(dat)

#Deseasonalization
to_deseasonalize <- c("gdp_nominal_qoq","imp_price_qoq", "exp_price_qoq", 
                      grep("cpi.*", colnames(dat), value=T))
dat_unseas <- dat
for(colname in to_deseasonalize){
  current_ts <- ts(dat[,colname], start=c(2005, 2), frequency=4)
  current_seas <- seas(current_ts, transform.function="none",
                       #regression.aictest=NULL, outlier=NULL,
                       #automdl=NULL,
                       seats.noadmiss="yes")
  dat_unseas[,colname] <- as.numeric(final(current_seas))
}

data_draw <- dat_unseas[,c("date", "oil_USD_qoq", "miacr_31", 
                           "neer_qoq", "gdp_real_SA_qoq", "cpi_all_qoq")]
colnames(data_draw) <- c("Date", "Oil price gr. rate", 
                         "MIACR 31-180 days", 
                         "NEER gr. rate", 
                         "Real GDP gr. rate",
                         "CPI gr. rate")
data_draw$Date <- as.Date(data_draw$Date)
data_draw_long <- gather(data_draw, "Variable", "Value", -Date)
months <- month(data_draw$Date)
ggplot(data_draw_long, aes(x = Date, y = Value), size = 0.1) + geom_hline(yintercept = 0, col="grey") +
  geom_line() + 
  facet_wrap(~ Variable, ncol=3, scales="free") + 
  theme_minimal() + #+theme_bw()
  scale_x_date(breaks=pretty_breaks(), date_labels="%y", 
               date_breaks="1 year") + 
ggsave(filename="time_series.eps",
       path="/Users/rutra/ВШЭ/Магистратура/Thesis/Text/figures/",
       device="eps",
       width=320, height=210, dpi=320, units = "mm", limitsize=FALSE)

data_draw_long_perr <- data_draw_long[data_draw_long$Variable %in% 
                                        c("NEER gr. rate", "CPI gr. rate"),]
ggplot(data_draw_long_perr, aes(x = Date, y = Value, linetype=Variable), size = 0.1) + 
  geom_hline(yintercept = 0, col="grey")+ geom_line() + theme_minimal() + #+theme_bw()
  #facet_wrap(~ Variable, nrow=2, scales="free") + 
  scale_x_date(breaks=pretty_breaks(), date_labels="%b %y", 
               date_breaks = "6 months", expand=c(0,0)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top")
ggsave(filename="neer_cpi.eps",
       path="/Users/rutra/ВШЭ/Магистратура/Thesis/Text/figures/",
       device="eps",
       width=200, height=100, dpi=320, units = "mm", limitsize=FALSE)

for(colname in colnames(data_draw)[-1]){
  cat(colname, "\n")
  print(PP.test(data_draw[,colname]))
  print(adf.test(data_draw[,colname]))
}

corstarsl <- function(x){ 
  require(Hmisc) 
  x <- as.matrix(x) 
  R <- rcorr(x)$r 
  p <- rcorr(x)$P 
  
  ## define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
  
  ## trunctuate the matrix that holds the correlations to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1] 
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x)) 
  diag(Rnew) <- paste(diag(R), " ", sep="") 
  rownames(Rnew) <- colnames(x) 
  colnames(Rnew) <- paste(colnames(x), "", sep="") 
  
  ## remove upper triangle
  Rnew <- as.matrix(Rnew)
  Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
  Rnew <- as.data.frame(Rnew) 
  
  ## remove last column and return the matrix (which is now a data frame)
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  return(Rnew) 
}

corstarsl(data_draw[,-1])

### IRF
irf_st_data <- read.csv("/Users/rutra/ВШЭ/Магистратура/Thesis/R workspaces/st_irf.csv", check.names = FALSE)
irf_st_long <- gather(irf_st_data, "Variable", "Value", -Period)
irf_st_long$Variable <- str_replace_all(irf_st_long$Variable, " shock ", ' shock\n')
irf_st_long_1 <- irf_st_long[irf_st_long$Variable %in% c(grep("Global persistent *", irf_st_long$Variable, value=T),
                                                         grep("Global demand *", irf_st_long$Variable, value=T),
                                                         grep("Monetary *", irf_st_long$Variable, value=T)),]
irf_st_long_2 <- irf_st_long[irf_st_long$Variable %in% c(grep("Exchange rate *", irf_st_long$Variable, value=T),
                                                         grep("Supply *", irf_st_long$Variable, value=T),
                                                         grep("Demand *", irf_st_long$Variable, value=T)),]
  
ggplot(irf_st_long_1, aes(x = Period, y = Value), size = 0.1) + geom_line() + 
  facet_wrap(~ Variable, nrow=3, scales="free") + geom_hline(yintercept = 0, col="grey") +
  theme_minimal() + scale_x_continuous(breaks=1:12, expand=c(0,0)) +
  theme(text = element_text(size=10 ), axis.text.x = element_text(size=8), 
        axis.text.y = element_text(size=8 ))
  #+theme_bw()
  #scale_x_date(breaks=pretty_breaks(), date_labels="%y", 
  #             date_breaks="1 year")
  
ggsave(filename="irf_st_1.eps",
       path="/Users/rutra/ВШЭ/Магистратура/Thesis/Text/figures/",
       device="eps",
       width=320, height=210, dpi=320, units = "mm", limitsize=FALSE)

ggplot(irf_st_long_2, aes(x = Period, y = Value), size = 0.1) + geom_line() + 
  facet_wrap(~ Variable, nrow=3, scales="free") + geom_hline(yintercept = 0, col="grey") +
  theme_minimal() + scale_x_continuous(breaks=1:12, expand=c(0,0))+
  theme(text = element_text(size=10 ), axis.text.x = element_text(size=8 ), 
        axis.text.y = element_text(size=8 ))#+theme_bw()
  #scale_x_date(breaks=pretty_breaks(), date_labels="%y", 
  #             date_breaks="1 year")
  
ggsave(filename="irf_st_2.eps",
         path="/Users/rutra/ВШЭ/Магистратура/Thesis/Text/figures/",
         device="eps",
         width=320, height=210, dpi=320, units = "mm", limitsize=FALSE)


#Historical decomposition
hd_st <- read.csv("/Users/rutra/ВШЭ/Магистратура/Thesis/R workspaces/st_hd.csv", check.names = FALSE)
hd_st_joined <- data.frame(
  `Global persistent shock` = hd_st$`Cumulative effect of  oil_USD shock on  cpi_all_mom`,
  `Global demand shock` = hd_st$`Cumulative effect of  imp_price_mom shock on  cpi_all_mom`,
  `Monetary shock` = hd_st$`Cumulative effect of  miacr_31 shock on  cpi_all_mom`,
  `Exchange rate shock` = hd_st$`Cumulative effect of  neer_mom shock on  cpi_all_mom`,
  `Supply shock` = hd_st$`Cumulative effect of  gdp_per_cap_SA shock on  cpi_all_mom`,
  `Demand shock` = hd_st$`Cumulative effect of  cpi_all_mom shock on  cpi_all_mom`,
  check.names = FALSE
  )
hd_st_joined$Date <- as.Date(seq(as.yearmon("2005-7"),as.yearmon("2020-8"), 1/12)) #model's order matters!
hd_st_joined_bar <- gather(hd_st_joined, "Variable", "Value", -Date)
ggplot(hd_st_joined_bar, aes(x = Date, y = Value, fill = Variable)) + 
  geom_bar(stat="identity") + #scale_x_date(breaks=hd_st_joined_bar$Date) + 
  scale_x_date(breaks=pretty_breaks(), date_labels="%b %y", 
      date_breaks = "6 months", expand=c(0,0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 8),
        axis.title.x = element_text(vjust = 0.5, size = 9),
        axis.title.y = element_text(vjust = 0.5, size = 9),
        legend.position="bottom") +
  xlab("Date") + ylab("Contribution to CPI gr. r.") +
  geom_line(aes(x=rep(hd_st_joined$Date, 6), y=rep(hd_st$`Demeaned series  cpi_all_mom`, 6)), col="black", size=0.3) +##383838
  scale_fill_brewer(palette="Dark2") + 

ggsave(filename="hd_st_full.eps",
       path="/Users/rutra/ВШЭ/Магистратура/Thesis/Text/figures/",
       device="eps",
       width=320, height=210, dpi=320, units = "mm", limitsize=FALSE)

hd_st_cut <- hd_st
hd_st_cut$Date <- as.Date(seq(as.yearmon("2005-7"),as.yearmon("2020-8"), 1/12))
hd_st_cut <- hd_st_cut[hd_st_cut$Date >= as.Date("2017-01-01"),]
hd_st_joined_cut <- data.frame(
  `Global persistent shock` = hd_st_cut$`Cumulative effect of  oil_USD shock on  cpi_all_mom`,
  `Global demand shock` = hd_st_cut$`Cumulative effect of  imp_price_mom shock on  cpi_all_mom`,
  `Monetary shock` = hd_st_cut$`Cumulative effect of  miacr_31 shock on  cpi_all_mom`,
  `Exchange rate shock` = hd_st_cut$`Cumulative effect of  neer_mom shock on  cpi_all_mom`,
  `Supply shock` = hd_st_cut$`Cumulative effect of  gdp_per_cap_SA shock on  cpi_all_mom`,
  `Demand shock` = hd_st_cut$`Cumulative effect of  cpi_all_mom shock on  cpi_all_mom`,
  check.names = FALSE
)
hd_st_joined_cut <- hd_st_joined[hd_st_joined$Date >= as.Date("2017-01-01"), ]
hd_st_joined_cut_bar <- gather(hd_st_joined_cut, "Variable", "Value", -Date)
ggplot(hd_st_joined_cut_bar, aes(x = Date, y = Value, fill = Variable)) + 
  geom_bar(stat="identity") +
  scale_x_date(breaks=pretty_breaks(), date_labels="%b %y", 
               date_breaks = "2 months", expand=c(0,0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 8),
        axis.title.x = element_text(vjust = 0.5, size = 9),
        axis.title.y = element_text(vjust = 0.5, size = 9),
        legend.position="bottom") +
  xlab("Date") + ylab("Contribution to CPI gr. r.") +
  geom_line(aes(x=rep(hd_st_joined_cut$Date, 6), 
                y=rep(hd_st_cut$`Demeaned series  cpi_all_mom`, 6)), col="black", size=0.3) +##383838
  scale_fill_brewer(palette="Dark2") 

ggsave(filename="hd_st_cut.eps",
       path="/Users/rutra/ВШЭ/Магистратура/Thesis/Text/figures/",
       device="eps",
       width=320, height=210, dpi=320, units = "mm", limitsize=FALSE)

#FEVD
