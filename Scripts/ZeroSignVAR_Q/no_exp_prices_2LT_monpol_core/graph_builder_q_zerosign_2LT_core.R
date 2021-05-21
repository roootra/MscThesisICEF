library(ggplot2)
library(rio)
library(lubridate)
library(seasonal)
library(reshape2)
library(tidyr)
library(scales)
library(tseries)
library(stringr)
library(zoo)

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

#STOP HERE
data_draw <- dat_unseas[,c("date", "oil_USD_qoq", "miacr_31", 
                           "neer_qoq", "gdp_real_SA_qoq", "cpi_all_qoq", "cpi_core_qoq")]
colnames(data_draw) <- c("Date", "Oil price gr. rate", 
                         "MIACR 31-180 days", 
                         "NEER gr. rate", 
                         "Real GDP gr. rate",
                         "CPI gr. rate")
data_draw$`NEER gr. rate` <- data_draw$`NEER gr. rate`*-1
data_draw$Date <- as.Date(data_draw$Date)
data_draw_long <- gather(data_draw, "Variable", "Value", -Date)
months <- month(data_draw$Date)
ggplot(data_draw_long, aes(x = Date, y = Value), size = 0.1) + geom_hline(yintercept = 0, col="grey") +
  geom_line() + xlab("Year") +
  facet_wrap(~ Variable, ncol=3, scales="free") + 
  theme_minimal() + #+theme_bw()
  scale_x_date(breaks=pretty_breaks(), date_labels="%y", 
               date_breaks="1 year") 
ggsave(filename="time_series.eps",
       path="/Users/rutra/ВШЭ/Магистратура/Thesis/Text/figures/",
       device="eps",
       width=320, height=210, dpi=320, units = "mm", limitsize=FALSE)

data_draw_long_perr <- data_draw_long[data_draw_long$Variable %in% 
                                        c("NEER gr. rate", "CPI gr. rate"),]
ggplot(data_draw_long_perr, aes(x = Date, y = Value, linetype=Variable), size = 0.1) + 
  geom_hline(yintercept = 0, col="grey") + geom_line() + theme_minimal() + #+theme_bw()
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
#START HERE
### IRF
irf_mat <- import("/Users/rutra/ВШЭ/Магистратура/Thesis/Scripts/ZeroSignVAR_Q/no_exp_prices_2LT_monpol_core/Model_PERR_Q_rc/tables/irfs.mat")
varnames <- c("real GDP", "core CPI", "int. rate", "NEER", "oil price")
shocknames <- c("Supply", "Demand", "Monetary", "Exchange rate", "Global persistent")
colnames_irf <- c()
for(varname in varnames){
  new_names <- paste0(shocknames, " shock\n to ", varname)
  colnames_irf <- c(colnames_irf, new_names)
}
irf_mat_ctm <- as.data.frame(irf_mat$irfCTM)
colnames(irf_mat_ctm) <- colnames_irf
irf_mat_ctm$Period <- 0:(NROW(irf_mat_ctm)-1)
irf_mat_ctm <- irf_mat_ctm[1:9,]
irf_mat_long <- gather(irf_mat_ctm, "Variable", "Value", -Period)
#irf_st_long$Variable <- str_replace_all(irf_st_long$Variable, " shock ", ' shock\n')
irf_mat_long_1 <- irf_mat_long[irf_mat_long$Variable %in% c(grep("Global persistent *", irf_mat_long$Variable, value=T),
                                                         grep("Monetary *", irf_mat_long$Variable, value=T),
                                                         grep("Exchange rate *", irf_mat_long$Variable, value=T)),]
irf_mat_long_2 <- irf_mat_long[irf_mat_long$Variable %in% c(grep("Supply *", irf_mat_long$Variable, value=T),
                                                         grep("Demand *", irf_mat_long$Variable, value=T)),]
  
ggplot(irf_mat_long_1, aes(x = Period, y = Value), size = 0.1) + geom_line() + 
  facet_wrap(~ Variable, nrow=3, scales="free") + geom_hline(yintercept = 0, col="grey") +
  theme_minimal() + scale_x_continuous(breaks=1:12, expand=c(0,0)) +
  theme(text = element_text(size=10 ), axis.text.x = element_text(size=8), 
        axis.text.y = element_text(size=8 ))
  #+theme_bw()
  #scale_x_date(breaks=pretty_breaks(), date_labels="%y", 
  #             date_breaks="1 year")
  
ggsave(filename="irf_1_LT_core.eps",
       path="/Users/rutra/ВШЭ/Магистратура/Thesis/Text/figures/",
       device="eps",
       width=320, height=210, dpi=320, units = "mm", limitsize=FALSE)

ggplot(irf_mat_long_2, aes(x = Period, y = Value), size = 0.1) + geom_line() + 
  facet_wrap(~ Variable, nrow=2, scales="free") + geom_hline(yintercept = 0, col="grey") +
  theme_minimal() + scale_x_continuous(breaks=1:12, expand=c(0,0))+
  theme(text = element_text(size=10 ), axis.text.x = element_text(size=8 ), 
        axis.text.y = element_text(size=8 ))#+theme_bw()
  #scale_x_date(breaks=pretty_breaks(), date_labels="%y", 
  #             date_breaks="1 year")
  
ggsave(filename="irf_2_LT_core.eps",
         path="/Users/rutra/ВШЭ/Магистратура/Thesis/Text/figures/",
         device="eps",
         width=320, height=210, dpi=320, units = "mm", limitsize=FALSE)


#Historical decomposition
hd_mat <- import("/Users/rutra/ВШЭ/Магистратура/Thesis/Scripts/ZeroSignVAR_Q/no_exp_prices_2LT_monpol_core/Model_PERR_Q_rc/tables/HistDecompCTM.mat")
hd_mat_ctm <- hd_mat$HistDecomp
hd_mat_ctm_cpi <- hd_mat_ctm[,6:10]
hd_mat_ctm_cpi_df <- as.data.frame(hd_mat_ctm_cpi)
colnames(hd_mat_ctm_cpi_df) <- paste0(shocknames, " shock")

hd_mat_ctm_cpi_df$Date <- as.Date(seq(as.yearqtr("2005-3"),as.yearqtr("2020-4"), 1/4)) #model's order matters!
hd_mat_ctm_cpi_df_bar <- gather(hd_mat_ctm_cpi_df, "Variable", "Value", -Date)
hd_mat_ctm_cpi_df_bar$Date <-  as.Date(hd_mat_ctm_cpi_df_bar$Date)
demeaned_cpi <- dat_unseas$cpi_core_qoq - mean(dat_unseas$cpi_all_qoq)

ggplot(hd_mat_ctm_cpi_df_bar, aes(x = Date, y = Value, fill = Variable)) + 
  geom_bar(stat="identity") +  
  scale_x_date(breaks=hd_mat_ctm_cpi_df_bar$Date, date_labels="%b %y", 
      date_breaks = "4 months", expand=c(0,0)) +
  theme_minimal() + geom_hline(yintercept = 0, col="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 8),
        axis.title.x = element_text(vjust = 0.5, size = 9),
        axis.title.y = element_text(vjust = 0.5, size = 9),
        legend.position="bottom") +
  xlab("Date") + ylab("Contribution to CPI gr. r.") +
  geom_line(aes(x=as.Date(rep(hd_mat_ctm_cpi_df$Date, 5)), y=rep(demeaned_cpi[-1], 5)), col="black", size=0.3) +##383838
  scale_fill_brewer(palette="Dark2")

ggsave(filename="hd_cpi_full_LT_core.eps",
       path="/Users/rutra/ВШЭ/Магистратура/Thesis/Text/figures/",
       device="eps",
       width=320, height=210, dpi=320, units = "mm", limitsize=FALSE)

ggplot(hd_mat_ctm_cpi_df_bar, aes(x = Date, y = Value, fill = Variable)) + 
  geom_bar(stat="identity") +  
  scale_x_date(breaks=hd_mat_ctm_cpi_df_bar$Date, date_labels="%b %y", 
               expand=c(0,0), 
               limits=c(as.Date("2015-01-02"), NA)) +
  ylim(c(-0.025, 0.02)) +
  theme_minimal() + geom_hline(yintercept = 0, col="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8),
        axis.text.y = element_text(vjust = 0.5, hjust=1, size = 8),
        axis.title.x = element_text(vjust = 0.5, size = 9),
        axis.title.y = element_text(vjust = 0.5, size = 9),
        legend.position="bottom") +
  xlab("Date") + ylab("Contribution to CPI gr. r.") + 
  geom_line(aes(x=as.Date(rep(hd_mat_ctm_cpi_df$Date, 5)), 
                y=rep(demeaned_cpi[-1], 5)), col="black", size=0.3) +
  scale_fill_brewer(palette="Dark2")

ggsave(filename="hd_cpi_cut_LT_core.eps",
       path="/Users/rutra/ВШЭ/Магистратура/Thesis/Text/figures/",
       device="eps",
       width=320, height=210, dpi=320, units = "mm", limitsize=FALSE)

#FEVD
