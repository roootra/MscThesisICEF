library(ggplot2)
library(rio)
library(lubridate)
library(seasonal)
library(reshape2)
library(tidyr)
library(scales)
library(tseries)

Sys.setlocale("LC_TIME", "C")
dat <- import("/Users/rutra/ВШЭ/Магистратура/Thesis/Data/Aggregated data/Data Monthly.xlsx", sheet=5)
dat <- dat[NROW(dat):1,]

#Deseasonalization
to_deseasonalize <- c("ind_output_yoy","imp_price_mom", "exp_price_mom", 
                      grep("cpi.*", colnames(dat), value=T))
dat_unseas <- dat
for(colname in to_deseasonalize){
  current_ts <- ts(dat[,colname], start=c(2005, 3), frequency=12)
  current_seas <- seas(current_ts, transform.function="none",
                       #regression.aictest=NULL, outlier=NULL,
                       #automdl=NULL,
                       seats.noadmiss="yes")
  dat_unseas[,colname] <- as.numeric(final(current_seas))
}

data_draw <- dat_unseas[,c("date", "oil_USD", "imp_price_mom", 
                           "miacr_31", "neer_mom", "gdp_per_cap_SA", "cpi_all_mom")]
colnames(data_draw) <- c("Date", "Oil price gr. rate", "World import price gr. rate", 
                         "MIACR 31-180 days", "NEER gr. rate", "GDP per capita gr. rate",
                         "CPI gr. rate")
data_draw$Date <- as.Date(data_draw$Date)
data_draw_long <- gather(data_draw, "Variable", "Value", -Date)
months <- month(data_draw$Date)
ggplot(data_draw_long, aes(x = Date, y = Value), size = 0.1) + geom_line() + 
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
  geom_line() + theme_minimal() + #+theme_bw()
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

