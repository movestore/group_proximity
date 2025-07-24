library(move2)
library(ggplot2)
library(dplyr)
library(units)
# Herd protection ICARUS MPIAB 2025 Apel - 6400315030 - 5cow + 2horse
# Herd protection ICARUS MPIAB FVA 2025 Marc Boehler - 6214896791 - 58 Capra hircus(19), Bos taurus(48)
# Herd protection ICARUS MPIAB FVA 2025 Winterhalter - 6278745305 - 86 cow

ap <- movebank_download_study(6400315030, sensor_type_id="gps") #, taxon_canonical_name="Bos taurus")
bo <- movebank_download_study(6214896791, sensor_type_id="gps")
hist(mt_speed(bo))
table(mt_track_data(bo)$taxon_canonical_name)
ggplot()+
  geom_sf(data=mt_track_lines(bo), aes(color=taxon_canonical_name))
wi <- movebank_download_study(6278745305, sensor_type_id="gps")
head(mt_time(wi))

####################

ap <- movebank_download_study(6400315030, sensor_type_id="gps")
# bo <- movebank_download_study(6214896791, sensor_type_id="gps")
# wi <- movebank_download_study(6278745305, sensor_type_id="gps")

cws <- wi

## arrange track
cws <- dplyr::arrange(cws, mt_track_id(cws), mt_time(cws)) 
## remove duplicated timestamps
cws <- cws %>% 
  mutate(n_na = rowSums(is.na(pick(everything())))) %>%
  dplyr::arrange(n_na) %>%
  mt_filter_unique(criterion='first') %>% 
  dplyr::arrange(mt_track_id()) %>% 
  dplyr::arrange(mt_track_id(),mt_time())

cw <- filter_track_data(cws, taxon_canonical_name=="Bos taurus")

cw_speeds <-  as.numeric(mt_speed(cw, "m/s"))
hist(cw_speeds)
quantile(cw_speeds, probs=seq(0.99,1,0.0001),na.rm=T)

max_speed <- set_units(1, m/s) 
cw_f <- cw %>% filter(mt_speed(.)<=max_speed | is.na(mt_speed(.)))

# cw_f <- cw
# while(any(mt_speed(cw_f) > max_speed, na.rm = TRUE)){
#   cw_f <- cw_f %>% filter(mt_speed(.)<=max_speed | is.na(mt_speed(.)))
# }

ggplot()+ geom_sf(data=mt_track_lines(cw), aes(color=individual_local_identifier))+theme_bw()#+ facet_wrap(~individual_local_identifier)
ggplot()+ geom_sf(data=mt_track_lines(cw_f), aes(color=individual_local_identifier), show.legend = F)+theme_bw()


# cw_tl <-  as.numeric(mt_time_lags(cw_f, "min"))
# hist(cw_tl)
# quantile(cw_tl, probs=seq(0.99,1,0.0001),na.rm=T)
# 
# cw_f$segments <- mt_segments(cw_f)
# cw_f$cw_tl <- cw_tl
# 
# ggplot()+ geom_sf(data=cw_f,aes(geometry = segments, color=as.numeric(location_error_numerical)))+theme_bw()
# 
mt_track_id_column(cw_f)
names(cw_f)
names(mt_track_data(cw_f))
cw_f <- mt_as_event_attribute(cw_f,individual_local_identifier)
df <- data.frame(ts=mt_time(cw_f), id=cw_f$individual_local_identifier)

length(unique(df$id))

plot(df$ts,df$id)

library(lubridate)
df16jun <- df[day(df$ts)==10,]
plot(df16jun$ts,df16jun$id)

df16jun12 <- df16jun[hour(df16jun$ts)==12,]
plot(df16jun12$ts,df16jun12$id)

ts_min <- min(df16jun12$ts, na.rm = TRUE)
ts_max <- max(df16jun12$ts, na.rm = TRUE)
ten_min_intervals <- seq(from = ts_min, to = ts_max, by = "10 min")
abline(v = ten_min_intervals, col = "red", lty = "dotted")

barplot((minute(df16jun12$ts)))

df16jun12$mins <- minute(df16jun12$ts)

ggplot(data=df16jun12[df16jun12$mins<30,])+geom_bar(aes(mins))

summary(mt_time_lags(cw_f,"min"))

