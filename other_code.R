library(move2)
library(ggplot2)
library(dplyr)
library(units)
library(sf)
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
bo <- movebank_download_study(6214896791, sensor_type_id="gps")
wi <- movebank_download_study(6278745305, sensor_type_id="gps")

cws <- bo

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

cw_f <- mt_as_event_attribute(cw_f,individual_local_identifier, .keep=T)
ggplot()+
  geom_sf(data=mt_track_lines(cw_f), aes(color=individual_local_identifier), show.legend = F)+facet_wrap(~individual_local_identifier)

sort(unique(cw_f$individual_local_identifier))
# Herd_2F4AWC
# Herd_YWWEWR
## 24.6. tag is off?

prind <- filter_track_data(cw_f, individual_local_identifier %in% c("Herd_YWWEWR"))#, "Herd_YWWEWR"))
ggplot()+
  geom_sf(data=mt_track_lines(prind), aes(color=individual_local_identifier))

prind$date <- lubridate::date(mt_time(prind))                    
ggplot()+
  geom_sf(data=prind, aes(color=as.factor(date)))#+facet_wrap(~individual_local_identifier)

mt_track_id(prind) <- as.factor(prind$date)

ggplot()+geom_sf(data=mt_track_lines(prind))+facet_wrap(~deployment_id)


##################
library(lubridate)
tr <- range(mt_time(cw_f))

lubridate::ceiling_date(tr[1],unit = "10 minutes")
lubridate::floor_date(tr[2],unit = "10 minutes")

# tseq <- seq.POSIXt(ceiling_date(tr[1],unit = "10 minutes"),floor_date(tr[2],unit = "10 minutes"),by="10 min")


cw_f_int <- mt_interpolate(cw_f, "10 min", omit = T)
if(!"individual_local_identifier"%in%names(cw_f_int)){
  cw_f_int <- mt_as_event_attribute(cw_f_int,individual_local_identifier, .keep=T)
}
if("group_id"%in%names(mt_track_data(cw_f_int))){
  cw_f_int <- mt_as_event_attribute(cw_f_int,group_id, .keep=T)
}

cw_red <- cw_f_int |> select(timestamp,geometry, individual_local_identifier, group_id)

aeqd_crs <- mt_aeqd_crs(cw_red, "centroid", "m")
cw_red_pj <- sf::st_transform(cw_red, aeqd_crs)

cw_red_sf <- cw_red_pj; class(cw_red_sf) <- class(cw_red_pj) %>% setdiff("move2")

cw_red_sf_L <- split(cw_red_sf,cw_red_sf$timestamp)

lgthL <- lapply(cw_red_sf_L, nrow)
lgth <- unlist(lgthL)
summary(lgth)
plot(unique(cw_red_sf$timestamp),lgth)
which(lgth==median(lgth))[1]

minNbAn <- 40 ## nb animal in proximity
proxdist <- 20 ## mtrs around cow

prox_L <- lapply(which(lgth>=minNbAn), function(x){

tt <- cw_red_sf_L[[x]]

mtrx <- st_coordinates(tt)
rownames(mtrx) <- tt$individual_local_identifier
dst <- round(dist(mtrx,upper=T))

# hc <- hclust(dst)
# plot(hc)
# dend <- as.dendrogram(hc)
# 
# car_type <- factor(tt$group_id)
# n_car_types <- length(unique(car_type))
# cols_4 <- colorspace::rainbow_hcl(n_car_types, c = 70, l  = 50)
# col_car_type <- cols_4[car_type]
# labels_colors(dend) <- col_car_type[order.dendrogram(dend)]
# plot(dend)
# legend("topleft", legend = levels(car_type), fill = cols_4)



dstm <- as.matrix(dst)
dist_l <- lapply(1:ncol(dstm), function(x){
  dind <- dstm[,x]
  length(dind[dind>0&dind<proxdist])
})

df <- data.frame(time=tt$timestamp, neigh=unlist(dist_l))

return(df)
})

prox_df <- do.call("rbind",prox_L)

sub <- prox_df[lubridate::date(prox_df$time)=="2025-05-29",]

library(ggplot2)
library(suncalc)
library(lubridate)

crds <- st_coordinates(cw_red)



lat <- mean(crds[,2], rm.na=T)
lon <- mean(crds[,1], rm.na=T)  
date <- as.Date("2025-05-29")

sun_times <- getSunlightTimes(date = date, lat = lat, lon = lon, keep = c("sunrise", "sunset"), tz = "UTC")
sunrise_utc <- sun_times$sunrise
sunset_utc  <- sun_times$sunset

shade <- data.frame(
  xmin = sunset_utc,
  xmax = sunrise_utc + 1,  # Add 1 day if sunrise is on the next day
  ymin = -Inf,
  ymax = Inf
)

summary_df <- sub %>%
  group_by(time) %>%
  summarize(
    mean = mean(neigh, na.rm = TRUE),
    sd = sd(neigh, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci = qt(0.975, df = n - 1) * se
  )


ggplot(summary_df, aes(x = time, y = mean)) +
  geom_line(color = "blue") +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.2, fill = "blue") +
  # geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "Mean Value with 95% Confidence Interval",
    y = "Mean value",
    x = "Time"
  )+
  annotate("rect",
           xmin = sunset_utc,
           xmax = sunrise_utc,
           ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "blue")
