library(move2)
library(ggplot2)
library(dplyr)
library(units)
library(sf)
library(suncalc)
library(lubridate)
library(data.table) # for rleid

# logger.fatal(), logger.error(), logger.warn(), logger.info(), logger.debug(), logger.trace()

# Herd protection ICARUS MPIAB 2025 Apel - 6400315030 - 5cow + 2horse
# Herd protection ICARUS MPIAB FVA 2025 Marc Boehler - 6214896791 - 58 Capra hircus(19), Bos taurus(48)
# Herd protection ICARUS MPIAB FVA 2025 Winterhalter - 6278745305 - 86 cow

ap <- movebank_download_study(6400315030, sensor_type_id="gps")
bo <- movebank_download_study(6214896791, sensor_type_id="gps", timestamp_start=as.POSIXct(Sys.time()-(5*24*60*60), tz="UTC"))
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

### apps in WF to clean data
cw <- filter_track_data(cws, taxon_canonical_name=="Bos taurus")
max_speed <- set_units(1, m/s) 
cw_f <- cw %>% filter(mt_speed(.)<=max_speed | is.na(mt_speed(.)))

data <- cw_f


### parameters
timesteps <- "10 min"
minNbAn <- 40 ## nb animal in proximity
proxdist <- 20 ## mtrs around cow

#seelasXdays <- 5 ## this should be regulated by the movebank download app
#dateToDisplayAboveAll <- Sys.Date()-1 ## this should be always last 24h of the data

rFunction = function(data, ) {
  ## interpolate data to the selected timestep
  data_int <- mt_interpolate(data, timesteps, omit = T)
  ## make indiv name the track id
  if("individual_local_identifier"!=mt_track_id_column(data_int)){
    data_int <- mt_as_event_attribute(data_int,individual_local_identifier, .keep=T)
    data_int <- mt_set_track_id_column(data_int,"individual_local_identifier") 
  }
  ## move columns of interest to event data table
  if("group_id"%in%names(mt_track_data(data_int))){
    data_int <- mt_as_event_attribute(data_int,group_id, .keep=F)
  }
  if(!"group_id"%in%names(data_int)){
    data_int$group_id <- NA
    logger.info("The animals of the dataset are not assigned to any groups/herds. If you have this information you can add it in movebank under the attribute 'animal group ID'. See documentation for details.")
    }
  ## reduce data set to needed columns
  data_int_red <- data_int |> select(timestamp,geometry, individual_local_identifier, group_id)
  
  ## add sunrise/sunset & day/night to each loc and local time for mean coords of data set
  crds <- st_coordinates(data_int_red)
  lat <- mean(crds[,2], rm.na=T)
  lon <- mean(crds[,1], rm.na=T)
  
  local_tz <- lutz::tz_lookup_coords(lat, lon, method = "accurate")
  data_int_red$local_time <- with_tz(mt_time(data_int_red), local_tz)
  
  
  srst <- data.frame(date = lubridate::date(data_int_red$local_time), lat = st_coordinates(data_int_red)[,2], lon = st_coordinates(data_int_red)[,1])
  
  data_int_red$sunrise <- getSunlightTimes(data=srst, keep = c("sunrise"), tz = local_tz)["sunrise"][,1]
  data_int_red$sunset <- getSunlightTimes(data=srst, keep = c("sunset"), tz = local_tz)["sunset"][,1]
  
  data_int_red <- data_int_red %>%
    mutate(day_night = if_else(mt_time(.) >= sunrise & mt_time(.) < sunset, "day", "night"))  %>%
    mutate(period_run = rleid(day_night)) %>%
    group_by(period_run) %>%
    mutate(day_night_date = paste0(as.Date(first(timestamp))," ", day_night)) %>%
    ungroup()
  data_int_red$period_run <- NULL
  
  ## create table to match indiv with group_id again in later steps
  indGroup <- as.data.frame(data_int_red) |> select(individual_local_identifier,group_id)
  indGroup <- indGroup[!duplicated(indGroup$individual_local_identifier),]
  gr_nb <- data.frame(table(indGroup$group_id))
  colnames(gr_nb) <- c("group_id","gr_nb")
  indGroup <- merge(indGroup,gr_nb, by="group_id")
  indGroup$group_id_nb <- paste0(indGroup$group_id," (",indGroup$gr_nb,")") ## add in() the tot bumber of indv in that group
  
  ## reprojecting to center of trajectory
  aeqd_crs <- mt_aeqd_crs(data_int_red, "centroid", "m")
  data_int_red_pj <- sf::st_transform(data_int_red, aeqd_crs)
  ## tranf to sf
  data_int_red_sf <- data_int_red_pj; class(data_int_red_sf) <- class(data_int_red_pj) %>% setdiff("move2")
  
  ## split data by timestamps
  data_int_red_sf_L <- split(data_int_red_sf,data_int_red_sf$timestamp)
  lgthGr <- unlist(lapply(data_int_red_sf_L, nrow))
  
  plot(unique(data_int_red_sf$timestamp),lgthGr, xlab="", ylab="Number of animals", ylim=c(min(lgthGr),max(lgthGr)+5))
  abline(h=40, col="red")
  legend("topright", legend="minimum number of animal tagged",
         col="red", lty=1, cex=0.8)
  
  ## only will calculate distance if min number of animals is present in the timeslot
  ## counting all indiv
  # x <- which(lgthGr >= minNbAn)[1]
  prox_L <- lapply(which(lgthGr >= minNbAn), function(x){
    tmsl <- data_int_red_sf_L[[x]]
    mtrx <- st_coordinates(tmsl)
    rownames(mtrx) <- tmsl$individual_local_identifier
    dstm <- as.matrix(round(dist(mtrx, upper = TRUE)))
    diag(dstm) <- NA  # Exclude self-distances
    
    # Count neighbors for all individuals at once
    neigh_counts <- colSums(dstm > 0 & dstm < proxdist, na.rm = TRUE)
    data.frame(
      time = tmsl$timestamp[1],
      local_time = tmsl$local_time[1],
      day_night = tmsl$day_night[1],
      sunrise = tmsl$sunrise[1],
      sunset = tmsl$sunset[1],
      individual_local_identifier = colnames(dstm),
      nb_neigh = neigh_counts
    )
  })
  prox_df <- do.call(rbind, prox_L)
  prox_df <- merge(prox_df,indGroup, by.x = "individual_local_identifier", by.y = "individual_local_identifier")
  # head(prox_df)
  
  
  # ## counting only those of the same herd
  # # x <- which(lgthGr >= minNbAn)[1]
  # prox_h_L <- lapply(which(lgthGr >= minNbAn), function(x){
  #   tmsl <- data_int_red_sf_L[[x]]
  #   mtrx <- st_coordinates(tmsl)
  #   rownames(mtrx) <- tmsl$individual_local_identifier
  #   dstm <- as.matrix(round(dist(mtrx, upper = TRUE)))
  #   diag(dstm) <- NA  # Exclude self-distances
  #   
  #   # Count neighbors for all individuals at once
  #   neigh_counts <- colSums(dstm > 0 & dstm < proxdist, na.rm = TRUE)
  #   data.frame(
  #     time = tmsl$timestamp[1],
  #     indiv = colnames(dstm),
  #     nb_neigh = neigh_counts
  #   )
  # })
  # prox_df_h <- do.call(rbind, prox_h_L)
  # prox_df_h <- merge(prox_df_h,indGroup, by.x = "indiv", by.y = "individual_local_identifier")
  # 
  
  
  
  
 #  ## display shaded day night,  display local time on axis where useful
   
  # crds <- st_coordinates(data_int_red)
  # lat <- mean(crds[,2], rm.na=T)
  # lon <- mean(crds[,1], rm.na=T)
  # 
  # local_tz <- lutz::tz_lookup_coords(lat, lon, method = "accurate")
  # prox_df$local_time <- with_tz(prox_df$time, local_tz)
  # 
 
 #  ## do 1 plot all indiv lumbipng all data to 0-12h + one day over (to see if it different for avarage) => something seems off with the one day....
  prox_df$time_lc <- format(prox_df$local_time, format = "%H:%M:%S")
  prox_df_summry<- prox_df %>%
    group_by(time_lc) %>%
    summarize(
      mean = mean(nb_neigh, na.rm = TRUE),
      sd = sd(nb_neigh, na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n),
      lower_ci = mean - qt(0.975, n - 1) * se,
      upper_ci = mean + qt(0.975, n - 1) * se,
      day_night =  unique(day_night),               
      sunrise = unique(sunrise),
      sunset = unique(sunset)
    )
  
  
  prox_df_summry$time_lc <- as.POSIXct(paste0(max(lubridate::date(prox_df$local_time))," ",prox_df_summry$time_lc))

  sub_day <- prox_df[lubridate::date(prox_df$time)%in%max(lubridate::date(prox_df$local_time)),]
  sub_day_summry<- sub_day %>%
    group_by(local_time) %>%
    summarize(
      mean = mean(nb_neigh, na.rm = TRUE),
      sd = sd(nb_neigh, na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n),
      lower_ci = mean - qt(0.975, n - 1) * se,
      upper_ci = mean + qt(0.975, n - 1) * se,
      day_night =  unique(day_night),               
      sunrise = unique(sunrise),
      sunset = unique(sunset)
    )

  head(sub_day_summry)


  sub_summry_ss <- sub_day_summry[!duplicated(lubridate::date(sub_day_summry$local_time)),]
  sr <- sub_summry_ss$sunrise
  ss <- sub_summry_ss$sunset
  
  ggplot(prox_df_summry, aes(x = time_lc, y = mean))+
    annotate("rect",
             xmin = sr,# + days(1),  # Adjust if necessary
             xmax = ss,
             ymin = -Inf, ymax = Inf,
             alpha = 0.2, fill = "yellow")+
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "blue") +
   geom_line(data=sub_day_summry, aes(x = local_time, y = mean),color = "red") +
    geom_ribbon(data=sub_day_summry, aes(x = local_time, y = mean,ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "red") +
    # geom_point(size = 2) +
    theme_minimal() +
    labs(
      title = "Mean Value with 95% Confidence Interval",
      y = "Mean value",
      x = "Time"
    )

  
  ## do 1 plot all indiv last 5 days 
  
    # subset to last 5 days
  # tslts <- gsub(" .*$", "", names(prox_L))
  # tslts <- tslts[!duplicated(tslts)]
  # 
  # sub <- prox_df#[lubridate::date(prox_df$time)%in%tail(tslts,seelasXdays),]
  sub_summry<- prox_df %>%
    group_by(local_time) %>%
    summarize(
      mean_nb_neigh = mean(nb_neigh, na.rm = TRUE),
      sd = sd(nb_neigh, na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n),
      lower_ci = mean_nb_neigh - qt(0.975, n - 1) * se,
      upper_ci = mean_nb_neigh + qt(0.975, n - 1) * se,
      day_night =  unique(day_night),               
      sunrise = unique(sunrise),
      sunset = unique(sunset)
    )
  
  sub_summry_ss <- sub_summry[!duplicated(lubridate::date(sub_summry$local_time)),]
  sr <- sub_summry_ss$sunrise
  ss <- sub_summry_ss$sunset
  
  ggplot(sub_summry, aes(x = local_time, y = mean_nb_neigh))+
    annotate("rect",
             xmin = sr,# + days(1),  # Adjust if necessary
             xmax = ss,
             ymin = -Inf, ymax = Inf,
             alpha = 0.2, fill = "gold")+
    # geom_rect(
    #   aes(xmin = sunrise, xmax = sunset, fill = day_night), 
    #   ymin = -Inf, ymax = Inf, alpha = 0.8, 
    #   data = sub_summry_ss
    # ) +
    geom_line(color = "blue") +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "blue") +
    # geom_point(size = 2) +
    theme_minimal() +
    labs(
      title = "Mean Value with 95% Confidence Interval",
      y = "Mean number of neighbours",
      x = ""
    )
  
  ## do 1 plot all indiv last 5 days + line of a selected individual - optional - or create pdf with one per indiv?
  
  ## do 1 plot. dist calc only for within group, one line per group. Add in ledgend tot nb ind in group - optional (or of group exisits?)
  
  
  # subset to last 5 days
  tslts <- gsub(" .*$", "", names(prox_L))
  tslts <- tslts[!duplicated(tslts)]
  
  sub <- prox_df[lubridate::date(prox_df$time)%in%tail(tslts,seelasXdays),]
  sub_summry<- sub %>%
    group_by(local_time,group_id_nb) %>%
    summarize(
      mean = mean(nb_neigh, na.rm = TRUE),
      sd = sd(nb_neigh, na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n),
      lower_ci = mean - qt(0.975, n - 1) * se,
      upper_ci = mean + qt(0.975, n - 1) * se
    )
  
  # timeSeq <- data.frame(local_time=seq.POSIXt(min(sub_summry$local_time),max(sub_summry$local_time), by=timesteps))
  sun_times <- getSunlightTimes(date = as.Date(tail(tslts,seelasXdays)), lat = lat, lon = lon, keep = c("sunrise", "sunset"), tz = local_tz)
  sunrise_utc <- sun_times$sunrise
  sunset_utc  <- sun_times$sunset
  
  ggplot(sub_summry, aes(x = local_time, y = mean, group=group_id_nb))+
    annotate("rect",
             xmin = sunrise_utc,# + days(1),  # Adjust if necessary
             xmax = sunset_utc,
             ymin = -Inf, ymax = Inf,
             alpha = 0.2, fill = "yellow")+
    geom_line(aes(color=group_id_nb)) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "blue") +
    # geom_point(size = 2) +
    theme_minimal() +
    labs(
      title = "Mean Value with 95% Confidence Interval",
      y = "Mean value",
      x = "Time"
    )+ facet_wrap(~group_id_nb)
  
  #####################
  ## plot last 5 days and today as dots?
  ggplot(prox_df)+geom_boxplot(aes(individual_local_identifier,nb_neigh, color=day_night))+facet_wrap(~group_id_nb, scale="free_x")
  
  prox_df_tst <- prox_df
  prox_df_tst$t_range <- "last 5 days"
  prox_df_tst$t_range[lubridate::date(prox_df$local_time)==max(lubridate::date(prox_df$local_time))] <- "last day"

  ggplot(prox_df_tst)+
    geom_boxplot(aes(individual_local_identifier,nb_neigh, color=t_range, fill=day_night))+
    scale_color_manual(values=c("last 5 days"="green4", "last day"="red"))+
    scale_fill_manual(values=c("day"="orange", "night"="cyan"))+
    facet_wrap(~group_id, scale="free_x")+
    theme_bw()+
    guides(
      color = guide_legend("time range"),   # Legend for boxplot colors
      fill  = guide_legend("time slot"),     # Legend for points fill
    )
  #####################
  ### nuber of other cows in buff X throught the day /night
  ###########
  sub_day <- prox_df[lubridate::date(prox_df$time)%in%dateToDisplayAboveAll,]
  sub_day_summry<- sub_day %>%
    group_by(local_time) %>%
    summarize(
      mean = mean(nb_neigh, na.rm = TRUE),
      sd = sd(nb_neigh, na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n),
      lower_ci = mean - qt(0.975, n - 1) * se,
      upper_ci = mean + qt(0.975, n - 1) * se
    )
  
  
  ################
  ## daily distance move. split up by night and day, show boxplox for all data and dot of last day??
  
  data_int_red_L <- split(data_int_red, as.character(mt_track_id(data_int_red)))
  
  # ind <- data_int_red_L[["Herd_Goat_1"]]
  dit_moved_L <- lapply(data_int_red_L, function(ind){
    print(unique(mt_track_id(ind)))
    ind_L <- split(ind, ind$day_night_date)
    cumDistDayL <- lapply(ind_L, function(x){sum(mt_distance(x, "km"),na.rm=T)})
    cumDistDay <- data.frame(do.call("rbind",cumDistDayL))
    colnames(cumDistDay) <- "cumDistDay"
    cumDistDay$day_night_date <- rownames(cumDistDay)
    # head(cumDistDay)

    maxNetDispL <- lapply(ind_L, function(x){max(st_distance(x=x[-nrow(x),],y=x[-1,], by_element=T),na.rm=T)})
    
    maxNetDisp <- data.frame(do.call("rbind",maxNetDispL))
    colnames(maxNetDisp) <- "maxNetDisp"
    maxNetDisp$day_night_date <- rownames(maxNetDisp)
    units(maxNetDisp$maxNetDisp) <- "km"
    # head(maxNetDisp)
    
    df <- data.frame(individual_local_identifier=ind$individual_local_identifier, group_id=ind$group_id , day_night_date=ind$day_night_date)
    distdf <- merge(maxNetDisp,cumDistDay, by="day_night_date")
df_dist <- merge(df,distdf,by="day_night_date")
df_dist <- df_dist[!duplicated(df_dist),]
    return(df_dist)
  })
  
  dist_moved_df <- data.frame(do.call("rbind",dit_moved_L))
  
  library(dplyr)
  library(tidyr)
  
  dist_moved_df <- dist_moved_df %>%
    separate(day_night_date, into = c("date", "day_night"), sep = " ")
  
  
  # ggplot(dist_moved_df)+geom_point(aes(date,maxNetDisp, color=day_night))+facet_wrap(~group_id)
  # ggplot(dist_moved_df)+geom_boxplot(aes(cumDistDay,individual_local_identifier, color=day_night))#+facet_wrap(~group_id)
  # ggplot(dist_moved_df)+geom_boxplot(aes(cumDistDay,group_id, color=day_night))#+facet_wrap(~group_id)
  
  dist_moved_df_lstday <- dist_moved_df[dist_moved_df$date==max(dist_moved_df$date),]
  
  ggplot(dist_moved_df)+
    geom_boxplot(aes(individual_local_identifier,cumDistDay, color=day_night))+
    geom_point(data=dist_moved_df_lstday, aes(individual_local_identifier,cumDistDay,fill=day_night, shape=day_night), size=2, color="black")+ 
    scale_shape_manual(values=c("day"=23, "night"=24))+
    scale_color_manual(values=c("day"="orange", "night"="cyan"))+
    scale_fill_manual(values=c("day"="orange", "night"="cyan"))+
    facet_wrap(~group_id, scale="free_x")+
    theme_bw()+
    guides(
      color = guide_legend("Last 5 days"),   # Legend for boxplot colors
     fill  = guide_legend("last day"),     # Legend for points fill
      shape = guide_legend("last day")     # Legend for points shape
    )
  
  
  
  ggplot(dist_moved_df)+
    geom_boxplot(aes(individual_local_identifier,maxNetDisp, color=day_night))+
    geom_point(data=dist_moved_df_lstday, aes(individual_local_identifier,maxNetDisp,fill=day_night, shape=day_night), size=2, color="black")+ 
    scale_shape_manual(values=c("day"=23, "night"=24))+
    scale_color_manual(values=c("day"="orange", "night"="cyan"))+
    scale_fill_manual(values=c("day"="orange", "night"="cyan"))+
    facet_wrap(~group_id, scale="free_x")+
    theme_bw()+
    guides(
      color = guide_legend("Last 5 days"),   # Legend for boxplot colors
      fill  = guide_legend("last day"),     # Legend for points fill
      shape = guide_legend("last day")     # Legend for points shape
    )
  
  
  ##############
  return(data)
}
