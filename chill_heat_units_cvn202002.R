

#  This code analyses hourly temperature data and produces time
#  series graphs of accumulated chill and heat units for the
#  Willcox Bench station of the Arizona Meteorological (AZMET) 
#  Network (https://cals.arizona.edu/azmet). Graphs are used in
#  the 2020 February issue of the Climate Viticulture Newsletter:
#  https://cals.arizona.edu/research/climategem/content/climate-viticulture-newsletter-2020-february

#  Hourly temperature units are degrees Celsius.

#  Chill units are calculated as chilling hours and chill 
#  portions. See details from package 'chillR'. Heat units are
#  calculated as growing degree days.

#  Author:
#  Jeremy Weiss, Climate and Geospatial Extension Scientist
#  School of Natural Resources and the Environment
#  University of Arizona
#  520-626-8063, jlweiss@email.arizona.edu


##################################################################
##  A. LOAD, FILTER, AND TRANSFORM AZMET STATION DATA
##################################################################


#####  SETUP ITEMS INTRODUCED IN THIS SECTION


#  Load required libraries.
library( "tidyverse" )

#  Data collection at the Willcox Bench station started in middle
#  2016 and continues at present. 'obs_yrs' format matches data
#  filename convention.
obs_yrs <- seq( 16,20 )


#####


#  Data column names are from:
#  https://cals.arizona.edu/azmet/raw2003.htm
col_names <- c( "Year","JDay","Hour","Temp","rh","vpd",
                "sr","prcp","4sm","20sm","wavg","wvm",
                "wvd","wdsd","mws","reto","avp","dp" )

#  Load data, which is in files for individual years.
for ( iyr in 1:length( obs_yrs ) ) {
  obs <- read.table( paste0( "39", obs_yrs[ iyr ], "rh.txt" ),
                     header = FALSE,
                     sep=',',
                     fill = TRUE )
  colnames( obs ) <- col_names
  obs <- select( obs, Year, JDay, Hour, Temp )
  obs_date <- as.Date( obs$JDay-1,
                       origin = paste0( "20", obs_yrs[ iyr ], "-01-01" ) )
  obs[ "Date" ] <- as.character( obs_date )
  obs[ "Month" ] <- as.numeric( format( obs_date, "%m" ) )
  obs[ "Day" ]  <- as.numeric( format( obs_date, "%d" ) )
  rm( obs_date )
  
  #  Concatenate station data from different years.
  if ( iyr == 1 ) {
    stn_data <- obs
  } else {
    stn_data <- rbind( stn_data, obs )
  }
  
  rm( obs )
}

rm( col_names, iyr )


##################################################################
##  B. CALCULATE ACCUMULATED HEAT UNITS
##################################################################


#####  SETUP ITEMS INTRODUCED IN THIS SECTION


#  Load required libraries.
library( "lubridate" )

#  Set base temperature (in degrees F) for growing degree day
#  calculation, following National Phenology Network approach to 
#  calculate growing degree days, described here:
#  https://www.usanpn.org/data/agdd_maps.
base_temp <- 50

#  Set the start and end days-of-year for the period of interest
#  for analysis by non-leap-year values. These values should match
#  the overall period of interest (pertaining to both chill and
#  heat units), in order to facilitate plotting below. Heat units
#  will start accumulating on December 1 for a given winter year.
#  This is accounted for in the code for this section.
doy_start <- 305  # November 1
doy_end <- 120  # April 30


#####


#  Aggregate hourly data to daily data.
stn_data_daily <- stn_data %>%
  group_by( Year, JDay, Date, Month, Day ) %>%
  summarize( avgTempC = mean( Temp, na.rm = TRUE ) )

#  Convert daily temperature units from degrees Celsius to degrees
#  Fahrenheit, as we will follow National Phenology Network 
#  approach to calculate growing degree days, described here:
#  https://www.usanpn.org/data/agdd_maps.
stn_data_daily[ "avgTempF" ] <- NA
stn_data_daily$avgTempF <- ( 1.8 * stn_data_daily$avgTempC ) + 32

#  Calculate growing degree days (GDDs), based on 'avgTempF'. If 
#  the Tmean value for a given day is present and less than or 
#  equal to 'base_temp', GDD equals 0. Otherwise, GDD equals the 
#  difference between the daily average temperature and 
#  'base_temp'.
stn_data_daily[ "GDDs" ] <- NA

for ( d in 1:nrow( stn_data_daily ) ) {
  if ( stn_data_daily$avgTempF[ d ] <= base_temp ) {
    stn_data_daily$GDDs[ d ] <- 0
  } else {
    stn_data_daily$GDDs[ d ] <- stn_data_daily$avgTempF[ d ] - base_temp
  }
}

rm( d )

#  Calculate accumulated growing degree days (AGDDs) and put into
#  new dataframe, accounting for leap years. Start date for 
#  calculations is December 1 of a given winter.
years <- seq( as.numeric( paste0( 20, obs_yrs[ 1 ] ) ),
              as.numeric( paste0( 20, obs_yrs[ length( obs_yrs ) ] ) ) )
col_names <- c( "Winter", "iDay", "JDay", "AGGDs" )

for ( iyr in 1:length( years )-1 ) {
  
  winter_yr1 <- years[ iyr ]
  winter_yr2 <- years[ iyr+1 ]
  
  #####  CONDITION 1
  
  #  Condition 1: when first calendar year of a winter year is a
  #  leap year
  if ( leap_year( winter_yr1 ) == TRUE & leap_year( winter_yr2 ) == FALSE ) {
    
    #  Filter data by given winter year.
    winter_data <- rbind( stn_data_daily %>%
                            filter( Year == winter_yr1 & JDay >= doy_start+1 ),
                          stn_data_daily %>%
                            filter( Year == winter_yr2 & JDay <= doy_end ) )
    
    #  Initialize and setup dataframe that will store AGDD values 
    #  for the given winter year.
    winter_agdds <- data.frame( matrix( data = NA, nrow = nrow( winter_data ), ncol = 4 ) )
    colnames( winter_agdds ) <- col_names
    winter_agdds$Winter <- paste( years[ iyr ], years[ iyr+1 ], sep = "-" )
    winter_agdds$iDay <- seq( 1, nrow( winter_data ) )
    winter_agdds$JDay <- winter_data$JDay
    
    #  Calculate AGDD values for the given winter year.
    for ( d in 1:nrow( winter_data ) ) {
      if ( winter_data$JDay[ d ] >= doy_start+1 & winter_data$JDay[ d ] <= doy_start+30 ) {
        winter_agdds$agdds[ d ] = NA  #  November values
      } else if ( winter_data$JDay[ d ] == doy_start+30+1 ) {
        winter_agdds$AGDDs[ d ] = winter_data$GDDs[ d ]  #  December 1 value
      } else {
        winter_agdds$AGDDs[ d ] = winter_data$GDDs[ d ] + winter_agdds$AGDDs[ d-1 ]  #  all other day values
      }
    }
    
  }
  
  
  
  
  
  
  
  
  
  
  
  #####  CONDITION 2
  
  #  Condition 2: when neither calendar years of a winter year are
  #  leap years
  else if ( leap_year( winter_yr1 ) == FALSE & leap_year( winter_yr2 ) == FALSE ) {
    
    #  Filter data by given winter year.
    winter_data <- rbind( stn_data_daily %>%
                            filter( Year == winter_yr1 & JDay >= doy_start ),
                          stn_data_daily %>%
                            filter( Year == winter_yr2 & JDay <= doy_end ) )
    
    #  Initialize and setup dataframe that will store AGDD values 
    #  for the given winter year.
    winter_agdds <- data.frame( matrix( data = NA, nrow = nrow( winter_data ), ncol = 4 ) )
    colnames( winter_agdds ) <- col_names
    winter_agdds$Winter <- paste( years[ iyr ], years[ iyr+1 ], sep = "-" )
    winter_agdds$iDay <- seq( 1, nrow( winter_data ) )
    winter_agdds$JDay <- winter_data$JDay
    
    #  Calculate AGDD values for the given winter year.
    
    
    
    #  Condition 2: when neither calendar years of a winter year is
    #  a leap year
    for ( d in 1:nrow( y ) ) {
      if ( y$JDay[ d ] >= 305 & y$JDay[ d ] < 335 ) {
        aggds1718$CGGDs[ d ] = NA
      } else if ( y$JDay[ d ] == 335 ) {
        aggds1718$CGGDs[ d ] = y$GDDs[ d ]
      } else {
        aggds1718$CGGDs[ d ] = y$GDDs[ d ] + aggds1718$CGGDs[ d-1 ]
      }
    }
    
    for ( d in 1:nrow( y ) ) {
      if ( y$JDay[ d ] >= 305 & y$JDay[ d ] < 335 ) {
        aggds1819$CGGDs[ d ] = NA
      } else if ( y$JDay[ d ] == 335 ) {
        aggds1819$CGGDs[ d ] = y$GDDs[ d ]
      } else {
        aggds1819$CGGDs[ d ] = y$GDDs[ d ] + aggds1819$CGGDs[ d-1 ]
      }
    }
    
  }
  
  
  
  
  #####  CONDITION 3
  
  #  Condition 3: when second calendar year of a winter year is a
  #  leap year
  else if ( leap_year( winter_yr1 ) == FALSE & leap_year( winter_yr2 ) == TRUE ) {
    winter_data <- rbind( stn_data_daily %>%
                            filter( Year == winter_yr1 & JDay >= doy_start ),
                          stn_data_daily %>%
                            filter( Year == winter_yr2 & JDay <= doy_end+1 ) )
    
    
    #  Initialize and setup dataframe that will store AGDD values 
    #  for the given winter year.
    winter_agdds <- data.frame( matrix( data = NA, nrow = nrow( winter_data ), ncol = 4 ) )
    colnames( winter_agdds ) <- col_names
    winter_agdds$Winter <- paste( years[ iyr ], years[ iyr+1 ], sep = "-" )
    winter_agdds$iDay <- seq( 1, nrow( winter_data ) )
    winter_agdds$JDay <- winter_data$JDay
    
    
    #  Condition 3: when second calendar year of a winter year is a
    #  leap year
    for ( d in 1:nrow( y ) ) {
      if ( y$JDay[ d ] >= 305 & y$JDay[ d ] < 335 ) {
        aggds1920$CGGDs[ d ] = NA
      } else if ( y$JDay[ d ] == 335 ) {
        aggds1920$CGGDs[ d ] = y$GDDs[ d ]
      } else {
        aggds1920$CGGDs[ d ] = y$GDDs[ d ] + aggds1920$CGGDs[ d-1 ]
      }
    }
    
    rm( d )
    
    
    
    
    
    
  }
  
  #####
  
  
  
  
  
  
  

  
 
  
  
  
  
  
  
  
  
  
  
  
  #  Concatenate AGDD data from different winters.
  if ( iyr == 1 ) {
    agdds <- winter_agdds
  } else {
    agdds <- rbind( agdds, winter_agdds )
  }
  
  rm( winter_agdds, winter_data, winter_yr1, winter_yr2 )
  
}
rm( col_names, iyr )


##################################################################
##  C. CALCULATE ACCUMULATED CHILLING UNITS
##################################################################


#####  SETUP ITEMS INTRODUCED IN THIS SECTION


#  Load required libraries.
library( "chillR" )





#####


#  calculate daily chill portions for winter 2016-2017
stn_data_x <- filter( stn_data, Year == 2016 | Year == 2017 )
x <- chilling( hourtemps = stn_data_x,
               Start_JDay = 306,  # leap year value for November 1
               End_JDay = 120 )

# initialize and label columns of dataframe
chill1617 <- data.frame( matrix( data = NA,
                              nrow = x$Data_days[ which( x$Season == "2016/2017" ) ],
                              ncol = 7 ) )
colnames( chill1617 ) <- c( "Winter",
                         "iDay",
                         "JDay",
                         "Chilling_Hours",
                         "Utah_Model",
                         "Chill_portions",
                         "GDH" )
chill1617$Winter <- "2016-2017"
rm( x )

#  calculate daily chill portions for 2016 and write to 
#  dataframe
d <- 1
for ( jday in 306:366 ) {
  chill1617$iDay[ d ] <- d
  chill1617$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 305,  # leap year value for October 31
                 End_JDay = jday )
  chill1617$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2016 ) ]
  chill1617$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2016 ) ]
  chill1617$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2016 ) ]
  chill1617$GDH[ d ] <- x$GDH[ which( x$End_year == 2016 ) ]
  d <- d+1
  rm( x )
}
rm( jday )

#  calculate daily chill portions for 2017 and write to 
#  dataframe
for ( jday in 1:120 ) {
  chill1617$iDay[ d ] <- d
  chill1617$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 305,  # leap year value for October 31
                 End_JDay = jday )
  chill1617$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2017 ) ]
  chill1617$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2017 ) ]
  chill1617$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2017 ) ]
  chill1617$GDH[ d ] <- x$GDH[ which( x$End_year == 2017 ) ]
  d <- d+1
  rm( x )
}
rm( d,jday,stn_data_x )

matplot( chill1617$iDay,
         cbind( chill1617$Chilling_Hours / 10, chill1617$Chill_portions ),
         col = c( "black","red" ), type = "l")


#####

#  calculate daily chill portions for winter 2017-2018
stn_data_x <- filter( stn_data, Year == 2017 | Year == 2018 )
x <- chilling( hourtemps = stn_data_x,
               Start_JDay = 305,  # non-leap year value for November 1
               End_JDay = 120 )

# initialize and label columns of dataframe
chill1718 <- data.frame( matrix( data = NA,
                              nrow = x$Data_days[ which( x$Season == "2017/2018" ) ],
                              ncol = 7 ) )
colnames( chill1718 ) <- c( "Winter",
                         "iDay",
                         "JDay",
                         "Chilling_Hours",
                         "Utah_Model",
                         "Chill_portions",
                         "GDH" )
chill1718$Winter <- "2017-2018"
rm( x )

#  calculate daily chill portions for 2017 and write to 
#  dataframe
d <- 1
for ( jday in 305:365 ) {
  chill1718$iDay[ d ] <- d
  chill1718$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1718$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2017 ) ]
  chill1718$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2017 ) ]
  chill1718$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2017 ) ]
  chill1718$GDH[ d ] <- x$GDH[ which( x$End_year == 2017 ) ]
  d <- d+1
  rm( x )
}
rm( jday )

#  calculate daily chill portions for 2018 and write to 
#  dataframe
for ( jday in 1:120 ) {
  chill1718$iDay[ d ] <- d
  chill1718$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1718$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2018 ) ]
  chill1718$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2018 ) ]
  chill1718$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2018 ) ]
  chill1718$GDH[ d ] <- x$GDH[ which( x$End_year == 2018 ) ]
  d <- d+1
  rm( x )
}
rm( d,jday,stn_data_x )

matplot( chill1718$iDay,
         cbind( chill1718$Chilling_Hours / 10, chill1718$Chill_portions ),
         col = c( "black","red" ), type = "l")



#####

#  calculate daily chill portions for winter 2018-2019
stn_data_x <- filter( stn_data, Year == 2018 | Year == 2019 )
x <- chilling( hourtemps = stn_data_x,
               Start_JDay = 305,  # non-leap year value for November 1
               End_JDay = 120 )

# initialize and label columns of dataframe
chill1819 <- data.frame( matrix( data = NA,
                              nrow = x$Data_days[ which( x$Season == "2018/2019" ) ],
                              ncol = 7 ) )
colnames( chill1819 ) <- c( "Winter",
                         "iDay",
                         "JDay",
                         "Chilling_Hours",
                         "Utah_Model",
                         "Chill_portions",
                         "GDH" )
chill1819$Winter <- "2018-2019"
rm( x )

#  calculate daily chill portions for 2018 and write to 
#  dataframe
d <- 1
for ( jday in 305:365 ) {
  chill1819$iDay[ d ] <- d
  chill1819$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1819$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2018 ) ]
  chill1819$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2018 ) ]
  chill1819$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2018 ) ]
  chill1819$GDH[ d ] <- x$GDH[ which( x$End_year == 2018 ) ]
  d <- d+1
  rm( x )
}
rm( jday )

#  calculate daily chill portions for 2019 and write to 
#  dataframe
for ( jday in 1:120 ) {
  chill1819$iDay[ d ] <- d
  chill1819$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1819$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2019 ) ]
  chill1819$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2019 ) ]
  chill1819$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2019 ) ]
  chill1819$GDH[ d ] <- x$GDH[ which( x$End_year == 2019 ) ]
  d <- d+1
  rm( x )
}
rm( d,jday,stn_data_x )

matplot( chill1819$iDay,
         cbind( chill1819$Chilling_Hours / 10, chill1819$Chill_portions ),
         col = c( "black","red" ), type = "l")



#####

#  calculate daily chill portions for winter 2019-2020
stn_data_x <- filter( stn_data, Year == 2019 | Year == 2020 )
x <- chilling( hourtemps = stn_data_x,
               Start_JDay = 305,  # non-leap year value for November 1
               End_JDay = max( stn_data20$JDay ) )  #  ***UPDATE***

# initialize and label columns of dataframe
chill1920 <- data.frame( matrix( data = NA,
                              nrow = x$Data_days[ which( x$Season == "2019/2020" ) ],
                              ncol = 7 ) )
colnames( chill1920 ) <- c( "Winter",
                         "iDay",
                         "JDay",
                         "Chilling_Hours",
                         "Utah_Model",
                         "Chill_portions",
                         "GDH" )
chill1920$Winter <- "2019-2020"
rm( x )

#  calculate daily chill portions for 2019 and write to 
#  dataframe
d <- 1
for ( jday in 305:365 ) {
  chill1920$iDay[ d ] <- d
  chill1920$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1920$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2019 ) ]
  chill1920$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2019 ) ]
  chill1920$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2019 ) ]
  chill1920$GDH[ d ] <- x$GDH[ which( x$End_year == 2019 ) ]
  d <- d+1
  rm( x )
}
rm( jday )

#  calculate daily chill portions for 2020 and write to 
#  dataframe
for ( jday in 1:max( stn_data20$JDay ) ) {  #  ***UPDATE***
  chill1920$iDay[ d ] <- d
  chill1920$JDay[ d ] <- jday
  x <- chilling( hourtemps = stn_data_x,
                 Start_JDay = 304, # non-leap year value for October 31
                 End_JDay = jday )
  chill1920$Chilling_Hours[ d ] <- x$Chilling_Hours[ which( x$End_year == 2020 ) ]
  chill1920$Utah_Model[ d ] <- x$Utah_Model[ which( x$End_year == 2020 ) ]
  chill1920$Chill_portions[ d ] <- x$Chill_portions[ which( x$End_year == 2020 ) ]
  chill1920$GDH[ d ] <- x$GDH[ which( x$End_year == 2020 ) ]
  d <- d+1
  rm( x )
}
rm( d,jday,stn_data_x )

matplot( chill1920$iDay,
         cbind( chill1920$Chilling_Hours / 10, chill1920$Chill_portions ),
         col = c( "black","red" ), type = "l")


#####

#  put all of the chilling portion data together and graph
chill_data <- rbind( chill1617, chill1718, chill1819, chill1920 )


#  FIGURE SETUP

#  Load required libraries.
library( "gridExtra" )


#  To help with setting axis breaks on the figures, create a
#  function that rounds numbers to the nearest specified base.
myround <- function( x,base ) {
  base*round( x/base )
}



#  Create a ggplot object for the chill portions data.
cp <- ggplot( data = chill_data ) +
  
  #  Add chill portions data.
  geom_line( aes( x = iDay,
                  y = Chill_portions,
                  color = Winter ),
             lineend = "round",
             linetype = "solid",
             size = 1.0 ) +
  
  #  Add the graph title.
  ggtitle( "accumulated chill portions" ) +
  
  #  Add the subtitle, x/y axis labels, and caption.
  labs( subtitle = "AZMET Willcox Bench station",
        x = "day",
        #y = "chilling hours",
        caption = paste0( "\ndata source: AZMET (cals.arizona.edu/azmet)",
                          "\nchill calculations: R package 'chillR'" ) ) +
  
  #  Set the line colors.
  scale_color_brewer( palette = "Dark2" ) +
  
  #  Specify the breaks, gridlines, and limits of both plot axes.
  #scale_x_continuous( breaks = c( 1, 15, 31, 45, 62, 76, 93, 107, 121, 135, 152, 166, 182 ),
  #                    labels = c( "Nov 1", "Nov 15", "Dec 1", "Dec 15", "Jan 1", "Jan 15",
  #                                "Feb 1", "Feb 15", "Mar 1", "Mar 15", "Apr 1", "Apr 15", "May 1" ),
  #                    limits = c( ( min( chill_data$iDay ) ),
  #                                ( max( chill_data$iDay+1 ) ) ) ) +
  scale_x_continuous( breaks = c( 1, 31, 62, 93, 121, 152, 182 ),
                      labels = c( "Nov 1", "Dec 1", "Jan 1", "Feb 1", "Mar 1", "Apr 1", "May 1" ),
                      limits = c( ( min( chill_data$iDay ) ),
                                  ( max( chill_data$iDay+1 ) ) ) ) +
  scale_y_continuous( breaks = seq( 0,
                                    myround( ceiling( max( chill_data$Chill_portions ) ),10 ),
                                    by = 10 ),
                      limits = c( 0, max( chill_data$Chill_portions ) ) ) +
  
  #  Specify the ggplot theme, or overall appearance, of the
  #  graph with the font set to 'mono'.
  theme_light( base_family = "mono" ) +
  
  #  Further customize the appearance of the graph.
  theme( axis.line = element_blank(),
         axis.text.x = element_text( color = "gray30",
                                     size = 8 ),
         axis.text.y = element_text( color="gray30",
                                     size = 8 ),
         axis.ticks.length = unit( 0.0,"mm" ),
         axis.title.x = element_text( color = "gray30",
                                      face = "bold",
                                      size = 10 ),
         axis.title.y = element_blank(),
         legend.title = element_text( color = "black",
                                      size = 8 ),
         legend.position = "right",
         panel.border = element_blank(),
         panel.grid.major.x = element_line( color = "gray80",
                                            size = 0.3 ),
         panel.grid.major.y = element_line( color = "gray80",
                                            size = 0.3 ),
         panel.grid.minor.x = element_blank(),
         panel.grid.minor.y = element_blank(),
         plot.caption = element_text( color = "white",
                                      hjust = 0,
                                      size = 6 ),
         plot.margin = unit( c( 2,2,2,5 ),"mm" ),
         plot.subtitle = ( element_text( size = 10 ) ), 
         plot.title = ( element_text( face = "bold",
                                      size = 12 ) ) )

#  Create a ggplot object for the chilling hours data.
ch <- ggplot( data = chill_data ) +
  
  #  Add chilling hours data.
  geom_line( aes( x = iDay,
                  y = Chilling_Hours,
                  color = Winter ),
             lineend = "round",
             linetype = "longdash",
             size = 1.0 ) +
  
  #  Add the graph title.
  ggtitle( "accumulated chilling hours" ) +
  
  #  Add the subtitle, x/y axis labels, and caption.
  labs( subtitle = "AZMET Willcox Bench station",
        x = "day",
        #y = "chilling hours",
        caption = paste0( "\ndata source: AZMET (cals.arizona.edu/azmet)",
                          "\nchill calculations: R package 'chillR'" ) ) +
  
  #  Set the line colors.
  scale_color_brewer( palette = "Dark2" ) +
  
  #  Specify the breaks, gridlines, and limits of both plot axes.
  #scale_x_continuous( breaks = c( 1, 15, 31, 45, 62, 76, 93, 107, 121, 135, 152, 166, 182 ),
  #                    labels = c( "Nov 1", "Nov 15", "Dec 1", "Dec 15", "Jan 1", "Jan 15",
  #                                "Feb 1", "Feb 15", "Mar 1", "Mar 15", "Apr 1", "Apr 15", "May 1" ),
  #                    limits = c( ( min( chill_data$iDay ) ),
  #                                ( max( chill_data$iDay+1 ) ) ) ) +
  scale_x_continuous( breaks = c( 1, 31, 62, 93, 121, 152, 182 ),
                      labels = c( "Nov 1", "Dec 1", "Jan 1", "Feb 1", "Mar 1", "Apr 1", "May 1" ),
                      limits = c( ( min( chill_data$iDay ) ),
                                  ( max( chill_data$iDay+1 ) ) ) ) +
  scale_y_continuous( breaks = seq( 0,
                                    myround( ceiling( max( chill_data$Chilling_Hours ) ),100 ),
                                    by = 200 ),
                      limits = c( 0, max( chill_data$Chilling_Hours ) ) ) +
  
  #  Specify the ggplot theme, or overall appearance, of the
  #  graph with the font set to 'mono'.
  theme_light( base_family = "mono" ) +
  
  #  Further customize the appearance of the graph.
  theme( axis.line = element_blank(),
         axis.text.x = element_text( color = "gray30",
                                     size = 8 ),
         axis.text.y = element_text( color="gray30",
                                     size = 8 ),
         axis.ticks.length = unit( 0.0,"mm" ),
         axis.title.x = element_text( color = "gray30",
                                      face = "bold",
                                      size = 10 ),
         axis.title.y = element_blank(),
         legend.title = element_text( color = "black",
                                      size = 8 ),
         legend.position = "right",
         panel.border = element_blank(),
         panel.grid.major.x = element_line( color = "gray80",
                                            size = 0.3 ),
         panel.grid.major.y = element_line( color = "gray80",
                                            size = 0.3 ),
         panel.grid.minor.x = element_blank(),
         panel.grid.minor.y = element_blank(),
         plot.caption = element_text( color = "black",
                                      hjust = 0,
                                      size = 6 ),
         plot.margin = unit( c( 2,2,2,5 ),"mm" ),
         plot.subtitle = ( element_text( size = 10 ) ), 
         plot.title = ( element_text( face = "bold",
                                      size = 12 ) ) )


#  PLOT THE TWO GRAPHS IN A VERTICAL STACK AND SAVE FIGURE

gA = ggplotGrob( cp )
gB = ggplotGrob( ch )
gB$widths <- gA$widths
grid.arrange( gA, gB )


#  Save the figure as a .png file.
ggsave( "chill-chill-time-series-stack-2020feb.png",
        #plot=multiplot( ch, cp, cols = 1 ),
        plot = grid.arrange( gA, gB ),
        device = "png",
        path = NULL,
        scale = 1,
        width = 6.0,
        height = 8.0,
        units = "in",
        dpi = 600 )


#####




#  Create a ggplot object for the heat units data.
hu <- ggplot( data = aggds_data ) +
  
  #  Add chilling hours data.
  geom_line( aes( x = iDay,
                  y = CGGDs,
                  color = Winter ),
             lineend = "round",
             linetype = "dotted",
             size = 1.0 ) +
  
  #  Add the graph title.
  ggtitle( "accumulated growing degree days" ) +
  
  #  Add the subtitle, x/y axis labels, and caption.
  labs( subtitle = "AZMET Willcox Bench station",
        x = "day",
        #y = "chilling hours",
        caption = paste0( "\ndata source: AZMET (cals.arizona.edu/azmet)",
                          "\nchill calculations: R package 'chillR'" ) ) +
  
  #  Set the line colors.
  scale_color_brewer( palette = "Dark2" ) +
  
  #  Specify the breaks, gridlines, and limits of both plot axes.
  #scale_x_continuous( breaks = c( 1, 15, 31, 45, 62, 76, 93, 107, 121, 135, 152, 166, 182 ),
  #                    labels = c( "Nov 1", "Nov 15", "Dec 1", "Dec 15", "Jan 1", "Jan 15",
  #                                "Feb 1", "Feb 15", "Mar 1", "Mar 15", "Apr 1", "Apr 15", "May 1" ),
  #                    limits = c( ( min( chill_data$iDay ) ),
  #                                ( max( chill_data$iDay+1 ) ) ) ) +
  scale_x_continuous( breaks = c( 1, 31, 62, 93, 121, 152, 182 ),
                      labels = c( "Nov 1", "Dec 1", "Jan 1", "Feb 1", "Mar 1", "Apr 1", "May 1" ),
                      limits = c( ( min( chill_data$iDay ) ),
                                  ( max( chill_data$iDay+1 ) ) ) ) +
  scale_y_continuous( breaks = seq( 0,
                                    myround( ceiling( max( aggds_data$CGGDs, na.rm = TRUE ) ),100 ),
                                    by = 100 ),
                      limits = c( 0, max( aggds_data$CGGDs, na.rm = TRUE ) ) ) +
  
  #  Specify the ggplot theme, or overall appearance, of the
  #  graph with the font set to 'mono'.
  theme_light( base_family = "mono" ) +
  
  #  Further customize the appearance of the graph.
  theme( axis.line = element_blank(),
         axis.text.x = element_text( color = "gray30",
                                     size = 8 ),
         axis.text.y = element_text( color="gray30",
                                     size = 8 ),
         axis.ticks.length = unit( 0.0,"mm" ),
         axis.title.x = element_text( color = "gray30",
                                      face = "bold",
                                      size = 10 ),
         axis.title.y = element_blank(),
         legend.title = element_text( color = "black",
                                      size = 8 ),
         legend.position = "right",
         panel.border = element_blank(),
         panel.grid.major.x = element_line( color = "gray80",
                                            size = 0.3 ),
         panel.grid.major.y = element_line( color = "gray80",
                                            size = 0.3 ),
         panel.grid.minor.x = element_blank(),
         panel.grid.minor.y = element_blank(),
         plot.caption = element_text( color = "black",
                                      hjust = 0,
                                      size = 6 ),
         plot.margin = unit( c( 2,2,2,5 ),"mm" ),
         plot.subtitle = ( element_text( size = 10 ) ), 
         plot.title = ( element_text( face = "bold",
                                      size = 12 ) ) )


#  PLOT THE TWO GRAPHS IN A VERTICAL STACK AND SAVE FIGURE

gA = ggplotGrob( cp )
gC = ggplotGrob( hu )
gC$widths <- gC$widths
grid.arrange( gA, gC )


#  Save the figure as a .png file.
ggsave( "chill-heat-time-series-stack-2020feb.png",
        #plot=multiplot( ch, cp, cols = 1 ),
        plot = grid.arrange( gA, gC ),
        device = "png",
        path = NULL,
        scale = 1,
        width = 6.0,
        height = 8.0,
        units = "in",
        dpi = 600 )





#####


