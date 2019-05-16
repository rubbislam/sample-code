# install.packages("devtools")
# devtools::install_github("ropensci/RSelenium")
library(RSelenium)
# tutorial:https://blog.gtwang.org/r/rselenium-r-selenium-browser-web-scraping-tutorial/
# shell('docker pull selenium/standalone-chrome')
# shell('docker run -d -p 4445:4444 selenium/standalone-chrome')
# remDr <- remoteDriver(remoteServerAddr="localhost",
#                       port=4444,
#                       browser="chrome")
# java -Dwebdriver.chrome.driver=D:\chromedriver.exe -Dwebdriver.gecko.driver=D:\geckodriver.exe -jar selenium-server-standalone-3.13.0.jar
# remDr$open()
# 
start="2016-01-01"
end="2018-12-31"
ICAO="TJBQ"
country.abbr="pr"
city="hatillo"


# climate_scraper <- function(start="2016-01-01", end=Sys.Date(), ICAO="FNLU", country.abbr="ao", city="Luanda") {
  if (is.null(ICAO)) {
    stop("Please enter location in ICAO code.")
  }
  if (is.null(country.abbr)) {
    stop("Please enter Country 2-digit abbreviation.")
  }
  if (is.null(city)) {
    stop("Please enter city name.")
  }
  # if (!require(RCurl)) {
  #   install.packages("RCurl", repos = "http://cran.us.r-project.org")
  # }
  if (!require(rvest)) {
    install.packages("rvest", repos = "http://cran.us.r-project.org")
  }
  
  if (!require(xml2)) {
    install.packages("xml2", repos = "http://cran.us.r-project.org")
  }
  if (!require(doParallel)) {
    install.packages("doParallel", repos = "http://cran.us.r-project.org")
  }
  if (!require(foreach)) {
    install.packages("foreach", repos = "http://cran.us.r-project.org")
  }
  if (!require(RSelenium)) {
    install.packages("foreach", repos = "http://cran.us.r-project.org")
  }
  
  
  # library(RCurl)
  library(rvest)
  library(xml2)
  library(doParallel)
  library(foreach)
  library(RSelenium)
  
  dates <- format(seq(as.Date(start), as.Date(end), by="months"),"%Y-%m-%d")
  # dates <- c( "2016-11-01", "2016-12-01", "2017-09-01", "2017-10-01", "2017-12-01", "2018-09-01")
  (cl <- (detectCores() - 1) %>%  makeCluster) %>% registerDoParallel
  # open a remoteDriver for each node on the cluster
  clusterEvalQ(cl, {
    library(RSelenium)
    remDr <- remoteDriver(
      remoteServerAddr="localhost",
      port=4444,
      browser="chrome"
    )
    remDr$open()
  })
  
  dat <- 
    foreach(i=1:length(dates), .combine=rbind, .packages=c("rvest","xml2","RSelenium")) %dopar% {
              web <- paste0("https://www.wunderground.com/history/monthly/",
                            country.abbr, "/", tolower(city), "/", ICAO, "/date/", substr(dates[i], 1, 7))
              
              # remDr <- remoteDriver(remoteServerAddr="localhost",
              #                       port=4444+i-1,
              #                       browser="chrome")
              # remDr$open()
              
              remDr$navigate(web)
              Sys.sleep(10)
              output <- remDr$getPageSource()
              # remDr$close()
              doc <- read_html(output[[1]])
              tabs <- doc %>% html_nodes(xpath="//table[@class='days']//tbody//tr") %>% html_text()
              # rm(doc, tabs)
              
              if (length(tabs)>0) {
                tt <- strsplit(tabs, split="\n ")
                tt <- lapply(tt, function(x) {
                  y <- gsub(" ","",x)
                  return(y[y!=""])
                })
                
                
                dd <- tt[2:length(tt)]
                
                # print(dates[i])
                # print(length(which(tt[[1]]=="Max")))
                # print(length(which(tt[[1]]=="Avg")))
                # print(length(which(tt[[1]]=="Min")))
                
                
                no.row <- which(tt[[1]]=="Max")[[1]]-1
                mat <- matrix(NA, ncol=19, nrow=no.row)
                
                ind <- seq(1, length(dd), no.row)
                col.ind <- 1
                for (f in ind) {
                  mat[1:no.row,col.ind:(col.ind+length(dd[[f]])-1)] <- matrix(unlist(dd[f:(f+no.row-1)]),nrow=no.row,byrow=TRUE)
                  col.ind <- col.ind + length(dd[[f]])
                }
                
                df <- as.data.frame(mat, stringsAsFactors=FALSE)
                df <- df[-1,]
                df$V1 <- paste0(substr(dates[i], 1, 7), "-", ifelse(nchar(df$V1)>1, df$V1, paste0(0,df$V1)))
                dat <- df
              } else {
                dat <- NULL
              }
              # dat2 <- rbind(dat, dat2)
    }
  
  ind <- which(substr(dates,1,7) %in% unique(substr(dat$V1,1,7)) == FALSE)
  dates[ind]
  
  clusterEvalQ(cl, {
    remDr$close()
  })
  stopImplicitCluster()
  # stop Selenium Server
  # selServ$stop()
  colnames(dat) <- c("date", "temp.max", "temp.avg", "temp.min", "dewpoint.max",
                     "dewpoint.avg", "dewpoint.min", "humidity.max", "humidity.avg",
                     "humidity.min", "windspeed.max", "windspeed.avg",
                     "windspeed.min", "pressure.max", "pressure.avg",
                     "pressure.min", "precipation.max", "precipation.avg",
                     "precipation.max")
  dat$date <- as.Date(dat$date)
  dat <- dat[order(dat$date, decreasing=FALSE),]
  # return(dat)
# }

# dat <- rbind(dat, dat3)

write.csv(dat, "Arecibo_201601-201812.csv", quote=FALSE, row.names=FALSE)

# test <- climate_scraper(start="2000-01-01")