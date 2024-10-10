library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

headers <- c("item","bepoch",	"ecode","label","onset","diff","dura","b_flags","a_flags","enable", "bin")

load_long_logfile <- function(path){
  require(tidyr)
  return_df <- data.frame()
  files <- list.files(path,full.names = T)
    for (i in 1 :length(files)){
      subno <- str_extract(files[i], "[0-9][0-9][0-9]")
      file <- read.table(files[i], header = T)
      file$subno <- subno
      return_df <- rbind(file, return_df)
    }
  return_df <- return_df %>%
    arrange(subno,block) %>%
    group_by(subno) %>%
    mutate(order = seq(1,320,1)) %>%
    filter(!subno == "001")
  return(return_df %>% arrange(subno,order,block))
}

load_eventcode_files <- function(path, headers){
  return_list <- list()
  dirs <- list.files(path = path)
  files <- list.files(path = paste(path,dirs[-1],sep = ""),pattern = "ppevent.txt",full.names = TRUE)
  list_files <- map(files, function(x){
    d <- read.table(x)
    colnames(d) <- headers
    d$subno <- str_extract(x, "[0-9][0-9][0-9]")
    d$subno <- as.factor(d$subno)
    return(d)
  })
  return(list_files)
}
  
  
event_files <- load_eventcode_files("data/", headers = headers)
  
logfiles <- load_long_logfile("Logfiles")


change_event_codes <- function(events, logfiles){
  events_stim <- map(events, function(x){
    return(x %>% filter(ecode == 42 | ecode == 41))
  })
  logfiles$new_codes <- ifelse(logfiles$cond == 1, 61,
                               ifelse(logfiles$cond == 2, 71,
                                      ifelse(logfiles$cond == 0, 81,NA)))
  
      for (i in 1:nrow(logfiles)) {
        logfiles[nrow(logfiles) + 1,] <- list(logfiles$order[i],
                                              logfiles$cond[i],
                                              logfiles$time[i],
                                              logfiles$button[i],
                                              logfiles$rt[i],
                                              logfiles$correct[i],
                                              logfiles$block[i],
                                              logfiles$counterbalance[i],
                                              logfiles$subno[i],
                                              logfiles$new_codes[i] + 1)
      }
  logfiles <- arrange(logfiles, subno,order)
  long_eventcodes <- map_dfr(events_stim,function(x){
      return(x)
    })$item
  
  bind_events <- cbind(logfiles,new_order = long_eventcodes)
  new_events <- map(events, function(x){
    final <- x %>%
      left_join(select(bind_events, new_codes,new_order), by = c("subno","item" = "new_order")) %>%
      mutate(ecode = coalesce(new_codes, ecode)) %>%
      select(-new_codes, -subno)
  })
  
  return(list(logfiles,new_events))
}

x <- change_event_codes(event_files, logfiles)


subs <- unique(logfiles$subno)

map2(x[[2]],subs, function(x,y){
  return(write.table(x, file = paste("new_codes/",y,"_codes.txt", sep = ""),sep = "\t"))
  })

