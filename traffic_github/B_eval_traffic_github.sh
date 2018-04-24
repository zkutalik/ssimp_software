#!/usr/bin/env Rscript




## ////////////////////
## (B) using the functions
## ////////////////////

args <- commandArgs(trailingOnly = TRUE)
folder.args <- (as.character(args[1])) ## your ssimp line

library(tidyverse) ## yep - this one you need: to read and write data 

## bundle import all data
dir(folder.args)
raw <- map(paste0(folder.args,"/", dir(folder.args)), function(x) read_csv(x) %>% mutate(filename=x))

    
## bind_row + unique
dat <- bind_rows(raw) %>% unique() %>% rename(uniquevisitors_cloners = `unique_visitors/cloners`)

## append filename as column > extract traffic, referrer, clone later
dat <- dat %>% mutate(type = case_when(
  grepl("traffic", filename) ~ "traffic",
  grepl("clone", filename)  ~ "clone",
  grepl("referrer", filename)  ~ "referrer"))


## clones
ggplot(data = dat %>% filter(type == "clone")) + facet_wrap(~type) + geom_path(aes(date, clones), color = "#4daf4a") + geom_path(aes(date, uniquevisitors_cloners), color = "#377eb8")

## traffic
ggplot(data = dat %>% filter(type == "traffic")) + facet_wrap(~type) + geom_path(aes(date, views), color = "#4daf4a") + geom_path(aes(date, uniquevisitors_cloners), color = "#377eb8")

## unique visitors
qplot(date, uniquevisitors_cloners, data = dat %>% filter(type %in% c("traffic","clone"))) + facet_wrap(~type) + geom_path()

## sum(views) is ok, but unique visitors should be ...
dat %>% group_by(site) %>% summarize(sum(views))
dat %>% group_by(site) %>% summarize(sum(uniquevisitors_cloners ))

## run this file with
## chmod a+x ssimp_chunks.sh
## ./B_traffic_github.sh csv
