library(ggtree)
library(dplyr)
library(grepl)

# Read in .txt
human <- read.table('alltime_human.txt') |>
  filter("2013", alltime_human.txt$text)
