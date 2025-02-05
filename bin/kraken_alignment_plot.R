library(ggplot2)
library(dplyr)

theme_set(theme_bw())

args <- commandArgs(trailingOnly = TRUE)

csv <- args[1]
workdir <- args[2]

setwd(workdir)

df <- read.csv(csv)

# trim any leading spaces from 3rd column in dataframe
df$percent_classified_reads <- sub("^\\s+", "", df$percent_classified_reads)

#convert 3rd column to double
df$percent_classified_reads <- as.double(df$percent_classified_reads)

# Round the 'Value' column to 2 decimal places
df <- df %>%
  mutate(percent_classified_reads = round(percent_classified_reads, 2))

# Pull the top n scores by category - will give us 2n rows because there are two algorithm categories
top_hits <- df %>%
  group_by(algorithm) %>%
  slice_max(order_by = percent_classified_reads, n = 10)

# Randomizing the order of the rows - this is so when we go to create this bar chart, the order of species isn't always highest to lowest or vice-versa - makes it look better in my opinion
top_hits <- top_hits[sample(nrow(top_hits)), ]

categories <- top_hits[, 1, drop = FALSE]
categories <- categories$Species.name
categories <- unique(categories)
categories <- rev(categories)


# Creating new df with kraken rows in top half in the same order as the "categories" var
top_hits_mod <- data.frame()

for (i in categories) {
  
  for (x in 1:nrow(top_hits)) {
   
    if (top_hits[x, "algorithm"] == "kraken" & top_hits[x, "Species.name"] == i) {
      
      top_hits_mod <- rbind(top_hits[x, ], top_hits_mod) 
      
    }
     
  }
  
}

# Appending alignment rows to bottom half of top_hits_mod in the same order as the "categories" var - this will give us ordered kraken rows at the top of the df and ordred alignment rows at the bottom half of the df
for (i in categories) {
  
  for (x in 1:nrow(top_hits)) {
    
    if (top_hits[x, "algorithm"] == "alignment" & top_hits[x, "Species.name"] == i) {
      
      top_hits_mod <- rbind(top_hits[x, ], top_hits_mod) 
      
    }
    
  }
  
}

#top_hits_mod <- top_hits_mod[order(top_hits_mod$Species.name), ]

# This for loop will calculate the position of each data label for each row
top_hits_mod$label_position <- NA

kraken_count <- 0
alignment_count <- 0

for (i in 1:nrow(top_hits_mod)) {
  
  if (top_hits_mod[i, "algorithm"] == "kraken") {
    
    if (kraken_count == 0) {
      
      top_hits_mod[i, "label_position"] <- top_hits_mod[i, "percent_classified_reads"] - 0.5 * top_hits_mod[i, "percent_classified_reads"]
      kraken_sum = top_hits_mod[i, "percent_classified_reads"]
      
    } else {
      
      top_hits_mod[i, "label_position"] <- kraken_sum + top_hits_mod[i, "percent_classified_reads"] - 0.5 * top_hits_mod[i, "percent_classified_reads"]
      kraken_sum = kraken_sum + top_hits_mod[i, "percent_classified_reads"]
      
    }
    
    kraken_count <- i
    
    
  } else {
    
    if (alignment_count == 0) {
      
      top_hits_mod[i, "label_position"] <- top_hits_mod[i, "percent_classified_reads"] - 0.5 * top_hits_mod[i, "percent_classified_reads"]
      alignment_sum = top_hits_mod[i, "percent_classified_reads"]
      
    } else {
      
      top_hits_mod[i, "label_position"] <- alignment_sum + top_hits_mod[i, "percent_classified_reads"] - 0.5 * top_hits_mod[i, "percent_classified_reads"]
      alignment_sum = alignment_sum + top_hits_mod[i, "percent_classified_reads"]
      
    }
    
    alignment_count <- i   
    
  }
  
}


# Adding column that will define the data label text displayed on the graph
top_hits_mod$label_percent <- NA
top_hits_mod$label_percent = paste0(sprintf("%.1f", top_hits_mod$percent_classified_reads), "%")


colors <- c(
  "#FF6666", "#33CC33", "#6666FF", "#CCCC33", # Dark Red, Dark Green, Dark Blue, Dark Yellow
  "#CC33CC", "#33CCCC", "#CC6633", "#CC9933", # Dark Magenta, Dark Cyan, Burnt Orange, Dark Gold
  "#66CC33", "#9933CC", "#33CC99", "#3366CC", # Forest Green, Deep Purple, Sea Green, Royal Blue
  "#CC66CC", "#3399CC", "#CC9966", "#99CC66", # Plum, Sky Blue, Tan, Moss Green
  "#CCCC66", "#66CCCC", "#CC3399", "#6699CC", # Olive, Teal, Raspberry, Steel Blue
  "#CC99CC", "#99CCCC", "#CC9933", "#CC3399"  # Mauve, Soft Aqua, Goldenrod, Fuchsia
)

plot_1 <- ggplot(top_hits_mod, aes(x = algorithm, fill = factor(Species.name, levels = categories), y = percent_classified_reads)) +
  geom_col(width = 0.6) +
  theme(legend.title = element_blank()) +
  geom_text(aes(y = label_position, label = label_percent), size = 3) +
  coord_flip() +
  scale_fill_manual(values=colors, 
                    labels=categories)

ggsave("plot.png", width = 1920, height = 1080, dpi = 200, units = "px")
