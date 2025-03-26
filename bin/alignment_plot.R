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
# df <- df %>%
#   mutate(percent_classified_reads = round(percent_classified_reads, 2))

# print(df)
# quit()

# Pull the top n scores by category - will give us 2n rows because there are two algorithm categories
if ( nrow(df) < 10 ) {

  slice_num <- nrow(df)

} else {

  slice_num <- 10

}

top_hits <- df %>%
  slice_max(n = slice_num, order_by = percent_classified_reads)

# Randomizing the order of the rows - this is so when we go to create this bar chart, the order of species isn't always highest to lowest or vice-versa - makes it look better in my opinion
#top_hits <- top_hits[sample(nrow(top_hits)), ]

df_temp <- top_hits %>%
  arrange(desc(percent_classified_reads))

# curating list of species to be used as level parameter in ggplot
categories <- df_temp[, 1, drop = FALSE]
categories <- categories$Species.name
categories <- unique(categories)

# this will put "unclassified" data at the end of the array
if ("unclassified" %in% categories) {
  categories <- c(categories[categories != "unclassified"], "unclassified")
}


top_hits_mod <- data.frame()

for (i in categories) {
  
  for (x in 1:nrow(top_hits)) {
    
    if (top_hits[x, "algorithm"] == "alignment" & top_hits[x, "Species.name"] == i) {
      
      top_hits_mod <- rbind(top_hits[x, ], top_hits_mod) 
      
    }
    
  }
  
}

### This block will calculate the number and percent of "other" species classified by the "alignment" algorithm and add them as a new row
total_percent_alignment <- sum(top_hits_mod$percent_classified_reads[top_hits_mod$algorithm == "alignment"])

if (total_percent_alignment < 99.99) {

  percent_other_alignment <- 100 - total_percent_alignment

  new_row_other <- data.frame(
  Species.name = "other",
    num_reads = round(sum(top_hits_mod$num_reads[top_hits_mod$algorithm == "alignment"]) / (total_percent_alignment / 100)) - sum(top_hits_mod$num_reads[top_hits_mod$algorithm == "alignment"]),
    percent_classified_reads = percent_other_alignment,
    algorithm = "alignment"
  )
  top_hits_mod <- rbind(top_hits_mod, new_row_other)

  categories <- c(categories, "other")

}
###

# This for loop will calculate the position of each data label for each row
top_hits_mod$label_position <- NA

alignment_count <- 0

for (i in 1:nrow(top_hits_mod)) {
    
  if (alignment_count == 0) {
    
    top_hits_mod[i, "label_position"] <- top_hits_mod[i, "percent_classified_reads"] - 0.5 * top_hits_mod[i, "percent_classified_reads"]
    alignment_sum = top_hits_mod[i, "percent_classified_reads"]
    
  } else {
    
    top_hits_mod[i, "label_position"] <- alignment_sum + top_hits_mod[i, "percent_classified_reads"] - 0.5 * top_hits_mod[i, "percent_classified_reads"]
    alignment_sum = alignment_sum + top_hits_mod[i, "percent_classified_reads"]
    
  }

  alignment_count <- i   

}


# Adding column that will define the data label text displayed on the graph
top_hits_mod$label_percent <- NA
top_hits_mod$label_percent = paste0(sprintf("%.1f", top_hits_mod$percent_classified_reads), "%") # round the percentages to the first decimal place so they fit better in the bar plot

if ("unclassified" %in% categories) {

  # because the previous rounding can cause 99.99% to turn into 100.0% (which can be the case for unclassified reads sometimes) we are going to round 'unclassified' label percentages that equal 100.0 to two decimal places to give us a less misleading value
  if (top_hits_mod$label_percent[top_hits_mod$Species.name == "unclassified"] == "100.0%" & top_hits_mod$percent_classified_reads[top_hits_mod$Species.name == "unclassified"] < 99.999999999) {

    top_hits_mod$label_percent[top_hits_mod$Species.name == "unclassified"] <- paste0(sprintf("%.2f", top_hits_mod$percent_classified_reads), "%")

  }

}

colors <- c(
  "#FF7666", "#5FD1B7",  # Lighter Red → Teal  
  "#FFA463", "#6FA9E3",  # Lighter Orange → Blue  
  "#B885E6", "#F8DB4A",  # Lighter Purple → Yellow  
  "#6FDA8B", "#EA7BA5",  # Lighter Green → Pink  
  "#8899B3", "#F4B843",  # Soft Slate Blue → Goldenrod  
  "#5D7A99", "#E08E5C",  # Soft Blue Gray → Burnt Orange  
  "#C7A0E6", "#E57373",  # Medium Purple → Soft Crimson  
  "#69D4C3", "#D49EE2",  # Strong Teal → Soft Violet  
  "#E57373", "#81C784",  # Soft Red → Fresh Green  
  "#FF7666", "#64B5F6",  # Lighter Bold Red → Vibrant Blue  
  "#AF7AC5", "#F9E79F",  # Soft Purple → Warm Yellow  
  "#566F8E", "#C2C2C2"   # Deep Gray-Blue → Light Gray  
)

# counting length of category array
# category_len <- length(categories)

if ("unclassified" %in% categories && "other" %in% categories) {

  colors[2] <- "#D3D3D3"
  colors[1] <-  "#A9A9A9"

} else if ("unclassified" %in% categories || "other" %in% categories) {

  colors[1] <- "#D3D3D3"

}

# converting percent values below 1% to 1% so they are more easily visible on the graph
top_hits_mod <- top_hits_mod %>%
  mutate(percent_classified_reads = ifelse(percent_classified_reads < 1, 1, percent_classified_reads))

plot_1 <- ggplot(top_hits_mod, aes(x = algorithm, fill = factor(Species.name, levels = rev(factor(categories))), y = percent_classified_reads)) +
  geom_col(width = 0.5) +
  theme(legend.title = element_blank(), legend.position = "top") +
  theme(axis.text.x = element_blank()) +
  geom_text(aes(label = label_percent), position = position_stack(vjust = 0.5), size = 2.5) +  
  scale_fill_manual(values=colors) +
  coord_flip()

ggsave("plot.png", width = 14, height = 6, dpi = 150)
