library(ggplot2)
library(readr)
library(dplyr)

raw <- read_csv("/Users/cmdb/Desktop/JOHNSON_rotation/dynamin_model/PythonCode/RAWDATA/Dynamin1AB_kinetics_RawData.csv") %>%
  dplyr::rename(AverageFluorescence = Aerage, Time_s = `...1`) %>%
  dplyr::select(AverageFluorescence, Time_s) %>%
  dplyr::filter(Time_s >= -10 & Time_s <= 10)


final <- ggplot(raw, aes(x = Time_s, y = AverageFluorescence)) +
  geom_point(color = 'green') + 
  geom_line(color = 'green') +
  geom_vline(xintercept = c(-4.00, 4.00), linetype = 'dashed', color = 'red') +
  labs(title = "Dynamin During Endocytosis",
       x = "Time (s)",
       y = "Average Fluorescence (A.U.)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 10)
  )

final


#ggsave(final, "/Users/cmdb/Desktop/JOHNSON_rotation/dynamin_model/PythonCode/FIGURES/.png", dpi=900)