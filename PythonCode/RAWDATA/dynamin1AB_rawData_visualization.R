library(ggplot2)
library(readr)
library(dplyr)

raw <- read_csv("/Users/cmdb/Desktop/JOHNSON_rotation/dynamin_model/PythonCode/RAWDATA/Dynamin1AB_kinetics_RawData.csv") %>%
  dplyr::rename(AverageFluorescence = Aerage, Time_s = `...1`) %>%
  dplyr::select(AverageFluorescence, Time_s) %>%
  dplyr::filter(Time_s >= -20 & Time_s <= 20)


ggplot(raw, aes(x = Time_s, y = AverageFluorescence)) +
  geom_point(color = 'green') +
  geom_vline(xintercept = c(-4.00, 4.00), linetype = 'dashed', color = 'red') +
  labs(title = "Dynamin dynamics during Endocytosis",
       x = "Time (s)",
       y = "Average Fluorescence (A.U.)")


print(colnames(raw))
