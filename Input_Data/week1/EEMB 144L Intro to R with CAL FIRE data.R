##### EEMB 144L Intro to R with CAL FIRE data #####
#Kerri Luttrell
#10/5/2020

#installing packages

#collection of useful packages
install.packages("tidyverse")
library(tidyverse)
#way to import data
install.packages("readxl")
library(readxl)


#### Load Data ####



# .xlsx files can have multiple pages

excel_sheets("2013_2019_CALFIRE_Redbook.xlsx")

calfire.metadata <- read_excel("2013_2019_CALFIRE_Redbook.xlsx", sheet="Metadata")

calfire.data <- read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet=2)

