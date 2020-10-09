##### EEMB 144L Intro to R with CAL FIRE data #####
#Kerri Luttrell
#10/5/2020

#installing packages

#collection of useful packages

library(tidyverse)

#way to import data

library(readxl)



#### Load Data ####



# .xlsx files can have multiple pages

excel_sheets("2013_2019_CALFIRE_Redbook.xlsx")
#we can store each sheet as a separate data frame
calfire.metadata <- read_excel("2013_2019_CALFIRE_Redbook.xlsx", sheet="Metadata")

calfire.data <- read_excel("2013_2019_CALFIRE_Redbook.xlsx", sheet="Data")

#### Intial Data Exploration ####

names(calfire.data)
# [1] "Incident"             "County_Unit"         
#[3] "Fire_Name"            "Start_Date"          
#[5] "Controlled_Date"      "Origin_DPA_Agency"   
#[7] "Total_Acres_Burned"   "Veg_Type"            
#[9] "Cause"                "Structures_Destroyed"
#[11] "Structures_Damaged"   "Fire_Fatalities"     
#[13] "Civil_Fatalities"

dim(calfire.data)
#378 rows (observations) and 13 columns(variables)

class(calfire.data)
#tells us what the object is (tbl_df=table dataframe)

head(calfire.data)
#shows first six lines of set
tail(calfire.data)
#last six lines

#single columns in dataset can be referred to using '$'
county<- calfire.data$County_Unit

names(calfire.data)
#look at variables and ask question, what was maximium amount of acres burned
max_acres <- max(calfire.data$Total_Acres_Burned, na.rm=TRUE)
#na.rm says if there are NA's in the column, dont pay attention to them when running analyses
# without na.rm, below command results in NA as max of structures destroyed
max(calfire.data$Structures_Destroyed, na.rm=TRUE)

##### Basic Data Wrangling (dplyr fucntions) #####

#subsetting data to look at only a few variables

df1 <- select(calfire.data, County_Unit:Controlled_Date, Total_Acres_Burned:Civil_Fatalities)
#maintains specified order
view(df1)

df2 <- filter(df1, County_Unit
              %in% c('SANTA BARBARA',"VENTURA","LOS ANGELES", 'ORANGE', 'SAN DIEGO') & Total_Acres_Burned >= 500 | Fire_Name=="THOMAS" )
df3 <- arrange(df2, desc(Start_Date), (Total_Acres_Burned))

df4 <- mutate_at(df3, vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0)
#change these specific columns in this way

#to calculate something and add a new column
df5 <- mutate(df4, Fatalies
              = Fire_Fatalities + Civil_Fatalities)

load(lubridate)

df6<- mutate(df5,
             interv = interval
             (Start_Date, Controlled_Date),
             dur= as.duration(interv),
             days= as.numeric(dur, "days"))

socal.fires <- calfire.data %>% 
  select(County_Unit:Controlled_Date, Total_Acres_Burned:Civil_Fatalities) %>% 
  filter(County_Unit
         %in% c('SANTA BARBARA',"VENTURA","LOS ANGELES", 'ORANGE', 'SAN DIEGO') & Total_Acres_Burned >= 500 | Fire_Name=="THOMAS" ) %>% 
  arrange(desc(Start_Date), (Total_Acres_Burned)) %>% 
  mutate_at(vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0) %>% 
  mutate(Fatalies
         = Fire_Fatalities + Civil_Fatalities) %>% 
  mutate(interv = interval
         (Start_Date, Controlled_Date),
         dur= as.duration(interv),
         days= as.numeric(dur, "days"))
ggplot(socal.fires, aes(x= Start_Date, y= Total_Acres_Burned))+
  geom_point(aes(color= County_Unit))+
   ggtitle("CA South Coast Major Fires \n2014 - 2018")+
  labs(x="", y="Total Acres Burned", color = "County")+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor= element_blank())+
  facet_grid(rows= "County_Unit", scales = "free")
  

plot.data <- socal.fires  %>% 
  rename(county = County_Unit,
         acres= Total_Acres_Burned,
         start= Start_Date,
         end= Controlled_Date) %>% 
  mutate (year = year(start),
         county= ifelse(county== "VENTURA/SANTA BARBARA", "VENTURA",
                        county)) 
incidents <- plot.data %>% 
  group_by(county,year) %>% 
  mutate(ave_acres = mean(acres, na.rm =T)) %>% 
  ungroup()

incidents.plot <- incidents %>% 
  ggplot(aes(x= year, y=n )) +
  geom_point(aes(color =county)) +
  geom_line(aes(color = county)) +
  labs(title = "CA South Coast Major Fire Incidents \n 2013-2018", x = "", y ="Incidents", color="County") +
  theme_bw()+
  facet_grid(rows= "county", scales= "free")




 
