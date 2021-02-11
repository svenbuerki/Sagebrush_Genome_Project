#1. Load library
require(googlesheets4)

#2. Authenticate your Google account

#Get url of Google spreadsheet
url <- "https://docs.google.com/spreadsheets/d/1ZVrQ8ZiX6d5ghVKgoY6jI1lXTErAsddI15k_sGXDWJ4/edit?usp=sharing"

#Obtain meta-data on document
metaDoc <- gs4_get(url)

#List names of sheets in url
sheetsID <- metaDoc$sheets$name

curDat <- read_sheet(url, sheet=sheetsID[1])
