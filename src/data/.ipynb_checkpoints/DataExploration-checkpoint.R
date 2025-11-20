library(readr)
library(rlang)
library(skimr)
library(data.table)



# Read datasets --------------------------------------------------------
directory_path = "/home/jsancheg/Descargas/CI5plus_Summary_Legacy"
file1_path = file.path(directory_path,"Europe.csv")
file2_path = file.path(directory_path,"Europe_Cases.csv")
file3_path = file.path(directory_path,"Europe_Pops.csv")

Europe <- fread(file1_path)
Europe_cases <- fread(file2_path)
Europe_Pops <- fread(file3_path)
str(Europe)
skim(Europe)

setDT(Europe)
setDT(Europe_cases)
setDT(Europe_Pops)


# Explore datasets --------------------------------------------------------

# Colnames of the dataset Europe
colnames(Europe)
head(Europe)

# Colnames of the dataset Europe cases
colnames(Europe_cases)
head(Europe_cases)


# Create a cancer dictionary ----------------------------------------------

cancer_dictionary <- data.frame( CANCER = c(1:29),
cancer_classification = c("All cancers but non-melanoma skin (C00-96, but C44)",
"Lip, oral cavity and pharynx (C00-14)",
"Oesophagus (C15)",
"Stomach (C16)",
"Colon (C18)",
"Rectum and anus (C19-21)",
"Liver and intrahepatic bile ducts (C22)",
"Gallbladder and extrahepatic ducts (C23-24)",
"Pancreas (C25)",
"Larynx (C32)",
"Trachea, bronchus and lung (C33-34)",
"Bone (C40-41)",
"Melanoma of skin (C43)",
"Connective and soft tissue (C47+C49)",
"Breast1 (C50)",
"Cervix uteri (C53)",
"Corpus uteri (C54)",
"Ovary (C56)",
"Prostate (C61)",
"Testis (C62)",
"Kidney and renal pelvis (C64-65)",
"Bladder (C67)",
"Eye (C69)",
"Brain, central nervous system (C70-72)",
"Thyroid (C73)",
"Hodgkin lymphoma (C81)",
"Non-Hodgkin lymphoma2 (C82-86, C96)",
"Multiple myeloma and immunoproliferative diseases (C88+C90)",
"Leukaemia (C91-95)"
),
additional_description = c("Includes HIV disease resuting in cancer", 
rep("",28))
)

setDT(cancer_dictionary)
head(cancer_dictionary)


# Transform datasets ------------------------------------------------------

# Add a column to represent sex as factor

Europe_cases[,SEX_CHARACTER := factor(SEX, levels = c(1,2), labels = c("M","F"))]
Europe_cases[, unique(CANCER)]
colnames(Europe_cases)

colnames(cancer_dictionary)
class(Europe_cases$SEX)
Europe_cases$SEX_CHARACTER


# Merge datasets to add the cancer description

aux_europe_cases_cancer_description <- merge(
  Europe_cases,
  cancer_dictionary,
  by = "CANCER",
  all.x = TRUE
)


head(Europe_cases_cancer_description)
TOTAL_N <- "TOTAL_N"
exclude_cols <- c("SEX", "SEX_CHARACTER",
                  "YEAR", "CANCER", "cancer_classification",
                  "additional_description", "TOTAL")

Europe_cases_cancer_description <- aux_europe_cases_cancer_description %>%
  select(-REGISTRY) %>%
  group_by(CANCER, cancer_classification, 
           additional_description,SEX, 
           SEX_CHARACTER, YEAR) %>%
  summarise(across(where(is.numeric), sum),
            .groups = "drop")

Europe_cases_cancer_description = data.table(Europe_cases_cancer_description)

colnames(Europe_cases_cancer_description)
head(Europe_cases_cancer_description)
nrow(aux_europe_cases_cancer_description)
nrow(Europe_cases_cancer_description)

# Check all the numeric columns sums the quantity in column TOTAL
Europe_cases_cancer_description[,sum(N_UNK)]
Europe_cases_cancer_description[, sum(TOTAL)]
Europe_cases_cancer_description[,(TOTAL_N) := rowSums(.SD), 
                                .SDcols = !exclude_cols]
Europe_cases_cancer_description[,list(SUM_TOTAL = sum(TOTAL), 
                                      SUM_TOTAL_N = sum(TOTAL_N) ) ]

# Check the period of time that the data set reports cancer cases
Europe_cases_cancer_description[,list(MIN_YEAR = min(YEAR), 
                                      MAX_YEAR = max(YEAR))]

# Number of cancer cases by year
Europe_cases_cancer_description[, list(TOTAL_CASES = sum(TOTAL) ), 
                                by = YEAR][order(YEAR)]

# Number of cancer cases by cancer type and year
cancer_cases_per_year <- Europe_cases_cancer_description[, 
                                list(TOTAL_CASES = sum(TOTAL)), 
          by = list(cancer_classification,YEAR)][order(cancer_classification,YEAR)]

colnames(Europe_cases_cancer_description)

# Cancer cases by sex and year
cancer_cases_per_sex_year <- Europe_cases_cancer_description[, list(TOTAL_CASES = sum(TOTAL)), 
          by = list(cancer_classification,SEX_CHARACTER,YEAR)][order(cancer_classification,SEX_CHARACTER,YEAR)]

head(cancer_cases_per_year,100)
table(cancer_cases_per_year$YEAR)
summary(cancer_cases_per_year$TOTAL_CASES)
sum(cancer_cases_per_year$TOTAL_CASES==0)

head(cancer_cases_per_sex_year,100)
table(cancer_cases_per_sex_year$YEAR)

skim(cancer_cases_per_year)

skim(cancer_cases_per_sex_year)
colnames(cancer_cases_per_sex_year)

cancer_cases_per_year[, list(INCIDENCES = sum(TOTAL_CASES) ) ]

cancer_cases_per_sex_year[, list(INCIDENCES = sum(TOTAL_CASES) ) ]

plot(density(cancer_cases_per_sex_year$TOTAL_CASES, bw = "sj"))

# cancer cases that are not prevalent in men
zero_cancer_cases_sex1 <- cancer_cases_per_sex_year[SEX_CHARACTER=="M" & TOTAL_CASES ==0]
zero_cancer_cases_sex2 <- cancer_cases_per_sex_year[SEX_CHARACTER=="F" & TOTAL_CASES ==0]

length(zero_cancer_cases_sex1)
length(zero_cancer_cases_sex1$cancer_classification)

nrow(zero_cancer_cases_sex2)
table(zero_cancer_cases_sex2$cancer_classification)

colnames(Europe_cases_cancer_description)

# 1. Define the columns to be kept as identifiers (id.vars)
id_cols <- c("CANCER", "SEX", "YEAR",
             "SEX_CHARACTER","cancer_classification",
              "TOTAL", "TOTAL_N")

# 2. Define the columns to be melted (measure.vars)
# These are the age-split columns

measure_cols <- c("N0_4", "N5_9", "N10_14", "N15_19", "N20_24",
                  "N25_29", "N30_34", "N35_39", "N40_44", "N45_49",
                  "N50_54", "N55_59", "N60_64", "N65_69", "N70_74",
                  "N75_79", "N80_84", "N85+", "N_UNK")

# 3. Apply 

europe_cases_cancer_long <- melt(
  Europe_cases_cancer_description,
  id.vars = id_cols,
  measure.vars = measure_cols,
  variable.name = "Age_Group", # New column for the original column names
  value.name = "Cases_by_Age" # New column for the values
)

head(europe_cases_cancer_long)

original_cols <- names(europe_cases_cancer_long)
new_order <- c()

for (col in original_cols){
  # Check if 'SEX' is next, insert 'SEX_CHARACTER" first
  
  if (col == "SEX") {
    new_order <- c(new_order, "SEX_CHARACTER")
  }
  
  if (col == "CANCER") {
    new_order <- c(new_order, "cancer_classification")
  }
  
  if (col == "TOTAL_N") {
    new_order <- c(new_order, "YEAR", "Age_Group", "Cases_by_Age",
                   "TOTAL")
  }
  # Only Include the column if it hasn't been moved yet
  if (!col %in% c("SEX_CHARACTER", "cancer_classification", 
                   "YEAR", 
                  "Age_Group", "Cases_by_Age", "TOTAL")) {
    new_order <- c(new_order, col)
  }
}

# Remove any dupplicated/unwanted columns that might have slipped through
new_order <- unique(new_order)

# Set the new order
setcolorder(europe_cases_cancer_long, new_order)
colnames(europe_cases_cancer_long)
head(europe_cases_cancer_long)

# write the transform data in the transformed_data folder

europe_cancer_cases_csv <- europe_cases_cancer_long[, 
                                                    -c("TOTAL_N", "CANCER", "SEX")]


colnames(europe_cancer_cases_csv)
head(europe_cancer_cases_csv)


aux_europe_by_cancer_cases_by_sex <- europe_cancer_cases_csv[, 
                                                          -c( "TOTAL")]
head(aux_europe_by_cancer_cases_by_sex)

# 1. Define the identifying columns (columns that remain as rows)
# Everything except the columns to spread (SEX_CHARACTER) and the valye columns
# id_cols <- c("cancer_classification", "YEAR", "Age_Group")

# 2. Define the columns containing the values to be spread
# value_cols <- c("TOTAL", "Cases_by_Age")

# 3. Apply dcast() to pivot the table

europe_by_cancer_cases_by_sex <- dcast(aux_europe_by_cancer_cases_by_sex,
                                          cancer_classification + YEAR + Age_Group ~ SEX_CHARACTER,
                                          value.var = "Cases_by_Age",
                                          fun.aggregate = sum)

europe_by_cancer_cases_by_sex <- europe_by_cancer_cases_by_sex[,TOTAL:= M + F]

head(europe_by_cancer_cases_by_sex)

processed_dir <-  "/home/jsancheg/git_environment/Europe_Cancer_Incidence/data/processed"

processed_file <- file.path(processed_dir,"europe_cancer_cases_long_format.csv")

processed_file_total_by_sex <- file.path(processed_dir,"total_europe_cancer_cases_by_sex.csv")
processed_file_by_age_group_sex <- file.path(processed_dir,"europe_by_cancer_cases_and_sex.csv")


write_csv(europe_cancer_cases_csv, processed_file)
write_csv(europe_total_cancer_cases_by_sex, processed_file_total_by_sex)
write_csv(europe_by_cancer_cases_by_sex, processed_file_by_age_group_sex)

