# ============================================================
# Clinical Genomics SQL Database: Real TCGA LUAD Data
# Author: Gokul Selvaraj
# GitHub: GokulSelvaraj-Scientist
# Description: Build a SQLite database from real TCGA clinical
#              data and run analytical SQL queries
# Data: Real TCGA-LUAD patient data (downloaded via TCGAbiolinks)
# Note: Dataset contains 118 deceased patients with matched
#       RNA-seq and clinical data from TCGA-LUAD cohort
# ============================================================

library(RSQLite)
library(DBI)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

cat("============================================================\n")
cat("TCGA LUAD Clinical Genomics SQL Database\n")
cat("Real Patient Data from The Cancer Genome Atlas\n")
cat("============================================================\n\n")

# ============================================================
# STEP 1: LOAD AND PREPARE REAL TCGA DATA
# ============================================================

cat("=== Loading real TCGA data ===\n")

clinical <- read.csv("tcga_luad_clinical_processed.csv", stringsAsFactors = FALSE)
subtypes  <- read.csv("tcga_luad_subtype_assignments.csv", stringsAsFactors = FALSE)

cat("Clinical data:", nrow(clinical), "real TCGA-LUAD patients\n")

clinical_clean <- clinical %>%
  filter(!is.na(OS_time), OS_time > 0) %>%
  mutate(
    os_years         = round(OS_time / 365.0, 2),
    patient_id       = bcr_patient_barcode,
    age_at_diagnosis = as.integer(age),
    tumor_stage      = stage,
    vital_status     = "Dead",  # All patients in this cohort are deceased
    gender           = toupper(gender),
    # Create age groups for analysis
    age_group = case_when(
      age < 60  ~ "Under 60",
      age < 70  ~ "60-69",
      age < 80  ~ "70-79",
      TRUE      ~ "80+"
    ),
    # Create survival groups based on OS_time
    survival_group = case_when(
      OS_time < 365  ~ "Short (<1 year)",
      OS_time < 730  ~ "Medium (1-2 years)",
      OS_time < 1825 ~ "Long (2-5 years)",
      TRUE           ~ "Very Long (5+ years)"
    ),
    smoking_status = case_when(
      smoking == "Never"   ~ "Never Smoker",
      smoking == "Former"  ~ "Former Smoker",
      smoking == "Current" ~ "Current Smoker",
      TRUE ~ "Unknown"
    )
  )

cat("Patients loaded:", nrow(clinical_clean), "\n")
cat("Stage I:", sum(clinical_clean$tumor_stage == "Stage I"), "\n")
cat("Stage IV:", sum(clinical_clean$tumor_stage == "Stage IV"), "\n")
cat("Age range:", min(clinical_clean$age_at_diagnosis, na.rm=TRUE), "-",
    max(clinical_clean$age_at_diagnosis, na.rm=TRUE), "\n")

# ============================================================
# STEP 2: BUILD SQLITE DATABASE
# ============================================================

cat("\n=== Building SQLite database ===\n")

if (file.exists("tcga_luad_genomics.db")) file.remove("tcga_luad_genomics.db")
con <- dbConnect(SQLite(), "tcga_luad_genomics.db")

# Create schema
dbExecute(con, "CREATE TABLE cancer_types (
  cancer_type_id INTEGER PRIMARY KEY,
  cancer_code TEXT, cancer_name TEXT,
  tissue_origin TEXT, icd10_code TEXT)")

dbExecute(con, "CREATE TABLE patients (
  patient_id TEXT PRIMARY KEY,
  cancer_type_id INTEGER,
  age_at_diagnosis INTEGER,
  age_group TEXT,
  gender TEXT,
  tumor_stage TEXT,
  smoking_status TEXT,
  vital_status TEXT,
  survival_group TEXT)")

dbExecute(con, "CREATE TABLE survival_outcomes (
  outcome_id INTEGER PRIMARY KEY,
  patient_id TEXT,
  overall_survival_days INTEGER,
  os_years REAL,
  vital_status TEXT)")

dbExecute(con, "CREATE TABLE molecular_subtypes (
  subtype_id INTEGER PRIMARY KEY,
  patient_id TEXT,
  subtype TEXT,
  pc1 REAL,
  pc2 REAL)")

# Insert data
dbExecute(con, "INSERT INTO cancer_types VALUES
  (1,'LUAD','Lung Adenocarcinoma','Lung','C34.1')")

patients_df <- clinical_clean %>%
  mutate(cancer_type_id = 1) %>%
  select(patient_id, cancer_type_id, age_at_diagnosis, age_group,
         gender, tumor_stage, smoking_status, vital_status, survival_group)
dbWriteTable(con, "patients", patients_df, append=TRUE, row.names=FALSE)

survival_df <- clinical_clean %>%
  mutate(outcome_id = row_number(),
         overall_survival_days = as.integer(OS_time)) %>%
  select(outcome_id, patient_id, overall_survival_days, os_years, vital_status)
dbWriteTable(con, "survival_outcomes", survival_df, append=TRUE, row.names=FALSE)

subtypes_df <- subtypes %>%
  filter(sample %in% clinical_clean$patient_id) %>%
  mutate(subtype_id = row_number()) %>%
  select(subtype_id, patient_id=sample, subtype=Subtype, pc1=PC1, pc2=PC2)
dbWriteTable(con, "molecular_subtypes", subtypes_df, append=TRUE, row.names=FALSE)

cat("Database built successfully\n")
cat("Patients:", dbGetQuery(con,"SELECT COUNT(*) FROM patients")[1,1], "\n")

# ============================================================
# STEP 3: ANALYTICAL SQL QUERIES
# ============================================================

cat("\n=== Running analytical queries on real TCGA data ===\n")

# Query 1: Cohort overview
cat("\n-- Query 1: Cohort demographics --\n")
q1 <- dbGetQuery(con, "
  SELECT
    COUNT(p.patient_id)                                    AS total_patients,
    ROUND(AVG(p.age_at_diagnosis), 1)                      AS mean_age,
    MIN(p.age_at_diagnosis)                                AS min_age,
    MAX(p.age_at_diagnosis)                                AS max_age,
    SUM(CASE WHEN p.gender='MALE' THEN 1 ELSE 0 END)       AS n_male,
    SUM(CASE WHEN p.gender='FEMALE' THEN 1 ELSE 0 END)     AS n_female,
    ROUND(AVG(so.os_years), 2)                             AS mean_os_years,
    MIN(so.os_years)                                       AS min_os_years,
    MAX(so.os_years)                                       AS max_os_years
  FROM patients p
  JOIN survival_outcomes so ON p.patient_id = so.patient_id")
print(q1)

# Query 2: Survival by tumor stage
cat("\n-- Query 2: Survival by tumor stage --\n")
q2 <- dbGetQuery(con, "
  SELECT
    p.tumor_stage,
    COUNT(p.patient_id)                    AS n_patients,
    ROUND(AVG(so.os_years), 2)             AS mean_os_years,
    ROUND(MIN(so.os_years), 2)             AS min_os_years,
    ROUND(MAX(so.os_years), 2)             AS max_os_years,
    ROUND(AVG(so.overall_survival_days),0) AS mean_os_days
  FROM patients p
  JOIN survival_outcomes so ON p.patient_id = so.patient_id
  GROUP BY p.tumor_stage
  ORDER BY mean_os_years DESC")
print(q2)

# Query 3: Survival by molecular subtype
cat("\n-- Query 3: Survival by molecular subtype --\n")
q3 <- dbGetQuery(con, "
  SELECT
    ms.subtype,
    COUNT(p.patient_id)                    AS n_patients,
    ROUND(AVG(so.os_years), 2)             AS mean_os_years,
    ROUND(MIN(so.os_years), 2)             AS min_os_years,
    ROUND(MAX(so.os_years), 2)             AS max_os_years,
    ROUND(AVG(p.age_at_diagnosis), 1)      AS mean_age
  FROM patients p
  JOIN survival_outcomes so ON p.patient_id = so.patient_id
  JOIN molecular_subtypes ms ON p.patient_id = ms.patient_id
  GROUP BY ms.subtype
  ORDER BY mean_os_years DESC")
print(q3)

# Query 4: Smoking and survival
cat("\n-- Query 4: Smoking history and survival --\n")
q4 <- dbGetQuery(con, "
  SELECT
    p.smoking_status,
    COUNT(p.patient_id)                    AS n_patients,
    ROUND(AVG(p.age_at_diagnosis), 1)      AS mean_age,
    ROUND(AVG(so.os_years), 2)             AS mean_os_years,
    ROUND(MIN(so.os_years), 2)             AS min_os_years,
    ROUND(MAX(so.os_years), 2)             AS max_os_years
  FROM patients p
  JOIN survival_outcomes so ON p.patient_id = so.patient_id
  WHERE p.smoking_status != 'Unknown'
  GROUP BY p.smoking_status
  ORDER BY mean_os_years DESC")
print(q4)

# Query 5: Survival by age group
cat("\n-- Query 5: Survival by age group --\n")
q5 <- dbGetQuery(con, "
  SELECT
    p.age_group,
    COUNT(p.patient_id)                    AS n_patients,
    ROUND(AVG(p.age_at_diagnosis), 1)      AS mean_age,
    ROUND(AVG(so.os_years), 2)             AS mean_os_years,
    ROUND(AVG(so.overall_survival_days),0) AS mean_os_days
  FROM patients p
  JOIN survival_outcomes so ON p.patient_id = so.patient_id
  GROUP BY p.age_group
  ORDER BY CASE p.age_group
    WHEN 'Under 60' THEN 1
    WHEN '60-69' THEN 2
    WHEN '70-79' THEN 3
    WHEN '80+' THEN 4 END")
print(q5)

# Query 6: Survival group distribution
cat("\n-- Query 6: Survival duration distribution --\n")
q6 <- dbGetQuery(con, "
  SELECT
    p.survival_group,
    COUNT(p.patient_id)                    AS n_patients,
    ROUND(COUNT(p.patient_id)*100.0/118,1) AS pct_patients,
    ROUND(AVG(p.age_at_diagnosis),1)       AS mean_age,
    p.tumor_stage
  FROM patients p
  JOIN survival_outcomes so ON p.patient_id = so.patient_id
  GROUP BY p.survival_group, p.tumor_stage
  ORDER BY CASE p.survival_group
    WHEN 'Short (<1 year)' THEN 1
    WHEN 'Medium (1-2 years)' THEN 2
    WHEN 'Long (2-5 years)' THEN 3
    WHEN 'Very Long (5+ years)' THEN 4 END")
print(q6)

# Query 7: Subtype by stage
cat("\n-- Query 7: Molecular subtype by tumor stage --\n")
q7 <- dbGetQuery(con, "
  SELECT
    p.tumor_stage, ms.subtype,
    COUNT(*) AS n_patients,
    ROUND(AVG(so.os_years), 2) AS mean_os_years
  FROM patients p
  JOIN molecular_subtypes ms ON p.patient_id = ms.patient_id
  JOIN survival_outcomes so ON p.patient_id = so.patient_id
  GROUP BY p.tumor_stage, ms.subtype
  ORDER BY p.tumor_stage, ms.subtype")
print(q7)

# Query 8: Top 10 longest survivors
cat("\n-- Query 8: Top 10 longest surviving patients --\n")
q8 <- dbGetQuery(con, "
  SELECT
    p.patient_id, p.age_at_diagnosis, p.gender,
    p.tumor_stage, p.smoking_status,
    ms.subtype AS molecular_subtype,
    so.os_years, so.overall_survival_days
  FROM patients p
  JOIN survival_outcomes so ON p.patient_id = so.patient_id
  LEFT JOIN molecular_subtypes ms ON p.patient_id = ms.patient_id
  ORDER BY so.os_years DESC
  LIMIT 10")
print(q8)

# ============================================================
# STEP 4: VISUALISATIONS
# ============================================================

cat("\n=== Generating visualisations ===\n")

# Plot 1: Survival by stage
q2$tumor_stage <- factor(q2$tumor_stage, levels=c("Stage I","Stage IV"))
p1 <- ggplot(q2, aes(x=tumor_stage, y=mean_os_years, fill=tumor_stage)) +
  geom_bar(stat="identity", alpha=0.85, width=0.5) +
  geom_text(aes(label=paste0(mean_os_years, " yrs\nn=", n_patients)),
            vjust=-0.4, size=4, color="grey30") +
  scale_fill_manual(values=c("Stage I"="#2A9D8F", "Stage IV"="#E76F51")) +
  labs(title="Mean Overall Survival by Tumor Stage",
       subtitle="Real TCGA-LUAD Data (n=118 patients, all deceased)",
       x="Tumor Stage", y="Mean Overall Survival (years)",
       caption="Source: The Cancer Genome Atlas (TCGA)") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold"), legend.position="none") +
  ylim(0, max(q2$mean_os_years)*1.35)
ggsave("01_survival_by_stage.png", p1, width=7, height=6, dpi=300)
cat("Saved: 01_survival_by_stage.png\n")

# Plot 2: Survival by molecular subtype
subtype_colors <- c("Subtype 1"="#2A9D8F","Subtype 2"="#E76F51","Subtype 3"="#457B9D")
p2 <- ggplot(q3, aes(x=reorder(subtype,-mean_os_years), y=mean_os_years, fill=subtype)) +
  geom_bar(stat="identity", alpha=0.85, width=0.5) +
  geom_text(aes(label=paste0(mean_os_years, " yrs\nn=", n_patients)),
            vjust=-0.4, size=4, color="grey30") +
  scale_fill_manual(values=subtype_colors) +
  labs(title="Mean Overall Survival by Molecular Subtype",
       subtitle="Real TCGA-LUAD Data — RNA-seq derived subtypes",
       x="Molecular Subtype", y="Mean Overall Survival (years)",
       caption="Subtypes derived from unsupervised PCA clustering of gene expression") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold"), legend.position="none") +
  ylim(0, max(q3$mean_os_years)*1.35)
ggsave("02_survival_by_subtype.png", p2, width=7, height=6, dpi=300)
cat("Saved: 02_survival_by_subtype.png\n")

# Plot 3: Survival by smoking
q4_plot <- q4 %>% filter(smoking_status != "Unknown")
p3 <- ggplot(q4_plot, aes(x=reorder(smoking_status,-mean_os_years),
                           y=mean_os_years, fill=smoking_status)) +
  geom_bar(stat="identity", alpha=0.85, width=0.5) +
  geom_text(aes(label=paste0(mean_os_years, " yrs\nn=", n_patients)),
            vjust=-0.4, size=4, color="grey30") +
  scale_fill_manual(values=c("Never Smoker"="#2A9D8F",
                              "Former Smoker"="#E9C46A",
                              "Current Smoker"="#E76F51")) +
  labs(title="Survival by Smoking History",
       subtitle="Real TCGA-LUAD Data",
       x="Smoking Status", y="Mean Overall Survival (years)",
       caption="Source: The Cancer Genome Atlas (TCGA)") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold"), legend.position="none") +
  ylim(0, max(q4_plot$mean_os_years)*1.35)
ggsave("03_survival_by_smoking.png", p3, width=7, height=6, dpi=300)
cat("Saved: 03_survival_by_smoking.png\n")

# Plot 4: Survival by age group
p4 <- ggplot(q5, aes(x=age_group, y=mean_os_years, fill=age_group)) +
  geom_bar(stat="identity", alpha=0.85, width=0.5) +
  geom_text(aes(label=paste0(mean_os_years, " yrs\nn=", n_patients)),
            vjust=-0.4, size=4, color="grey30") +
  scale_fill_brewer(palette="Blues", direction=-1) +
  labs(title="Survival by Age Group at Diagnosis",
       subtitle="Real TCGA-LUAD Data",
       x="Age Group", y="Mean Overall Survival (years)",
       caption="Source: The Cancer Genome Atlas (TCGA)") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold"), legend.position="none") +
  ylim(0, max(q5$mean_os_years)*1.35)
ggsave("04_survival_by_age.png", p4, width=8, height=6, dpi=300)
cat("Saved: 04_survival_by_age.png\n")

# Plot 5: OS time distribution
os_data <- dbGetQuery(con, "
  SELECT p.tumor_stage, p.smoking_status, p.age_group,
         so.os_years, ms.subtype
  FROM patients p
  JOIN survival_outcomes so ON p.patient_id = so.patient_id
  LEFT JOIN molecular_subtypes ms ON p.patient_id = ms.patient_id
  WHERE p.smoking_status != 'Unknown'")

p5 <- ggplot(os_data, aes(x=tumor_stage, y=os_years, fill=tumor_stage)) +
  geom_violin(alpha=0.7, trim=FALSE) +
  geom_boxplot(width=0.12, fill="white", outlier.shape=NA) +
  geom_jitter(width=0.08, size=1.5, alpha=0.4) +
  scale_fill_manual(values=c("Stage I"="#2A9D8F", "Stage IV"="#E76F51")) +
  labs(title="Overall Survival Distribution by Stage",
       subtitle="Real TCGA-LUAD Data",
       x="Tumor Stage", y="Overall Survival (years)",
       caption="Source: The Cancer Genome Atlas (TCGA)") +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(face="bold"), legend.position="none")
ggsave("05_os_distribution.png", p5, width=7, height=6, dpi=300)
cat("Saved: 05_os_distribution.png\n")

# Save results
write.csv(q1, "query1_cohort_overview.csv",     row.names=FALSE)
write.csv(q2, "query2_survival_by_stage.csv",   row.names=FALSE)
write.csv(q3, "query3_survival_by_subtype.csv", row.names=FALSE)
write.csv(q4, "query4_smoking_survival.csv",    row.names=FALSE)
write.csv(q5, "query5_age_survival.csv",        row.names=FALSE)

dbDisconnect(con)

cat("\n=== Analysis Complete ===\n")
cat("Real TCGA patients analyzed:", nrow(clinical_clean), "\n")
cat("SQL queries executed: 8\n")
cat("Visualisations generated: 5\n")
