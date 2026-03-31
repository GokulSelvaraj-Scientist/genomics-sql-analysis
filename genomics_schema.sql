-- ============================================================
-- Clinical Genomics SQL Database: Real TCGA LUAD Data
-- Author: Gokul Selvaraj
-- GitHub: GokulSelvaraj-Scientist
-- Description: Build and query a clinical genomics database
--              using real patient data from The Cancer Genome
--              Atlas Lung Adenocarcinoma (TCGA-LUAD) cohort
-- Database: SQLite
-- Data Source: TCGA via TCGAbiolinks (https://portal.gdc.cancer.gov/)
-- ============================================================

-- ============================================================
-- PART 1: SCHEMA CREATION
-- ============================================================

DROP TABLE IF EXISTS survival_outcomes;
DROP TABLE IF EXISTS molecular_subtypes;
DROP TABLE IF EXISTS patients;
DROP TABLE IF EXISTS cancer_types;

-- Cancer type reference
CREATE TABLE cancer_types (
    cancer_type_id   INTEGER PRIMARY KEY,
    cancer_code      TEXT NOT NULL UNIQUE,
    cancer_name      TEXT NOT NULL,
    tissue_origin    TEXT NOT NULL,
    icd10_code       TEXT
);

-- Patient clinical data
CREATE TABLE patients (
    patient_id            TEXT PRIMARY KEY,
    cancer_type_id        INTEGER NOT NULL,
    age_at_diagnosis      INTEGER,
    gender                TEXT,
    tumor_stage           TEXT,
    smoking_status        TEXT,
    vital_status          TEXT,
    FOREIGN KEY (cancer_type_id) REFERENCES cancer_types(cancer_type_id)
);

-- Molecular subtype assignments from PCA clustering
CREATE TABLE molecular_subtypes (
    subtype_id    INTEGER PRIMARY KEY,
    patient_id    TEXT NOT NULL UNIQUE,
    subtype       TEXT NOT NULL,
    pc1           REAL,
    pc2           REAL,
    FOREIGN KEY (patient_id) REFERENCES patients(patient_id)
);

-- Survival outcomes
CREATE TABLE survival_outcomes (
    outcome_id                     INTEGER PRIMARY KEY,
    patient_id                     TEXT NOT NULL UNIQUE,
    overall_survival_days          INTEGER,
    progression_free_survival_days INTEGER,
    vital_status                   TEXT,
    os_years                       REAL,
    FOREIGN KEY (patient_id) REFERENCES patients(patient_id)
);

-- Indexes
CREATE INDEX idx_patients_cancer  ON patients(cancer_type_id);
CREATE INDEX idx_patients_stage   ON patients(tumor_stage);
CREATE INDEX idx_survival_patient ON survival_outcomes(patient_id);
CREATE INDEX idx_subtype_patient  ON molecular_subtypes(patient_id);
