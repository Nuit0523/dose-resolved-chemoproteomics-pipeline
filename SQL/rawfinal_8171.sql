-- =========================================================
-- 8171 Chemoproteomics SQL Pipeline
-- Purpose:
--   1. Import CurvCurator concentration-response output.
--   2. Identify significant dose-response hit sites.
--   3. Count hit sites and unique genes.
--   4. Select one best site per gene for enrichment analysis.
--   5. Generate summary tables for thesis writing.
--
-- Database: PostgreSQL
-- Input file: E:/gradthesis/curvecurator/8171/8171_conc_curve.tsv
-- Import method: DBeaver TSV import into rawfinal_8171
-- =========================================================


-- =========================================================
-- Step 0. Clean Previous Objects
-- =========================================================

DROP TABLE IF EXISTS rawfinal_8171 CASCADE;
DROP VIEW IF EXISTS gene_list;
DROP VIEW IF EXISTS best_site_per_gene;
DROP VIEW IF EXISTS hit_sites;
DROP TABLE IF EXISTS rawfinal_8171;

-- =========================================================
-- 8171 Chemoproteomics SQL Pipeline
-- Purpose:
--   1. Import CurvCurator concentration-response output.
--   2. Identify significant dose-response hit sites.
--   3. Count hit sites and unique genes.
--   4. Select one best site per gene for enrichment analysis.
--   5. Generate summary tables for thesis writing.
--
-- Database: PostgreSQL
-- Input file: E:/gradthesis/curvecurator/8171/8171_conc_curve.tsv
-- Import method: DBeaver TSV import into rawfinal_8171
-- =========================================================


-- =========================================================
-- Step 0. Clean Previous Objects
-- =========================================================

DROP VIEW IF EXISTS gene_list;
DROP VIEW IF EXISTS best_site_per_gene;
DROP VIEW IF EXISTS hit_sites;
DROP TABLE IF EXISTS rawfinal_8171;


-- =========================================================
-- Step 1. Create Raw Table
-- =========================================================

CREATE TABLE rawfinal_8171 (
    name text,
    "Raw 1" double precision,
    "Raw 2" double precision,
    "Raw 3" double precision,
    "Raw 4" double precision,
    "Raw 5" double precision,
    "Raw 6" double precision,
    "Imputation N" bigint,
    "Imputation Position" text,
    "Normalized 1" double precision,
    "Normalized 2" double precision,
    "Normalized 3" double precision,
    "Normalized 4" double precision,
    "Normalized 5" double precision,
    "Normalized 6" double precision,
    "Ratio 1" double precision,
    "Ratio 2" double precision,
    "Ratio 3" double precision,
    "Ratio 4" double precision,
    "Ratio 5" double precision,
    "Ratio 6" double precision,
    "Signal Quality" double precision,
    pec50 double precision,
    "Curve Slope" double precision,
    "Curve Front" double precision,
    "Curve Back" double precision,
    "Curve Fold Change" double precision,
    "Curve AUC" double precision,
    "Curve RMSE" double precision,
    "Curve R2" double precision,
    "pEC50 Error" double precision,
    "Curve Slope Error" double precision,
    "Curve Front Error" double precision,
    "Curve Back Error" double precision,
    "Null Model" double precision,
    "Null RMSE" double precision,
    "Curve F_Value" double precision,
    "Curve P_Value" double precision,
    "Curve Log P_Value" double precision,
    "Curve P_Value adjusted" double precision,
    "Curve Log P_Value adjusted" double precision,
    "Curve Regulation" text
);
--check empty table
SELECT COUNT(*) AS row_count
FROM rawfinal_8171;
SELECT COUNT(*) AS column_count
FROM information_schema.columns
WHERE table_name = 'rawfinal_8171';
--import the raw data
COPY rawfinal_8171
FROM 'E:/gradthesis/curvecurator/8171/8171_conc_curve.tsv'
WITH (
    FORMAT csv,
    HEADER true,
    DELIMITER E'\t',
    NULL '',
    QUOTE E'\b'
);
--check the imported table
SELECT COUNT(*) AS total_rows
FROM rawfinal_8171;
--check if it has missing value
SELECT
    COUNT(*) AS missing_key_values
FROM rawfinal_8171
WHERE name IS NULL
   OR pec50 IS NULL
   OR "Curve R2" IS NULL
   OR "Curve Fold Change" IS NULL
   OR "Curve Log P_Value adjusted" IS NULL;
--check the regulation
SELECT
    COUNT(*) AS missing_curve_regulation
FROM rawfinal_8171
WHERE "Curve Regulation" IS NULL
   OR "Curve Regulation" = '';
--check the deviation of the regulation
SELECT
    COALESCE(NULLIF("Curve Regulation", ''), 'blank') AS curve_regulation,
    COUNT(*) AS row_count
FROM rawfinal_8171
GROUP BY COALESCE(NULLIF("Curve Regulation", ''), 'blank')
ORDER BY row_count DESC;
----------------------------------build up the hit_sites view---------------------------
CREATE OR REPLACE VIEW hit_sites AS
SELECT *
FROM rawfinal_8171
WHERE "Curve Log P_Value adjusted" >= 1.30103
  AND pec50 >= 4
  AND "Curve R2" >= 0.99
  AND ABS("Curve Fold Change") >= 2;
--check the hit ammount
SELECT COUNT(*) AS hit_site_count
FROM hit_sites;
--check the top20 hits
SELECT
    split_part(name, '_', 1) AS gene,
    name,
    pec50,
    "Curve Fold Change",
    "Curve R2",
    "Curve Log P_Value adjusted",
    "Curve Regulation"
FROM hit_sites
ORDER BY pec50 DESC
LIMIT 20;
-------------------------count the unique genes and multu-site genes-----------------------
--count unique genes
SELECT COUNT(DISTINCT split_part(name, '_', 1)) AS unique_gene_count
FROM hit_sites;
--count the hits dits of each genes (multiple hits genes)
SELECT
    split_part(name, '_', 1) AS gene,
    COUNT(*) AS site_count
FROM hit_sites
GROUP BY split_part(name, '_', 1)
HAVING COUNT(*) > 1
ORDER BY site_count DESC, gene;
--count the single-stie genes and multi-site genes
SELECT COUNT(*) AS genes_with_one_site
FROM (
    SELECT split_part(name, '_', 1) AS gene
    FROM hit_sites
    GROUP BY split_part(name, '_', 1)
    HAVING COUNT(*) = 1
) t;

SELECT COUNT(*) AS genes_with_multiple_sites
FROM (
    SELECT split_part(name, '_', 1) AS gene
    FROM hit_sites
    GROUP BY split_part(name, '_', 1)
    HAVING COUNT(*) > 1
) t;
------------------------------build best_site_per_gene--------------------------
CREATE OR REPLACE VIEW best_site_per_gene AS
SELECT DISTINCT ON (split_part(name, '_', 1))
       *
FROM hit_sites
ORDER BY
    split_part(name, '_', 1),
    pec50 DESC,
    "Curve Log P_Value adjusted" DESC,
    "Curve R2" DESC,
    ABS("Curve Fold Change") DESC,
    name;
--check the amount
SELECT COUNT(*) AS best_site_gene_count
FROM best_site_per_gene;
--check the top 20
SELECT
    split_part(name, '_', 1) AS gene,
    name,
    pec50,
    "Curve Fold Change",
    "Curve R2",
    "Curve Log P_Value adjusted",
    "Curve Regulation"
FROM best_site_per_gene
ORDER BY gene
LIMIT 20;
---------------------------build gene_list------------------------
CREATE OR REPLACE VIEW gene_list AS
SELECT
    split_part(name, '_', 1) AS gene
FROM best_site_per_gene
ORDER BY gene;
--check the amount of gene
SELECT COUNT(*) AS gene_count
FROM gene_list;
--check the top30
SELECT *
FROM gene_list
LIMIT 30;
--check if the gene is repeated
SELECT
    gene,
    COUNT(*) AS n
FROM gene_list
GROUP BY gene
HAVING COUNT(*) > 1;
------------------------------summary stats------------------------------
SELECT
    COUNT(*) AS hit_sites,
    COUNT(DISTINCT split_part(name, '_', 1)) AS unique_genes,
    AVG(pec50) AS mean_pec50,
    MIN(pec50) AS min_pec50,
    MAX(pec50) AS max_pec50,
    AVG("Curve Fold Change") AS mean_fold_change,
    MIN("Curve Fold Change") AS min_fold_change,
    MAX("Curve Fold Change") AS max_fold_change,
    AVG("Curve R2") AS mean_r2,
    MIN("Curve R2") AS min_r2,
    MAX("Curve R2") AS max_r2
FROM hit_sites;
--regulation: up/down/not
SELECT
    "Curve Regulation",
    COUNT(*) AS site_count
FROM hit_sites
GROUP BY "Curve Regulation"
ORDER BY site_count DESC, "Curve Regulation";
--------------------------------export all hit sites------------------------
SELECT
    split_part(name, '_', 1) AS gene,
    name,
    pec50,
    "Curve Fold Change",
    "Curve R2",
    "Curve Log P_Value adjusted",
    "Curve Regulation"
FROM hit_sites
ORDER BY
    pec50 DESC,
    "Curve Log P_Value adjusted" DESC,
    "Curve R2" DESC,
    ABS("Curve Fold Change") DESC,
    name;