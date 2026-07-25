select table_name
from information_schema.tables
where table_schema = 'public';

select table_schema, table_name
from information_schema."tables" t 
where table_name = 'protein_sites';

select*
from final_8171_analysis
limit 5;

select *
from "8171_conc_curve"
limit 5;

select count (*)
from "8171_conc_curve";

select 
	MIN("pEC50"),
	max("pEC50"),
	avg("pEC50")
from "8171_conc_curve";

select count (*)
from "8171_conc_curve"
where "pEC50" is null;

select 
	MIN("Curve R2"),
	max("Curve R2"),
	avg("Curve R2")
from "8171_conc_curve";

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE "Curve R2" < 0.8;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE "Curve R2" < 0.5;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE "Curve R2" < 0.2;

SELECT
    MIN("Curve Fold Change"),
    MAX("Curve Fold Change"),
    AVG("Curve Fold Change")
FROM "8171_conc_curve";

SELECT
    MIN("Signal Quality"),
    MAX("Signal Quality"),
    AVG("Signal Quality")
FROM "8171_conc_curve";

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE "pEC50" > 4;

SELECT
    "Name",
    "pEC50",
    "Curve R2",
    "Curve Fold Change"
FROM "8171_conc_curve"
ORDER BY "pEC50"
LIMIT 20;

SELECT
    "Name",
    "Curve Fold Change",
    "pEC50",
    "Curve R2"
FROM "8171_conc_curve"
ORDER BY "Curve Fold Change"
LIMIT 20;

SELECT
    "Name",
    "pEC50",
    "Curve R2",
    "Signal Quality",
    "Curve Fold Change"
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
  AND "Curve R2" < 0.5
ORDER BY "Curve R2";、

SELECT
    COUNT(*) OVER() AS total_count,
    "Name",
    "pEC50",
    "Curve R2",
    "Signal Quality",
    "Curve Fold Change"
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
  AND "Curve R2" < 0.5
ORDER BY "Curve R2";

SELECT
    AVG("Signal Quality")
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
AND "Curve R2" < 0.5;

SELECT
    AVG("Signal Quality")
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
AND "Curve R2" >= 0.5;

SELECT
    AVG("Curve Fold Change")
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
AND "Curve R2" < 0.5;

SELECT
AVG(ABS("Curve Fold Change"))
FROM "8171_conc_curve"
WHERE "pEC50">=4
AND "Curve R2"<0.5;

SELECT
AVG(ABS("Curve Fold Change"))
FROM "8171_conc_curve"
WHERE "pEC50">=4
AND "Curve R2">=0.5;

SELECT
    MIN("Curve RMSE"),
    MAX("Curve RMSE"),
    AVG("Curve RMSE")
FROM "8171_conc_curve";

SELECT
    AVG("Curve RMSE")
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
AND "Curve R2" < 0.5;

SELECT
    AVG("Curve RMSE")
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
AND "Curve R2" >= 0.5;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE "Curve RMSE" > 0.2;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE "Curve RMSE" > 0.2
  AND "Curve R2" < 0.5;

SELECT
    MIN("Curve Slope"),
    MAX("Curve Slope"),
    AVG("Curve Slope")
FROM "8171_conc_curve";

SELECT AVG("Curve Slope")
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
  AND "Curve R2" < 0.5;

SELECT AVG("Curve Slope")
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
  AND "Curve R2" >= 0.5;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE "Curve Slope" = 10;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE "Curve Slope" = 10
AND "Curve R2" < 0.5;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
AND "Curve R2" < 0.5
AND "Curve Slope" = 10;

SELECT
    "Name",
    "pEC50",
    "Curve R2",
    "Curve Slope",
    "Curve Fold Change",
    "Curve RMSE",
    "Signal Quality"
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
  AND "Curve R2" < 0.5
  AND "Curve Slope" = 10;

SELECT
    "Curve Regulation",
    COUNT(*) AS protein_count
FROM "8171_conc_curve"
WHERE
    "pEC50" >= 4
AND "Curve R2" < 0.5
AND "Curve Slope" = 10
GROUP BY "Curve Regulation"
ORDER BY protein_count DESC;

SELECT
    MIN("Curve P_Value"),
    MAX("Curve P_Value"),
    AVG("Curve P_Value")
FROM "8171_conc_curve";

SELECT AVG("Curve P_Value")
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
AND "Curve R2" < 0.5;

SELECT AVG("Curve P_Value")
FROM "8171_conc_curve"
WHERE "pEC50" >= 4
AND "Curve R2" >= 0.5;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE
    "pEC50" >= 4
AND "Curve R2" < 0.5
AND "Curve P_Value" < 0.05;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE
    "pEC50" >= 4
AND "Curve R2" >= 0.5
AND "Curve P_Value" < 0.05;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE
    "pEC50" >= 4
AND "Curve R2" >= 0.5
AND "Curve P_Value" >= 0.05;

SELECT
    "Name",
    "Curve P_Value",
    "Curve R2",
    "pEC50"
FROM "8171_conc_curve"
WHERE
    "pEC50" >= 4
ORDER BY "Curve P_Value"
LIMIT 20;

SELECT
    MIN("Curve Log P_Value adjusted"),
    MAX("Curve Log P_Value adjusted"),
    AVG("Curve Log P_Value adjusted")
FROM "8171_conc_curve";

SELECT
    AVG("Curve Log P_Value adjusted")
FROM "8171_conc_curve"
WHERE
    "pEC50" >= 4
AND "Curve R2" < 0.5;

SELECT
    AVG("Curve Log P_Value adjusted")
FROM "8171_conc_curve"
WHERE
    "pEC50" >= 4
AND "Curve R2" >= 0.5;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE
    "pEC50" >=4
AND "Curve Log P_Value adjusted" >=1.30103;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE
    "pEC50" >=4
AND "Curve R2"<0.5
AND "Curve Log P_Value adjusted">=1.30103;

SELECT COUNT(*)
FROM "8171_conc_curve"
WHERE
    "pEC50" >=4
AND "Curve R2">=0.5
AND "Curve Log P_Value adjusted">=1.30103;

SELECT *
FROM "8171_conc_curve"
WHERE "Curve Log P_Value adjusted">=1.30103
  AND "pEC50" >= 4.0
  AND "Curve R2" >= 0.99
  AND ABS("Curve Fold Change") >= 2;

SELECT COUNT(*) AS hit_count
FROM "8171_conc_curve"
WHERE "Curve Log P_Value adjusted" >= 1.30103
  AND "pEC50" >= 4.0
  AND "Curve R2" >= 0.99
  AND ABS("Curve Fold Change") >= 2;

SELECT
    COUNT(DISTINCT split_part("Name",'_',1)) AS gene_number
FROM "8171_conc_curve"
WHERE "Curve Log P_Value adjusted" >= 1.30103
  AND "pEC50" >= 4
  AND "Curve R2" >= 0.99
  AND ABS("Curve Fold Change") >= 2;

SELECT
    split_part("Name",'_',1) AS gene,
    COUNT(*) AS site_count
FROM "8171_conc_curve"
WHERE "Curve Log P_Value adjusted" >= 1.30103
  AND "pEC50" >= 4
  AND "Curve R2" >= 0.99
  AND ABS("Curve Fold Change") >= 2
GROUP BY gene
HAVING COUNT(*) > 1
ORDER BY site_count DESC;

ALTER TABLE public.newtable
ALTER COLUMN "default~quantified_peptides" TYPE text;

SELECT
    column_name,
    data_type,
    character_maximum_length
FROM information_schema.columns
WHERE table_name = 'newtable'
  AND column_name = 'default~quantified_peptides';

ALTER TABLE public.newtable
ALTER COLUMN "Protein Id" TYPE text,
ALTER COLUMN gene_symbol TYPE text,
ALTER COLUMN prot_description TYPE text,
ALTER COLUMN motif TYPE text,
ALTER COLUMN best_search_name TYPE text,
ALTER COLUMN file TYPE text,
ALTER COLUMN redundancy TYPE text,
ALTER COLUMN "sequence" TYPE text,
ALTER COLUMN "default~quantified_peptides" TYPE text;
TRUNCATE TABLE public.newtable;
