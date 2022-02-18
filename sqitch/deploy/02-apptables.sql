-- Deploy covid-pipeline:02-apptables to pg
-- requires: 01-appschema

BEGIN;

  -- Create tables, enums, et al for the SARS COV2 Application

  SET LOCAL search_path = sars_cov_2;

  CREATE TABLE IF NOT EXISTS "sample" (
    "sample_id"      SERIAL PRIMARY KEY
   ,"lab_id"         VARCHAR
   ,"date_collected" TIMESTAMP WITHOUT TIME ZONE
   ,"data_sequenced" TIMESTAMP WITHOUT TIME ZONE
   ,"sequencing_run" VARCHAR
   ,"sample_number"  VARCHAR
   ,"pangolin_lineage" VARCHAR
   ,"genbank_id"     INTEGER
   ,"gisaid_id"      VARCHAR
   ,"genome_length"  INTEGER
  );

  COMMENT ON TABLE "sample" IS
    'Sample data';
    COMMENT ON COLUMN "sample"."sample_id" IS
      'Primary key, serial sequence generated in the database';
    COMMENT ON COLUMN "sample"."lab_id" IS
      'Lab sample identifier';
    COMMENT ON COLUMN "sample"."date_collected" IS
      'timestamp of sample collection';
    COMMENT ON COLUMN "sample"."data_sequenced" IS
      'timestamp of sample sequencing';
    COMMENT ON COLUMN "sample"."sequencing_run" IS
      'identifier of sequencing run';
    COMMENT ON COLUMN "sample"."sample_number" IS
      'Lab sample identifier (may be redundant with lab_id)';
    COMMENT ON COLUMN "sample"."pangolin_lineage" IS
      'Viral phylogenetic lineage';
    COMMENT ON COLUMN "sample"."genbank_id" IS
      'GeneBank accession of virus genome';
    COMMENT ON COLUMN "sample"."gisaid_id" IS
      'GISAID identifier (GISAID is a global science initiative and primary source that provides open-access to genomic data of the novel coronavirus responsible for COVID-19)';
    COMMENT ON COLUMN "sample"."genome_length" IS
      'Number of base pair in the virus genome';


  CREATE TABLE IF NOT EXISTS "sample_qc" (
    "sample_id" INTEGER
      REFERENCES "sample" ("sample_id")
      ON DELETE NO ACTION ON UPDATE NO ACTION DEFERRABLE INITIALLY DEFERRED
   ,"pct_N_bases" DOUBLE PRECISION
   ,"pct_covered_bases"	DOUBLE PRECISION
   ,"longest_no_N_run" INTEGER
   ,"num_aligned_reads"	DOUBLE PRECISION
   ,"qc_pass" BOOLEAN
   ,"qc_plot" BYTEA
   ,"pipeline_version" VARCHAR
  );

  COMMENT ON TABLE "sample_qc" IS
    'Sample Quality Control, QC data for a sample';
    COMMENT ON COLUMN "sample_qc"."sample_id" IS
      'Primary key, and foreign key to sample table';
    COMMENT ON COLUMN "sample_qc"."pct_N_bases" IS
      'percentage of N bases in sequenced genome';
    COMMENT ON COLUMN "sample_qc"."pct_covered_bases" IS
      'percentage of covered bases in sequenced genome';
    COMMENT ON COLUMN "sample_qc"."longest_no_N_run" IS
      'longest sequence run without N base';
    COMMENT ON COLUMN "sample_qc"."num_aligned_reads" IS
      'Number of reads succesfully aligned and mapped';
    COMMENT ON COLUMN "sample_qc"."qc_pass" IS
      'Has sample passed QC';
    COMMENT ON COLUMN "sample_qc"."qc_plot" IS
      'QC plot image';
    COMMENT ON COLUMN "sample_qc"."pipeline_version" IS
      'mapping pipeline version';

COMMIT;
