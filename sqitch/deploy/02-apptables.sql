-- Deploy ps-bahrain-covid:02-apptables to pg
-- requires: 01-appschema

BEGIN;

  -- Create tables, enums, et al for the SARS COV2 (Bahrain) Application

  SET LOCAL search_path = sars_cov_2;

  CREATE TYPE "gender" AS ENUM (
    'M'
   ,'F'
   ,'U'
  );

  CREATE TYPE "governerate" AS ENUM (
    'Capital'
   ,'Northern'
   ,'Southern'
   ,'Muharraq'
  );

  CREATE TYPE "hospital_admittance" AS ENUM (
    'No'
   ,'Standard'
   ,'ICU'
  );

  -- NOTE: If geo coords are used often, replace lat/long with PostGIS
  -- https://postgis.net/install/
  -- Initially PostGIS extension not used for simplicity
  CREATE TABLE IF NOT EXISTS "area" (
    "name"      VARCHAR PRIMARY KEY
   ,"latitude"  DOUBLE PRECISION
   ,"longitude"	DOUBLE PRECISION
   ,"colour"  	BYTEA
  );

  COMMENT ON TABLE "area" IS
    'Geographic "Area" location';
    COMMENT ON COLUMN "area"."name" IS
      'Primary key, the area name';
    COMMENT ON COLUMN "area"."latitude" IS
      'The angular distance of the location north or south of the equator';
    COMMENT ON COLUMN "area"."longitude" IS
      'the angular distance of the location east or west of the Greenwich meridian';
    COMMENT ON COLUMN "area"."colour" IS
      '';

  CREATE TABLE IF NOT EXISTS "block" (
    "number"    VARCHAR PRIMARY KEY
   ,"latitude"  DOUBLE PRECISION
   ,"longitude" DOUBLE PRECISION
   ,"colour"    BYTEA
  );

  COMMENT ON TABLE "block" IS
    'Geographic "Block" location';
    COMMENT ON COLUMN "block"."number" IS
      'Primary key, the block number';
    COMMENT ON COLUMN "block"."latitude" IS
      'The angular distance of the location north or south of the equator';
    COMMENT ON COLUMN "block"."longitude" IS
      'the angular distance of the location east or west of the Greenwich meridian';
    COMMENT ON COLUMN "block"."colour" IS
      '';

  CREATE TABLE IF NOT EXISTS "comorbidity" (
    "comorbidity_id" INTEGER PRIMARY KEY
   ,"description"    TEXT
  );

  COMMENT ON TABLE "comorbidity" IS
    'Comorbidity details';
    COMMENT ON COLUMN "comorbidity"."comorbidity_id" IS
      'Primary key';
    COMMENT ON COLUMN "comorbidity"."description" IS
      'Description of the comorbidity';

  CREATE TABLE IF NOT EXISTS "sample" (
    "sample_id"      SERIAL PRIMARY KEY
   ,"lab_id"         VARCHAR
   ,"date_collected" TIMESTAMP WITHOUT TIME ZONE
   ,"data_sequenced" TIMESTAMP WITHOUT TIME ZONE
   ,"sequencing_run" VARCHAR
   ,"gender"	     "gender" NOT NULL DEFAULT 'U'
   ,"age"            INTEGER -- Why not DOB?
   ,"nationality"    VARCHAR
   ,"governerate"    "governerate"
   ,"area_name"      VARCHAR
      REFERENCES "area" ("name")
      ON DELETE NO ACTION ON UPDATE NO ACTION DEFERRABLE INITIALLY DEFERRED
   ,"block_number"   VARCHAR
      REFERENCES     "block" ("number")
      ON DELETE NO ACTION ON UPDATE NO ACTION DEFERRABLE INITIALLY DEFERRED
   ,"sample_number"  INTEGER
   ,"ct_value"       DOUBLE PRECISION
   ,"symptoms"       TEXT
   ,"travel_exposure"  VARCHAR
   ,"hospital_admittance" "hospital_admittance"
   ,"pangolin_lineage" VARCHAR
   ,"genbank_id"     INTEGER
   ,"gisaid_id"      VARCHAR
   ,"genome_length"  INTEGER
   ,"mrn"            INTEGER
  );

  COMMENT ON TABLE "sample" IS
    'Sample data';
    COMMENT ON COLUMN "sample"."sample_id" IS
      'Primary key, serial sequence generated in the database';
    COMMENT ON COLUMN "sample"."lab_id" IS
      '';
    COMMENT ON COLUMN "sample"."date_collected" IS
      'timestamp of sample collection';
    COMMENT ON COLUMN "sample"."data_sequenced" IS
      'timestamp of sample sequencing';
    COMMENT ON COLUMN "sample"."sequencing_run" IS
      '';
    COMMENT ON COLUMN "sample"."gender" IS
      'Gender of person';
    COMMENT ON COLUMN "sample"."age" IS
      'Age of person in years';
    COMMENT ON COLUMN "sample"."nationality" IS
      'Nationality of person';
    COMMENT ON COLUMN "sample"."governerate" IS
      '';
    COMMENT ON COLUMN "sample"."area_name" IS
      'foreign key to area table';
    COMMENT ON COLUMN "sample"."block_number" IS
      'foreign key to block table';
    COMMENT ON COLUMN "sample"."sample_number" IS
      '';
    COMMENT ON COLUMN "sample"."ct_value" IS
      '';
    COMMENT ON COLUMN "sample"."symptoms" IS
      '';
    COMMENT ON COLUMN "sample"."travel_exposure" IS
      '';
    COMMENT ON COLUMN "sample"."hospital_admittance" IS
      '';
    COMMENT ON COLUMN "sample"."pangolin_lineage" IS
      '';
    COMMENT ON COLUMN "sample"."genbank_id" IS
      '';
    COMMENT ON COLUMN "sample"."gisaid_id" IS
      '';
    COMMENT ON COLUMN "sample"."genome_length" IS
      '';
    COMMENT ON COLUMN "sample"."mrn" IS
      'Individual National Identity number';


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
      '';
    COMMENT ON COLUMN "sample_qc"."pct_covered_bases" IS
      '';
    COMMENT ON COLUMN "sample_qc"."longest_no_N_run" IS
      '';
    COMMENT ON COLUMN "sample_qc"."num_aligned_reads" IS
      '';
    COMMENT ON COLUMN "sample_qc"."qc_pass" IS
      '';
    COMMENT ON COLUMN "sample_qc"."qc_plot" IS
      '';
    COMMENT ON COLUMN "sample_qc"."pipeline_version" IS
      '';


  CREATE TABLE IF NOT EXISTS "sample_comorbidity" (
    "sample_id" INTEGER
      REFERENCES "sample" ("sample_id")
      ON DELETE NO ACTION ON UPDATE NO ACTION DEFERRABLE INITIALLY DEFERRED
   ,"comorbidity_id" INTEGER
      REFERENCES "comorbidity" ("comorbidity_id")
      ON DELETE NO ACTION ON UPDATE NO ACTION DEFERRABLE INITIALLY DEFERRED
   ,PRIMARY KEY ("sample_id", "comorbidity_id")
  );

  COMMENT ON TABLE "sample_comorbidity" IS
    'Sample comorbitiy, a linking table between the sample an comorbidity record';
    COMMENT ON COLUMN "sample_comorbidity"."sample_id" IS
      'Primary key with comorbidity_id, and foreign key to sample table';
    COMMENT ON COLUMN "sample_comorbidity"."comorbidity_id" IS
      'Primary key with sample_id, and foreign key to comorbidity table';


COMMIT;
