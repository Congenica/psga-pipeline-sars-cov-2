-- Deploy covid-pipeline:02-apptables to pg
-- requires: 01-appschema

BEGIN;

  -- Create tables, enums, et al for the SARS COV2 Application

  SET LOCAL search_path = sars_cov_2;

  CREATE TYPE "gender" AS ENUM (
    'M'
   ,'F'
   ,'U'
  );

  CREATE TYPE "governorate_name" AS ENUM (
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
      'Primary key, the area name, low level geolocation from where the sample was collected';
    COMMENT ON COLUMN "area"."latitude" IS
      'The angular distance of the location north or south of the equator';
    COMMENT ON COLUMN "area"."longitude" IS
      'the angular distance of the location east or west of the Greenwich meridian';
    COMMENT ON COLUMN "area"."colour" IS
      'The colour of this block in reports';

  CREATE TABLE IF NOT EXISTS "governorate" (
    "name"      "governorate_name" PRIMARY KEY
   ,"latitude"  DOUBLE PRECISION
   ,"longitude" DOUBLE PRECISION
   ,"colour"    BYTEA
  );

  COMMENT ON TABLE "governorate" IS
    'Geographic "Governorate" location';
    COMMENT ON COLUMN "area"."name" IS
      'Primary key, the area name, capitol level geolocation from where the sample was collected';
    COMMENT ON COLUMN "area"."latitude" IS
      'The angular distance of the location north or south of the equator';
    COMMENT ON COLUMN "area"."longitude" IS
      'the angular distance of the location east or west of the Greenwich meridian';
    COMMENT ON COLUMN "area"."colour" IS
      'The colour of this block in reports';

  CREATE TABLE IF NOT EXISTS "block" (
    "number"    VARCHAR PRIMARY KEY
   ,"latitude"  DOUBLE PRECISION
   ,"longitude" DOUBLE PRECISION
   ,"colour"    BYTEA
  );

  COMMENT ON TABLE "block" IS
    'Geographic "Block" location';
    COMMENT ON COLUMN "block"."number" IS
      'Primary key, the block number, low level geolocation from where the sample was collected';
    COMMENT ON COLUMN "block"."latitude" IS
      'The angular distance of the location north or south of the equator';
    COMMENT ON COLUMN "block"."longitude" IS
      'the angular distance of the location east or west of the Greenwich meridian';
    COMMENT ON COLUMN "block"."colour" IS
      'The colour of this block in reports';

  CREATE TABLE IF NOT EXISTS "comorbidity" (
    "comorbidity_id" INTEGER PRIMARY KEY
   ,"description"    TEXT
  );

  COMMENT ON TABLE "comorbidity" IS
    'Coded diseases presented in patient (eg Diabetes)';
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
   ,"governorate_name"    "governorate_name"
      REFERENCES "governorate" ("name")
      ON DELETE NO ACTION ON UPDATE NO ACTION DEFERRABLE INITIALLY DEFERRED
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
      'Lab sample identifier';
    COMMENT ON COLUMN "sample"."date_collected" IS
      'timestamp of sample collection';
    COMMENT ON COLUMN "sample"."data_sequenced" IS
      'timestamp of sample sequencing';
    COMMENT ON COLUMN "sample"."sequencing_run" IS
      'identifier of sequencing run';
    COMMENT ON COLUMN "sample"."gender" IS
      'Gender of host person';
    COMMENT ON COLUMN "sample"."age" IS
      'Age of host person in years';
    COMMENT ON COLUMN "sample"."nationality" IS
      'Nationality of host person';
    COMMENT ON COLUMN "sample"."governorate_name" IS
      'high level geolocation from where sample was collected';
    COMMENT ON COLUMN "sample"."area_name" IS
      'foreign key to area table, mid level geolocation from where sample was collected';
    COMMENT ON COLUMN "sample"."block_number" IS
      'foreign key to block table, low level geolocation from where sample was collected';
    COMMENT ON COLUMN "sample"."sample_number" IS
      'Lab sample identifier (may be redundant with lab_id)';
    COMMENT ON COLUMN "sample"."ct_value" IS
      'RNA quality measure';
    COMMENT ON COLUMN "sample"."symptoms" IS
      'Symptoms observed';
    COMMENT ON COLUMN "sample"."travel_exposure" IS
      'Country host visited recently';
    COMMENT ON COLUMN "sample"."hospital_admittance" IS
      'Hospital type host admited to';
    COMMENT ON COLUMN "sample"."pangolin_lineage" IS
      'Viral phylogenetic lineage';
    COMMENT ON COLUMN "sample"."genbank_id" IS
      'GeneBank accession of virus genome';
    COMMENT ON COLUMN "sample"."gisaid_id" IS
      'GISAID identifier (GISAID is a global science initiative and primary source that provides open-access to genomic data of the novel coronavirus responsible for COVID-19)';
    COMMENT ON COLUMN "sample"."genome_length" IS
      'Number of base pair in the virus genome';
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
