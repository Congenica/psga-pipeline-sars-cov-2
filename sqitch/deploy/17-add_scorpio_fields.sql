-- add scorpio fields as calculated by pangolin
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "scorpio_call" VARCHAR;
    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "scorpio_support" DOUBLE PRECISION;
    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "scorpio_conflict" DOUBLE PRECISION;
    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "note" VARCHAR;

    COMMENT ON COLUMN "sample"."scorpio_call"
        IS 'If a query is assigned a constellation by scorpio this call is output in this column. The full set of constellations searched by default can be found at the constellations repository';
    COMMENT ON COLUMN "sample"."scorpio_support"
        IS 'The support score is the proportion of defining variants which have the alternative allele in the sequence.';
    COMMENT ON COLUMN "sample"."scorpio_conflict"
        IS 'The conflict score is the proportion of defining variants which have the reference allele in the sequence. Ambiguous/other non-ref/alt bases at each of the variant positions contribute only to the denominators of these scores';
    COMMENT ON COLUMN "sample"."note"
        IS 'If any conflicts from the decision tree, this field will output the alternative assignments. If the sequence failed QC this field will describe why. If the sequence met the SNP thresholds for scorpio to call a constellation, it’ll describe the exact SNP counts of Alt, Ref and Amb (Alternative, reference and ambiguous) alleles for that call';

COMMIT;

