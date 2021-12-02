-- Remove pangolin_probability and add pangolin conflict / ambiguity_score to sample table
-- requires: 02-apptables

BEGIN;

    SET LOCAL search_path = sars_cov_2;

    ALTER TABLE "sample" DROP COLUMN IF EXISTS "pangolin_probability";

    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "pangolin_conflict" DOUBLE PRECISION;
    ALTER TABLE "sample" ADD COLUMN IF NOT EXISTS "pangolin_ambiguity_score" DOUBLE PRECISION;

    COMMENT ON COLUMN "sample"."pangolin_conflict"
        IS 'In the pangoLEARN decision tree model, a given sequence gets assigned to the most likely category based on known diversity. If a sequence can fit into more than one category, the conflict score will be greater than 0 and reflect the number of categories the sequence could fit into. If the conflict score is 0, this means that within the current decision tree there is only one category that the sequence could be assigned to.';
    COMMENT ON COLUMN "sample"."pangolin_ambiguity_score"
        IS 'This score is a function of the quantity of missing data in a sequence. It represents the proportion of relevant sites in a sequnece which were imputed to the reference values. A score of 1 indicates that no sites were imputed, while a score of 0 indicates that more sites were imputed than were not imputed. This score only includes sites which are used by the decision tree to classify a sequence.';

COMMIT;
