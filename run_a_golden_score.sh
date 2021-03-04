# TODO: replace with your dataset prefix
PLINK_PREFIX="/home/mvn373/fast/prs-repo/pre-calculated/genetic/TARGET_final"
GOLDEN_SCORE="/home/mvn373/fast/prs-repo/golden/01-weights/BMI-standardized.weights.Inter99-derived.LDpred_p3.0000e-01.txt"
OUTPUT_PREFIX="BMIPRS"

# [variant ID col.] [allele col.] [score col.]
plink --bfile "${PLINK_PREFIX}" --score "${GOLDEN_SCORE}" 1 4 6 sum header --out "${OUTPUT_PREFIX}"
