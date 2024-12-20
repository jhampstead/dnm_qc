#include <htslib/hts.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/vcf.h"

#define NBASES 4
#define BASES "ACGT"
#define UNICODE_SIZE 128

#define MATCH 3
#define MISMATCH -1
#define GAP_SCORE -2
#define NDIRECTIONS 3

#define NSAMPLES 3 // Trio consists of three samples
#define CHILD 0
#define MOTHER 1
#define FATHER 2
#define MAX_STORED 10000 // Some number that should be too high to ever get overflowed

#define ALIGN_SCORE_TAG "ALIGN_SCORE"
#define ALIGN_SCORE_DESC "##INFO=<ID=ALIGN_SCORE,Number=1,Type=Integer,Description=\"Score of alignment\">"
#define ALIGN_THRESH_TAG "ALIGN_THRESH"
#define ALIGN_THRESH_DESC "##INFO=<ID=ALIGN_THRESH,Number=1,Type=Integer,Description=\"Threshold used for alignment\">"

typedef enum direction {
    DIAGONAL = 0,
    VERTICAL,
    HORIZONTAL
} Direction;

typedef struct {
    int is_valid;
    int score;
} Move;

typedef struct {
    uint64_t start, end;
} Range;

typedef enum {
    MIXED = 0,
    PROBAND,
    PARENT
} VariantOrigin;

typedef struct {
    char *ref;
    size_t ref_len;
    char *alt;
    size_t alt_len;
} Variant;

Range get_max_variant_range(bcf1_t *v) {
    bcf_unpack(v, BCF_UN_ALL);
    Range range = { v->pos, v->pos + v->rlen };
    for (int i = 1; i < v->n_allele; i++) {
        const int len = strlen(v->d.allele[i]);
        int new_end = range.start + len;
        if (new_end <= range.end) continue;
        range.end = new_end;
    }

    return range;
}

static inline int ranges_overlap(Range r1, Range r2) {
    return ((r1.start < r2.end) & (r2.start < r1.end)); // Open-ended, so < instead of <=
}

void parse_variant_sample_genotypes(int *gt_array, int n_alts, int gt_samples[]) {
    for (int i = 0; i < NSAMPLES; i++) {
        gt_samples[i] = bcf_gt_allele(gt_array[i * 2]) & bcf_gt_allele(gt_array[i * 2 + 1]);
    }
}

VariantOrigin get_origin_from_genotypes(int gt_samples[]) {
    int proband = gt_samples[0];
    int parent = gt_samples[1] | gt_samples[2];

    if (proband & parent) return MIXED;
    if (proband) return PROBAND;
    return PARENT;
}

Variant normalise_variant(char *ref, size_t ref_len, char *alt, size_t alt_len) {
    // Remove prefix
    while ((ref[0] == alt[0]) & (ref_len > 1) & (alt_len > 1)) {
        ref++;
        ref_len--;

        alt++;
        alt_len--;
    }

    // Remove suffix
    while ((ref[ref_len-1] == alt[alt_len-1]) & (ref_len > 1) & (alt_len > 1)) {
        ref_len--;
        alt_len--;
    }

    char *new_ref = malloc(ref_len + 1);
    strncpy(new_ref, ref, ref_len);
    new_ref[ref_len] = '\0';

    char *new_alt = malloc(alt_len + 1);
    strncpy(new_alt, alt, alt_len);
    new_alt[alt_len] = '\0';

    return (Variant) {ref, ref_len, alt, alt_len};
}

Move get_score_direction(char a1, char a2, int original_scores[NDIRECTIONS], int substitution_matrix[UNICODE_SIZE][UNICODE_SIZE], int gap_score) {
    Move moves[NDIRECTIONS];
    int max_move_score = 0;
    for (int i = 0; i < NDIRECTIONS; i++) {
        int move_score = original_scores[i];

        if (i == DIAGONAL) move_score += substitution_matrix[(int) a1][(int) a2];
        else move_score += gap_score;

        // Keep track of what the best move is
        if (move_score > max_move_score) max_move_score = move_score;

        moves[i] = (Move) { move_score > 0, move_score};
        // fprintf(stdout, "%i:%i ", moves[i].is_valid, moves[i].score);
    }
    // fprintf(stdout, "\n");

    // Only let the highest move be valid;
    for (int i = 0; i < NDIRECTIONS; i++) {
        if (!moves[i].is_valid) continue;
        if (moves[i].score == max_move_score) continue;

        moves[i].is_valid = 0;
    }

    Move best_move = {0}; // Used in case of no valid moves
    for (int i = 0; i < NDIRECTIONS; i++) {
        if (!moves[i].is_valid) continue;

        best_move = moves[i];
        break;
    }

    return best_move;
}

int smith_waterman_score(char *s1, size_t s1_len, char *s2, size_t s2_len, int substitution_matrix[UNICODE_SIZE][UNICODE_SIZE], int gap_score) {
    // Pad with zeros on top & left side
    int n_rows = s1_len+1, n_cols = s2_len+1;
    int scoring_matrix[n_rows][n_cols];

    // Initialize first row to zeros
    for (int i = 0; i < n_rows; i++) scoring_matrix[i][0] = 0;

    // Initialize first column to zeros, skip 0 because it overlaps first row
    for (int i = 1; i < n_cols; i++) scoring_matrix[0][i] = 0;

    // Populate matrix with scores
    int best_score = 0; // Tracks the best alignment score
    for (int i = 1; i < n_rows; i++) {
        char c1 = s1[i-1];
        for (int j = 1; j < n_cols; j++) {
            char c2 = s2[j-1];

            int scores[NDIRECTIONS];
            scores[0] = scoring_matrix[i-1][j-1];
            scores[1] = scoring_matrix[i-1][j];
            scores[2] = scoring_matrix[i][j-1];

            Move move = get_score_direction(c1, c2, scores, substitution_matrix, gap_score);
            // fprintf(stdout, "%i:%i\t%i:%i\n", i, j, move.is_valid, move.score);
            if (!move.is_valid) {
                scoring_matrix[i][j] = 0;
                continue;
            }
            scoring_matrix[i][j] = move.score;

            if (move.score <= best_score) continue;
            best_score = move.score;
        }
    }

    // As we are not interested in the alignment itself, just the score, no traceback necessary!

    return best_score;
}

static inline int score_threshold(size_t ref_len, size_t alt_len) {
    return (ref_len < alt_len ? ref_len : alt_len); // min len
}

int main(int argc, char** argv) {
    if (argc != 3) {
        fprintf(stderr, "Incorrect amount of arguments passed.\n");
        fprintf(stdout, "Usage: ./infer_pileup <vcf_in> <vcf_out>\n");
        exit(1);
    }

    // Initialize substitution matrix, indexing is done based on int value of the characters (hash behaviour)
    int substitution_matrix[UNICODE_SIZE][UNICODE_SIZE];
    for (int i = 0; i < NBASES; i++) for (int j = 0; j < NBASES; j++) substitution_matrix[(int) BASES[i]][(int) BASES[j]] = i == j ? MATCH : MISMATCH;

    // Read in input VCF
    const char *vcf_in_path = argv[1];
    htsFile *fp_in = hts_open(vcf_in_path, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp_in);

    // Add info tags
    int ret = bcf_hdr_append(hdr, ALIGN_SCORE_DESC);
    if (ret < 0) {
        fprintf(stderr, "Error adding new INFO annotation for %s\n", ALIGN_SCORE_TAG);
        bcf_hdr_destroy(hdr);
        hts_close(fp_in);
        return 1;
    }

    ret = bcf_hdr_append(hdr, ALIGN_THRESH_DESC);
    if (ret < 0) {
        fprintf(stderr, "Error adding new INFO annotation for %s\n", ALIGN_THRESH_TAG);
        bcf_hdr_destroy(hdr);
        hts_close(fp_in);
        return 1;
    }

    // Check n of samples to be consistent with trio vcf
    const int nsamples = bcf_hdr_nsamples(hdr); 
    if (nsamples != NSAMPLES) {
        fprintf(stderr, "Incorrect number of samples. This is not a trio vcf.\n");
        exit(1);
    }

    // Prepare output VCF
    const char *vcf_out_path = argv[2];
    htsFile *fp_out = hts_open(vcf_out_path, "w");
    // TODO: Add tags here
    bcf_hdr_write(fp_out, hdr);

    // Process input vcf
    int n_stored = 0;
    bcf1_t *vs[MAX_STORED];

    bcf1_t *v = bcf_init();
    int ngt_arr = 0, *gt_arr = NULL;
    int stored_ngt_arr = 0, *stored_gt_arr = NULL;
    while (bcf_read(fp_in, hdr, v) >= 0) {
        bcf_unpack(v, BCF_UN_ALL);
        Range curr_range = get_max_variant_range(v);
        // fprintf(stdout, "%i:%i\n", v->rid, v->pos);

        // Get genotype information from samples
        bcf_get_genotypes(hdr, v, &gt_arr, &ngt_arr);
        int gt_samples[NSAMPLES];
        parse_variant_sample_genotypes(gt_arr, v->n_allele-1, gt_samples);
        VariantOrigin curr_origin = get_origin_from_genotypes(gt_samples);

        for (int i = n_stored-1; i >= 0; i--) {
            bcf1_t *stored_v = vs[i];
            Range stored_range = get_max_variant_range(stored_v);
            if (ranges_overlap(curr_range, stored_range)) break; // If it still overlaps, don't write it out just yet

            // If it doesn't overlap with anything current, write and remove
            bcf_write(fp_out, hdr, stored_v);
            bcf_destroy(vs[i]);
            vs[n_stored--] = NULL;
        }

        // If nothing is stored and variant should not be considered in a denovo check, just write it out
        // Order of variants prevents us from doing this if there are still stored variants
        if ((curr_origin == MIXED) & (v->n_allele == 2)) {
            if (n_stored == 0) {
                bcf_write(fp_out, hdr, v);
                continue;
            }

            vs[n_stored] = bcf_init();
            bcf_copy(vs[n_stored++], v);
            continue;
        }

        // Split if multiallelic
        Variant variants[v->n_allele-1];
        char *ref = v->d.allele[0];
        size_t ref_len = v->rlen;

        for (int i = 1; i < v->n_allele; i++) {
            char *alt = v->d.allele[i];
            size_t alt_len = strlen(alt);
            variants[i-1] = normalise_variant(ref, ref_len, alt, alt_len);
        }

        // Align multiallelic sites within variant
        if ((v->n_allele > 2) & !((gt_samples[CHILD] == gt_samples[FATHER]) | (gt_samples[CHILD] == gt_samples[MOTHER]))) { // Multiallelic denovo called!
            for (int i = 0; i < v->n_allele-1; i++) {
                Variant source = variants[i];
                for (int j = i+1; j < v->n_allele-1; j++) {
                    Variant target = variants[j];
                    if ((source.ref_len > 1) | (target.ref_len > 1)) {
                        int alignment_score = smith_waterman_score(source.ref, source.ref_len, target.ref, target.ref_len, substitution_matrix, GAP_SCORE);
                        int threshold = score_threshold(source.ref_len, target.ref_len);
                        if (alignment_score > threshold) {
                            if (gt_samples[CHILD] == 0) {
                                gt_arr[1] = 4;
                                bcf_update_genotypes(hdr, v, gt_arr, ngt_arr);
                            }

                            bcf_update_info_int32(hdr, v, ALIGN_SCORE_TAG, &alignment_score, 1);
                            bcf_update_info_int32(hdr, v, ALIGN_THRESH_TAG, &threshold, 1);
                            fprintf(stdout, "Score: %i>%i\t%s>%s\n", alignment_score, threshold, source.ref, target.ref);
                        }
                    }

                    if ((source.alt_len > 1) | (target.alt_len > 1)) {
                        int alignment_score = smith_waterman_score(source.alt, source.alt_len, target.alt, target.alt_len, substitution_matrix, GAP_SCORE);
                        int threshold = score_threshold(source.alt_len, target.alt_len);
                        if (alignment_score > threshold) {
                            if (gt_samples[CHILD] == 0) {
                                gt_arr[1] = 4;
                                bcf_update_genotypes(hdr, v, gt_arr, ngt_arr);
                            }

                            bcf_update_info_int32(hdr, v, ALIGN_SCORE_TAG, &alignment_score, 1);
                            bcf_update_info_int32(hdr, v, ALIGN_THRESH_TAG, &threshold, 1);
                            fprintf(stdout, "Score: %i>%i\t%s>%s\n", alignment_score, threshold, source.alt, target.alt);
                        }
                    }
                }
            }
        }

        // Check for possible merge between variants in probands and parents, to reduce denovos
        // WARN: This may falsely add a variant to a child if both parents have it called differently
        // WARN: This is not a problem for a denovo analysis as these "corrupted" variants will get removed anyway
        for (int i = n_stored-1; i >= 0; i--) {
            bcf1_t *stored_v = vs[i];
            bcf_unpack(stored_v, BCF_UN_ALL);

            // Get genotype information
            bcf_get_genotypes(hdr, v, &stored_gt_arr, &stored_ngt_arr);
            int stored_gt_samples[NSAMPLES];
            parse_variant_sample_genotypes(gt_arr, v->n_allele-1, stored_gt_samples);

            // If GTs child == father or child == mother, it is not denovo in child
            if (
                !((gt_samples[CHILD] == gt_samples[FATHER]) | (gt_samples[CHILD] == gt_samples[MOTHER])) |
                !((stored_gt_samples[CHILD] == stored_gt_samples[FATHER]) | (gt_samples[CHILD] == gt_samples[MOTHER]))
            ) {
                Variant stored_variants[stored_v->n_allele-1];
                char *ref = stored_v->d.allele[0];
                size_t ref_len = stored_v->rlen;

                for (int i = 1; i < stored_v->n_allele; i++) {
                    char *alt = stored_v->d.allele[i];
                    size_t alt_len = strlen(alt);
                    stored_variants[i-1] = normalise_variant(ref, ref_len, alt, alt_len);
                }

                for (int i = 0; i < v->n_allele-1; i++) {
                    Variant source = variants[i];
                    for (int j = 0; j < stored_v->n_allele-1; j++) {
                        Variant target = stored_variants[j];

                        if ((source.ref_len > 1) & (target.ref_len > 1)) {
                            int alignment_score = smith_waterman_score(source.ref, source.ref_len, target.ref, target.ref_len, substitution_matrix, GAP_SCORE);
                            int threshold = score_threshold(source.ref_len, target.ref_len);
                            if (alignment_score > threshold) {
                                if (gt_samples[CHILD] == 0) {
                                    gt_arr[1] = 4;
                                    bcf_update_genotypes(hdr, v, gt_arr, ngt_arr);
                                    bcf_update_info_int32(hdr, v, ALIGN_SCORE_TAG, &alignment_score, 1);
                                    bcf_update_info_int32(hdr, v, ALIGN_THRESH_TAG, &threshold, 1);
                                }

                                if (stored_gt_samples[CHILD] == 0) {
                                    gt_arr[1] = 4;
                                    bcf_update_genotypes(hdr, stored_v, gt_arr, ngt_arr);
                                    bcf_update_info_int32(hdr, stored_v, ALIGN_SCORE_TAG, &alignment_score, 1);
                                    bcf_update_info_int32(hdr, stored_v, ALIGN_THRESH_TAG, &threshold, 1);

                                }
                                fprintf(stdout, "Score: %i>%i\t%s>%s\n", alignment_score, threshold, source.ref, target.ref);
                            }
                        }

                        if ((source.alt_len > 1) & (target.alt_len > 1)) {
                            int alignment_score = smith_waterman_score(source.alt, source.alt_len, target.alt, target.alt_len, substitution_matrix, GAP_SCORE);
                            int threshold = score_threshold(source.alt_len, target.alt_len);
                            if (alignment_score > threshold) {
                                if (gt_samples[CHILD] == 0) {
                                    gt_arr[1] = 4;
                                    bcf_update_genotypes(hdr, v, gt_arr, ngt_arr);
                                    bcf_update_info_int32(hdr, v, ALIGN_SCORE_TAG, &alignment_score, 1);
                                    bcf_update_info_int32(hdr, v, ALIGN_THRESH_TAG, &threshold, 1);

                                }

                                if (stored_gt_samples[CHILD] == 0) {
                                    gt_arr[1] = 4;
                                    bcf_update_genotypes(hdr, stored_v, gt_arr, ngt_arr);
                                    bcf_update_info_int32(hdr, stored_v, ALIGN_SCORE_TAG, &alignment_score, 1);
                                    bcf_update_info_int32(hdr, stored_v, ALIGN_THRESH_TAG, &threshold, 1);

                                }
                                fprintf(stdout, "Score: %i>%i\t%s>%s\n", alignment_score, threshold, source.alt, target.alt);
                            }
                        }
                    }
                }
            }
        }

        // Add to storage
        if (n_stored) memmove(&vs[1], &vs[0], n_stored * sizeof(bcf1_t *));
        vs[0] = bcf_init();
        bcf_copy(vs[0], v);
        n_stored++;
    }

    bcf_destroy(v);

    // Write out remaining stored records
    for (int i = 0; i < n_stored; i++) {
        bcf_write(fp_out, hdr, vs[i]);
        bcf_destroy(vs[i]);
    }

    hts_close(fp_in);
    hts_close(fp_out);
}

// int main(int argc, char **argv) {
//     int substitution_matrix[UNICODE_SIZE][UNICODE_SIZE];
//     for (int i = 0; i < NBASES; i++) for (int j = 0; j < NBASES; j++) substitution_matrix[(int) BASES[i]][(int) BASES[j]] = i == j ? MATCH : MISMATCH;
//
//     char *s1 = "TGCA";
//     char *s2 = "CGAATA";
//
//     int alignment_score = smith_waterman_score(s1, 4, s2, 6, substitution_matrix, GAP_SCORE);
//     
//     fprintf(stdout, "%i\n", alignment_score);
// }

