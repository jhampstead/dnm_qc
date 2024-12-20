#include <stdio.h>
#include <stdlib.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

// Function to add INFO field to VCF header
void add_info_field(bcf_hdr_t *vcf_hdr, const char *info_tag) {
    int info_id = bcf_hdr_id2int(vcf_hdr, BCF_DT_ID, info_tag);
    if (info_id < 0) {
        int buffer_size = snprintf(NULL, 0, "##INFO=<ID=%s,Number=1,Type=Integer,Description=\"Flag indicating if variant is within BED region specified in tag\">\n", info_tag) + 1; // +1 for the null terminator
        char *hdr_str = (char *)malloc(buffer_size);
        sprintf(hdr_str, "##INFO=<ID=%s,Number=1,Type=Integer,Description=\"Flag indicating if variant is within BED region specified in tag\">\n", info_tag);
        bcf_hdr_append(vcf_hdr, hdr_str);
        free(hdr_str);
    }
}

void annotate_vcf(const char *vcf_file, const char *bed_file, tbx_t *bed_idx, const char *info_tag, const char *output_file) {
    htsFile *vcf_fp = bcf_open(vcf_file, "r");
    bcf_hdr_t *vcf_hdr = bcf_hdr_read(vcf_fp);

    htsFile *bed_fp = hts_open(bed_file, "r");

    // Add the INFO field to VCF header
    add_info_field(vcf_hdr, info_tag);

    htsFile *out_fp = bcf_open(output_file, "w");
    bcf_hdr_write(out_fp, vcf_hdr);

    bcf1_t *vcf_rec = bcf_init();
    hts_itr_t *itr;
    while (bcf_read(vcf_fp, vcf_hdr, vcf_rec) == 0) {
        const char *chrom = bcf_hdr_id2name(vcf_hdr, vcf_rec->rid);
        int pos = vcf_rec->pos + 1;  // VCF positions are 0-based, BED positions are 1-based
        itr = tbx_itr_queryi(bed_idx, bcf_hdr_name2id(vcf_hdr, chrom), pos - 1, pos);
        if (itr == NULL) {
            fprintf(stderr, "Failed to query region %s:%d-%d\n", chrom, pos - 1, pos);
            exit(1);
        }

        kstring_t s = {0};
        // Annotate VCF record
        int annot = 0;
        if (tbx_itr_next(bed_fp, bed_idx, itr, &s) >= 0) {
            annot = 1;
            bcf_update_info_int32(vcf_hdr, vcf_rec, info_tag, &annot, 1);
        } else {
            bcf_update_info_int32(vcf_hdr, vcf_rec, info_tag, &annot, 1);
        }
        bcf_write(out_fp, vcf_hdr, vcf_rec);

        free(s.s);
    }
    hts_itr_destroy(itr);

    tbx_destroy(bed_idx);
    bcf_destroy(vcf_rec);
    bcf_hdr_destroy(vcf_hdr);
    bcf_close(vcf_fp);
    hts_close(bed_fp);
    bcf_close(out_fp);

}

int main(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <INFO_TAG> <input.bed> <input.vcf> <output.vcf>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char *info_tag = argv[1];
    const char *bed_file = argv[2];
    const char *vcf_file = argv[3];
    const char *output_file = argv[4];

    const char *tbx_ext = ".tbi";
    size_t len = strlen(bed_file) + strlen(tbx_ext) + 1;
    char *idx = (char*)malloc(len);
    strcpy(idx, bed_file);
    strcat(idx, tbx_ext);

    tbx_t *bed_idx = tbx_index_load(idx);

    annotate_vcf(vcf_file, bed_file, bed_idx, info_tag, output_file);

    return EXIT_SUCCESS;
}

