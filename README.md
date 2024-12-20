# dnm_qc
Executables and resource files for DNM meta-analysis filtering and QC

## Compiling from source
The C code in `src/` can be compiled from source using `cmake --build build`. Htslib is required to compile and execute. If htslib isn't installed on your system, you can download it [here](https://www.htslib.org/download/) or [install from source](https://github.com/samtools/htslib). By default, `cmake` looks for htslib in standard directories `/usr/include` and `/usr/local/include`. If htslib is installed elsewhere on your system, you can tell `cmake` where to look for it by specifying `include_directories(/path/to/external/lib/include)` within `CMakeLists.txt`.

Building should generate the following executables in `build/`:
* infer_pileup
* find_adjacent_variants
* annotate_bed_regions

## infer_pileup
`infer_pileup` is designed to address two primary concerns: 
1. DNMs that occur at multiallelic sites/in regions with many overlapping indel calls in proband, and
2. Indels that are aligned differently in proband and parent because trios were not joint called.

The program performs local Smith-Waterman alignment on overlapping indels and generates an alignment score. If the alignment score is greater than a threshold (set to the minimum length of the reference/alternate allele by default), the genotype of the putative parent-of-origin is altered to reflect a heterozygous variant call so the variant is no longer called *de novo* by basic genotype filtering strategies. There are some issues with this approach generally, but as a quick and dirty way to get rid of these dodgy indels it does the job.

![Depiction of genotype modification in proband](https://i.imgur.com/R6IZds7.png)

`infer_pileup` will print variants with a score greater than the threshold to stdout in the following format during runtime:
```
Score: 87>29	TATATATATATATATATATATATATATAT>CATATATATATATATATATATATATATATATAT
Score: 27>9	TACACACAC>TACACACACAC
Score: 15>5	CACAC>TACACACACACACAC
Score: 3>1	TTATA>A
```

`infer_pileup` expects two command line arguments: an input VCF `<vcf_in>` and an output VCF `<vcf_out>`. `<vcf_in>` is expected to be a merged trio VCF containing the genotypes of the proband and both parents, where the proband is expected to be GT[0]. These VCFs can be compressed or uncompressed, indexed or unindexed. Multiallelic sites can be split or unsplit. For testing, merged trio VCFs were generated using bcftools:

```
bcftools merge --missing-to-ref proband.vcf mother.vcf father.vcf -Oz -o output.vcf]
```

## find_adjacent_variants

`find_adjacent_variants` is designed to determine whether multiple a) indels, b) SNVs, or c) both exist within a specified radius. The program annotates these variants with a FORMAT/MNV tag (although many are not actually MNVs...) = 1, so they can be easily filtered using `bcftools`:

```
bcftools view -e 'FORMAT/MNV=1' input.vcf -o output.vcf
bcftools filter -e 'FORMAT/MNV=1' input.vcf -o output.vcf
```
Variants that do not have other variants co-occuring within the specified radius will be given a FORMAT/MNV tag = 0.

`find_adjacent_variants` expects the following command line arguments:
* An input VCF file,
* A variant type to consider (one of 'snv', 'indel', or 'both'),
* A radius in which to look for adjacent variants,
* An output VCF file

Similar to previous, VCFs can be compressed or uncompressed, indexed or unindexed. Multiallelic sites can be split or unsplit. `<vcf_in>` can be a single- or multi-sample VCF file, but note the program will consider co-occurring variants across all samples in the VCF file (so if the mother has an indel called at position 5 and the proband has an indel called at position 10, these calls will both have FORMAT/MNV = 1).

## annotate_bed_regions

Petr has provided a variety of region files we can consider excluding calls from because they tend to be noisy or otherwise poor in short read data:
* repeat-masker.b37.30pct-div.bed.gz
* paralogs.single.sim70.b37.txt.gz
* segdups.b37.txt.gz

While `bcftools view` and `bcftools filter` can both filter a VCF to .bed regions, sometimes I want to see the impact this filtering would have without explicitly removing variants from an input VCF file. To facilitate this, I wrote `annotate_bed_regions`. This program takes an input VCF file, .bed file, and INFO tag name, and annotates the output vcf with INFO/<TAG_NAME> specifying whether the variant occurs within a .bed region (= 1) or does not (= 0). These can then be easily filtered later using `bcftools`:

```
bcftools view -e 'INFO/<TAG_NAME>=1' input.vcf -o output.vcf
bcftools filter -e 'INFO/<TAG_NAME>=1' input.vcf -o output.vcf
```
As previous, VCF files can be compressed or uncompressed, indexed or unindexed. 


