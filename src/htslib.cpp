#include <Rcpp.h>

#include "htslib/hts.h"
#include "htslib/hts_defs.h"
#include "htslib/regidx.h"
#include "htslib/tbx.h"
#include "htslib/kseq.h"
#include <errno.h>

using namespace Rcpp;

typedef struct
{
    char *regions_fname, *targets_fname;
    int print_header, header_only, cache_megs, download_index, separate_regs;
}
args_t;

static void error(const char *format, ...) {
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  exit(EXIT_FAILURE);
}

static void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) HTS_NORETURN
error_errno(const char *format, ...)
{
    va_list ap;
    int eno = errno;
    fflush(stdout);
    if (format) {
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }
    if (eno) {
        fprintf(stderr, "%s%s\n", format ? ": " : "", strerror(eno));
    } else {
        fprintf(stderr, "\n");
    }
    fflush(stderr);
    exit(EXIT_FAILURE);
}

#define IS_GFF (1 << 0)
#define IS_BED (1 << 1)
#define IS_SAM (1 << 2)
#define IS_VCF (1 << 3)
#define IS_BCF (1 << 4)
#define IS_BAM (1 << 5)
#define IS_CRAM (1 << 6)
#define IS_TXT (IS_GFF | IS_BED | IS_SAM | IS_VCF)


static int query_regions(args_t *args, tbx_conf_t *conf, char *fname, const char *fnidx,
                         char **regs, int nregs, FILE *out_fp) {
  int i;

  htsFile *fp = hts_open(fname, "r");
  if (!fp)
    error_errno("Could not open \"%s\"", fname);
  enum htsExactFormat format = hts_get_format(fp)->format;

  if (args->cache_megs)
    hts_set_cache_size(fp, args->cache_megs * 1048576);

  regidx_t *reg_idx = NULL;
  if (args->targets_fname) {
    reg_idx = regidx_init(args->targets_fname, NULL, NULL, 0, NULL);
    if (!reg_idx)
      error_errno("Could not build region list for \"%s\"",
                  args->targets_fname);
  }

  if (format == vcf || format == sam || format == bed ||
      format == text_format || format == unknown_format) {
    tbx_t *tbx = tbx_index_load3(
        fname, fnidx, args->download_index ? HTS_IDX_SAVE_REMOTE : 0);
    if (!tbx)
      error_errno("Could not load .tbi/.csi index of %s", fname);
    kstring_t str = {0, 0, 0};
    if (args->print_header) {
      int ret;
      while ((ret = hts_getline(fp, KS_SEP_LINE, &str)) >= 0) {
        if (!str.l || str.s[0] != tbx->conf.meta_char)
          break;
        if (fprintf(out_fp, "%s\n", str.s) < 0)
        // if (puts(str.s) < 0)
          error_errno("Error writing to file");
      }
      if (ret < -1)
        error_errno("Reading \"%s\" failed", fname);
    }
    if (!args->header_only) {
      int nseq;
      const char **seq = NULL;
      if (reg_idx) {
        seq = tbx_seqnames(tbx, &nseq);
        if (!seq)
          error_errno("Failed to get sequence names list");
      }
      for (i = 0; i < nregs; i++) {
        int ret, found = 0;
        hts_itr_t *itr = tbx_itr_querys(tbx, regs[i]);
        if (!itr)
          continue;
        while ((ret = tbx_itr_next(fp, tbx, itr, &str)) >= 0) {
          if (reg_idx &&
              !regidx_overlap(reg_idx, seq[itr->curr_tid], itr->curr_beg,
                              itr->curr_end - 1, NULL))
            continue;
          if (!found) {
            if (args->separate_regs)
              printf("%c%s\n", conf->meta_char, regs[i]);
            found = 1;
          }
          if (fprintf(out_fp, "%s\n", str.s) < 0)
          // if (puts(str.s) < 0)
            error_errno("Failed to write to file");
        }
        if (ret < -1)
          error_errno("Reading \"%s\" failed", fname);
        tbx_itr_destroy(itr);
      }
      free(seq);
    }
    free(str.s);
    tbx_destroy(tbx);
  } else if (format == bam)
    error("Please use \"samtools view\" for querying BAM files.\n");

  if (reg_idx)
    regidx_destroy(reg_idx);
  if (hts_close(fp))
    error_errno("hts_close returned non-zero status: %s", fname);

  for (i = 0; i < nregs; i++) {
    free(regs[i]);
  }
    
  free(regs);
  return 0;
}

// [[Rcpp::export]]
void c_read_tabix_table(std::string file_path, StringVector regions, std::string output_file, std::string index_path, bool download_index = false) {
  args_t args;
  memset(&args, 0, sizeof(args_t));
  args.cache_megs = 10;
  args.download_index = download_index ? 1 : 0;
  args.print_header = 1;

  tbx_conf_t conf = tbx_conf_bed;

  // char *reg[10];
  char **reg = (char**)malloc(regions.size() * sizeof(char*));
  for (int i = 0; i < regions.size(); ++i) {
    reg[i] = (char*)malloc((strlen(regions(i)) + 1) * sizeof(char));
    strcpy(reg[i], regions(i));
  }

  char *c_file_path = (char*)malloc(strlen(file_path.c_str()) * sizeof(char));
  strcpy(c_file_path, file_path.c_str());

  FILE *fp = fopen(output_file.c_str(), "wb");
  query_regions(&args, &conf, c_file_path, index_path.empty() ? NULL : index_path.c_str(), reg, regions.size(), fp);
  fclose(fp);
  free(c_file_path);
}
