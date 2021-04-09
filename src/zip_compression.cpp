#include <Rcpp.h>
#include <stdio.h>
#include <errno.h>
#include <zlib.h>
#include <fcntl.h>
#include "htslib/bgzf.h"
#include "htslib/tbx.h"

void _zip_error(const char *txt, const char *err, int infd, int outfd)
{
    close(infd);
    close(outfd);
    err ? Rf_error(txt, err) : Rf_error(txt);
}

void _zip_open(std::string input_file_path, std::string output_file_path, int *infd, int *outfd, bool append)
{
    int iflag = O_RDONLY, oflag = O_WRONLY | O_CREAT;
#ifdef _WIN32
    iflag |= O_BINARY;
    oflag |= O_BINARY;
#endif

    if (append)
        oflag |= O_APPEND;
    else
        oflag |= O_TRUNC;

    *infd = open(input_file_path.c_str(), iflag);
    if (0 > *infd)
        Rf_error("opening 'input_file_path': %s", strerror(errno));

    /* we overwrite existing files here */
    *outfd = open(output_file_path.c_str(), oflag, 0666);
    if (0 > *outfd) {
        close(*infd);
        Rf_error("opening 'output_file_path': %s", strerror(errno));
    }
}

// [[Rcpp::export]]
void bgzip(std::string input_file_path, std::string output_file_path, bool append = false)
{
    static const int BUF_SIZE = 64 * 1024;
    void *buffer;
    int infd, outfd, cnt;
    gzFile in;
    BGZF *outp;

    buffer = R_alloc(BUF_SIZE, sizeof(void *));

    _zip_open(input_file_path, output_file_path, &infd, &outfd, append);
    in = gzdopen(infd, "rb");
    if (NULL == in)
        _zip_error("opening input 'input_file_path'", NULL, infd, outfd);
    outp = bgzf_dopen(outfd, "w");
    if (NULL == outp)
        _zip_error("opening output 'output_file_path'", NULL, infd, outfd);

    while (0 < (cnt = gzread(in, buffer, BUF_SIZE)))
        if (0 > bgzf_write(outp, buffer, cnt))
            _zip_error("writing compressed output", NULL, infd, outfd);
    if (0 > cnt)
        _zip_error("reading compressed input: %s",
                   strerror(errno), infd, outfd);

    if (0 > bgzf_close(outp))
        Rf_error("closing compressed output");
    if (gzclose(in) != Z_OK)
        _zip_error("closing input after compression", NULL, infd, outfd);
}


// [[Rcpp::export]]
void build_tabix_index(std::string file_path) {
    tbx_conf_t conf = tbx_conf_bed;

    if (bgzf_is_bgzf(file_path.c_str()) != 1)
        Rf_error("file does not appear to be bgzip'd");
    if (tbx_index_build(file_path.c_str(), 0, &conf) == -1)
        Rf_error("index build failed");
}