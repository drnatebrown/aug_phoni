#ifndef _SUBSAMPLING_SA_HH
#define _SUBSAMPLING_SA_HH

#include "common.hpp"
#include "sdsl/int_vector.hpp"
#include <malloc_count.h>
#include <r_index.hpp>
#include <string>

void read_samples(std::string filename, sdsl::int_vector<> &samples);

void read_samples(std::string filename, sdsl::int_vector<> &samples) {
  struct stat filestat_samples;
  struct stat filestat_bwt;
  struct stat filestat_bwt_heads;
  FILE *fd;

  std::string base_filename = filename.substr(0, filename.find_last_of("."));

  if ((fd = fopen(filename.c_str(), "r")) == nullptr)
    error("open() file " + filename + " failed");

  int fn = fileno(fd);
  if (fstat(fn, &filestat_samples) < 0)
    error("stat() file " + filename + " failed");

  if (filestat_samples.st_size % SSABYTES != 0)
    error("invilid file " + filename);

  FILE *fd_bwt;

  if ((fd_bwt = fopen((base_filename + ".bwt").c_str(), "r")) == nullptr)
    error("open() file " + base_filename + ".bwt failed");

  int fn_bwt = fileno(fd_bwt);
  if (fstat(fn_bwt, &filestat_bwt) < 0)
    error("stat() file " + base_filename + ".bwt failed");

  auto n = filestat_bwt.st_size;
  int log_n = bitsize(uint64_t(n));

  fclose(fd_bwt);

  FILE *fd_bwt_heads;

  if ((fd_bwt_heads = fopen((base_filename + ".bwt.heads").c_str(), "r")) ==
      nullptr)
    error("open() file " + base_filename + ".bwt.heads failed");

  int fn_bwt_heads = fileno(fd_bwt_heads);
  if (fstat(fn_bwt_heads, &filestat_bwt_heads) < 0)
    error("stat() file " + base_filename + ".bwt.heads failed");

  ulint r = filestat_bwt_heads.st_size;

  fclose(fd_bwt_heads);

  ulint length = filestat_samples.st_size / (2 * SSABYTES);
  // Check that the length of the file is 2*r elements of 5 bytes
  assert(length == r);

  // Create the vector
  samples = int_vector<>(r, 0, log_n);

  // Read the vector
  uint64_t left = 0;
  uint64_t right = 0;
  size_t i = 0;
  while (fread((char *)&left, SSABYTES, 1, fd) &&
         fread((char *)&right, SSABYTES, 1, fd)) {
    ulint val = (right ? right - 1 : n - 1);
    assert(bitsize(uint64_t(val)) <= log_n);
    samples[i++] = val;
  }

  fclose(fd);
}

#endif /* end of include guard: _SUBSAMPLING_SA_HH */
