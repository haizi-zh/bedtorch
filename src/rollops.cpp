#include <Rcpp.h>
#include <string.h>
#include <assert.h>
#include <math.h>
using namespace Rcpp;
using namespace std;

#define ALIGN_CENTER 1
#define ALIGN_LEFT 2
#define ALIGN_RIGHT 3

#define ROLLSUM 1
#define ROLLMEAN 2

NumericVector rollapply(NumericVector x, int k, double (*func)(NumericVector, const int, const int, const int, const bool, double*, int*), bool na_pad = false, bool na_rm = false, int align = 1) {
  int mode = align;
  assert(mode > 0);

  // Initialize output vector
  int n = x.size() - k + 1;
  if (n <= 0) {
    return NULL;
  }
  if (k >= x.size()) return NULL;

  int out_len = n;
  if (na_pad) {
    out_len = x.size();
  }

  NumericVector out(out_len);
  
  double *buf = (double *)malloc(k * sizeof(double));
  int buf_ptr = -1;

  // Based on mode, determine the start and end index
  int start = -1;
  int end = -1;
  // Sliding window start offset
  int offset = -1;
  // int shift_r = -1;

  if (mode == ALIGN_LEFT) {
    // align: left
    start = 0;
  } else if (mode == ALIGN_CENTER) {
    // align: center
    start = (int)floor((k - 1) / 2);
  } else if (mode == ALIGN_RIGHT) {
    // align: right
    start = k - 1;
  } else {
    assert(false);
  }
  end = start + n;
  offset = -start;

  for (int i = start; i < end; ++i) {
    // double v = func(x, i + offset, k, na_rm);
    double v = func(x, i - start, i + offset, k, na_rm, buf, &buf_ptr);
    if (na_pad)
      out[i] = v;
    else
      out[i - start] = v;
  }

  free(buf);

  if (na_pad) {
    // Set padding areas to NA
    for (int i = 0; i < start; ++i) {
      out[i] = NA_REAL;
    }
    for (int i = end; i < out_len; ++i) {
      out[i] = NA_REAL;
    }
  }

  return out;
}

void print_buffer(double *buf, int *buf_ptr, int buf_len) {
  printf("Buffer pointer: %d\n", *buf_ptr);
    for (int i = 0; i < buf_len; ++i) {
      printf("%f\t", buf[i]);
    }
  printf("\n");
}

double calc_sum(NumericVector x, const int idx, const int pos, const int k, const bool na_rm, 
                double *buf, int *buf_ptr)
{
  // mode: ROLLSUM or ROLLMEAN
  // buf_ptr: point to the current available slot
  static double total;
  // static int count;
  // How many NAs in the window?
  static int na_count;

  if (idx == 0)
  {
    // Initialization
    total = 0;
    na_count = 0;
    // count = 0;

    // Mimic the state right before a standard loop
    // First element is the dummy. Will be poped soon
    buf[0] = 0;
    // count += 1;
    for (int i = 1; i < k; ++i)
    {
      const double v = x[pos + i - 1];
      if (NumericVector::is_na(v))
      {
        na_count += 1;
      }
      else
      {
        // count += 1;
        total += v;
      }
      buf[i] = v;
    }

    // Now in buf:
    // 0    x[pos],   x[pos+1], ...   x[pos + k - 2]
    // First element: 0 (dummy)
    // Last element points to pos + k - 2. Now we're reay to update x[pos + k -1] at buf[0]
    *buf_ptr = 0;

    // printf("Initialization");
    // print_buffer(buf, buf_ptr, k);
  }

  // Updating. pop/push
  // For the current sliding window, fetch the last element, and append to the end of the buffer
  // Then pop the first element of the buffer
  double popped = buf[*buf_ptr];
  const double pushed = x[pos + k - 1];
  buf[*buf_ptr] = pushed;
  *buf_ptr += 1;
  // Cycle back
  if (*buf_ptr >= k)
    *buf_ptr = 0;

  if (NumericVector::is_na(pushed)) {
    na_count += 1;
  } 
  
  if (NumericVector::is_na(popped)) {
    na_count -= 1;
  }

  // Remove NAs and do the calculation
  if (!NumericVector::is_na(pushed))
    total += pushed;

  if (!NumericVector::is_na(popped))
    total -= popped;

  // printf("#%d Pushed: %f, popped: %f, total: %f, # of NAs: %d, count: %d\n", idx, pushed, popped, total, na_count, count);
  // print_buffer(buf, buf_ptr, k);

  if (na_count > 0 && !na_rm)
    return NA_REAL;
  else {
    return total;
  } 
}


double calc_mean(NumericVector x, const int idx, const int pos, const int k, const bool na_rm, 
                double *buf, int *buf_ptr)
{
  // mode: ROLLSUM or ROLLMEAN
  // buf_ptr: point to the current available slot
  static double total;
  static int count;
  // How many NAs in the window?
  static int na_count;

  if (idx == 0)
  {
    // Initialization
    total = 0;
    na_count = 0;
    count = 0;

    // Mimic the state right before a standard loop
    // First element is the dummy. Will be poped soon
    buf[0] = 0;
    count += 1;
    for (int i = 1; i < k; ++i)
    {
      const double v = x[pos + i - 1];
      if (NumericVector::is_na(v))
      {
        na_count += 1;
      }
      else
      {
        count += 1;
        total += v;
      }
      buf[i] = v;
    }

    // Now in buf:
    // 0    x[pos],   x[pos+1], ...   x[pos + k - 2]
    // First element: 0 (dummy)
    // Last element points to pos + k - 2. Now we're reay to update x[pos + k -1] at buf[0]
    *buf_ptr = 0;

    // printf("Initialization");
    // print_buffer(buf, buf_ptr, k);
  }

  // Updating. pop/push
  // For the current sliding window, fetch the last element, and append to the end of the buffer
  // Then pop the first element of the buffer
  double popped = buf[*buf_ptr];
  const double pushed = x[pos + k - 1];
  buf[*buf_ptr] = pushed;
  *buf_ptr += 1;
  // Cycle back
  if (*buf_ptr >= k)
    *buf_ptr = 0;

  if (NumericVector::is_na(pushed)) {
    na_count += 1;
  } else {
    count += 1;
  }
  if (NumericVector::is_na(popped)) {
    na_count -= 1;
  } else {
    count -= 1;
  }

  // Remove NAs and do the calculation
  if (!NumericVector::is_na(pushed))
    total += pushed;

  if (!NumericVector::is_na(popped))
    total -= popped;

  // printf("#%d Pushed: %f, popped: %f, total: %f, # of NAs: %d, count: %d\n", idx, pushed, popped, total, na_count, count);
  // print_buffer(buf, buf_ptr, k);

  if (na_count > 0 && !na_rm)
    return NA_REAL;
  else {
    return total / count;
  } 
}


// double calc_sum(NumericVector x, const int idx, const int pos, const int k, const bool na_rm,
//                 double *buf, int *buf_ptr)
// {
//   return calc_sum_mean(x, idx, pos, k, na_rm, buf, buf_ptr, ROLLSUM);
// }

// double calc_mean(NumericVector x, const int idx, const int pos, const int k, const bool na_rm,
//                 double *buf, int *buf_ptr)
// {
//   return calc_sum_mean(x, idx, pos, k, na_rm, buf, buf_ptr, ROLLMEAN);
// }

// [[Rcpp::export]]
NumericVector c_rollmean(NumericVector x, int k, bool na_pad = false, bool na_rm = false, int align = 1) {
  return rollapply(x, k, &calc_mean, na_pad, na_rm, align);
}


// [[Rcpp::export]]
NumericVector c_rollsum(NumericVector x, int k, bool na_pad = false, bool na_rm = false, int align = 1) {
  return rollapply(x, k, &calc_sum, na_pad, na_rm, align);
}

