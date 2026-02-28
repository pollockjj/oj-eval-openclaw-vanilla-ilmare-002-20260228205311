#include "int2048.h"

namespace sjtu {

// Naive multiplication
static std::vector<long long> mul_naive(const std::vector<long long>& a, const std::vector<long long>& b) {
  if (a.empty() || b.empty()) return {0};
  if (a.size() == 1 && a[0] == 0) return {0};
  if (b.size() == 1 && b[0] == 0) return {0};
  if (a.size() == 1) return mul_small(b, a[0]);
  if (b.size() == 1) return mul_small(a, b[0]);
  
  std::vector<long long> res(a.size() + b.size() + 1, 0);
  for (size_t i = 0; i < a.size(); i++) {
    long long carry = 0;
    for (size_t j = 0; j < b.size(); j++) {
      __int128 prod = (__int128)a[i] * b[j] + res[i + j] + carry;
      res[i + j] = prod % BASE;
      carry = prod / BASE;
    }
    if (carry) res[i + b.size()] += carry;
  }
  trim_vec(res);
  return res;
}

// Karatsuba multiplication
static std::vector<long long> karatsuba(const std::vector<long long>& a, const std::vector<long long>& b) {
  size_t n = std::max(a.size(), b.size());
  
  if (n <= 64) return mul_naive(a, b);
  
  size_t m = n / 2;
  
  std::vector<long long> a0(a.begin(), a.begin() + std::min(m, a.size()));
  std::vector<long long> a1;
  if (a.size() > m) a1.assign(a.begin() + m, a.end());
  if (a0.empty()) a0 = {0};
  if (a1.empty()) a1 = {0};
  trim_vec(a0);
  trim_vec(a1);
  
  std::vector<long long> b0(b.begin(), b.begin() + std::min(m, b.size()));
  std::vector<long long> b1;
  if (b.size() > m) b1.assign(b.begin() + m, b.end());
  if (b0.empty()) b0 = {0};
  if (b1.empty()) b1 = {0};
  trim_vec(b0);
  trim_vec(b1);
  
  std::vector<long long> z0 = karatsuba(a0, b0);
  std::vector<long long> z2 = karatsuba(a1, b1);
  
  std::vector<long long> sum_a = add_abs(a0, a1);
  std::vector<long long> sum_b = add_abs(b0, b1);
  std::vector<long long> z1 = karatsuba(sum_a, sum_b);
  z1 = sub_abs(z1, z0);
  z1 = sub_abs(z1, z2);
  
  std::vector<long long> res(std::max({z0.size() + 0UL, z1.size() + m, z2.size() + 2*m}) + 1, 0);
  
  for (size_t i = 0; i < z0.size(); i++) res[i] += z0[i];
  
  long long carry = 0;
  for (size_t i = 0; i <= z1.size() || carry; i++) {
    __int128 sum = carry;
    if (m + i < res.size()) {
      sum += res[m + i];
      if (i < z1.size()) sum += z1[i];
      res[m + i] = sum % BASE;
      carry = sum / BASE;
    } else {
      break;
    }
  }
  
  carry = 0;
  for (size_t i = 0; i <= z2.size() || carry; i++) {
    __int128 sum = carry;
    if (2*m + i < res.size()) {
      sum += res[2*m + i];
      if (i < z2.size()) sum += z2[i];
      res[2*m + i] = sum % BASE;
      carry = sum / BASE;
    } else {
      break;
    }
  }
  
  trim_vec(res);
  return res;
}

// Multiply absolute values
std::vector<long long> mul_abs(const std::vector<long long>& a, const std::vector<long long>& b) {
  return karatsuba(a, b);
}

// Divide absolute values, returns quotient and remainder
std::vector<long long> div_abs(const std::vector<long long>& a, const std::vector<long long>& b, std::vector<long long>& rem) {
  if (compare_abs(a, b) < 0) {
    rem = a;
    return {0};
  }
  
  std::vector<long long> q(a.size() - b.size() + 1, 0);
  rem = a;
  
  for (int i = (int)q.size() - 1; i >= 0; i--) {
    long long lo = 0, hi = BASE - 1;
    while (lo < hi) {
      long long mid = (lo + hi + 1) / 2;
      std::vector<long long> prod = mul_small(b, mid);
      prod.insert(prod.begin(), i, 0);
      
      if (compare_abs(prod, rem) <= 0) {
        lo = mid;
      } else {
        hi = mid - 1;
      }
    }
    q[i] = lo;
    
    if (lo > 0) {
      std::vector<long long> prod = mul_small(b, lo);
      prod.insert(prod.begin(), i, 0);
      rem = sub_abs(rem, prod);
    }
  }
  
  trim_vec(q);
  trim_vec(rem);
  return q;
}

} // namespace sjtu
