// Integer 1:
// Implement a signed big integer class that only needs to support simple addition and subtraction

// Integer 2:
// Implement a signed big integer class that supports addition, subtraction, multiplication, and division, and overload related operators

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>

// Do not use "using namespace std;"

namespace sjtu {

// Constants
static const long long BASE = 1000000000LL;  // 10^9
static const int BASE_DIGITS = 9;

// Helper functions at namespace level
inline void trim_vec(std::vector<long long>& v) {
  while (v.size() > 1 && v.back() == 0) v.pop_back();
}

inline int compare_abs(const std::vector<long long>& a, const std::vector<long long>& b) {
  if (a.size() != b.size()) return (a.size() > b.size()) ? 1 : -1;
  for (int i = (int)a.size() - 1; i >= 0; i--) {
    if (a[i] != b[i]) return (a[i] > b[i]) ? 1 : -1;
  }
  return 0;
}

inline std::vector<long long> add_abs(const std::vector<long long>& a, const std::vector<long long>& b) {
  std::vector<long long> res(std::max(a.size(), b.size()) + 1, 0);
  long long carry = 0;
  for (size_t i = 0; i < res.size(); i++) {
    long long sum = carry;
    if (i < a.size()) sum += a[i];
    if (i < b.size()) sum += b[i];
    res[i] = sum % BASE;
    carry = sum / BASE;
  }
  trim_vec(res);
  return res;
}

inline std::vector<long long> sub_abs(const std::vector<long long>& a, const std::vector<long long>& b) {
  std::vector<long long> res(a.size(), 0);
  long long borrow = 0;
  for (size_t i = 0; i < a.size(); i++) {
    long long diff = a[i] - borrow - (i < b.size() ? b[i] : 0);
    if (diff < 0) {
      diff += BASE;
      borrow = 1;
    } else {
      borrow = 0;
    }
    res[i] = diff;
  }
  trim_vec(res);
  return res;
}

inline std::vector<long long> mul_small(const std::vector<long long>& a, long long b) {
  if (b == 0) return {0};
  std::vector<long long> res(a.size() + 2, 0);
  long long carry = 0;
  for (size_t i = 0; i < a.size(); i++) {
    __int128 prod = (__int128)a[i] * b + carry;
    res[i] = prod % BASE;
    carry = prod / BASE;
  }
  size_t idx = a.size();
  while (carry) {
    res[idx++] = carry % BASE;
    carry /= BASE;
  }
  trim_vec(res);
  return res;
}

// Naive multiplication for small numbers
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
  
  if (n <= 32) return mul_naive(a, b);
  
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

class int2048 {
private:
  std::vector<long long> digits_;  // Little-endian (least significant first)
  bool negative_;  // true if negative

  void trim_() {
    trim_vec(digits_);
    if (digits_.size() == 1 && digits_[0] == 0) negative_ = false;
  }

public:
  // Constructors
  int2048() : digits_(1, 0), negative_(false) {}
  
  int2048(long long x) : negative_(x < 0) {
    if (x < 0) x = -x;
    if (x == 0) {
      digits_ = {0};
    } else {
      digits_.clear();
      while (x > 0) {
        digits_.push_back(x % BASE);
        x /= BASE;
      }
    }
  }
  
  int2048(const std::string& s) : digits_(1, 0), negative_(false) {
    read(s);
  }
  
  int2048(const int2048& other) : digits_(other.digits_), negative_(other.negative_) {}

  // The parameter types of the following functions are for reference only, you can choose to use constant references or not
  // If needed, you can add other required functions yourself
  // ===================================
  // Integer1
  // ===================================

  // Read a big integer
  void read(const std::string& s) {
    digits_.clear();
    negative_ = false;
    
    size_t start = 0;
    if (!s.empty() && s[0] == '-') {
      negative_ = true;
      start = 1;
    } else if (!s.empty() && s[0] == '+') {
      start = 1;
    }
    
    while (start < s.size() - 1 && s[start] == '0') start++;
    
    int len = s.size() - start;
    digits_.reserve((len + BASE_DIGITS - 1) / BASE_DIGITS);
    
    for (int i = len; i > 0; i -= BASE_DIGITS) {
      int st = std::max(0, i - BASE_DIGITS);
      long long num = 0;
      for (int j = st; j < i; j++) {
        num = num * 10 + (s[start + j] - '0');
      }
      digits_.push_back(num);
    }
    
    trim_();
  }
  
  // Output the stored big integer, no need for newline
  void print() const {
    if (negative_ && !(digits_.size() == 1 && digits_[0] == 0)) {
      std::putchar('-');
    }
    std::printf("%lld", digits_.back());
    for (int i = (int)digits_.size() - 2; i >= 0; i--) {
      std::printf("%09lld", digits_[i]);
    }
  }

  // Add a big integer
  int2048& add(const int2048& other) {
    if (negative_ == other.negative_) {
      digits_ = add_abs(digits_, other.digits_);
    } else {
      int cmp = compare_abs(digits_, other.digits_);
      if (cmp >= 0) {
        digits_ = sub_abs(digits_, other.digits_);
      } else {
        digits_ = sub_abs(other.digits_, digits_);
        negative_ = other.negative_;
      }
    }
    trim_();
    return *this;
  }
  
  // Return the sum of two big integers
  friend int2048 add(int2048 a, const int2048& b) {
    return a.add(b);
  }

  // Subtract a big integer
  int2048& minus(const int2048& other) {
    if (negative_ != other.negative_) {
      digits_ = add_abs(digits_, other.digits_);
    } else {
      int cmp = compare_abs(digits_, other.digits_);
      if (cmp >= 0) {
        digits_ = sub_abs(digits_, other.digits_);
      } else {
        digits_ = sub_abs(other.digits_, digits_);
        negative_ = !negative_;
      }
    }
    trim_();
    return *this;
  }
  
  // Return the difference of two big integers
  friend int2048 minus(int2048 a, const int2048& b) {
    return a.minus(b);
  }

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const {
    return *this;
  }
  
  int2048 operator-() const {
    int2048 res(*this);
    if (!(res.digits_.size() == 1 && res.digits_[0] == 0)) {
      res.negative_ = !res.negative_;
    }
    return res;
  }

  int2048& operator=(const int2048& other) {
    if (this != &other) {
      digits_ = other.digits_;
      negative_ = other.negative_;
    }
    return *this;
  }

  int2048& operator+=(const int2048& other) {
    return add(other);
  }
  
  friend int2048 operator+(int2048 a, const int2048& b) {
    return a += b;
  }

  int2048& operator-=(const int2048& other) {
    return minus(other);
  }
  
  friend int2048 operator-(int2048 a, const int2048& b) {
    return a -= b;
  }

  int2048& operator*=(const int2048& other) {
    bool result_negative = negative_ != other.negative_;
    digits_ = mul_abs(digits_, other.digits_);
    negative_ = result_negative;
    trim_();
    return *this;
  }
  
  friend int2048 operator*(int2048 a, const int2048& b) {
    return a *= b;
  }

  int2048& operator/=(const int2048& other) {
    std::vector<long long> rem;
    digits_ = div_abs(digits_, other.digits_, rem);
    
    bool result_negative = negative_ != other.negative_;
    
    // Floor division: if signs differ and has remainder, add 1 to absolute value
    if (result_negative && !(rem.size() == 1 && rem[0] == 0)) {
      std::vector<long long> one = {1};
      digits_ = add_abs(digits_, one);
    }
    
    negative_ = result_negative;
    trim_();
    return *this;
  }
  
  friend int2048 operator/(int2048 a, const int2048& b) {
    return a /= b;
  }

  int2048& operator%=(const int2048& other) {
    int2048 q = *this / other;
    q *= other;
    *this -= q;
    return *this;
  }
  
  friend int2048 operator%(int2048 a, const int2048& b) {
    return a %= b;
  }

  friend std::istream& operator>>(std::istream& is, int2048& x) {
    std::string s;
    is >> s;
    x.read(s);
    return is;
  }
  
  friend std::ostream& operator<<(std::ostream& os, const int2048& x) {
    if (x.negative_ && !(x.digits_.size() == 1 && x.digits_[0] == 0)) {
      os << '-';
    }
    os << x.digits_.back();
    for (int i = (int)x.digits_.size() - 2; i >= 0; i--) {
      os << std::setfill('0') << std::setw(BASE_DIGITS) << x.digits_[i];
    }
    return os;
  }

  friend bool operator==(const int2048& a, const int2048& b) {
    return a.negative_ == b.negative_ && a.digits_ == b.digits_;
  }
  
  friend bool operator!=(const int2048& a, const int2048& b) {
    return !(a == b);
  }
  
  friend bool operator<(const int2048& a, const int2048& b) {
    if (a.negative_ != b.negative_) {
      return a.negative_;
    }
    int cmp = compare_abs(a.digits_, b.digits_);
    if (a.negative_) {
      return cmp > 0;
    }
    return cmp < 0;
  }
  
  friend bool operator>(const int2048& a, const int2048& b) {
    return b < a;
  }
  
  friend bool operator<=(const int2048& a, const int2048& b) {
    return !(a > b);
  }
  
  friend bool operator>=(const int2048& a, const int2048& b) {
    return !(a < b);
  }
};

} // namespace sjtu
