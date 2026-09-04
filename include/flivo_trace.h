// Determinism trace (diagnostic only; enabled by setting FLIVO_TRACE=<path>).
// Dumps full-precision state and data hashes at pipeline boundaries so two runs
// can be diffed to find the first stage that diverges. When FLIVO_TRACE is
// unset every entry point returns immediately and behaviour is unchanged.
#ifndef FLIVO_TRACE_H_
#define FLIVO_TRACE_H_

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>

namespace flivo_trace
{

inline std::ofstream &stream()
{
  static std::ofstream f;
  static bool opened = false;
  if (!opened)
  {
    const char *p = getenv("FLIVO_TRACE");
    if (p) f.open(p, std::ios::out);
    opened = true;
  }
  return f;
}

inline bool on() { return stream().is_open(); }

// FNV-1a over raw bytes: order-sensitive, so reordering is caught too.
inline uint64_t hashBytes(const void *data, size_t n, uint64_t h = 1469598103934665603ULL)
{
  const unsigned char *p = static_cast<const unsigned char *>(data);
  for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= 1099511628211ULL; }
  return h;
}

inline void raw(const char *fmt_tag, int a, int b, double v1, double v2, uint64_t h)
{
  if (!on()) return;
  std::ofstream &f = stream();
  f << a << " " << fmt_tag << " i=" << b
    << std::scientific << std::setprecision(17)
    << " v1=" << v1 << " v2=" << v2 << " h=" << h << "\n";
  f.unsetf(std::ios::scientific);
}

} // namespace flivo_trace

#endif // FLIVO_TRACE_H_
