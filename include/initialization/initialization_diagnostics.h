#pragma once
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

// Observations only: no thresholds, admission decisions or estimator changes.
struct InitializationDiagnostics
{
  bool enabled = true;
  std::ofstream log;

  void initializeLog(const std::string &path)
  {
    if (!enabled) return;
    log.open(path);
    if (!log) throw std::runtime_error("Cannot create initialization diagnostics log");
    log << "timestamp,stage,matches,points,match_ratio,residual_m,translation_correction_m,"
           "rotation_correction_rad,velocity_correction_mps,gravity_correction_mps2\n";
  }

  void record(double time, const char *stage, size_t matches, size_t points,
              double residual, double translation, double rotation,
              double velocity, double gravity)
  {
    if (!enabled) return;
    // Replayed timestamps can precede the initialization observations. Missing
    // matches and non-finite diagnostics are recorded without rejecting data.
    log << std::setprecision(17) << time << ',' << stage << ',' << matches << ',' << points << ','
        << (points ? double(matches) / points : 0) << ',' << residual << ',' << translation << ','
        << rotation << ',' << velocity << ',' << gravity << '\n' << std::flush;
  }
};
