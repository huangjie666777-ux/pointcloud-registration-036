#include "icp3d/types.h"

namespace icp3d {

const char *toString(TerminationReason reason) {
  switch (reason) {
  case TerminationReason::Converged:
    return "Converged";
  case TerminationReason::MaxIterationsReached:
    return "MaxIterationsReached";
  case TerminationReason::InsufficientCorrespondences:
    return "InsufficientCorrespondences";
  case TerminationReason::DegenerateGeometry:
    return "DegenerateGeometry";
  }
  return "Unknown";
}

} // namespace icp3d
