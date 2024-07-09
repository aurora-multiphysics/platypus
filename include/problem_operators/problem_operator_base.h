#pragma once

namespace platypus
{

/**
 * Base ProblemOperator class. Required to solve issues with inheritance.
 */
class ProblemOperatorBase
{
protected:
  // NB: - protected constructor to only allow derived classes to be constructed.
  ProblemOperatorBase() = default;

  // Virtual destructor to prevent resource leaks.
  virtual ~ProblemOperatorBase() = default;
};

} // platypus