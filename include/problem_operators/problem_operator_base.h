#pragma once

namespace platypus
{

/**
 * Base ProblemOperator class. Required to solve issues with inheritance.
 *
 * TODO: - move the problem operator interface code into this class. Base this
 * on the hepahestus ProblemOperatorBase class (see Hephaestus PR #108).
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