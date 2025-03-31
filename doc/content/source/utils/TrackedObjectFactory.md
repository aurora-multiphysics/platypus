# TrackedObjectFactory

## Summary

`TrackedObjectFactory` is a factory class used to create objects to in
Platypus and retaining a shared pointer to them.

## Overview

`TrackedObjectFactory` is used to create and manage tracked objects in Platypus. Iterator methods are
provided to allow access to all of the objects that have been created.

`TrackedObjectFactory` objects are used to create and store `mfem::Coefficient`, `mfem::VectorCoefficient`,
and `mfem::MatrixCoefficient` derived objects added to the MFEM problem.

End users should not usually need to interact with the `TrackedObjectFactory` directly. Developers wanting
to add new objects tracked by the `TrackedObjectFactory` should do so using the `TrackedObjectFactory::make`
method.
