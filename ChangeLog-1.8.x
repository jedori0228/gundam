Fixes relative to 1.8.1

Issue #499 : Update to a new version of the cpp-generic-toolkit submodule. The cpp-generic-toolkit fix sets the branch status to enable branches that are used for the data.

Issue #494 : Rearrange the application of the Fritschle-Carlson criteria so that it can be tested.  This also applies a fix to Fritschle-Carlson so that the end of the spline are handled in a more reasonable fashion.

Issue #493 : Add tests for ComputeGeneralSplines to the validation suite.  As a special bonus it also adds tests for CalculateUniformSpline.h, and resolves unused variable warnings for both those functions.

Fixes relative to 1.8.0

Issue #485 : A job will correctly continue when the likelihood returns an infinite value. The job will stop if the likelihood returns a NaN.

Issue #486 : Fix the Catmull-Rom splines so that they have symmetric behavior. This changes the extrapolation behavior for both Catmull-Rom and Catmull-Rom,monotonic splines. The Catmull-Rom monotonic splines are updated to use the full Fritsche-Carlson criteria so that the interpolation is smoother.