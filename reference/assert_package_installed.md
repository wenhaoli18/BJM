# Assert that an optional (Suggests-only) package is installed

Shared input-validation helper for functions that depend on a package
listed in `Suggests` rather than `Imports` (so that most users – and
CRAN's own checks – are unaffected by a large, optional dependency they
never call). Every call site into that package must still be fully
namespaced (`pkg::fun()`), never `@importFrom`; this helper just turns a
missing package into a clear, actionable error instead of a cryptic
"could not find function" failure deep inside the calling code.

## Usage

``` r
assert_package_installed(pkg, context)
```
