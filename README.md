# lapack

Dense linear algebra for sysl — solving, least squares, inverse, determinant, LU, Cholesky, QR,
eigenvalues and the singular value decomposition — bound to [LAPACK](https://www.netlib.org/lapack/)
through its Fortran interface, for all four of its element types.

*An independent sysl package for LAPACK. It is not affiliated with or endorsed by LAPACK's authors
(the University of Tennessee, the University of California Berkeley, the University of Colorado
Denver and NAG), Netlib, Apple (Accelerate) or the OpenBLAS project.*

```sysl
import sh.sysl.lapack.{eigh, mat, solve, svdvals, vector}

var a = [4.0, 1.0, 1.0, 3.0]                  // column-major: [[4, 1], [1, 3]]
val m = mat(a, 2, 2)

print(solve(m, vector([1.0, 2.0])).unwrap().data)   // [0.0909…, 0.6363…]
print(eigh(m).unwrap().values)                      // [2.381…, 4.618…]: (7 ∓ √5) / 2
print(svdvals(m).unwrap())                          // the same two, descending
```

```
dependencies {
  lapack { git = "github.com/sysl-lang/lapack", version = "0.1.0" }
}
```

The coordinate names the package, `lapack`; the module a program imports is `sh.sysl.lapack`.

## Which LAPACK

**macOS: Accelerate, with nothing to install.** `@link("lapack")` becomes `-llapack`, which the macOS
SDK answers with Accelerate's `libLAPACK.dylib` — its classic interface (`dgesv_`, 32-bit integers).

**OpenBLAS, on macOS or anywhere:** name its directory at the build.

```
brew install openblas
sysl test . --link-path /opt/homebrew/opt/openblas/lib
```

Homebrew's `liblapack.dylib` there is OpenBLAS (a link to `libopenblas.dylib`), and `--link-path`
also writes the run-time path. **On Linux**, `-llapack` is the distribution's LAPACK — install
`liblapack-dev` or `libopenblas-dev`. `otool -L` / `ldd` on a binary says which one it got.

The suite passes against both Accelerate and OpenBLAS.

## What it is

Two layers. **`sh.sysl.lapack.c`** is LAPACK's Fortran interface, one declaration per routine per
type, every argument an address and each `character` argument's hidden length explicit.
**`sh.sysl.lapack`** is what a program imports:

| | |
|---|---|
| **`Mat[T]`** | a column-major view: `data: []T`, `rows`, `cols`, `ld` (the distance between columns). `mat(data, m, n)`, `strided(data, m, n, ld)`, `from_ptr(p, m, n, ld)` (storage C or an interpreter owns), `vector(xs)`, `zeros(m, n)`; `at`, `set`, `column`. **Nothing is copied to make one**, so a block of a larger matrix, a `Buf`'s storage or a nex array is handed over where it lies |
| solving | `solve(A, B)` (`?gesv`), `lstsq(A, B)` (`?gelsd`, minimum norm, any rank), `solve_triangular(A, B, upper)` (`?trtrs`), `lu_solve(lu, B)` (`?getrs`), `chol_solve(L, B, upper)` (`?potrs`) |
| matrices | `inv(A)` (`?getrf` + `?getri`), `det(A)` (`?getrf`; exactly 0 for a zero pivot), `rcond(A)` (`?gecon`, 1-norm) |
| factorizations | `lu(A)` → `Lu { l, u, perm, factors, ipiv }` with `A[perm[i], :] == (L U)[i, :]`; `chol(A, upper)`; `qr(A)` → thin `Qr { q, r }` (`?geqrf` + `?orgqr`/`?ungqr`) |
| eigenvalues | `eigvals(A)`, `eig(A)` → `Eig { values, vectors }`, always complex (`?geev`); `eigh(A)` → `Eigh { values, vectors }`, real ascending values (`?syevd`/`?heevd`) |
| SVD | `svd(A)` → thin `Svd { u, s, v }` with A = U diag(s) Vᴴ; `svdvals(A)` (`?gesdd`) |
| in place | `solve_inplace(A, ipiv, B)`, `lu_inplace(A, ipiv)`, `chol_inplace(A, upper)`, `qr_inplace(A, tau)`: work in the caller's storage, answer LAPACK's `info` |
| `Scalar` | the trait over `f32`, `real`, `Complex[f32]`, `Complex[real]` that every function is generic over; `T.gesv(…)`, `T.gesvd(…)` … reach each routine at the right letter |
| `lamch[F](c)` | LAPACK's machine parameters |

**These are the functions of nex's linear-algebra prelude, and each makes the call nex's compiler
emits** — the same routine, its input copied into the same dense layout, the workspace sized by the
same query or formula — so the nex interpreter, which calls this package, answers bit for bit what
compiled nex answers on the same LAPACK.

**Every failure is a `Result`.** A value form answers `Err(Error)` naming the function, the typed
routine and the cause, from LAPACK's `info`: *"solve: dgesv: matrix is singular (U(2,2) = 0)"*,
*"chol: dpotrf: matrix is not positive definite (leading minor 2)"*, *"eigh: dsyevd: failed to
converge (info 3)"*, and a shape that does not fit before LAPACK is called (*"solve: the matrix is
not square (2 × 3)"*). `Error.info()` gives LAPACK's number back. The value forms never write their
arguments; workspace queries are made internally.

## What is deliberately absent

- **LAPACKE.** The Fortran routines are what nex calls, and bit-identity with it is the point; see
  `SURVEY.md`.
- **`slamch`, `slange`, `clange`.** A Fortran `REAL` function returns a `float` from gfortran
  (OpenBLAS) and a `double` from f2c (Accelerate's classic interface), and the link — not the source —
  decides which library answers, so no single declaration is right. `lamch[f32]` answers `slamch`'s
  values from IEEE 754 constants; `dlamch`, `dlange`, `zlange` return `DOUBLE PRECISION` everywhere
  and are bound.
- **Pleasant forms of `?gels` and `?gesvd`.** `lstsq` and `svd` use `?gelsd` and `?gesdd`, as nex
  does; the other two are reachable as `T.gels` / `T.gesvd` and tested there.
- **Banded, tridiagonal, packed and generalized problems** — not in nex's first set either.
- **The ILP64 interface** (`$NEWLAPACK$ILP64`, OpenBLAS's `64_` builds): matrices past 2³¹ − 1 rows
  are refused as a shape error instead.

## Testing

```
sysl test .
sysl test . --link-path /opt/homebrew/opt/openblas/lib
SYSL_EXTRA_CFLAGS="-fsanitize=address -g" sysl test .
```

The oracles are closed forms, never LAPACK's own answer read back: systems built from a chosen
solution, the 4 × 4 Hilbert matrix's integer inverse and its determinant 1/6048000, a Vandermonde
determinant, a rotation's eigenvalues cos θ ± i sin θ, the second-difference matrix's spectrum
2 − 2 cos(kπ/(n+1)), a rank-one matrix's single singular value, a line through collinear points; and
residuals recomputed in plain sysl — ‖A X − B‖, ‖A − Q R‖, ‖A V − V Λ‖, ‖A − U Σ Vᴴ‖, QᴴQ = I.

LAPACK is BSD-licensed (modified BSD); this binding is ISC.
