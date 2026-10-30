# Survey of the LAPACK routines bound

Every routine this package declares, classified before the binding was written: **expressible**
(declares directly), **shim-able** (needs C of its own), or **needs the language to move**.

## Which interface: Fortran, not LAPACKE

LAPACK is Fortran. Its Fortran symbols (`dgesv_`) take every argument by address — `integer` as a
32-bit `int *` under LP64, `character` as a `char *` with its length passed as a hidden trailing
`size_t` (gfortran) — and store matrices column-major with a leading dimension. LAPACKE is a C layer
over those routines: a `matrix_layout` flag, values by value, and on its row-major path a transposing
copy before and after the same Fortran routine.

**This package binds the Fortran routines, because nex's compiler does** (`nex/compiler/emit_lapack.sysl`
calls `dgesv_` and its siblings directly, by reference, hidden lengths of 1 appended, `lwork = -1`
workspace queries). The nex interpreter calls this package so that interpreted and compiled nex run
the *same routine on the same bytes* and answer the same bits; going through LAPACKE would put a
different code path — and, row-major, different data — between the two. The Fortran symbols are also
the one interface Accelerate and OpenBLAS both export under the same names.

**Accelerate: the classic interface, as nex links it.** `liblapack.tbd` in the macOS SDK exports both
`_dgesv_` (classic, LP64, f2c conventions) and `_dgesv$NEWLAPACK` (the `ACCELERATE_NEW_LAPACK`
interface, a different symbol). nex calls `dgesv_`, so this does: `@link("lapack")` resolves to
`/System/Library/Frameworks/Accelerate.framework/…/vecLib.framework/…/libLAPACK.dylib`, the library
`-framework Accelerate` reaches. `otool -L` on a test binary shows it.

## Result

| bucket | count | what |
|---|---|---|
| expressible | 64 routines + 3 functions | `?gesv ?getrf ?getrs ?getri ?gecon ?potrf ?potrs ?geqrf ?gels ?gelsd ?geev ?gesdd ?gesvd ?trtrs` (× s, d, c, z), `sorgqr dorgqr cungqr zungqr`, `ssyevd dsyevd cheevd zheevd`; `dlamch`, `dlange`, `zlange` |
| shim-able | 0 | — |
| needs the language to move | 0 | — |
| not declared, by the library's ABI | 3 | `slamch`, `slange`, `clange` — see below |

`nm -gU libopenblas.dylib` and the SDK's `liblapack.tbd` both export all 64 routines (the count was
taken by grep over each) and all six `?lamch` / `?lange`.

## The shapes, and how each declares

| shape | how |
|---|---|
| `INTEGER` scalar argument (n, lda, info, lwork, …) | `*i32`; a value is passed as `&int(n)`, storage made at the `&` |
| `CHARACTER*1` argument (trans, uplo, jobz, …) | `*u8`, and a trailing `usize` length per character argument, passed as 1 |
| `REAL` / `DOUBLE PRECISION` array or scalar | `*f32` / `*real` |
| `COMPLEX` / `COMPLEX*16` array | `*Complex[f32]` / `*Complex[real]` — `sysl.math.complex.Complex` is two fields of its width, Fortran's layout |
| `INTEGER` array (ipiv, iwork) | `*i32` |
| `DOUBLE PRECISION FUNCTION` (`dlamch`, `dlange`, `zlange`) | `-> real`: a `double` in both gfortran's and f2c's convention |
| workspace query (`lwork = -1`) | the routine called twice; the size read from `work(1)` (`Width.to_count`) |

## What is not declared, and why

**`slamch`, `slange`, `clange` — a `REAL FUNCTION` has two ABIs.** gfortran (OpenBLAS) returns it as
a C `float`; f2c (Accelerate's classic interface) returns it as a `double` — `clapack.h` declares
`__CLPK_doublereal slamch_(char *)`. Which library answers is decided at the link, not in the source,
so no one declaration is right for both, and a wrong one reads a register the callee did not write.
This is the library's ABI split, not a language gap. `lamch[f32]` answers `slamch`'s values from
`f32`'s IEEE 754 constants (checked against `f32.epsilon()` and `f32.max_value()` in the tests); a
single-precision norm is computed in sysl where one is needed (`rcond`).

## Language findings

None blocked the binding. Two were met and are reported to the sysl session:

1. **A `[]const`-to-`[]` mismatch on an argument whose callee's type parameter is solved from the
   expected type is not reported, and the enclosing test's body is dropped** — clang then fails on
   `use of undefined value '@…$<test name>'`. The test helper here takes `[]const real`, which is
   what it should have taken anyway. Reduction (sysl 0.8.2), as a `sysl run` file:

   ```sysl
   struct Box[T]
       v: T

   mk[T: Zero](xs: []real) -> Box[T] = Box(T.zero())

   same[T](a: Box[T], b: Box[T]) -> bool = true

   val TABLE: [2]real = [16.0, -120.0]

   print(same(Box(1.0), mk(TABLE)))
   ```

   answers `ld: symbol(s) not found … "_main"`; the same call inside a function compiles and runs
   with no diagnostic; inside a `@test` it is the clang error above. `first(1, TABLE)` against
   `first[T](n: T, xs: []real)` — `T` solved from a sibling argument — gets the proper refusal (*"a
   '[]const real' views elements it may not write …"*), which is what `arrays.md § []const T` says
   every one of these should get.
2. **`T::Real.zero()` cannot be written** in a generic body — `T::X` in expression position is
   reflection — but nothing needs it: a literal (`var w: []T::Real = [0.0; n]`) or inference from the
   declared type reaches the projection, so no contortion was written.
