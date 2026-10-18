picogl.backend.gl.enums.target.renderbuffer

You've found the cause: there are two different `bravais.h` headers with incompatible function signatures. The `index` version expects five arguments, while the `util` version expects six.

Your failing call passes six arguments:

```
bravais(&om, cl_om, crit, sys_id, output_type, error);
```

But `index/bravais.h` declares:

```
void bravais(
    struct matrix_3x3 *,
    struct matrix_3x3 *,
    float crit[45],
    char sys_id[45][3],
    int output_type
);
```

Meanwhile, `util/bravais.h` and `src/dps/bravais/bravais.h` declare the six-argument version.

## Recommended fix

The definition in `index/bravais.c` is the decisive evidence. Inspect it alongside the includes and Makefile:

```
sed -n '1,100p' index/bravais.c
sed -n '1,70p' index/bravais_fortran.c
sed -n '1,60p' index/bravais.h
grep -nE 'bravais_fortran|bravais\.o|CPPFLAGS|CFLAGS|INCLUDES' index/Makefile Makefile 2>/dev/null
```

If `index/bravais.c` defines the six-argument function, then `index/bravais.h` is stale or incorrect. Update its declaration to match the definition exactly.

If `index/bravais.c` defines the five-argument function, then the call in `bravais_fortran.c` is wrong: remove the final argument only if the caller does not need error reporting and the function's intended contract confirms that change.

Do not simply add the missing prototype to `bravais_fortran.c`. That would suppress the compiler error without resolving the incompatible declarations.

## Check the include paths, too

Because `bravais_fortran.c` includes `"bravais.h"`, the compiler generally searches the source file's own directory first. Therefore, in `index/bravais_fortran.c`, it is likely picking up `index/bravais.h`, not `util/bravais.h`.

You can verify which header is selected with:

```
clang -H -c index/bravais_fortran.c -o /tmp/bravais_fortran.o
```

The broader lesson for this port is to avoid maintaining multiple headers with the same basename and different interfaces. A shared header under a uniquely named path, such as `dps/bravais/bravais.h`, would make this class of error less likely. Before consolidating the headers, however, verify whether the five-argument and six-argument variants represent intentionally different APIs.