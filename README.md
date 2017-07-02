# jay

A simple utility that, given some initial terms of an unknown generating
function `F`, conjectures a functional equation that `F` satisfies. Inspired by
[this talk](http://jaypantone.com/talks/2017-02-24-GA-Tech-Experimental.pdf) by
[Jay Pantone](http://jaypantone.com/).

## Synopsis

```
Usage: jay [OPTION]...
Given initial terms of an unknown generating function F,
conjecture a functional equation that F satisfies.

  -D,  --diff        the highest derivative of F to appear (default: 0)
  -fd, --fdeg        the highest power of F to appear (default: 1)
  -d,  --deg         the highest powers of variables to appear
  -t,  --type        the type of the generating function: ogf, egf (default: ogf)
       --hops        generate output that HOPS understands
  -h,  --help        display this help and exit
  -q,  --quiet       don't log any unnecessary output

The program expects the initial terms to be specified on standard
input in the following format. On the first line an integer 'n',
denoting the number of variables in the generating function. On
each of the following lines, until EOF, one initial term must be
specified in the format 'i_1 i_2 ... i_n c_{i_1, i_2, ..., i_n}',
meaning that 'c_{i_1, i_2, ..., i_n}' is the coefficient of
'x_1^{i_1} x_2^{i_2} ... x_n^{i_n}' in F. The list of initial
terms must be exhaustive, meaning that if a coefficient is listed,
then all coefficients with smaller indices must be listed as
well.
```

## Installing

```
$ ./autogen.sh
$ ./configure
$ make
# make install
```

Or, using the Nix package manager:
```
$ nix-env -f ./default.nix -i
```

## License

MIT: see the [LICENSE](https://github.com/SuprDewd/jay/blob/master/LICENSE)
file.

