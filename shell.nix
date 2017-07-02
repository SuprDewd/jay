with import <nixpkgs> {}; {
  env = stdenv.mkDerivation {
    name = "jay_env";
    buildInputs = [ stdenv gcc gmp gmpxx ];
    C_INCLUDE_PATH="${gmp.dev}/include:${gmpxx.dev}/include";
  };
}
