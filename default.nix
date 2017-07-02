with import <nixpkgs> {};
stdenv.mkDerivation {
  name = "jay";
  version = "0.1";
  src = ./.;
  buildInputs = [ stdenv gcc gmp gmpxx autoreconfHook ];
  C_INCLUDE_PATH="${gmp.dev}/include:${gmpxx.dev}/include";
}
