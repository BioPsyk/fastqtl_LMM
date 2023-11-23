{ pkgs }:

with pkgs; mkShell {
  nativeBuildInputs = [
    cmake
    pkg-config
  ];

  buildInputs = [
    bzip2 zlib lzma
  ];
}
