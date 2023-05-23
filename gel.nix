{ nixpkgs, ... }:
let
  pkgs = import nixpkgs {
    system = "x86_64-linux";
  };
in
  with pkgs;
  stdenv.mkDerivation {
    name = "gel";
    version = "0.0.2";

  # https://nixos.org/nix/manual/#builtin-filterSource
  src = builtins.filterSource
  (path: type: lib.cleanSourceFilter path type
  && baseNameOf path != "doc/*"
  && baseNameOf path != "data/*"
  && baseNameOf path != ".github/*") ./.;

  cmakeFlags =["-DUse_GLGraphics=OFF"];

  enableParallelBuilding = true;
  nativeBuildInputs = [ cmake pkg-config ];

  # for optional opengl related dependencies:
  #buildInputs = [
    #libGL glfw3 x11 libGLU xorg.libXdmcp
  #];

}
