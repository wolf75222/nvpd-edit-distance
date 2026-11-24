#!/usr/bin/env bash
# Compare generated instructions with the pinned sources at identical compiler paths.
set -euo pipefail
build=$(cd "${1:?provide a configured CMake build}" && pwd)
upstream=$build/_deps/wfa2-src
tmp=$(mktemp -d "${TMPDIR:-/tmp}/wfa-objects.XXXXXXXX")
trap 'rm -rf -- "$tmp"' EXIT
cd "$tmp"
units=(wavefront_backtrace wavefront_pcigar)
case "$(uname -s)" in
    Darwin) section=__TEXT,__text; target=apple-arm64; flags=(-arch arm64) ;;
    Linux) section=.text; target=epyc-x86_64; flags=(-march=native); units+=(wavefront_extend_kernels_avx) ;;
    *) exit 2 ;;
esac
hash() {
    if command -v sha256sum >/dev/null; then sha256sum "$1" | cut -d' ' -f1
    else shasum -a 256 "$1" | cut -d' ' -f1; fi
}
echo "target,unit,variant,text_sha256,disassembly_relocations_sha256"
for unit in "${units[@]}"; do
    for variant in upstream portable; do
        if [[ $variant == upstream ]]; then source=$upstream/wavefront
        else source=$build/wfa2-portable; fi
        cp "$source/$unit.c" "$unit.c"
        "${CC:-cc}" -O3 -DNDEBUG "${flags[@]}" -I"$upstream" -I"$upstream/wavefront" \
            -I"$upstream/system" -c "$unit.c" -o "$unit.o"
        "${OBJCOPY:-objcopy}" --dump-section "$section=bytes" "$unit.o" copy.o
        "${OBJDUMP:-objdump}" -dr "$unit.o" > asm
        echo "$target,$unit,$variant,$(hash bytes),$(hash asm)"
    done
done
