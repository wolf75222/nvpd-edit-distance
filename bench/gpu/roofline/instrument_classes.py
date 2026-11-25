"""Build the isolated two-class A/B translation from the production single-class kernel."""

from pathlib import Path
from shutil import copyfile
import sys


def replace_one(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"unexpected source around {old}")
    return text.replace(old, new)


def main():
    source, out = map(Path, sys.argv[1:])
    header = (source / "local_tiles.cuh").read_text()
    header = replace_one(header, "void local_row_tile(", "void class_row_tile(")
    header = replace_one(
        header,
        "const int begin = static_cast<int>(blockIdx.x) * columns;",
        "const int begin = starts[blockIdx.x * local_tile_rows + local_tile_rows - 1];",
    )
    (out / "local_tiles.cuh").write_text(header)
    # The packed kernel also includes local_tiles.cuh. Both must resolve
    # the same instrumented header, even when packing is disabled.
    copyfile(source / "packed_match.cuh", out / "packed_match.cuh")
    text = (source / "cuda.cu").read_text()
    text = replace_one(
        text, '#include "local_tiles.cuh"', '#include "local_tiles.cuh"\n#include "width_classes.cuh"'
    )
    text = replace_one(
        text,
        "    auto& indices = indices_;",
        "    const auto groups = detail::group_widths(plan,n,columns);\n    auto& indices = indices_;",
    )
    old = """                const auto bytes = 2 * plan.pitch * sizeof(int);
                detail::local_row_tile<<<tile_blocks, threads, bytes>>>(prev, curr, indices.data,
                    boundaries_.data, symbols, i, count, n, columns, plan.pitch);"""
    new = """                const auto small_bytes=2*groups.small_pitch*sizeof(int);
                detail::class_row_tile<<<groups.small_count,threads,small_bytes>>>(prev,curr,indices.data,
                    boundaries_.data,symbols,i,count,n,columns,groups.small_pitch);
                check(cudaGetLastError(), "launch small halo class");
                if (tile_blocks > groups.small_count) {
                    const auto large_bytes=2*groups.large_pitch*sizeof(int);
                    detail::class_row_tile<<<tile_blocks-groups.small_count,threads,large_bytes>>>(prev,curr,indices.data,
                        boundaries_.data+groups.small_count*tile_rows,symbols,i,count,n,columns,groups.large_pitch);
                }"""
    (out / "class.cu").write_text(replace_one(text, old, new))


if __name__ == "__main__":
    main()
