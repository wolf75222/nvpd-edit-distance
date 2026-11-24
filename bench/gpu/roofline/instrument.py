"""Insert a comparison dispatcher into the exact archived CUDA baseline."""

import argparse
import hashlib
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    source = args.source.read_bytes()
    if (
        hashlib.sha256(source).hexdigest()
        != "32b153293cb01d2c503ac947e0a1a8fd817561dabc77d92732abfe0aafe4e090"
    ):
        raise ValueError("expected b5d08d1:src/cuda/cuda.cu")
    text = '#include "paired_rows.cuh"\nint roof_variant = 0;\n' + source.decode()
    marker = "    } else {\n        for (int i = 1; i <= m; ++i) {"
    dispatch = """    } else if (!trace && m >= 2 && roof_variant != 0) {
        int* curr = table.data + stride;
        int i = 1;
        for (; i < m; i += 2) {
            const int* a1 = indices.data + matches.symbol[static_cast<unsigned char>(a[i-1])] * stride;
            const int* a2 = indices.data + matches.symbol[static_cast<unsigned char>(a[i])] * stride;
            if (roof_variant == 1)
                detail::pair_rows<false><<<blocks, threads>>>(prev, curr, a1, a2, i, n);
            else
                detail::pair_rows<true><<<blocks, threads>>>(prev, curr, a1, a2, i, n);
            check(cudaGetLastError(), "launch paired rows");
            std::swap(prev, curr);
        }
        if (i <= m) {
            const int s = matches.symbol[static_cast<unsigned char>(a[i-1])];
            next_row<<<blocks, threads>>>(prev, curr, indices.data + s * stride, i, n);
            check(cudaGetLastError(), "launch odd final row");
            prev = curr;
        }
"""
    if text.count(marker) != 1:
        raise ValueError("missing comparison insertion point")
    args.out.write_text(text.replace(marker, dispatch + marker))


if __name__ == "__main__":
    main()
