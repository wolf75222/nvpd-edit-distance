"""Report file-page residency as seen by Linux mincore, without touching FASTA data."""
import argparse
import csv
import ctypes
import mmap
import os
from pathlib import Path

libc = ctypes.CDLL(None, use_errno=True)
libc.mmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_long]
libc.mmap.restype = ctypes.c_void_p
libc.mincore.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_ubyte)]
libc.munmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t]


def resident_pages(path, page_size):
    with path.open("rb") as stream:
        length = os.fstat(stream.fileno()).st_size
        pages = (length + page_size - 1) // page_size
        address = libc.mmap(None, length, 0, mmap.MAP_PRIVATE, stream.fileno(), 0)
        if address == ctypes.c_void_p(-1).value:
            raise OSError(ctypes.get_errno(), "mmap", path)
        try:
            state = (ctypes.c_ubyte * pages)()
            if libc.mincore(address, length, state) == -1:
                raise OSError(ctypes.get_errno(), "mincore", path)
            return length, pages, sum(value & 1 for value in state)
        finally:
            libc.munmap(address, length)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("phase", choices=("initial", "after_hash", "after_trials"))
    args = parser.parse_args()
    page_size = os.sysconf("SC_PAGESIZE")
    with args.manifest.open() as inputs, args.output.open("a", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        if output.tell() == 0:
            writer.writerow(("phase", "file", "bytes", "page_bytes", "pages", "resident_pages", "status"))
        for row in csv.reader(inputs, delimiter="\t"):
            path = Path(row[0])
            # Linux can mask page residency when the caller neither owns nor can write the file.
            if path.stat().st_uid != os.geteuid() and not os.access(path, os.W_OK):
                writer.writerow((args.phase, path.name, path.stat().st_size, page_size, "", "", "access-limited"))
                continue
            try:
                length, pages, resident = resident_pages(path, page_size)
            except OSError as error:
                writer.writerow((args.phase, path.name, path.stat().st_size, page_size, "", "", f"errno={error.errno}"))
            else:
                writer.writerow((args.phase, path.name, length, page_size, pages, resident, "ok"))


if __name__ == "__main__":
    main()
