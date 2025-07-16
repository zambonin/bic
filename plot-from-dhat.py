#!/usr/bin/env uv run
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "matplotlib",
# ]
# ///
from argparse import ArgumentParser
from io import BytesIO
from json import loads, JSONDecodeError
from re import findall
from sys import stdin, exit, stdout

from matplotlib import pyplot as plt
from matplotlib import colors


class CacheAccessMatrix:
    def __init__(self, input_lines: list[str]):
        if not input_lines:
            exit(0)

        header_line = input_lines[0]
        json_string = "".join(input_lines[1:])

        self._parse_config(header_line)

        try:
            self.data = loads(json_string)
        except JSONDecodeError:
            exit(1)

        self.matrix = self._process_data()

    def _parse_config(self, header_line: str) -> None:
        params = {
            k.strip(): v for k, v in findall(r"(\w+)\s*=\s*([\d.]+)", header_line)
        }
        try:
            n_base = int(params["n"])
            k_base = int(params["k"])
            self.bit_length = int(params["b"]) // 8
            cache_type = int(params["c"])
        except (KeyError, ValueError):
            exit(1)

        if cache_type == 0:
            exit(0)
        elif cache_type == 1:
            self.target_function = "bin_build_cache"
            self.n, self.k = n_base + k_base + 1, k_base
        elif cache_type == 2:
            self.target_function = "comb_build_cache"
            self.n, self.k = n_base + 1, k_base + 1
        else:
            exit(1)

    def _decode_acc(self, acc_list: list[int]) -> list[int]:
        decoded = []
        i = 0
        while i < len(acc_list):
            val = acc_list[i]
            if val < 0:
                run_count = -val
                i += 1
                if i >= len(acc_list):
                    break
                value_to_repeat = acc_list[i]
                if self.bit_length > 0 and run_count % self.bit_length == 0:
                    new_run_count = run_count // self.bit_length
                    decoded.extend([value_to_repeat] * new_run_count)
                else:
                    decoded.extend([value_to_repeat] * run_count)
            else:
                decoded.append(val)
            i += 1
        return decoded

    def _process_data(self) -> list[list[int]]:
        all_accesses = []
        frame_table = self.data.get("ftbl", [])

        for pp in self.data.get("pps", []):
            if self._check_relevance(pp, frame_table):
                if "acc" in pp:
                    all_accesses = self._decode_acc(pp["acc"])
                    break

        if not all_accesses:
            return []

        total_elements = self.n * self.k
        final_list = (all_accesses + [0] * total_elements)[:total_elements]
        return [final_list[i * self.k : (i + 1) * self.k] for i in range(self.n)]

    def _check_relevance(self, pp: dict[str, any], frame_table: list[str]) -> bool:
        return any(
            self.target_function in frame_table[frame_index]
            for frame_index in pp.get("fs", [])
            if 0 <= frame_index < len(frame_table)
        )

    def plot(self) -> None:
        if not self.matrix:
            return

        buffer = BytesIO()
        plot_data = [[(val if val > 0 else 1e-9) for val in row] for row in self.matrix]
        cmap = plt.get_cmap("hot").copy()
        cmap.set_under(color="black")
        aspect_ratio = (self.k / self.n) if self.n > 0 else "auto"

        plt.imshow(
            plot_data,
            cmap=cmap,
            interpolation="nearest",
            norm=colors.LogNorm(vmin=1),
            aspect=aspect_ratio,
        )
        plt.axis("off")
        plt.savefig(buffer, format="png", bbox_inches="tight", pad_inches=0, dpi=150)

        buffer.seek(0)
        stdout.buffer.write(buffer.getvalue())

    def __repr__(self) -> str:
        if not self.matrix:
            return ""
        return "\n".join(" ".join(f"{num:5d}" for num in row) for row in self.matrix)


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--print", action="store_true")
    args = parser.parse_args()

    try:
        processor = CacheAccessMatrix(stdin.readlines())
    except Exception:
        exit(1)

    if args.print:
        print(processor)
    else:
        processor.plot()


if __name__ == "__main__":
    main()
