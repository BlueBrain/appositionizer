#!/usr/bin/env python3

import sys
import argparse
import anytree
import numpy as np

import zirkel


def print_tree(tree):
    for pre, _, node in anytree.RenderTree(tree):
        indent = f"{pre}{node.region_name}"
        print(indent.ljust(50))


def print_header(desc_width, file=sys.stdout):
    print(
        " ".ljust(desc_width),
        *[s.rjust(10) for s in ["min [s]", "max [s]", "avg [s]", "avg [%]"]],
        file=file,
    )

def print_node(indent, node, key, desc_width, file=sys.stdout):
    values = [
        node.data_point(key + "_min"),
        node.data_point(key + "_max"),
        node.data_point(key + "_avg"),
        100 * node.data_point(key + "_avg_rel"),
    ]

    fmt = "{:10.3f}"
    print(indent.ljust(desc_width), *[fmt.format(v) for v in values], file=file)

def print_stats(tree, key, regions, file=sys.stdout):
    if regions is None:
        print_tree_stats(tree, key, file=file)
    else:
        print_selected_stats(tree, key, regions, file=file)

def print_tree_stats(tree, key, file=sys.stdout):
    print_header(40, file=file)
    for pre, _, node in anytree.RenderTree(tree):
        indent = f"{pre}{node.region_name}"
        print_node(indent, node, key, 40, file=file)


def print_selected_stats(tree, key, regions, file=sys.stdout):
    print_header(20, file=file)

    for region in regions:
        node = tree[region]
        print_node(region.split("/")[-1], node, key, 20, file=file)


def compute_simple_stats(node, total_time):
    def abs_rel(key, func):
        abs = func(time)
        node._data[f"{cat}_{key}"] = abs
        node._data[f"{cat}_{key}_rel"] = abs / total_time

    for cat in ["incl_time", "excl_time"]:
        time = node.data_point(cat)
        abs_rel("min", np.min)
        abs_rel("max", np.max)
        abs_rel("avg", np.mean)


def main():
    parser = argparse.ArgumentParser(description="Parse multi-threaded hatchet files.")
    parser.add_argument("--incl-time", action="store_true")
    parser.add_argument("--excl-time", action="store_true")
    parser.add_argument("--regions", nargs="*", default=None)
    parser.add_argument("profile_name")

    args = parser.parse_args()

    tree = zirkel.load_tree(args.profile_name, format="hatchet")
    tree.scan(zirkel.InclusiveTimeScan(excl_key="excl_time", incl_key="incl_time"))
    total_time = np.max(tree.data_point("incl_time"))
    regions = args.regions


    for node in anytree.PreOrderIter(tree):
        compute_simple_stats(node, total_time=total_time)

    if args.incl_time:
        print("Inclusive time:")
        print_stats(tree, "incl_time", regions)
        
    if args.excl_time:
        print("Exclusive time:")
        print_stats(tree, "excl_time", regions)


if __name__ == "__main__":
    main()
