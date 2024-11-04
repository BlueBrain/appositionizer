#!/usr/bin/env python3

import anytree
import numpy as np

import zirkel


def print_tree(tree):
    for pre, _, node in anytree.RenderTree(tree):
        indent = f"{pre}{node.region_name}"
        print(indent.ljust(50))


tree = zirkel.load_tree("caliper.json", format="mpi_report")
zirkel.print_data(tree, "inclusive_time_rank_avg")
