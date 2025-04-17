import matplotlib.pyplot as plt
from matplotlib import ticker
from Bio import Phylo
import seaborn as sns
import numpy as np


class TreePlot:
    def plotTree(self,
                 figsize=None,
                 tree_context={'lines.linewidth': 0.75},
                 leaf_ext_params={'c': 'gray', 'linestyle': 'dotted'},
                 strain_sheet=None,
                 lineage_lw=2):
        with plt.rc_context(tree_context):
            # generate the axis
            self.fig, ax = plt.subplots(figsize=figsize)
            # plot the tree
            Phylo.draw(self.tree,
                       axes=ax, show_confidence=False,
                       label_func=lambda x: None, do_show=False)
            # rescale x limits
            ax.set_xlim(ax.get_xlim()[0], np.max(self.x_ends))
            # draw dotted lines to line up with x max
            for x, y in zip(self.x_ends, self.y_ends):
                ax.plot([ax.get_xlim()[1], x], [y, y], **leaf_ext_params)
            # draw lineage lines if any are provided
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            if strain_sheet is not None:
                for leaf, y in zip(self.labels, self.y_ends):
                    ax.plot(
                        [xlim[1], xlim[1]],
                        [y - 0.505, y + 0.505],
                        c=strain_sheet.loc[leaf, 'lineage_color'],
                        lw=lineage_lw,
                        solid_capstyle='butt')
            ax.set_xlim(xlim[0], xlim[1] + (xlim[1] / 100))
            ax.set_ylim(ylim)
            # clean up the axes
            sns.despine(ax=ax, left=True, bottom=True)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            # save local axis as tree ax, also add to list
            self.tree_ax = ax
            self.axes.append(ax)

    def addParasiteAxis(self, width, visible=False):
        ax_pos = self.axes[-1].get_position()
        if len(self.axes) == 1:
            offset = self.offset_first
        else:
            offset = self.offset
        # add axis
        parasite_ax = self.fig.add_axes((ax_pos.x0 + ax_pos.width + offset,
                                    ax_pos.y0,
                                    width,
                                    ax_pos.height))
        if not visible:
            sns.despine(ax=parasite_ax, left=True, bottom=True)
            parasite_ax.set_xticks([])
            parasite_ax.set_yticks([])
        self.axes.append(parasite_ax)
        return parasite_ax

    def addHeatmap(self, width, input_data, y_offsets=(0.25, -0.3), **kwargs):
        # add parasite
        ax = self.addParasiteAxis(width)
        # plot the heatmap
        sns.heatmap(input_data.loc[self.labels, :], ax=ax, **kwargs)
        # add spacing
        ax.set_ylim(ax.get_ylim()[0] + y_offsets[0],
                    ax.get_ylim()[1] + y_offsets[1])
        # remove ticks
        ax.set_yticks([])
        return ax  # return axis for additional modification

    def addMeanPlot(
            self, width, input_data,
            ax=None,
            point_kwargs={'c': 'k', 's': 2, 'zorder': 2},
            line_origin=0,
            line_kwargs={'c': 'k', 'lw': 0.75, 'zorder': 1},
            show_grid=False,
            grid_lw=0.5):
        """
        input data: pandas Series with index as strain names (matching tree).
        """
        # add parasite
        if ax is None:
            ax = self.addParasiteAxis(width, visible=True)
        # plot datapoints
        if point_kwargs is not None:
            point_y = []
            for label in input_data.keys():
                point_y.append(self.leaf_xy_dict[label][1])
            ax.scatter(input_data.values, point_y, **point_kwargs)
        # plot line
        if line_kwargs is not None:
            for label in np.unique(input_data.index):
                datapoint = input_data.loc[label]
                y_value = self.leaf_xy_dict[label][1]
                ax.plot(
                    [line_origin, datapoint],
                    [y_value, y_value], **line_kwargs)
        # axis visualization options
        ax.set_yticks(np.arange(len(self.labels)))
        if show_grid:
            ax.grid(True, which='major', axis='y', lw=grid_lw, color='silver', linestyle='dotted')
            ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
            ax.grid(True, which='minor', axis='y', lw=grid_lw, color='gray', linestyle='dotted')
            ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
        ax.axes.yaxis.set_ticklabels([], minor=False)
        ax.axes.yaxis.set_ticklabels([], minor=True)
        ax.tick_params(axis='y', length=0, which='both')
        ax.xaxis.set_label_position('top')
        sns.despine(ax=ax, left=True, bottom=True, top=False)
        ax.tick_params(axis='x', top=True)
        ax.set_ylim(self.tree_ax.get_ylim())

    def get_x_positions(self):
        depths = self.tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = self.tree.depths(unit_branch_lengths=True)
        return depths

    def get_y_positions(self):
        maxheight = self.tree.count_terminals()
        # Rows are defined by the tips
        heights = {
            tip: maxheight - i for i, tip in enumerate(reversed(self.tree.get_terminals()))}

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (
                heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2

        if self.tree.root.clades:
            calc_row(self.tree.root)
        return heights

    def __init__(self, phylo_obj=None, root_name=None, offset_first=0.18, offset=0.04):
        self.tree = phylo_obj
        # get x and y information for all leaves
        leaf_x_dict = {}
        for clade, xpos in self.get_x_positions().items():
            if clade.name is not None:
                leaf_x_dict[clade.name] = xpos
        leaf_y_dict = {}
        for clade, ypos in self.get_y_positions().items():
            if clade.name is not None:
                leaf_y_dict[clade.name] = ypos
        self.leaf_xy_dict = {}
        for name, x in leaf_x_dict.items():
            self.leaf_xy_dict[name] = (x, leaf_y_dict[name])
        # remove root from dict so it is not plotted
        if root_name is not None:
            try:
                self.leaf_xy_dict.pop(root_name)
            except KeyError:
                print('warning, root not found')
        # also store raw label/x/y values
        self.labels = list(self.leaf_xy_dict.keys())
        self.x_ends, self.y_ends = np.asarray(
            [[x, y] for x, y in self.leaf_xy_dict.values()]).T
        # dictionary for data axes
        self.axes = []
        # defaults
        self.offset_first = offset_first
        self.offset = offset


def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format

    from: https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick